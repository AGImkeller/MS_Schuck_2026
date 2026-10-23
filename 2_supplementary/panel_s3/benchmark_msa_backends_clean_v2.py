#!/usr/bin/env python3
"""
Run a LongAIRR group-wise MSA backend benchmark from a LongAIRR
2_subsampled_groups.fasta file to final consensus FASTA.

This script keeps the LongAIRR/pRESTO downstream consensus step unchanged and
only swaps the per-group MSA backend:

  muscle38  -> MUSCLE v3.8-style stdin/stdout alignment
  muscle5   -> MUSCLE v5 -align/-super5 via /dev/stdin and /dev/stdout
  mafft     -> MAFFT via /dev/stdin, aligned FASTA on stdout
  abpoa     -> abPOA row-column MSA output (-r 1) via /dev/stdin

Workflow:
  input FASTA -> longairr_split_align.py -> chunk_*.fasta
  chunk_*.fasta -> group-wise MSA -> *_align-pass.fasta
  concatenate aligned chunks -> seqkit rmdup -n
  BuildConsensus.py -> longairr_header_metadata.py

The same chunking, grouping, header restoration, uppercase normalization,
parallelization and BuildConsensus command are used for all backends.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from collections import OrderedDict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from io import StringIO
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class BenchmarkError(RuntimeError):
    """Controlled benchmark error with a user-facing message."""


def log(msg: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


def run_command(
    cmd: Sequence[str],
    log_file: Path,
    stdout_path: Optional[Path] = None,
    cwd: Optional[Path] = None,
) -> None:
    """Run a command and append command/stdout/stderr to benchmark.log."""
    with log_file.open("a") as handle:
        handle.write("\n$ " + " ".join(map(shlex.quote, cmd)) + "\n")
        handle.flush()

        if stdout_path is None:
            result = subprocess.run(
                list(cmd),
                cwd=str(cwd) if cwd else None,
                stdout=handle,
                stderr=subprocess.STDOUT,
                text=True,
            )
        else:
            with stdout_path.open("w") as out_handle:
                result = subprocess.run(
                    list(cmd),
                    cwd=str(cwd) if cwd else None,
                    stdout=out_handle,
                    stderr=handle,
                    text=True,
                )

    if result.returncode != 0:
        raise BenchmarkError(
            f"Command failed with exit code {result.returncode}: {' '.join(cmd)}\n"
            f"See log: {log_file}"
        )


def which_or_error(executable: str, label: str) -> str:
    path = shutil.which(executable)
    if path is None:
        raise BenchmarkError(f"Could not find {label} executable on PATH: {executable}")
    return path


def parse_chunk_index(path_or_name: str | Path) -> int:
    """Parse numeric chunk index from raw or aligned chunk FASTA names.

    Supported names:
      chunk_001.fasta
      chunk_001_align-pass.fasta
    """
    name = Path(path_or_name).name
    match = re.search(r"chunk_(\d+)(?:_align-pass)?\.fasta$", name)
    if not match:
        raise BenchmarkError(f"Could not parse chunk index from: {name}")
    return int(match.group(1))


def detect_available_cpus() -> int:
    """Return CPUs available to this process, respecting Linux CPU affinity."""
    try:
        affinity = os.sched_getaffinity(0)  # type: ignore[attr-defined]
        if affinity:
            return len(affinity)
    except Exception:
        pass
    return os.cpu_count() or 1


def resolve_jobs(value: int, label: str) -> int:
    if value < 0:
        raise BenchmarkError(f"--{label} must be >= 0; use 0 for automatic CPU-count mode.")
    if value == 0:
        return max(1, detect_available_cpus())
    return max(1, value)


# -----------------------------------------------------------------------------
# FASTA/group helpers
# -----------------------------------------------------------------------------


def extract_field(description: str, field: str) -> Optional[str]:
    prefix = f"{field}="
    for token in description.replace("|", " ").split():
        if token.startswith(prefix):
            return token.split("=", 1)[1]
    return None


def group_records_by_field(
    records: Iterable[SeqRecord],
    field: str,
) -> "OrderedDict[str, List[SeqRecord]]":
    groups: "OrderedDict[str, List[SeqRecord]]" = OrderedDict()
    for rec in records:
        group_id = extract_field(rec.description, field)
        if group_id is None:
            raise BenchmarkError(f"Record is missing required field {field}: {rec.description}")
        groups.setdefault(group_id, []).append(rec)
    return groups


def records_to_fasta(records: Sequence[SeqRecord]) -> str:
    handle = StringIO()
    SeqIO.write(list(records), handle, "fasta-2line")
    return handle.getvalue()


def write_fasta_2line(records: Sequence[SeqRecord], path: Path) -> None:
    with path.open("w") as handle:
        SeqIO.write(list(records), handle, "fasta-2line")


def count_fasta_records(path: Path) -> int:
    count = 0
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count


def restore_headers_and_uppercase(
    alignment: MultipleSeqAlignment,
    original_records: Sequence[SeqRecord],
) -> List[SeqRecord]:
    """Restore original LongAIRR/pRESTO headers and uppercase aligned bases.

    MAFFT may emit lowercase bases; pRESTO BuildConsensus expects uppercase
    nucleotide symbols/gaps. If an aligner rewrites record identifiers, we fall
    back to restoring headers by row order, which is valid for per-group calls
    when output has the same number of rows as input.
    """
    desc_by_id = {rec.id: rec.description for rec in original_records}
    rec_by_id = {rec.id: rec for rec in original_records}
    restored: List[SeqRecord] = []

    for idx, aln in enumerate(alignment):
        if aln.id in desc_by_id:
            original = rec_by_id[aln.id]
        elif len(alignment) == len(original_records):
            original = original_records[idx]
        else:
            original = None

        if original is not None:
            aln.id = original.id
            aln.name = original.name
            aln.description = original.description

        aln.seq = Seq(str(aln.seq).upper())
        aln.letter_annotations = {}
        restored.append(aln)

    return restored


# -----------------------------------------------------------------------------
# Aligner command handling
# -----------------------------------------------------------------------------


DEFAULT_ALIGNER_ARGS = {
    "muscle38": "-diags -maxiters 2",
    "muscle5": "-perm none -perturb 0",
    "mafft": "--auto --quiet",
    "abpoa": "-r 1",
}


DEFAULT_EXECUTABLES = {
    "muscle38": "muscle",
    "muscle5": "muscle",
    "mafft": "mafft",
    "abpoa": "abpoa",
}


def get_tool_version(aligner: str, executable: str) -> str:
    """Return a compact version string. Version failures are non-fatal."""
    commands = {
        "muscle38": [executable, "-version"],
        "muscle5": [executable, "-version"],
        "mafft": [executable, "--version"],
        "abpoa": [executable, "--version"],
    }
    try:
        result = subprocess.run(
            commands[aligner],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=10,
        )
    except Exception as exc:
        return f"unknown ({exc})"

    if result.returncode != 0:
        return "unknown"

    text = result.stdout.strip()
    if not text:
        return "unknown"

    if aligner in {"muscle38", "muscle5"}:
        words = text.split()
        if len(words) >= 2:
            version = words[1]
        else:
            version = words[0]
        version = re.sub(r"^v", "", version)
        version = re.sub(r"\.linux.*$", "", version)
        return version

    return text.splitlines()[0].strip()


def validate_muscle_version(aligner: str, version: str) -> None:
    if aligner not in {"muscle38", "muscle5"}:
        return
    match = re.search(r"(\d+)", version)
    if not match:
        log(f"Warning: could not parse MUSCLE version from '{version}'.")
        return
    major = int(match.group(1))
    if aligner == "muscle38" and major >= 5:
        raise BenchmarkError(
            f"--aligner muscle38 requested, but executable looks like MUSCLE {version}."
        )
    if aligner == "muscle5" and major < 5:
        raise BenchmarkError(
            f"--aligner muscle5 requested, but executable looks like MUSCLE {version}."
        )


def add_muscle5_defaults(args: List[str], threads: int) -> List[str]:
    out = list(args)
    if "-nt" not in out and "-amino" not in out:
        out = ["-nt", *out]
    if threads > 0 and "-threads" not in out:
        out = [*out, "-threads", str(threads)]
    return out


def build_aligner_command(
    aligner: str,
    executable: str,
    aligner_args: Sequence[str],
    n_records: int,
    muscle5_super5_threshold: int,
) -> Tuple[List[str], str]:
    """Build a stdin/stdout command and return (cmd, msa_mode)."""
    if aligner == "muscle38":
        return [executable, *aligner_args], "muscle38"

    if aligner == "mafft":
        return [executable, *aligner_args, "/dev/stdin"], "mafft"

    if aligner == "abpoa":
        return [executable, *aligner_args, "/dev/stdin"], "abpoa"

    if aligner == "muscle5":
        mode = "super5" if n_records >= muscle5_super5_threshold else "align"
        return [executable, f"-{mode}", "/dev/stdin", "-output", "/dev/stdout", *aligner_args], mode

    raise BenchmarkError(f"Unsupported aligner: {aligner}")


def run_group_alignment(
    records: List[SeqRecord],
    group_id: str,
    aligner: str,
    executable: str,
    aligner_args: Sequence[str],
    muscle5_super5_threshold: int,
    group_timeout_sec: int,
) -> Tuple[List[SeqRecord], str]:
    """Run one per-group MSA and return aligned records plus MSA mode."""
    if len(records) < 2:
        return records, "none"

    cmd, mode = build_aligner_command(
        aligner=aligner,
        executable=executable,
        aligner_args=aligner_args,
        n_records=len(records),
        muscle5_super5_threshold=muscle5_super5_threshold,
    )
    timeout = group_timeout_sec if group_timeout_sec and group_timeout_sec > 0 else None

    try:
        child = subprocess.run(
            cmd,
            input=records_to_fasta(records),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired as exc:
        raise BenchmarkError(
            f"{aligner} timed out after {group_timeout_sec}s for group {group_id} "
            f"with {len(records)} records. Command: {' '.join(cmd)}"
        ) from exc

    if child.returncode != 0:
        raise BenchmarkError(
            f"{aligner} failed for group {group_id} with {len(records)} records.\n"
            f"Command: {' '.join(cmd)}\n"
            f"stderr/stdout:\n{child.stderr}\n{child.stdout[:1000]}"
        )

    if not child.stdout.strip():
        raise BenchmarkError(
            f"{aligner} produced empty stdout for group {group_id}.\n"
            f"Command: {' '.join(cmd)}\n"
            f"stderr:\n{child.stderr}"
        )

    try:
        alignment = AlignIO.read(StringIO(child.stdout), "fasta")
    except Exception as exc:
        raise BenchmarkError(
            f"{aligner} produced unreadable aligned FASTA for group {group_id}.\n"
            f"Command: {' '.join(cmd)}\n"
            f"stderr:\n{child.stderr}\n"
            f"stdout head:\n{child.stdout[:1000]}"
        ) from exc

    return restore_headers_and_uppercase(alignment, records), mode


# -----------------------------------------------------------------------------
# Chunk alignment
# -----------------------------------------------------------------------------


def progress_line(chunk_name: str, done: int, total: int, start: float, passed: int, failed: int) -> str:
    pct = int(round(100 * done / total)) if total else 100
    filled = int(round(pct / 5))
    bar = "#" * filled + " " * (20 - filled)
    elapsed = (time.time() - start) / 60
    return (
        f"PROGRESS> {chunk_name} {time.strftime('%H:%M:%S')} |{bar}| "
        f"{pct:3d}% ({done}/{total}) {elapsed:.1f} min PASS={passed} FAIL={failed}"
    )


def align_one_group(
    idx: int,
    chunk_name: str,
    group_id: str,
    group_records: List[SeqRecord],
    aligner: str,
    executable: str,
    aligner_args: Sequence[str],
    muscle5_super5_threshold: int,
    group_timeout_sec: int,
) -> Tuple[int, List[SeqRecord], Dict[str, object]]:
    start = time.time()
    aligned, mode = run_group_alignment(
        records=group_records,
        group_id=group_id,
        aligner=aligner,
        executable=executable,
        aligner_args=aligner_args,
        muscle5_super5_threshold=muscle5_super5_threshold,
        group_timeout_sec=group_timeout_sec,
    )

    lengths = [len(rec.seq) for rec in group_records]
    row = {
        "chunk": chunk_name,
        "group_index": idx,
        "group_id": group_id,
        "aligner": aligner,
        "msa_mode": mode,
        "n_reads": len(group_records),
        "input_min_len": min(lengths) if lengths else None,
        "input_max_len": max(lengths) if lengths else None,
        "aligned_len": len(aligned[0].seq) if aligned else None,
        "runtime_sec": round(time.time() - start, 6),
        "status": "pass",
        "error": "",
    }
    return idx, aligned, row


def align_chunk(
    chunk_path: str,
    aligned_dir: str,
    group_field: str,
    aligner: str,
    executable: str,
    aligner_args: Sequence[str],
    muscle5_super5_threshold: int,
    group_timeout_sec: int,
    group_jobs: int,
    sort_groups_by_size: bool,
    summary_min_reads: int,
    progress_step_pct: int,
    fast_singleton_chunk: bool,
) -> Tuple[str, int, int, List[Dict[str, object]]]:
    """Align all groups in one chunk and write *_align-pass.fasta."""
    chunk = Path(chunk_path)
    aligned_out = Path(aligned_dir) / f"{chunk.stem}_align-pass.fasta"
    start = time.time()

    if fast_singleton_chunk and chunk.name == "chunk_001.fasta":
        shutil.copyfile(chunk, aligned_out)
        n_records = count_fasta_records(chunk)
        print(
            f"CHUNK_END> {chunk.name} OUTPUT={aligned_out.name} "
            f"SEQUENCES={n_records} SETS={n_records} PASS={n_records} FAIL=0 "
            f"ELAPSED_MIN={(time.time() - start) / 60:.1f}",
            flush=True,
        )
        return str(aligned_out), n_records, n_records, []

    records = list(SeqIO.parse(str(chunk), "fasta"))
    groups = list(group_records_by_field(records, group_field).items())
    if sort_groups_by_size:
        groups.sort(key=lambda item: len(item[1]), reverse=True)

    total = len(groups)
    jobs = min(max(1, group_jobs), max(1, total))
    print(
        f"\nCHUNK_START> {chunk.name} SEQUENCES={len(records)} SETS={total} "
        f"GROUP_JOBS={jobs}",
        flush=True,
    )
    if total:
        print(progress_line(chunk.name, 0, total, start, 0, 0), flush=True)

    aligned_by_idx: Dict[int, List[SeqRecord]] = {}
    rows_by_idx: Dict[int, Dict[str, object]] = {}
    passed = 0
    failed = 0
    next_pct = max(1, progress_step_pct)

    def maybe_progress() -> None:
        nonlocal next_pct
        done = passed + failed
        while total and done >= total * next_pct / 100:
            print(progress_line(chunk.name, done, total, start, passed, failed), flush=True)
            next_pct += max(1, progress_step_pct)

    if jobs == 1:
        for idx, (group_id, group_records) in enumerate(groups):
            try:
                out_idx, aligned, row = align_one_group(
                    idx,
                    chunk.name,
                    group_id,
                    group_records,
                    aligner,
                    executable,
                    aligner_args,
                    muscle5_super5_threshold,
                    group_timeout_sec,
                )
                aligned_by_idx[out_idx] = aligned
                rows_by_idx[out_idx] = row
                passed += 1
            except Exception:
                failed += 1
                raise
            finally:
                maybe_progress()
    else:
        # Thread workers are sufficient here because each task launches an
        # external aligner subprocess. This avoids expensive pickling of
        # SeqRecord groups inside every chunk.
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            futures = {
                pool.submit(
                    align_one_group,
                    idx,
                    chunk.name,
                    group_id,
                    group_records,
                    aligner,
                    executable,
                    aligner_args,
                    muscle5_super5_threshold,
                    group_timeout_sec,
                ): (idx, group_id)
                for idx, (group_id, group_records) in enumerate(groups)
            }
            for future in as_completed(futures):
                idx, group_id = futures[future]
                try:
                    out_idx, aligned, row = future.result()
                    aligned_by_idx[out_idx] = aligned
                    rows_by_idx[out_idx] = row
                    passed += 1
                except Exception as exc:
                    failed += 1
                    raise BenchmarkError(f"Failed processing group {group_id} in {chunk.name}: {exc}") from exc
                finally:
                    maybe_progress()

    aligned_records: List[SeqRecord] = []
    summary_rows: List[Dict[str, object]] = []
    for idx, (_group_id, original_group_records) in enumerate(groups):
        aligned_records.extend(aligned_by_idx[idx])
        row = rows_by_idx[idx]
        if len(original_group_records) >= summary_min_reads or row.get("status") == "fail":
            summary_rows.append(row)

    write_fasta_2line(aligned_records, aligned_out)
    print(
        f"CHUNK_END> {chunk.name} OUTPUT={aligned_out.name} SEQUENCES={len(records)} "
        f"SETS={total} PASS={passed} FAIL={failed} ELAPSED_MIN={(time.time() - start) / 60:.1f}\n",
        flush=True,
    )
    return str(aligned_out), len(records), total, summary_rows


# -----------------------------------------------------------------------------
# LongAIRR/pRESTO commands
# -----------------------------------------------------------------------------


def build_consensus_command(
    library: str,
    aligned_fasta: Path,
    group_field: str,
    n_cons: int,
    cons_gap: float,
    cons_error: float,
    outdir: Path,
) -> List[str]:
    cmd = [
        "BuildConsensus.py",
        "-n",
        str(n_cons),
        "-s",
        str(aligned_fasta),
        "--bf",
        group_field,
        "--failed",
    ]

    if library == "bulk":
        cmd += ["--cf", "N_ORIG", "N_KEEP", "UMI", "--act", "majority", "majority", "majority"]
    elif library == "visium":
        cmd += [
            "--cf",
            "N_ORIG",
            "N_KEEP",
            "UMI",
            "X",
            "Y",
            "--act",
            "majority",
            "majority",
            "majority",
            "majority",
            "majority",
            "--pf",
            "SPBC",
        ]
    elif library == "visiumhd":
        cmd += [
            "--cf",
            "N_ORIG",
            "N_KEEP",
            "UMI",
            "BC1",
            "BC2",
            "X",
            "Y",
            "SPBC10X",
            "--act",
            "majority",
            "majority",
            "majority",
            "majority",
            "majority",
            "majority",
            "majority",
            "majority",
            "--pf",
            "SPBC",
        ]
    else:
        raise BenchmarkError(f"Unsupported library: {library}")

    cmd += [
        "--maxgap",
        str(cons_gap),
        "--maxerror",
        str(cons_error),
        "--outname",
        "groups_msa",
        "--outdir",
        str(outdir),
    ]
    return cmd


def write_group_summary(path: Path, rows: Sequence[Dict[str, object]]) -> None:
    fields = [
        "chunk",
        "group_index",
        "group_id",
        "aligner",
        "msa_mode",
        "n_reads",
        "input_min_len",
        "input_max_len",
        "aligned_len",
        "runtime_sec",
        "status",
        "error",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_run_summary(path: Path, rows: Dict[str, object]) -> None:
    with path.open("w") as handle:
        handle.write("key\tvalue\n")
        for key, value in rows.items():
            handle.write(f"{key}\t{value}\n")


# -----------------------------------------------------------------------------
# Main workflow
# -----------------------------------------------------------------------------


def prepare_outdir(outdir: Path, force: bool) -> None:
    if outdir.exists() and any(outdir.iterdir()):
        if not force:
            raise BenchmarkError(f"Output directory exists and is not empty: {outdir}. Use --force to overwrite.")
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)


def run_benchmark(args: argparse.Namespace) -> None:
    input_fasta = Path(args.input).resolve()
    outdir = Path(args.outdir).resolve()
    prepare_outdir(outdir, args.force)

    if not input_fasta.exists():
        raise BenchmarkError(f"Input FASTA does not exist: {input_fasta}")

    chunks_dir = outdir / "1_chunks"
    aligned_dir = outdir / "2_aligned_chunks"
    chunks_dir.mkdir(parents=True, exist_ok=True)
    aligned_dir.mkdir(parents=True, exist_ok=True)
    log_file = outdir / "benchmark.log"
    log_file.write_text("# LongAIRR MSA backend benchmark\n")

    aligner_exec = which_or_error(args.aligner_exec or DEFAULT_EXECUTABLES[args.aligner], args.aligner)
    aligner_args = shlex.split(args.aligner_args if args.aligner_args is not None else DEFAULT_ALIGNER_ARGS[args.aligner])
    if args.aligner == "muscle5":
        aligner_args = add_muscle5_defaults(aligner_args, args.muscle5_threads)

    aligner_version = get_tool_version(args.aligner, aligner_exec)
    validate_muscle_version(args.aligner, aligner_version)

    chunk_jobs = resolve_jobs(args.chunk_jobs, "chunk-jobs")
    group_jobs = resolve_jobs(args.group_jobs, "group-jobs")

    log(f"Input FASTA: {input_fasta}")
    log(f"Output dir:  {outdir}")
    log(f"Library:     {args.library}")
    log(f"Group field: {args.group_field}")
    log(f"Aligner:     {args.aligner} ({aligner_exec}, version {aligner_version})")
    log(f"Aligner args: {' '.join(aligner_args)}")
    log(f"Chunk jobs:  {chunk_jobs}; group jobs per chunk: {group_jobs}")

    split_script = which_or_error(args.split_script, "longairr_split_align.py")
    run_command(
        [
            split_script,
            str(input_fasta),
            "--group-field",
            args.group_field,
            "--max_umi_per_chunk",
            str(args.max_umi_per_chunk),
            "--outdir",
            str(chunks_dir),
        ],
        log_file=log_file,
    )

    chunk_paths = sorted(chunks_dir.glob("chunk_*.fasta"), key=parse_chunk_index)
    if not chunk_paths:
        raise BenchmarkError(f"No chunk_*.fasta files were created in {chunks_dir}")
    log(f"Created/found {len(chunk_paths)} chunks")

    all_rows: List[Dict[str, object]] = []
    aligned_paths: List[Path] = []
    n_records_total = 0
    n_groups_total = 0

    log("Running group-wise MSA")
    with ProcessPoolExecutor(max_workers=chunk_jobs) as pool:
        futures = {
            pool.submit(
                align_chunk,
                str(chunk),
                str(aligned_dir),
                args.group_field,
                args.aligner,
                aligner_exec,
                aligner_args,
                args.muscle5_super5_threshold,
                args.group_timeout_sec,
                group_jobs,
                args.sort_groups_by_size,
                args.summary_min_reads,
                args.progress_step_pct,
                args.fast_singleton_chunk,
            ): chunk
            for chunk in chunk_paths
        }
        for future in as_completed(futures):
            chunk = futures[future]
            try:
                aligned_out, n_records, n_groups, rows = future.result()
            except Exception as exc:
                raise BenchmarkError(f"Failed processing chunk {chunk.name}: {exc}") from exc

            aligned_paths.append(Path(aligned_out))
            all_rows.extend(rows)
            n_records_total += n_records
            n_groups_total += n_groups
            log(f"Aligned {chunk.name}: {n_records} reads, {n_groups} groups")

            if not args.keep_chunks:
                try:
                    chunk.unlink()
                except OSError:
                    pass

    write_group_summary(outdir / "group_msa_summary.tsv", all_rows)

    log("Concatenating aligned chunks")
    combined = outdir / "groups_aligned_combined.fasta"
    with combined.open("w") as out_handle:
        for path in sorted(aligned_paths, key=parse_chunk_index):
            with path.open() as in_handle:
                shutil.copyfileobj(in_handle, out_handle)

    if not combined.exists() or combined.stat().st_size == 0:
        raise BenchmarkError(f"Combined aligned FASTA is empty: {combined}")

    log("Removing duplicate sequence names with seqkit rmdup -n")
    dedup = outdir / "3_groups_aligned_combined_dedup.fasta"
    run_command(["seqkit", "rmdup", "-n", str(combined)], log_file=log_file, stdout_path=dedup)

    if not args.keep_combined:
        combined.unlink(missing_ok=True)

    log("Building consensus with BuildConsensus.py")
    run_command(
        build_consensus_command(
            library=args.library,
            aligned_fasta=dedup,
            group_field=args.group_field,
            n_cons=args.n_cons,
            cons_gap=args.cons_gap,
            cons_error=args.cons_error,
            outdir=outdir,
        ),
        log_file=log_file,
    )

    consensus_pass = outdir / "groups_msa_consensus-pass.fasta"
    if not consensus_pass.exists() or consensus_pass.stat().st_size == 0:
        raise BenchmarkError(f"BuildConsensus output missing or empty: {consensus_pass}")

    log("Rewriting consensus headers with longairr_header_metadata.py")
    header_script = which_or_error(args.header_script, "longairr_header_metadata.py")
    final_out = outdir / args.final_name
    run_command(
        [
            header_script,
            "--input",
            str(consensus_pass),
            "-o",
            str(final_out),
            "--library",
            args.library,
        ],
        log_file=log_file,
    )

    if not final_out.exists() or final_out.stat().st_size == 0:
        raise BenchmarkError(f"Final consensus output missing or empty: {final_out}")

    if not args.keep_aligned_chunks:
        log("Removing per-chunk aligned FASTA files")
        for path in aligned_paths:
            path.unlink(missing_ok=True)

    write_run_summary(
        outdir / "benchmark_summary.tsv",
        {
            "input": input_fasta,
            "outdir": outdir,
            "library": args.library,
            "group_field": args.group_field,
            "aligner": args.aligner,
            "aligner_exec": aligner_exec,
            "aligner_version": aligner_version,
            "aligner_args": " ".join(aligner_args),
            "max_umi_per_chunk": args.max_umi_per_chunk,
            "chunk_jobs_requested": args.chunk_jobs,
            "chunk_jobs_effective": chunk_jobs,
            "group_jobs_requested": args.group_jobs,
            "group_jobs_effective": group_jobs,
            "muscle5_super5_threshold": args.muscle5_super5_threshold,
            "muscle5_threads": args.muscle5_threads,
            "n_chunks": len(chunk_paths),
            "n_records_total": n_records_total,
            "n_groups_total": n_groups_total,
            "n_group_summary_rows": len(all_rows),
            "n_cons": args.n_cons,
            "cons_gap": args.cons_gap,
            "cons_error": args.cons_error,
            "final_output": final_out,
        },
    )

    log(f"Done. Final output: {final_out}")
    log(f"Run summary: {outdir / 'benchmark_summary.tsv'}")


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run a LongAIRR group-wise MSA backend benchmark."
    )
    parser.add_argument("--input", required=True, help="Input FASTA, usually 2_subsampled_groups.fasta")
    parser.add_argument("--outdir", required=True, help="Output directory for this benchmark run")
    parser.add_argument("--library", required=True, choices=["bulk", "visium", "visiumhd"])
    parser.add_argument("--group-field", required=True, choices=["UMI", "SPBCUMI", "UMISPBC"])

    parser.add_argument("--aligner", required=True, choices=["muscle38", "muscle5", "mafft", "abpoa"])
    parser.add_argument("--aligner-exec", help="Executable path/name. Defaults to muscle/mafft/abpoa depending on --aligner")
    parser.add_argument(
        "--aligner-args",
        help=(
            "Extra aligner arguments. If omitted, defaults are used: "
            "muscle38='-diags -maxiters 2', muscle5='-perm none -perturb 0', "
            "mafft='--auto --quiet', abpoa='-r 1'."
        ),
    )
    parser.add_argument(
        "--muscle5-super5-threshold",
        type=int,
        default=100,
        help="Use MUSCLE5 -super5 for groups with at least this many reads. Default: 100.",
    )
    parser.add_argument(
        "--muscle5-threads",
        "--muscle-threads",
        dest="muscle5_threads",
        type=int,
        default=1,
        help="Threads passed to MUSCLE5 via -threads. Use 0 to let MUSCLE5 choose. Default: 1.",
    )

    parser.add_argument("--max-umi-per-chunk", type=int, default=50, help="Passed to longairr_split_align.py")
    parser.add_argument(
        "--chunk-jobs",
        "--parallel-cjobs",
        dest="chunk_jobs",
        type=int,
        default=1,
        help="Chunks processed in parallel. Use 0 for all available CPUs. Default: 1.",
    )
    parser.add_argument(
        "--group-jobs",
        "--group-nproc",
        dest="group_jobs",
        type=int,
        default=1,
        help="Groups processed in parallel within each chunk. Use 0 for all available CPUs. Default: 1.",
    )
    parser.add_argument("--group-timeout-sec", type=int, default=3600, help="Timeout per group MSA; 0 disables")
    parser.add_argument(
        "--no-sort-groups-by-size",
        dest="sort_groups_by_size",
        action="store_false",
        help="Disable largest-first group scheduling inside each chunk.",
    )
    parser.set_defaults(sort_groups_by_size=True)
    parser.add_argument("--progress-step-pct", type=int, default=10, help="Per-chunk progress interval in percent")

    parser.add_argument("--n-cons", type=int, default=1, help="BuildConsensus.py -n value")
    parser.add_argument("--cons-gap", type=float, default=0.5, help="BuildConsensus.py --maxgap")
    parser.add_argument("--cons-error", type=float, default=0.2, help="BuildConsensus.py --maxerror")

    parser.add_argument("--split-script", default="longairr_split_align.py")
    parser.add_argument("--header-script", default="longairr_header_metadata.py")
    parser.add_argument("--final-name", default="4_longairr_consensus.fasta")

    parser.add_argument("--summary-min-reads", type=int, default=2)
    parser.add_argument("--no-fast-singleton-chunk", dest="fast_singleton_chunk", action="store_false")
    parser.set_defaults(fast_singleton_chunk=True)
    parser.add_argument("--keep-chunks", action="store_true", help="Keep chunk_*.fasta files after alignment")
    parser.add_argument("--keep-aligned-chunks", action="store_true", help="Keep per-chunk aligned FASTA files after finalization")
    parser.add_argument("--keep-combined", action="store_true", help="Keep non-deduplicated combined aligned FASTA")
    parser.add_argument("--force", action="store_true", help="Overwrite an existing non-empty output directory")

    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        run_benchmark(args)
    except BenchmarkError as exc:
        print(f"\n[ERROR] {exc}\n", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("\n[ERROR] Interrupted by user.\n", file=sys.stderr)
        return 130
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
