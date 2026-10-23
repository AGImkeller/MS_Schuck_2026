set -euo pipefail

#######
#
# Script subsamples input fastqs to specific %-subsamples listed below
#
# Can be used within longairr environment
#
#######



#FASTQ="/home/jonas/mnt_AG_Imkeller/group_members/schuck/nanoairr/nanoairr_data/nanoairr_sp12/basecalling/simplex_duplex.fastq"     
#OUTDIR="/home/jonas/Desktop/Airr_workflow/repo/8_rarefaction/subsamples_tail/"          

#FASTQ="/home/jonas/mnt_AG_Imkeller/group_members/schuck/nanoairr/nanoairr_data/pacbio/pacbio_s3_e52.fastq"    
#OUTDIR="/home/jonas/mnt_AG_Imkeller/group_members/schuck/nanoairr/manuscript/data/rarefaction/sp12_raw_subsamples/subsamples_pb/"          

FASTQ="/home/jonas/mnt_AG_Imkeller/group_members/schuck/nanoairr/nanoairr_data/longairr_sp28/basecall/simplex_duplex.fastq"     
OUTDIR="/home/jonas/mnt_AG_Imkeller/group_members/schuck/manuscripts/longairr/1_samples_longairr/1_longairr_revision/manuscript_analysis/rarefaction/subsamples_raw/subsamples_sp28/"          

SEED=42                                   

# fractions to regenerate (exactly those that broke)
#TARGETS=("0.215443" "0.464159")
#TARGETS=("0.600000" "0.750000" "0.900000")
#TARGETS=("0.001" "0.002154" "0.004642" "0.01" "0.021544" "0.046416" "0.1" "0.215443" "0.464159" "0.6" "0.75")

TARGETS=("0.001000" "0.002154" "0.004642" "0.010000" "0.021544" "0.046416" "0.100000" "0.215443" "0.464159" "0.600000" "0.750000" "0.900000" "1.000000")


# ---- helper: strict fastq validator ----
validate_fastq () {
  python3 - "$1" <<'PY'
import sys
fn=sys.argv[1]
i=0
def has_ws(s): return any(c.isspace() for c in s)
with open(fn,"rb") as fh:
    while True:
        h=fh.readline()
        if not h: break
        s=fh.readline(); p=fh.readline(); q=fh.readline()
        if not (s and p and q):
            print(f"TRUNCATED at record {i+1}", file=sys.stderr)
            sys.exit(2)
        i+=1
        h2=h.rstrip(b"\n\r").decode("utf-8","replace")
        s2=s.rstrip(b"\n\r").decode("utf-8","replace")
        p2=p.rstrip(b"\n\r").decode("utf-8","replace")
        q2=q.rstrip(b"\n\r").decode("utf-8","replace")
        if (not h2.startswith("@")) or (not p2.startswith("+")) or has_ws(s2) or has_ws(q2) or (len(s2)!=len(q2)):
            print(f"BAD record {i}: seq_len={len(s2)} qual_len={len(q2)} header={h2}", file=sys.stderr)
            sys.exit(1)
print(f"OK: validated {i} records", file=sys.stderr)
sys.exit(0)
PY
}

# ---- helper: fingerprint based on header stream (order-sensitive) ----
fingerprint_headers () {
  python3 - "$1" <<'PY'
import sys, hashlib
fn=sys.argv[1]
h=hashlib.sha256()
n=0
with open(fn,"rb") as f:
    while True:
        hdr=f.readline()
        if not hdr: break
        seq=f.readline(); plus=f.readline(); qual=f.readline()
        if not (seq and plus and qual):
            break
        hdr=hdr.rstrip(b"\r\n")
        h.update(hdr); h.update(b"\n")
        n+=1
print(f"{fn}\trecords={n}\tsha256(headers)={h.hexdigest()}")
PY
}

# basename handling (must match your script)
base="$(basename "${FASTQ}")"
base="${base%.gz}"
base="${base%.fastq}"
base="${base%.fq}"

for frac in "${TARGETS[@]}"; do
  frac_tag="${frac/./_}"
  out="${OUTDIR}/${base}_f${frac_tag}.fastq"
  log="${OUTDIR}/${base}_f${frac_tag}.log"

  echo "=== Regenerating $out (seed=$SEED, frac=$frac) ==="

  # show old fingerprint (even if file is broken, this often still works)
  if [ -f "$out" ]; then
    echo "[OLD]" | tee -a "$log"
    fingerprint_headers "$out" | tee -a "$log" || true
  fi

  tmp="${TMPDIR:-/tmp}/$(basename "$out").tmp.$$"

  # regenerate to tmp
  seqkit sample --proportion "$frac" -s "$SEED" -o "$tmp" "$FASTQ" >"$log" 2>&1

  # strict validate tmp
  validate_fastq "$tmp" >>"$log" 2>&1

  # atomic replace
  mv -f "$tmp" "$out"

  echo "[NEW]" | tee -a "$log"
  fingerprint_headers "$out" | tee -a "$log"
  echo "OK -> $out"
done
