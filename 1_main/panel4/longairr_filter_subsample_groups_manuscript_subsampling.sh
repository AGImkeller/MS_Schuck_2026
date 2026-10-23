#!/bin/bash


# Tested on LongAIRR version 7.1 | 1.1.0
# 03.02.2026 | 31.08.2026

#---------------------
# log message function

#######################################
# Log messages with timestamps
# 
# Formats and writes log messages to both console and a log file.
#
# Arguments:
#   $1 - Log level (INFO, WARNING, ERROR)
#   $2 - Nanoairr step/module name
#   $3 - Log message
#   $4 - Path to log file
#
# Output:
#   Logs messages in the format: [TIMESTAMP] [LEVEL] [STEP] MESSAGE
#######################################
log_message() {
  local level="$1"      # Log level: INFO, WARNING, ERROR
  local step="$2"
  local message="$3"    # Log message
  local log_file="$4"   # Path to log file

  # Get current timestamp
  local timestamp
  timestamp=$(date +"%Y-%m-%d %H:%M:%S")

  # Format the log message with timestamp and level
  local formatted_message="\n[${timestamp}] [${level}] [${step}] ${message}"

  # Append the message to log file
  if [[ -n "${log_file}" ]]; then
    echo -e "${formatted_message}" >> "${log_file}"
  fi
  echo -e "${formatted_message}"
}


#-------input parameters--------------------------------------------------------

#sizes=(5000 2500 2000 1500 1000 750 500 250 100 50 25 10 5 3 2 1) # List of subsample sizes
#sizes=(10000) # List of subsample sizes
#sizes=(1000 750 500 250 100 50 10 5 4 3 2 1) # List of subsample sizes
#sizes=(100 50 10 5 4 3 2 1) # List of subsample sizes


#s1 ont and pacbio
#sizes=(10000 100 50 10 5 4 3 2 1) # List of subsample sizes
#seeds=(1 2 3 4 5) # List of seeds

#s3_large
#sizes=(10000) # List of subsample sizes
sizes=(5000 2500 2000 1500 1000 750 500 250 100 50 10 5 4 3 2 1) # List of subsample sizes
#sizes=(1000 750 500 250 100 50 10 5 4 3 2 1) # List of subsample sizes
seeds=(1 2 3 4 5) # List of seeds
#seeds=(1) # List of seeds


### revision

#s1_pacbio
#input_fasta="/home/jonas/mnt_AG_Imkeller/group_members/schuck/manuscripts/longairr/1_samples_longairr/1_longairr_revision/manuscript_analysis/large_groups_subsampling/sp12_pb/1_large_groups/group_CGAGGAGCTTCCATATCTAAAGTCAGTA.fasta"
#base_output_dir="/home/jonas/mnt_AG_Imkeller/group_members/schuck/manuscripts/longairr/1_samples_longairr/1_longairr_revision/manuscript_analysis/large_groups_subsampling/sp12_pb/2_results/CGAGGAGCTTCCATATCTAAAGTCAGTA/"

#s1_ont
#input_fasta="/home/jonas/mnt_AG_Imkeller/group_members/schuck/manuscripts/longairr/1_samples_longairr/1_longairr_revision/manuscript_analysis/large_groups_subsampling/sp12_ont/1_large_groups/group_CGAGGAGCTTCCATATCTAAAGTCAGTA.fasta"
#base_output_dir="/home/jonas/mnt_AG_Imkeller/group_members/schuck/manuscripts/longairr/1_samples_longairr/1_longairr_revision/manuscript_analysis/large_groups_subsampling/sp12_ont/2_results/CGAGGAGCTTCCATATCTAAAGTCAGTA/"


#s3_large
input_fasta="/home/jonas/mnt_AG_Imkeller/group_members/schuck/manuscripts/longairr/1_samples_longairr/1_longairr_revision/manuscript_analysis/large_groups_subsampling/sp8/1_large_groups/group_ACGTAGATTGCTGATGTACCCCAGGAAG.fasta"
base_output_dir="/home/jonas/mnt_AG_Imkeller/group_members/schuck/manuscripts/longairr/1_samples_longairr/1_longairr_revision/manuscript_analysis/large_groups_subsampling/sp8/2_results/ACGTAGATTGCTGATGTACCCCAGGAAG_2/"

##s2 g1
#input_fasta="/home/jonas/mnt_AG_Imkeller/group_members/schuck/manuscripts/longairr/1_samples_longairr/1_longairr_revision/manuscript_analysis/large_groups_subsampling/sp7/1_large_groups/group1_CCGATCTCAACCTTATGAGTTTGAAACG.fasta"
#base_output_dir="/home/jonas/mnt_AG_Imkeller/group_members/schuck/manuscripts/longairr/1_samples_longairr/1_longairr_revision/manuscript_analysis/large_groups_subsampling/sp7/2_results/group1_CCGATCTCAACCTTATGAGTTTGAAACG/"

##s2 g2
#input_fasta="/home/jonas/mnt_AG_Imkeller/group_members/schuck/manuscripts/longairr/1_samples_longairr/1_longairr_revision/manuscript_analysis/large_groups_subsampling/sp7/1_large_groups/group2_CAGACCTGTAAGTGTTACGCTACCTCAC.fasta"
#base_output_dir="/home/jonas/mnt_AG_Imkeller/group_members/schuck/manuscripts/longairr/1_samples_longairr/1_longairr_revision/manuscript_analysis/large_groups_subsampling/sp7/2_results/group2_CAGACCTGTAAGTGTTACGCTACCTCAC/"

##s2 g3
#input_fasta="/home/jonas/mnt_AG_Imkeller/group_members/schuck/manuscripts/longairr/1_samples_longairr/1_longairr_revision/manuscript_analysis/large_groups_subsampling/sp8/1_large_groups/group3_AAGTGAGTCGGGTTTATTACTGGGCTGG.fasta"
#base_output_dir="/home/jonas/mnt_AG_Imkeller/group_members/schuck/manuscripts/longairr/1_samples_longairr/1_longairr_revision/manuscript_analysis/large_groups_subsampling/sp8/2_results/group3_AAGTGAGTCGGGTTTATTACTGGGCTGG/"

##s2 g4
#input_fasta="/home/jonas/mnt_AG_Imkeller/group_members/schuck/manuscripts/longairr/1_samples_longairr/1_longairr_revision/manuscript_analysis/large_groups_subsampling/sp8/1_large_groups/group4_GGGACCCGTATATCTTTTAATAGGAGTT.fasta"
#base_output_dir="/home/jonas/mnt_AG_Imkeller/group_members/schuck/manuscripts/longairr/1_samples_longairr/1_longairr_revision/manuscript_analysis/large_groups_subsampling/sp8/2_results/group4_GGGACCCGTATATCTTTTAATAGGAGTT/"



#-------

# Example input/output paths (adjust these)
# Sequence input
#input_fasta="/home/jonas/mnt_AG_Imkeller/group_members/schuck/nanoairr/manuscript/data/subsampling/new/sp8/1_large_groups/group_CCTTCAGTTAAAGTGATGTCGTGATTTG.fasta"
##last-one##input_fasta="/home/jonas/mnt_AG_Imkeller/group_members/schuck/nanoairr/manuscript/data/subsampling/new/sp8/1_large_groups/group_ACGTAGATTGCTGATGTACCCCAGGAAG.fasta"
#input_fasta="/home/jonas/mnt_AG_Imkeller/group_members/schuck/nanoairr/manuscript/data/subsampling/new/pacbio/1_large_groups/group_ATGGTCGCGTGGTTTCTATAAGGCATGA.fasta"
#input_fasta="/home/jonas/mnt_AG_Imkeller/group_members/schuck/nanoairr/manuscript/data/subsampling/new/sp12/1_large_groups/group_ATGGTCGCGTGGTTTCTATAAGGCATGA.fasta"
#input_fasta="/home/jonas/mnt_AG_Imkeller/group_members/schuck/nanoairr/manuscript/data/subsampling/new/pacbio/1_large_groups/group_CGCAGATCTTCACCCGGCGTTATTATAG.fasta"
#input_fasta="/home/jonas/mnt_AG_Imkeller/group_members/schuck/nanoairr/manuscript/data/subsampling/new/pacbio/1_large_groups/group_ATGGTCGCGTGGTTTCGTGTGGTCGGAT.fasta"
#input_fasta="/home/jonas/mnt_AG_Imkeller/group_members/schuck/nanoairr/manuscript/data/subsampling/new/sp12/1_large_groups/group_ATGGTCGCGTGGTTTCGTGTGGTCGGAT.fasta"
#input_fasta="/home/jonas/mnt_AG_Imkeller/group_members/schuck/nanoairr/manuscript/data/subsampling/new/sp12/1_large_groups/group_CGCAGATCTTCACCCGGCGTTATTATAG.fasta"
#input_fasta="/home/jonas/mnt_AG_Imkeller/group_members/schuck/nanoairr/manuscript/data/subsampling/new/sp12/1_large_groups/group_CGAGGAGCTTCCATATCTAAAGTCAGTA.fasta"
#input_fasta="/home/jonas/mnt_AG_Imkeller/group_members/schuck/nanoairr/manuscript/data/subsampling/new/pacbio/1_large_groups/group_CGAGGAGCTTCCATATCTAAAGTCAGTA.fasta"
#input_fasta="/home/jonas/mnt_AG_Imkeller/group_members/schuck/nanoairr/manuscript/data/subsampling/new/pacbio/1_large_groups/group_GGTCTTGAGCGCTCTTAGAACGGGCATT.fasta"

#g1-g4 subsample
#input_fasta="/home/jonas/mnt_AG_Imkeller/group_members/schuck/nanoairr/manuscript/data/subsampling/new/sp7/1_large_groups/g1_group_CCGATCTCAACCTTATGAGTTTGAAACG.fasta"
#input_fasta="/home/jonas/mnt_AG_Imkeller/group_members/schuck/nanoairr/manuscript/data/subsampling/new/sp7/1_large_groups/g2_group_CAGACCTGTAAGTGTTACGCTACCTCAC.fasta"
#input_fasta="/home/jonas/mnt_AG_Imkeller/group_members/schuck/nanoairr/manuscript/data/subsampling/new/sp8/1_large_groups/g3_group_AAGTGAGTCGGGTTTATTACTGGGCTGG.fasta"
#input_fasta="/home/jonas/mnt_AG_Imkeller/group_members/schuck/nanoairr/manuscript/data/subsampling/new/sp8/1_large_groups/g4_group_GGGACCCGTATATCTTTTAATAGGAGTT.fasta"



# Output
##last-one###base_output_dir="/home/jonas/mnt_AG_Imkeller/group_members/schuck/nanoairr/manuscript/data/subsampling/new/sp8/2_longairr_out/igl/"
#base_output_dir="/home/jonas/mnt_AG_Imkeller/group_members/schuck/nanoairr/manuscript/data/subsampling/new/pacbio/2_longairr_out/group_ATGGTCGCGTGGTTTCGTGTGGTCGGAT_005/"
#base_output_dir="/home/jonas/mnt_AG_Imkeller/group_members/schuck/nanoairr/manuscript/data/subsampling/new/sp12/2_longairr_out/CGCAGATCTTCACCCGGCGTTATTATAG_005/"
#base_output_dir="/home/jonas/mnt_AG_Imkeller/group_members/schuck/nanoairr/manuscript/data/subsampling/new/pacbio/2_longairr_out/GGTCTTGAGCGCTCTTAGAACGGGCATT_005/"



#g1-g4
#base_output_dir="/home/jonas/mnt_AG_Imkeller/group_members/schuck/nanoairr/manuscript/data/subsampling/new/sp7/2_longairr_out/g1/"
#base_output_dir="/home/jonas/mnt_AG_Imkeller/group_members/schuck/nanoairr/manuscript/data/subsampling/new/sp7/2_longairr_out/g2/"
#base_output_dir="/home/jonas/mnt_AG_Imkeller/group_members/schuck/nanoairr/manuscript/data/subsampling/new/sp8/2_longairr_out/g3/"
#base_output_dir="/home/jonas/mnt_AG_Imkeller/group_members/schuck/nanoairr/manuscript/data/subsampling/new/sp8/2_longairr_out/g4/"


#base_output_dir="/home/jonas/mnt_AG_Imkeller/group_members/schuck/nanoairr/manuscript/data/subsampling/new/sp12/2_longairr_out/1/"


log_file="${base_output_dir}/subsample.log"

# Create output directory if it doesn't exist
mkdir -p "$base_output_dir"

# this runs once per sample / input group
#----Subsampling of large umi groups and chunking / MSA-parallelization----#
trim_out="${input_fasta}"
trim_trimmed="${base_output_dir}/trimmed.fasta"


seqkit seq -w 0 -m 800 -M 3800 "${trim_out}" -o "${trim_trimmed}" >> "${log_file}" 2>&1 || {
  log_message "ERROR" "longairr collapse" \
    "Failed to filter sequences while executing: \
    seqkit seq -m \"${min_length}\" -M \"${max_length}\" -Q \"${min_qual}\" \"${trim_out}\" -o \"${trim_trimmed}\"" \
    "${log_file}"
  exit 1
}

#NanoPlot --fasta "${trim_trimmed}" -o "${base_output_dir}/nanoplot/after_trim/"

# Record start time before AlignSets
align_start=$(date +%s)
echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') Starting AlignSets for size=${size}, seed=${seed}" >> "${log_file}"

# Loop over sizes and seeds
for size in "${sizes[@]}"; do
    for seed in "${seeds[@]}"; do

        # Record start time before AlignSets
        align_start=$(date +%s)
        echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') Starting AlignSets for size=${size}, seed=${seed}" >> "${log_file}"

        #runtime info
        echo "Running subsample size: $size, seed: $seed"
        
        # Create unique output directory for this run
        run_output_dir="${base_output_dir}/n_${size}/s_${seed}/"
        mkdir -p "$run_output_dir"
        #intermediate file
        subsampled_trimmed="${run_output_dir}/subsampled_groups.fasta"
        failed_filtered="${run_output_dir}/failed_adaptive_filtered.fasta"
        #intermediate dirs
        chunked_reads="${run_output_dir}/chunks/"
        aligned_chunks="${run_output_dir}/aligned_chunks/"
        
        mkdir -p "${chunked_reads}" "${aligned_chunks}"
        mkdir -p "${chunked_reads}/failed_chunks/"


        # NOTE longairr_filter_subsample_groups.py is an internal script saved in `longairr_scripts/`
        longairr_filter_subsample_groups.py -i "${trim_trimmed}" -o "${subsampled_trimmed}" -n ${size} --adaptive-filter TRUE --group-field SPBCUMI --min-group-size 5 --bin-size 100 --peak-margin 0.05 --seed ${seed} --failed "${failed_filtered}" >> "${log_file}" 2>&1 || {
          log_message "ERROR" "longairr collapse" \
            "Failed to subsample sequences while executing: \
            longairr_filter_subsample_groups.py" \
            "${log_file}"
          exit 1
        }

        # NOTE longairr_split_align.py is an internal script saved in `longairr_scripts/`
        longairr_split_align.py "${subsampled_trimmed}" --group-field SPBCUMI --max_umi_per_chunk 35 --outdir "${chunked_reads}" >> "${log_file}" 2>&1 || {
          log_message "ERROR" "longairr collapse" \
            "Failed to create chunks while executing: \
            longairr_split_align.py" \
            "${log_file}"
          exit 1
        }

        parallel_cjobs=5

        # Parallel execution of AlignSets.py
        for chunk in ${chunked_reads}/chunk_*.fasta; do
          chunk_name=$(basename "$chunk" .fasta)
          (
            AlignSets.py muscle -s "${chunk}" --bf SPBCUMI --outname "${chunk_name}" --outdir "${aligned_chunks}" >> "${log_file}" 2>&1

            align_out="${aligned_chunks}/${chunk_name}_align-pass.fasta"
            
            if [[ -s "$align_out" ]]; then
              rm -f "$chunk" "${chunked_reads}/${chunk_name}.fasta"
            else
              log_message "ERROR" "longairr collapse" \
                "AlignSets on $chunk failed or output missing." \
                "${log_file}"
              mv "$chunk" "${chunked_reads}/failed_chunks/"
            fi
          ) &

          # Wait if too many background jobs are running
          while (( $(jobs -rp | wc -l) >= ${parallel_cjobs} )); do
            sleep 1
          done
        done

        wait # wait for all parallel processes to finish

        # Combine all MSA-chunks after all parallel jobs have completed
        for aligned_chunk in ${aligned_chunks}/*_align-pass.fasta; do
          cat "$aligned_chunk" >> "${run_output_dir}/groups_aligned_combined.fasta"
        done

        # Remove duplicates before consensus (delete(?)) - very defensive, but fast
        seqkit rmdup -n "${run_output_dir}/groups_aligned_combined.fasta" > "${run_output_dir}/3_groups_aligned_combined_dedup.fasta"
        rm "${run_output_dir}/groups_aligned_combined.fasta"

        BuildConsensus.py -n 1 -s "${run_output_dir}/3_groups_aligned_combined_dedup.fasta" --bf SPBCUMI --cf N_ORIG N_KEEP UMI X Y --act majority majority majority majority majority --pf SPBC --maxgap 0.5 --maxerror 0.2 --outname groups_msa --outdir "${run_output_dir}" >> "${log_file}" 2>&1 || {
          log_message "ERROR" "longairr collapse" \
            "BuildConsensus spatial failed." \
            "${log_file}"
          exit 1
        }
        
        consensus_out="${run_output_dir}groups_msa_consensus-pass.fasta"
        transformed_out="${run_output_dir}4_longairr_consensus.fasta"

        # Write consensus info (GROUP + CONSCOUNT + N_ORIG + N_KEEP+ UMI + SPBC + PRCOUNT) into sequence header (+ enumerate seqs with unique ID)
        # NOTE longairr_header_metadata.py is an internal script saved in 'longairr_scripts/'
        longairr_header_metadata.py --input "${consensus_out}" -o "${transformed_out}" --library visium >> "${log_file}" 2>&1 || {
          log_message "ERROR" "longairr collapse" \
            "longairr_header_metadata spatial failed." \
            "${log_file}"
          exit 1
        }
        
        # Record end time and compute elapsed duration
        align_end=$(date +%s)
        elapsed=$((align_end - align_start))
        printf "[INFO] %s Finished BuildConsensus for size=%s, seed=%s (elapsed: %02d:%02d:%02d)\n" \
          "$(date '+%Y-%m-%d %H:%M:%S')" "$size" "$seed" \
          $((elapsed/3600)) $((elapsed%3600/60)) $((elapsed%60)) >> "${log_file}"


    done
done