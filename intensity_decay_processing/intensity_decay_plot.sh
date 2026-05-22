#!/bin/bash

# Define the source and destination directories
#images=(/dls/mx/data/nr27313/nr27313-464/BIns/20250814_5um_9*.nxs)
images=(/dls/mx/data/nr27313/nr27313-464/BIns/*_*_*_*.nxs)
src_dir=$(pwd)
OUTPUT="doseVsIntensity_dmin.csv"
DOSE_PER_BATCH=5

COMBINED_OUTPUT="combined_doseVsIntensity_dmin.csv"

declare -A initialised_groups

echo 'Dose,Mean_Intensity,Group-ID' > "$src_dir/$COMBINED_OUTPUT"
echo "Combined output file initialised at $src_dir/$COMBINED_OUTPUT"

for wedge in "${images[@]}"; do
  echo "Processing file: $wedge"

  base_name=$(basename "$wedge" .nxs)
  IFS='_' read -r date_part size index batch <<< "$base_name"
  output_dir="$src_dir/dmin/$date_part/$size/$index"
  group_id="${date_part}_${size}_${index}"

  if [[ -z "${initialised_groups[$group_id]}" ]]; then
    echo "New group detected: $group_id. Initialising output file."
    mkdir -p "$output_dir"
    echo '"Dose, Mean_Intensity"' > "$output_dir/$OUTPUT"
    initialised_groups[$group_id]=1 #mark as initialised
    #output_dir="$src_dir/$date_part/$size/$index"
    #echo "Final output will be saved in: $output_dir"
  fi
  #directory="$src_dir/$date_part/$size/$index/$batch"
  directory="$output_dir/$batch"

  if [ -d $directory ]
  then
    echo "$directory exists: skipping"
    continue
  fi

  mkdir -p $directory
  cd $directory

  # Create the 'run.sh' file inside the new directory
  cat > run.sh << EOF
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --error=j_%J.err
#SBATCH --output=j_%J.out
#SBATCH -p cs04r

module load dials
#module load dials/latest

dials.import "$wedge"
#dials.find_spots "$src_dir/spots.phil" ice_rings.unit_cell=3.615,3.615,3.615,90,90,90 ice_rings.space_group=fm-3m ice_rings.width=0.01 ice_rings.filter=True imported.expt
dials.find_spots "$src_dir/spots.phil" d_min=2.1 mask=$src_dir/pixels.mask imported.expt
dials.index strong.refl imported.expt unit_cell=78,78,78,90,90,90 space_group=I23
dials.integrate d_min=2.1 mask=$src_dir/pixels.mask indexed.*
dials.python "$src_dir/get_avg_integrated_prf.py" integrated.refl > mean_prf_mask.txt
tail -n 1 mean_prf_mask.txt > tmp && mv tmp mean_prf_mask.txt 
EOF

    chmod +x run.sh
    # Submit the job to Slurm
    sbatch "run.sh"
    cd "$src_dir"
    echo "Submitted Slurm job for $base_name"

done

echo "Monitoring dials submissions to slurm for user $USER"

#if [ -n "$output_dir" ]; then
#    echo "Initialising output file: $output_dir/$OUTPUT"
#    echo '"Dose_(MGy)","Mean_Intensity"' > "$output_dir/$OUTPUT"
#else
#    echo "No files were processed, exiting."
#    exit 0
#fi

while true; do
    JOB_COUNT=$(squeue -u "$USER" -h | wc -l)
    if [ "$JOB_COUNT" -eq 0 ]; then
        echo "All jobs have finished executing".
        find "$src_dir" -type f -name "mean_prf.txt" | while read -r FILE; do
            PARENT_DIR=$(dirname "$FILE")
            INDEX_DIR=$(dirname "$PARENT_DIR")
            BATCH_NUM=$(basename "$PARENT_DIR")

            DOSE=$(( BATCH_NUM * DOSE_PER_BATCH ))
            OUTPUT_CSV="$INDEX_DIR/$OUTPUT"
            VALUE=$(head -n 1 "$FILE")

            #echo "$DOSE,$VALUE" >> "$output_dir/$OUTPUT"
            echo "$DOSE,$VALUE" >> "$OUTPUT_CSV"

            #Reconstruct the group_ID from the directory path
            RELATIVE_GROUP_PATH="${INDEX_DIR#"$src_dir/"}"
            RELATIVE_GROUP_ID=$(echo "$RELATIVE_GROUP_PATH" || tr '/' '_')

            #Append to combined file
            echo "$DOSE,$VALUE,$RELATIVE_GROUP_ID" >> "$src_dir/$COMBINED_OUTPUT"
        done
        echo "All results have been compiled! - location: $src_dir/$COMBINED_OUTPUT"
        break
    else
        echo "Still running"
        sleep 60
    fi
done

