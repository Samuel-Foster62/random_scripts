#!/bin/bash


# Define the source and destination directories
src_dir=$(pwd)
BASE=/PATH/TO/DATA

for dataset in group_1 group_2
do
    echo "Processing dataset: $dataset"
    images=$(ls -1v ${BASE}/${dataset}/*.nxs)
    for wedge in $images
    do
        echo $wedge
        name=`basename $wedge | cut -d '.' -f 1`
        date=`echo $wedge | cut -d '/' -f 7`
        echo $date
        directory=$src_dir/reproc/${dataset}/$name
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

#xia2 image=$wedge \
#  remove_blanks=True \
#  failover=True \
#  ice_rings.unit_cell=3.615,3.615,3.615,90,90,90 \
#  ice_rings.space_group=fm-3m \
#  ice_rings.width=0.01 \
#  ice_rings.filter=True \
#  find_spots.phil_file=/dls/science/groups/i02-1/scripts/processing/spots.phil
dials.import $wedge geometry.detector.distance=313.5
dials.find_spots /dls/science/groups/i02-1/scripts/processing/spots.phil imported.expt d_min=2.2 d_max=20
dials.index strong.refl imported.expt
dials.refine indexed.* scan_varying=True
dials.integrate refined.*

EOF

    chmod +x $directory/run.sh
        # Submit the job to Slurm
        cd $directory

        sbatch "run.sh"

        echo "Submitted Slurm job for $dataset_name"
    done
done

