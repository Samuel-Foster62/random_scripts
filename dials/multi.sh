#!/bin/bash

# Define the source and destination directories
src_dir=$(pwd)

i=1
while [[ -d "`printf %03d $i`" ]]; do
  let i++
done
subd="`printf %03d $i`"
newdir="$src_dir/$subd"
echo $subd
mkdir -p $newdir
cd $newdir || exit 1

  # Create the 'run.sh' file inside the new directory
  cat << EOF > "$newdir/run.sh"
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=24G
#SBATCH --time=23:00:00
#SBATCH --error=j_%J.err
#SBATCH --output=j_%J.out
#SBATCH --partition=cs04r

module load dials/latest
src_dir=$1
xia2.multiplex \\
  $(find "$src_dir" -type f \( -name "integrated.expt" -o -name "integrated.refl" \) -print | tr '\n' ' ') \\
  filtering.method=deltacchalf \\
  deltacchalf.stdcutoff=3 \\
  deltacchalf.mode=image_group \\
  deltacchalf.group_size=20 \\
  clustering.method=hierarchical \\
  clustering.output_clusters=False \\
  max_cluster_height=0.01 \\
  min_cluster_size=3
EOF

chmod +x $newdir/run.sh

  # Submit the job to Slurm
  cd $newdir

  sbatch "$newdir/run.sh"

  echo "Submitted Slurm job for $subd"
