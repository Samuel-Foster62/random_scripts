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

# File containing paths/templates to exclude
exclude_list="$src_dir/exclude_templates.txt"

cat << 'EOL' > "$exclude_list"
/PATH/TO/DATA/x_23.nxs
/PATH/TO/DATA/x_73.nxs
/PATH/TO/DATA/x_7.nxs
/PATH/TO/DATA/x_20.nxs
/PATH/TO/DATA/x_8.nxs
/PATH/TO/DATA/x_19.nxs
.
.
.
EOL

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
exclude_list="\$src_dir/exclude_templates.txt"
exclude_patterns="\$PWD/exclude_patterns.txt"

# Use Python to generate the list of files to include
python3 << 'PYEOL'
import os
from pathlib import Path

src_dir = "PATH/TO/DATA/PROCESSING"
exclude_file = "PATH/TO/DATA/PROCESSING/exclude_templates.txt"

# 1. Build a list of specific directory fragments to exclude
# Converts "/.../e1d1/x_5.nxs" into "e1d1/x_5"
exclude_patterns = []
with open(exclude_file, 'r') as f:
    for line in f:
        path_obj = Path(line.strip())
        if path_obj.suffix == '.nxs':
            family = path_obj.parent.name
            dataset = path_obj.stem
            # We add slashes to ensure we don't match x_5 against x_55
            exclude_patterns.append(f"/{family}/{dataset}/")

print(f"Excluding patterns: {exclude_patterns}")

# 2. Find files and filter
included_files = []
for root, dirs, files in os.walk(src_dir):
    for file in files:
        if file in ["integrated.expt", "integrated.refl"]:
            full_path = os.path.join(root, file)
            
            # If ANY pattern (like /e1d1/x_5/) is in the path, skip it
            if any(p in full_path for p in exclude_patterns):
                continue
                
            included_files.append(full_path)

# 3. Save to file
with open("files_for_xia2.txt", "w") as f:
    f.write(" ".join(included_files))
PYEOL


xia2.multiplex \\
  \$(cat files_for_xia2.txt) \\
  filtering.method=deltacchalf \\
  deltacchalf.stdcutoff=0.15 \\
  deltacchalf.mode=image_group \\
  deltacchalf.group_size=10 \\
  clustering.method=hierarchical \\
  clustering.output_clusters=True \\
  max_cluster_height=0.01 \\
  min_cluster_size=3
  # cosym.best_monoclinic_beta=True \\
  # resolution.d_min=2.5 \\
  # symmetry.cosym.relative_length_tolerance=0.1 \\
EOF

chmod +x $newdir/run.sh

  # Submit the job to Slurm
  cd $newdir

  sbatch "$newdir/run.sh"

  echo "Submitted Slurm job for $subd"
