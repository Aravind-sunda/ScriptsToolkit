# run this script to export the conda environments to the yaml files in this current directory
module load mamba
mamba activate

conda env list | awk 'NR>2 && $1 !~ /^\// && NF>0 {print $1, $NF}' | \
while read -r name path; do mamba env export -p "$path" --no-builds > "${name}.yaml"; done
