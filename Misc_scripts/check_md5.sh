
INPUT="/home/tmhaxs421/brannanlab/tmhaxs421/riboSTAMP_mouse/RawData"
for d in ${INPUT}/*; do
    if [ -f "$d/MD5.txt" ]; then
        echo "🔍 Checking: $d"
        (cd "$d" && md5sum -c MD5.txt)
    fi
done