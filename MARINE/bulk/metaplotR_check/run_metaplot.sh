/home/tmhaxs421/brannanlab/tmhaxs421/scripts/ScriptsToolkit/MARINE/bulk/helper_metaplotdist.py \
--genePred /home/tmhaxs421/brannanlab/tmhaxs421/scripts/ScriptsToolkit/MARINE/bulk/refseq/wgEncodeGencodeBasicV19.txt.gz \
--bed /home/tmhaxs421/brannanlab/tmhaxs421/scripts/ScriptsToolkit/MARINE/bulk/metaplotR_check/m6a.sorted.bed \
--out /home/tmhaxs421/brannanlab/tmhaxs421/scripts/ScriptsToolkit/MARINE/bulk/metaplotR_check/m6a.dist.measures.txt


/home/tmhaxs421/brannanlab/tmhaxs421/scripts/ScriptsToolkit/MARINE/bulk/helper_plot_metagene.py \
/home/tmhaxs421/brannanlab/tmhaxs421/scripts/ScriptsToolkit/MARINE/bulk/metaplotR_check/m6a.dist.measures.txt \
--labels m6A \
--out /home/tmhaxs421/brannanlab/tmhaxs421/scripts/ScriptsToolkit/MARINE/bulk/metaplotR_check/m6a_metagene_plot.png \
--dpi 300 \
--title "m6A metagene distribution" \
--figsize 6,4 
--rescale \

