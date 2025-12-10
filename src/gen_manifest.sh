
unzip -l -M raw_data/Stover_250204_Pietrakiak_demultiplxed-20251030T183257Z-1-001.zip | tail -n +4 | sed '$d' | sed '$d' | tr -s ' ' '\t' | cut -f 2,5 > raw_data/manifest.tsv
unzip -l -M raw_data/Stover_250204_Pietrakiak_demultiplxed-20251030T183257Z-1-002.zip | tail -n +4 | sed '$d' | sed '$d' | tr -s ' ' '\t' | cut -f 2,5 >> raw_data/manifest.tsv
unzip -l -M raw_data/Stover_250204_Pietrakiak_demultiplxed-20251030T183257Z-1-003.zip | tail -n +4 | sed '$d' | sed '$d' | tr -s ' ' '\t' | cut -f 2,5 >> raw_data/manifest.tsv

echo "manifest generated, file count: "
echo $(wc -l raw_data/manifest.tsv)
