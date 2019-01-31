#! /bin/bash

CHOCOPHLAN_PATH=/shares/CIBIO-Storage/CM/scratch/users/francesco.beghini/hg/chocophlan
bt2_path=/shares/CIBIO-Storage/CM/mir/tools/bowtie2-2.3.4.3
OUTFOLDER=${CHOCOPHLAN_PATH}/eval
gnu_parallel=/shares/CIBIO-Storage/CM/mir/tools/parallel-20171022/bin/parallel
markers_chunks=/shares/CIBIO-Storage/CM/scratch/users/francesco.beghini/hg/chocophlan/export/metaphlan2/mark/chunks_merged.fna
mpa2_markers=/shares/CIBIO-Storage/CM/scratch/users/francesco.beghini/hg/chocophlan/export/metaphlan2/mpa_v24_CHOCOPhlAn_0.2.1/mpa_v24_CHOCOPhlAn_0.2.1.fna

mkdir ${OUTFOLDER}
mkdir ${OUTFOLDER}/bin_chunks
mkdir ${OUTFOLDER}/bt2_idx
mkdir ${OUTFOLDER}/bt2_out

python3 ${CHOCOPHLAN_PATH}/src/split_markers.py ${mpa2_markers} ${OUTFOLDER}/mpa_v24_CHOCOPhlAn_0.2.1.chunks.fna

find ${CHOCOPHLAN_PATH}/data/ncbi -name *fna.gz | shuf | split -l 1000 - ${OUTFOLDER}/bin_chunks/bin

cat ${OUTFOLDER}/bin_chunks/* | xargs -I '{}'  bash -c "fn={};zcat {} | grep '>'  | cut -b 2- | cut -f1 -d ' '| awk -v fn={} '{print $1; print fn;}'" | paste - - > ${OUTFOLDER}/contig2genome.tsv &

${gnu_parallel} -j 5 "cat {} | xargs zcat > ${OUTFOLDER}/{/}.fna && ${bt2_path}/bowtie2-build --threads 20 ${OUTFOLDER}/{/}.fna ${OUTFOLDER}/bt2_idx/{/} && rm ${OUTFOLDER}/{/}.fna" ::: ${OUTFOLDER}/bin_chunks/*

${gnu_parallel} -j 20 "${bt2_path}/bowtie2 -f ${markers_chunks} \
                                          -x ${OUTFOLDER}/bt2_idx/{/} \
                                          -a --very-sensitive \
                                          --no-unal --threads 20 \
                                          -S ${OUTFOLDER}/bt2_out/{/}.sam 2> ${OUTFOLDER}/bt2_out/{/}.stats" :::  ${OUTFOLDER}/bin_chunks/* &

wait

awk '/^@/ {next}; {if($6 == "150M"){print $3, $1}}' ${OUTFOLDER}/bt2_out/*.sam 