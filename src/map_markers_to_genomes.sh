#! /bin/bash
export PATH=/shares/CIBIO-Storage/CM/mir/tools/bowtie2-2.3.4.3:$PATH
CHOCOPHLAN_PATH=/shares/CIBIO-Storage/CM/scratch/users/francesco.beghini/hg/chocophlan
bt2_path=/shares/CIBIO-Storage/CM/mir/tools/bowtie2-2.3.4.3
OUTFOLDER=${CHOCOPHLAN_PATH}/eval
gnu_parallel=/shares/CIBIO-Storage/CM/mir/tools/parallel-20171022/bin/parallel
mpa2_markers=/shares/CIBIO-Storage/CM/scratch/users/francesco.beghini/hg/chocophlan/export/metaphlan2/mpa_v25_CHOCOPhlAn_0.2/mpa_v25_CHOCOPhlAn_0.2.fna

mkdir ${OUTFOLDER}
mkdir ${OUTFOLDER}/bin_chunks
mkdir ${OUTFOLDER}/bt2_idx
mkdir ${OUTFOLDER}/bt2_out


find ${CHOCOPHLAN_PATH}/data/ncbi -name *fna.gz | shuf | split -l 1000 - ${OUTFOLDER}/bin_chunks/bin

cat ${OUTFOLDER}/bin_chunks/* | xargs -I '{}'  bash -c "fn={};zcat {} | grep '>'  | cut -b 2- | cut -f1 -d ' '| awk -v fn={} '{print $1; print fn;}'" | paste - - > ${OUTFOLDER}/contig2genome.tsv &

${gnu_parallel} -j 5 "cat {} | xargs zcat > ${OUTFOLDER}/{/}.fna && ${bt2_path}/bowtie2-build --threads 20 ${OUTFOLDER}/{/}.fna ${OUTFOLDER}/bt2_idx/{/} && rm ${OUTFOLDER}/{/}.fna" ::: ${OUTFOLDER}/bin_chunks/*

python3 ${CHOCOPHLAN_PATH}/src/split_markers.py ${mpa2_markers} ${mpa2_markers/.fna/.chunks.fna}
${gnu_parallel} -j 20 "${bt2_path}/bowtie2 -f ${mpa2_markers/.fna/.chunks.fna} \
                                          -x ${OUTFOLDER}/bt2_idx/{/} \
                                          -a --very-sensitive \
                                          --no-hd --no-sq \
                                          --no-unal --threads 4 \
                                          -S ${OUTFOLDER}/bt2_out/{/}.sam" :::  ${OUTFOLDER}/bin_chunks/*