#!/bin/bash
for var in $(ls ./alignment_ht2) #$(ls *_1.fq.gz)
do
echo $var
chmod +x /software/RNA/circRNA/call_cirRNA/circRNA_finder/starCirclesToBed.pl
chmod +x /software/RNA/circRNA/call_cirRNA/circRNA_finder/filterSpliceSiteCircles.pl
chmod +x /software/RNA/circRNA/call_cirRNA/circRNA_finder/nrForwardSplicedReads.pl
chmod +x /software/RNA/circRNA/call_cirRNA/circRNA_finder/postProcessStarAlignment.pl
sample=${var}  #%_1.fq.gz*};
dir_alignment='/work/circRNA/alignment_ht2'
dir_circfinder='/work/circRNA/CircFinder'
mkdir -p ${dir_circfinder}/${sample}
fa='/software/RNA/circRNA/call_cirRNA/ref/fa/gencode_v28_genome.fa'
gtf='/software/RNA/circRNA/call_cirRNA/ref/fa/gencode_v28_annotation.gtf'
star='/software/miniconda3/envs/RNAseq_py37/bin/STAR'
gdir='/software/RNA/circRNA/call_cirRNA/ref/fa/star_index'
unmap_fq=${dir_alignment}/${sample}/${sample}'_unmapped.fastq'
${star} --genomeDir ${gdir} --readFilesIn ${unmap_fq} --runThreadN 40 --chimSegmentMin 20 --chimScoreMin 1 --alignIntronMax 100000 \
--chimOutType Junctions SeparateSAMold --outFilterMismatchNmax 4 --alignTranscriptsPerReadNmax 100000 --outFilterMultimapNmax 2 --outFileNamePrefix ${dir_circfinder}/${sample}/${sample}
pstalign='/software/RNA/circRNA/call_cirRNA/circRNA_finder/postProcessStarAlignment.pl'
star_align=${dir_circfinder}/${sample}
circ_finder_out=${dir_circfinder}/${sample}
cd ${star_align}
perl ${pstalign} ${star_align} ${circ_finder_out}
done
