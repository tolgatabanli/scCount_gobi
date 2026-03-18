#!/usr/bin/bash

java -Xmx90g -jar /mnt/cip/home/t/tabanli/Desktop/scCount/out/artifacts/miniqut3_jar/miniqut3.jar count \
	-o /home/t/tabanli/Desktop/scCount/out/count_with_barcodes2 \
	-idx /mnt/cip/home/t/tabanli/Desktop/scCount/benchmark/idx/sccount_k31_g8.idx \
	-r1 /mnt/raidbio2/extdata/omics/SeqReads/pig_aorta/visium/read1_Vis_66-10.fastq.gz \
	-r2 /mnt/raidbio2/extdata/omics/SeqReads/pig_aorta/visium/read2_Vis_66-10.fastq.gz \
	-details \
	-batchSize 1000000 \
	-threads 30
