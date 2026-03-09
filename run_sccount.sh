#!/usr/bin/bash

java -jar /home/t/tabanli/Desktop/scCount/out/artifacts/scCount_jar/scCount.jar count \
	-o /home/t/tabanli/Desktop/scCount/out/count_out \
	-idx /home/t/tabanli/Desktop/scCount/out/sccount.idx \
	-r2 /mnt/raidbio2/extdata/omics/SeqReads/pig_aorta/visium/read2_Vis_36-15.fastq.gz \
	-batchSize 1000000 \
	-threads 24
