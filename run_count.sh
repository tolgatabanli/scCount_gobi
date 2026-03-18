#!/usr/bin/bash

java -Xmx90g -jar /mnt/cip/home/t/tabanli/Desktop/scCount/out/artifacts/miniqut3_jar/miniqut3.jar count \
	-o /home/t/tabanli/Desktop/scCount/out/pig_counts \
	-idx /home/t/tabanli/Desktop/scCount/benchmark/idx/sccount_k23_g11.idx \
	-r2 /mnt/raidbio2/extdata/omics/SeqReads/pig_aorta/visium/read2_Vis_36-15.fastq.gz \
	-batchSize 1000000 \
	-threads 32
