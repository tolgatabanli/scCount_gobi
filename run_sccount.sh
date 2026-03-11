#!/usr/bin/bash

java -jar /home/t/tabanli/Desktop/scCount/out/artifacts/scCount_jar/scCount.jar count \
	-o /home/t/tabanli/Desktop/scCount/out/count_out_simul_noMut \
	-idx /home/t/tabanli/Desktop/scCount/out/sccount_utrplus.idx \
	-r2 /mnt/cip/home/t/tabanli/Desktop/scCount/simulation/out_noMut/read2.fastq \
	-batchSize 500000 \
	-threads 8
