#!/usr/bin/bash

java -jar /home/t/tabanli/Desktop/scCount/out/artifacts/scCount_jar/scCount.jar count \
	-o /home/t/tabanli/Desktop/scCount/out/count_out_simul_noMut_gamma_protcoding \
	-idx /home/t/tabanli/Desktop/scCount/out/sccount_utrplus.idx \
	-r2 /mnt/cip/home/t/tabanli/Desktop/scCount/simulation/out_noMut_gamma_protcoding/read2.fastq \
	-batchSize 500000 \
	-threads 30
