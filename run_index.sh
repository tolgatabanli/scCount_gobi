#!/usr/bin/bash

java -jar /mnt/cip/home/t/tabanli/Desktop/scCount/out/artifacts/miniqut3_jar/miniqut3.jar index \
	-f /mnt/raidbio2/extdata/praktikum/genprakt/genprakt-ws25/Block/pig-genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz \
	-g /mnt/raidbio2/extdata/praktikum/genprakt/genprakt-ws25/Block/pig-genome/Sus_scrofa.Sscrofa11.1.115.chr.gtf.gz \
	-k 30 -minimLength 12 \
	-idx /mnt/cip/home/t/tabanli/Desktop/scCount/out/miniqut3_k30_g12.idx

