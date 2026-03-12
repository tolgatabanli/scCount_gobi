#!/usr/bin/bash

java -jar /home/t/tabanli/Desktop/scCount/out/artifacts/scCount_jar/scCount.jar index \
	-f /mnt/raidbio2/extdata/praktikum/genprakt/genprakt-ws25/Block/pig-genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz \
	-g /mnt/raidbio2/extdata/praktikum/genprakt/genprakt-ws25/Block/pig-genome/Sus_scrofa.Sscrofa11.1.115.chr.gtf.gz \
	-k 20 \
	-idx /mnt/cip/home/t/tabanli/Desktop/scCount/out/sccount_utrplus.idx

