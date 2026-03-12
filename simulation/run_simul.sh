#!/usr/bin/bash

java -jar scReadSimulator.jar \
	-readcounts random_counts.tsv \
	-fasta /mnt/raidbio2/extstud/praktikum/genprakt-ws25/gruppe_a/data/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
	-fidx /mnt/raidbio2/extstud/praktikum/genprakt-ws25/gruppe_a/data/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.fai \
	-gtf /mnt/raidbio2/extstud/praktikum/genprakt-ws25/gruppe_a/data/Sus_scrofa.Sscrofa11.1.gtf \
	-od out_noMut_gamma_protcoding \
	-length 75 -tailLength 500 -mutationrate 0.0
