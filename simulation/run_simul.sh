#!/usr/bin/bash

INPUT=$1
MUTRATE="${2:-0.0}"
OUTPUT="${INPUT%.*}_$MUTRATE"


java -jar scReadSimulator.jar \
	-readcounts $INPUT \
	-fasta /mnt/raidbio2/extstud/praktikum/genprakt-ws25/gruppe_a/data/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
	-fidx /mnt/raidbio2/extstud/praktikum/genprakt-ws25/gruppe_a/data/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.fai \
	-gtf /mnt/raidbio2/extstud/praktikum/genprakt-ws25/gruppe_a/data/Sus_scrofa.Sscrofa11.1.gtf \
	-od $OUTPUT \
	-length 75 -tailLength 500 -mutationrate $MUTRATE
