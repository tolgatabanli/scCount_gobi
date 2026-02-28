# Welcome to Progress Log

## 1st day

- Exploring the literature
  - Current tools use k-mer hashing and graphs to efficiently represent and map reads to the genome
  - edf
- Understanding input and Visium
  - Barcodes, UMI (read1), and sequence (read2)
  - Three-prime-end sequencing
- Idea: Use only a certain portion of the transcript sequences at the end to create hash indices

## 2nd day

- Prototyping:
  - Used minimap2, kallisto and cellranger(?) to check the assumption
  - Almost all reads map primarily to the transcript ends up to ca. 400-500 bps before end
- Determining a threshold
  - Only using last exon/UTR is not safe, transcript might have up to many short UTR exons such as 59!
- Program:
  - started with basic utilities (file readers and sequence fetchers), using htsjdk.samtools
  - For hashing, will either rely on Java hashing, or use BBHash or Sux4j tools

## 3rd day
- prototyping:
  - Making a robust comparison by filtering SAM alignments from above tools
- Program
  - Added unit tests to check positive/negative strand transcript seq fetches
  - 
## 4th day 
- prototyping: 
  - startet java implementation and finished gtf reader and fasta reader + adding utils that will be needed later 
- Determining how long the subset for the transcript sequence should be 
  - created bedfiles with the end region of the transcripts 
  - ran samtools bedcov but numbers are weird --> used bedtools multicov instead

## 5th day
- program
  - kmer encoding:
    - abandoned any type of hashing, encoding nucs as A: 00, C: 01, T: 10, G: 11
  - Counting idea:
    - Written the fundamental pseudocode for counting: from kmers create minimizers and map minimizers to transcripts.
  