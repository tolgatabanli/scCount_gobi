# Welcome to Progress Log

## 1st day
- Exploring the literature
  - Current tools use k-mer hashing, minimizers and/or deBruijn graphs to efficiently index genome and map reads
- Understanding input and Visium
  - Barcodes, UMI (read1), and sequence (read2)
  - Three-prime-end sequencing -> no need to consider strandedness
- Idea: Use only a certain portion of the transcript sequences at the end to create the index/mapping

## 2nd day
- Prototyping:
  - Used minimap2, kallisto and cellranger(?) to check the assumptions
  - Almost all reads map primarily to the transcript ends up to ca. 400-500 bps before the tail end
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
  - started Java implementation: GTF reader, transcriptomic seq fetcher (uses htsjdk.samtools), and other utils for KmerIterator 
- Determining how long tail sequence for transcripts should be: 
  - created BED files with the end region of the transcripts 
  - ran 'samtools bedcov' but numbers are weird --> used bedtools multicov -> ~ 500 bps is the elbow point

## 5th day
- For efficient mapping, we could use deBruijn graphs or Minimizer logic
- program
  - kmer encoding:
    - abandoned any type of hashing, encoding nucs as A: 00, C: 01, T: 10, G: 11 (least operations, based on char byte)
    - will use 'fastutils' for primitive type Sets, Iterators and Maps
- prototyping:
  - checked the distribution of the uniquely identifiable transcript ratio per gene based on last X bps of transcript sequences to see if the heuristic is theoretically promising

## 6th day
- program
  - Counting idea:
    - Written the fundamental pseudocode for counting: from kmers create minimizers and map minimizers to transcripts.
  - how to break ties and ambiguities in the mapping?
    - For a read, try to minimize the set of mapped transcripts of minimizers
    - Report the ambiguities in a summary file 
