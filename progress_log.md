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
## 4th day 
- prototyping: 
  - started Java implementation: GTF reader, transcriptomic seq fetcher (uses htsjdk.samtools), and other utils for KmerIterator 
- Determining how long tail sequence for transcripts should be: 
  - created BED files with the end region of the transcripts 
  - ran 'samtools bedcov' but numbers are weird --> used bedtools multicov -> ~ 500 bps is the elbow point
  - try to apply an exon filter for reads mapped by star and came to conclusion that the a only a few % were left 

## 5th day
- For efficient mapping, we could use deBruijn graphs or Minimizer logic
- program
  - kmer encoding:
    - abandoned any type of hashing, encoding nucs as A: 00, C: 01, T: 10, G: 11 (least operations, based on char byte)
    - will use 'fastutils' for primitive type Sets, Iterators and Maps
- prototyping:
  - checked the distribution of the uniquely identifiable transcript ratio per gene based on last X bps of transcript sequences to see if the heuristic is theoretically promising
- looked at average intron/junction size at the end of transcript to determine what the filter for junctions should be and trying to finalize the program logic (pseudo code and todo)

## 6th day
- program
  - Counting idea:
    - Written the fundamental pseudocode for counting: from kmers create minimizers and map minimizers to transcripts.
  - how to break ties and ambiguities in the mapping?
    - For a read, try to minimize the set of mapped transcripts of minimizers
    - Report the ambiguities in a summary file 
- looked at bifrost for the effcient de buijn graph construction 
- looked more at the implementation of kallisto to get a better idea of how to implement more efficiently 
## 7th day 
- tried to find out how space and cell ranger calls the alignment in the pipeline 
- looking through the mapping of space ranger (hella weird)
- looking through literature about proof of concept 
- found that the coverage plot is not counting uniquely for reads but for each transcript 
    - Report the ambiguities in a summary file

## 8th day 
- fixed coverage plot for the with merged bedfiles 
- downloaded human 5k_pbmc data from 10x genomics mapped it with star to look for coverage 
## 9th day 
- looked though the counts of the bam files using feature Counts with the default parameters 
## 10th day 
- looked through literature trying to explain the count differences, found some similar comparisons 
- tried the nextflow pipeline, failed because docker is not installed on the server 
## 11th day 
- trying different parameters for feature counts to see how they affect the counts (minoverlap fracoverlap)
  - spaceranger and star counts are more simliar with 70% overlap threshold 
  - the higher counts of spaceranger are eliminated by 70% threshold since the it mapps less strict to exon boundaries 
## 12th day 
- started read tracker program 

## 14th day
- program
  - got first count results!
  - 90% mapping rate (4635725 ambig, 47713064 total)
  - 0.95 correlation against STAR, but data points are not scatter well: the developed tool is biased towards capturing more counts (possibly an artifact of greedy rescue strategy)
  - processing 47M reads took ~ 20 minutes! -> gotta

## 15th day
- program
  - noticed that the genes with counts of ~ 4000x ours, where ours had reported < 100, are those with long UTRs 
  - and those that had 20x - 30x counts where our counts were > 100 have mapped reads with conserved mismatches, especially where A mutated -> changes minimizers
