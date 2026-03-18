#!/usr/bin/Rscript

library(DiagrammeR)
library(htmlwidgets)

petri_net_viz <- grViz("
digraph petri_net {
    # Decreased nodesep (between nodes) and ranksep (between columns)
    graph [rankdir = LR, nodesep = 0.2, ranksep = 0.4]

    # ============ PLACES (Circles: Data/States) ============
    node [shape = circle, fixedsize = true, width = 1.5, penwidth = 2, fontsize = 18]
    
    # Inputs & Params
    node [color = '#7E57C2'] 
    GTF; FASTA; FASTQ;
    K; MinimLength;
    
    # Intermediate States
    TranscriptTailSeqs [label = 'Transcript\\nTail Seqs']; 
    ReadMinims [label = 'Read\\nMinims'];
    Phase1Candidates [label = 'Phase 1\\nCandidates'];
    
    # Logic Results
    IndexFile [label = 'Index\\nFile']; 
    Unambiguous; 
    Ambiguous; 
    
    # Final Outputs
    CountsFile [label = 'Counts\\nFile']; 
    SummaryFile [label = 'Summary\\nFile'];

    # ============ TRANSITIONS (Rectangles: Processes) ============
    node [shape = box, fixedsize = false, width = 1.5, height = 1, 
          style = filled, fillcolor = '#E3F2FD', color = '#1976D2', fontsize = 18]
    
    FetchSeq [label = 'Fetch\\nSeq']; 
    ExtractMin1 [label = 'Extract Minims\\nWrite Index']; 
    ExtractMin2 [label = 'Extract\\nMinims']; 
    Phase1 [label = 'Phase 1']; 
    Phase2 [label = 'Phase 2']; 
    WriteOut [label = 'Write\\nOutput'];

    # ============ EDGES ============
    # Indexing Workflow
    GTF -> FetchSeq
    FASTA -> FetchSeq
    FetchSeq -> TranscriptTailSeqs
    
    TranscriptTailSeqs -> ExtractMin1
    K -> ExtractMin1
    MinimLength -> ExtractMin1
    ExtractMin1 -> IndexFile
    
    # Counting Workflow
    IndexFile -> Phase1
    FASTQ -> ExtractMin2
    K -> ExtractMin2
    MinimLength -> ExtractMin2
    ExtractMin2 -> ReadMinims
    ReadMinims -> Phase1
    
    Phase1 -> Unambiguous
    Phase1 -> Phase1Candidates
    
    Phase1Candidates -> Phase2
    Phase2 -> Unambiguous
    Phase2 -> Ambiguous
    
    Unambiguous -> WriteOut
    Ambiguous -> WriteOut
    WriteOut -> CountsFile
    WriteOut -> SummaryFile

    # Formatting
    label = '';
    labelloc = 't';
    fontsize = 16;
}
")

saveWidget(petri_net_viz, "petri_net.html", selfcontained = TRUE)