package org.gobiws26.Indexing;

import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import org.gobiws26.Querying.IndexGraph;

public class IndexData {
    public final IndexGraph graph;
    public final Int2ObjectOpenHashMap<String> int2TxString;
    public final Int2ObjectOpenHashMap<String> int2GeneString;
    public final Int2IntOpenHashMap txInt2GeneInt;

    public IndexData(IndexGraph graph,
                     Int2ObjectOpenHashMap<String> int2TxString,
                     Int2ObjectOpenHashMap<String> int2GeneString,
                     Int2IntOpenHashMap txInt2GeneInt) {
        this.graph = graph;
        this.int2TxString = int2TxString;
        this.int2GeneString = int2GeneString;
        this.txInt2GeneInt = txInt2GeneInt;
    }
}
