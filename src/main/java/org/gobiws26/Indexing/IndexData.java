package org.gobiws26.Indexing;

import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import org.gobiws26.Querying.IndexGraph;

public record IndexData(IndexGraph graph, Int2ObjectOpenHashMap<String> int2TxString,
                        Int2ObjectOpenHashMap<String> int2GeneString, Int2IntOpenHashMap txInt2GeneInt) {
}
