package org.gobiws26.Indexing;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import org.junit.jupiter.api.Test;

import static org.gobiws26.Test.*;
import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;

public class IndexerTests {
    @Test
    public void basicMinimizersFoundCorrectly() {
        // Test disabled: Old implementation stored tx->minimizers; new implementation stores triplet->txIds
        // This is now checked indirectly through query tests which verify triplet-based retrieval works.
        // Keeping test method signature for compatibility.
    }
}
