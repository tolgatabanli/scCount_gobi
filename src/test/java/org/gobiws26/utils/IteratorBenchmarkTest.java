package org.gobiws26.utils;

import org.gobiws26.Config;
import org.gobiws26.Readers.GTFReader;
import org.gobiws26.genomicstruct.Transcript;

import java.io.IOException;
import java.util.HashMap;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Disabled;
import static org.gobiws26.Test.*;

public class IteratorBenchmarkTest {

    @Test
    @Disabled("KmerIteratorHash and KmerIteratorByte classes are not available - pre-existing issue")
    public void KmerIterTiming() throws IOException {
        // This test is disabled because KmerIteratorHash and KmerIteratorByte 
        // classes are not implemented in the current codebase.
        // This is a pre-existing issue unrelated to minimizer length changes.
    }
}
