package org.gobiws26;

public class Config {
    public static int K;
    public static int minimLength = 8; // Default minimizer length

    /**
     * Validates that minimLength is set correctly relative to K.
     * Must be called after both K and minimLength are set.
     * @throws IllegalArgumentException if minimLength >= K
     */
    public static void validateMinimizerLength() {
        if (minimLength >= K) {
            throw new IllegalArgumentException(
                String.format("minimLength (%d) must be less than kmerLength (K=%d)", minimLength, K)
            );
        }
    }
}
