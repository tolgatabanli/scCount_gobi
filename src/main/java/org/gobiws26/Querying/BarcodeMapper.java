package org.gobiws26.Querying;

import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.Map;

/**
 * Thread-safe mapper for spatial barcode strings to integer identifiers.
 */
public class BarcodeMapper {
    private final Map<String, Integer> barcodeToInt = new ConcurrentHashMap<>();
    private final Map<Integer, String> intToBarcode = new ConcurrentHashMap<>();
    private final AtomicInteger nextId = new AtomicInteger(0);

    /**
     * Get or create an integer ID for a barcode string (thread-safe).
     * @param barcode the barcode string
     * @return the integer ID for this barcode
     */
    public int getId(String barcode) {
        return barcodeToInt.computeIfAbsent(barcode, k -> {
            int id = nextId.getAndIncrement();
            intToBarcode.put(id, barcode);
            return id;
        });
    }

    /**
     * Get the barcode string for an integer ID.
     * @param id the integer ID
     * @return the barcode string, or null if not found
     */
    public String getBarcode(int id) {
        return intToBarcode.get(id);
    }

    /**
     * Get the total number of unique barcodes seen.
     * @return the number of unique barcodes
     */
    public int size() {
        return nextId.get();
    }

    /**
     * Get all barcode IDs.
     * @return array of all barcode IDs in order
     */
    public int[] getAllIds() {
        int count = nextId.get();
        int[] ids = new int[count];
        for (int i = 0; i < count; i++) {
            ids[i] = i;
        }
        return ids;
    }
}
