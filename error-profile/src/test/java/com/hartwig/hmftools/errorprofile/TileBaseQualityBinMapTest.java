package com.hartwig.hmftools.errorprofile;

import static org.junit.Assert.assertEquals;

import java.util.Map;

import org.junit.Test;

import junit.framework.TestCase;

public class TileBaseQualityBinMapTest
{
    @Test
    public void testTileBaseQualityBinMap()
    {
        TileBaseQualityBinMap map = new TileBaseQualityBinMap();

        boolean firstOfPair = false;
        TileBaseQualityBin bin = new TileBaseQualityBin("A90924:7:HKKAHXX:1", 1, 1101, firstOfPair, 25, (byte)11);
        map.incrementCount(bin, false);
        map.incrementCount(bin, false);
        map.incrementCount(bin, true);

        // check values correct
        BaseQualityBinCounter.Count count = map.toTileBaseQualityMap().get(bin);
        assertEquals(2, count.getTotalCount());
        assertEquals(1, count.getErrorCount());

        firstOfPair = true;
        bin = new TileBaseQualityBin("A90924:7:HKKAHXX:1", 1, 1101, firstOfPair, 25, (byte)11);
        map.incrementCount(bin, false);
        map.incrementCount(bin, false);
        map.incrementCount(bin, false);
        map.incrementCount(bin, true);
        map.incrementCount(bin, true);

        // check values correct
        count = map.toTileBaseQualityMap().get(bin);
        assertEquals(3, count.getTotalCount());
        assertEquals(2, count.getErrorCount());
    }
}
