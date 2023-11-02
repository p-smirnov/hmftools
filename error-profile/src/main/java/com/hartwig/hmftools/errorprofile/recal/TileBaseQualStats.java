package com.hartwig.hmftools.errorprofile.recal;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.errorprofile.TileAdjustmentUnit;

class TileBaseQualStats
{
    public static class Count
    {
        int errorCount = 0;
        int totalCount = 0;

        double empiricalBaseQuality()
        {
            double p = ((double)errorCount) / totalCount;
            return -10 * Math.log10(p);
        }
    }

    public final TileAdjustmentUnit tileAdjustmentUnit;
    List<Count> positionCounts = new ArrayList<>();

    public TileBaseQualAdjustment adjustmentFunction = null;

    TileBaseQualStats(final TileAdjustmentUnit tileAdjustmentUnit)
    {
        this.tileAdjustmentUnit = tileAdjustmentUnit;
    }

    void addToCount(int position, int errorCount, int totalCount)
    {
        while(positionCounts.size() <= position)
        {
            positionCounts.add(new Count());
        }

        Count c = positionCounts.get(position);
        c.errorCount += errorCount;
        c.totalCount += totalCount;
    }
}
