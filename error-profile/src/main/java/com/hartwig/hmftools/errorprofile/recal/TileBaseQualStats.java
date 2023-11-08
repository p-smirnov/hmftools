package com.hartwig.hmftools.errorprofile.recal;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.errorprofile.BaseQualCalcs;
import com.hartwig.hmftools.errorprofile.TileAdjustmentKey;

class TileBaseQualStats
{
    public static class Count
    {
        int errorCount = 0;
        int totalCount = 0;
    }

    public final TileAdjustmentKey tileAdjustmentKey;
    List<Count> positionCounts = new ArrayList<>();

    public TileBaseQualAdjustment adjustmentFunction = null;

    public List<TileBaseQualOutlier> outliers = null;

    public double empiricalBaseQuality = Double.NaN;

    TileBaseQualStats(final TileAdjustmentKey tileAdjustmentKey)
    {
        this.tileAdjustmentKey = tileAdjustmentKey;
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
