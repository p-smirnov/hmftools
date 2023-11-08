package com.hartwig.hmftools.errorprofile.recal;

import com.hartwig.hmftools.errorprofile.BaseQualCalcs;
import com.hartwig.hmftools.errorprofile.TileAdjustmentKey;

public class TileBaseQualOutlier
{
    public TileAdjustmentKey tileAdjustmentKey;
    public int readPosition;
    public int errorCount;
    public int totalCount;

    public double empiricalBaseQuality()
    {
        return BaseQualCalcs.empiricalBaseQuality(errorCount, totalCount);
    }

    public TileBaseQualOutlier(final TileAdjustmentKey tileAdjustmentKey, final int readPosition, final int errorCount,
            final int totalCount)
    {
        this.tileAdjustmentKey = tileAdjustmentKey;
        this.readPosition = readPosition;
        this.errorCount = errorCount;
        this.totalCount = totalCount;
    }
}
