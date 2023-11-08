package com.hartwig.hmftools.errorprofile;

import org.jetbrains.annotations.NotNull;

public class TileAdjustmentKey
{
    public final @NotNull String flowcell;

    public final short lane;
    public final short tile;

    public final boolean firstOfPair;

    public final byte rawBaseQuality;

    public TileAdjustmentKey(@NotNull final String flowcell, final short lane, final short tile, final boolean firstOfPair, final byte rawBaseQuality)
    {
        this.flowcell = flowcell;
        this.lane = lane;
        this.tile = tile;
        this.firstOfPair = firstOfPair;
        this.rawBaseQuality = rawBaseQuality;
    }

    @Override
    public boolean equals(final Object o)
    {
        if(this == o)
        {
            return true;
        }
        if(!(o instanceof TileAdjustmentKey))
        {
            return false;
        }

        final TileAdjustmentKey that = (TileAdjustmentKey) o;

        if(lane != that.lane)
        {
            return false;
        }
        if(tile != that.tile)
        {
            return false;
        }
        if(rawBaseQuality != that.rawBaseQuality)
        {
            return false;
        }
        if(firstOfPair != that.firstOfPair)
        {
            return false;
        }
        return flowcell.equals(that.flowcell);
    }

    @Override
    public int hashCode()
    {
        int result = flowcell.hashCode();
        result = 31 * result + (int) lane;
        result = 31 * result + (int) tile;
        result = 31 * result + (int) rawBaseQuality;
        result = 31 * result + (firstOfPair ? 1 : 0);
        return result;
    }
}
