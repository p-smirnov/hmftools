package com.hartwig.hmftools.errorprofile;

// we get from read id
// @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>
public class TileBaseQualityBin
{
    public String flowcell;

    public int lane;
    public int tile;

    public boolean firstOfPair;

    public int readPosition;

    public byte ref;
    public byte alt;

    public int rawBaseQuality;

    public TileBaseQualityBin(final String flowcell, final int lane, final int tile, final boolean firstOfPair, final int readPosition,
            final byte ref, final byte alt, final int rawBaseQuality)
    {
        this.flowcell = flowcell;
        this.lane = lane;
        this.tile = tile;
        this.firstOfPair = firstOfPair;
        this.readPosition = readPosition;
        this.ref = ref;
        this.alt = alt;
        this.rawBaseQuality = rawBaseQuality;
    }

    @Override
    public boolean equals(final Object o)
    {
        if(this == o)
        {
            return true;
        }
        if(!(o instanceof TileBaseQualityBin))
        {
            return false;
        }

        final TileBaseQualityBin that = (TileBaseQualityBin) o;

        if(lane != that.lane)
        {
            return false;
        }
        if(tile != that.tile)
        {
            return false;
        }
        if(firstOfPair != that.firstOfPair)
        {
            return false;
        }
        if(readPosition != that.readPosition)
        {
            return false;
        }
        if(ref != that.ref)
        {
            return false;
        }
        if(alt != that.alt)
        {
            return false;
        }
        if(rawBaseQuality != that.rawBaseQuality)
        {
            return false;
        }
        return flowcell.equals(that.flowcell);
    }

    @Override
    public int hashCode()
    {
        // deliberately omit flowcell from hashcode, since it is inefficient to calculate
        int result = 1;
        result = 31 * result + lane;
        result = 31 * result + tile;
        result = 31 * result + (firstOfPair ? 1 : 0);
        result = 31 * result + readPosition;
        result = 31 * result + (int) ref;
        result = 31 * result + (int) alt;
        result = 31 * result + rawBaseQuality;
        return result;
    }
}
