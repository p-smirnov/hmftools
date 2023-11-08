package com.hartwig.hmftools.errorprofile;

import static htsjdk.samtools.util.SequenceUtil.A;
import static htsjdk.samtools.util.SequenceUtil.T;
import static htsjdk.samtools.util.SequenceUtil.C;
import static htsjdk.samtools.util.SequenceUtil.G;

// we get from read id
// @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>
// TODO: 1. use byte array for flowcell
//       2. improve the hash function
//       3. cache the hash code
public class TileBaseQualityBin
{
    public final String flowcell;

    public final short lane;
    public final short tile;

    public final boolean firstOfPair;

    public final short readPosition;

    //public final byte ref;
    //public final byte alt;

    public final byte rawBaseQuality;

    private final int hashCodeCache;

    public TileBaseQualityBin(final String flowcell, final int lane, final int tile, final boolean firstOfPair, final int readPosition,
            // final byte ref, final byte alt,
            final byte rawBaseQuality)
    {
        this.flowcell = flowcell;
        this.lane = (short)lane;
        this.tile = (short)tile;
        this.firstOfPair = firstOfPair;
        this.readPosition = (short)readPosition;
        //this.ref = ref;
        //this.alt = alt;
        this.rawBaseQuality = rawBaseQuality;
        hashCodeCache = calcHashCode();
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
        /*if(ref != that.ref)
        {
            return false;
        }
        if(alt != that.alt)
        {
            return false;
        }*/
        if(rawBaseQuality != that.rawBaseQuality)
        {
            return false;
        }
        return flowcell.equals(that.flowcell); // use intern string therefore this is ok
    }

    @Override
    public int hashCode()
    {
        return hashCodeCache;
    }

    private int calcHashCode()
    {
        final int PRIME = 109;
        int result = flowcell.hashCode();
        result = PRIME * result + (int) lane;
        result = PRIME * result + (int) tile;
        result = PRIME * result + (firstOfPair ? 1 : 0);
        result = PRIME * result + (int) readPosition;
        result = PRIME * result + (int) rawBaseQuality;
        return result;

        /*
        // usually 4 lanes
        int result = lane;

        // tile number first digit says whether it is top half or bottom half
        int tileFirstDigit = tile / 1000;
        result = 2 * result + tileFirstDigit;
        result = 2 * result + (firstOfPair ? 1 : 0);
        result = 151 * result + readPosition;
        //result = 4 * result + baseToHashCode(ref);
        //result = 4 * result + baseToHashCode(alt);
        result = 991 * result + tile % 1000;

        // previous parts produce a linear hash code, we want to spread it out
        result = 61 * result + rawBaseQuality;
        result = 8191 * result + flowcell.hashCode();
        return result;
         */
    }

    private static int baseToHashCode(byte base)
    {
        switch(base)
        {
            case A: return 0;
            case T: return 1;
            case C: return 2;
            case G: return 3;
            default: throw new IllegalArgumentException();
        }
    }
}
