package com.hartwig.hmftools.errorprofile;

// we get from read id
// @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>
public class BaseQualityBin
{
    public boolean firstOfPair;

    public int readPosition;

    public byte ref;
    public byte alt;

    // trinucleotide context is always the reference bases
    public byte trinucleotideContext0;
    public byte trinucleotideContext1;
    public byte trinucleotideContext2;

    public int rawBaseQuality;

    public byte[] getTrinucleotideContext()
    {
        return new byte[] { trinucleotideContext0, trinucleotideContext1, trinucleotideContext2 };
    }

    public String getTrinucleotideContextString()
    {
        return new String(getTrinucleotideContext());
    }

    public BaseQualityBin(final boolean firstOfPair, final int readPosition,
            final byte ref, final byte alt,
            final byte trinucleotideContext0, final byte trinucleotideContext1, final byte trinucleotideContext2,
            final int rawBaseQuality)
    {
        this.firstOfPair = firstOfPair;
        this.readPosition = readPosition;
        this.ref = ref;
        this.alt = alt;
        this.trinucleotideContext0 = trinucleotideContext0;
        this.trinucleotideContext1 = trinucleotideContext1;
        this.trinucleotideContext2 = trinucleotideContext2;
        this.rawBaseQuality = rawBaseQuality;
    }

    @Override
    public boolean equals(final Object o)
    {
        if(this == o)
        {
            return true;
        }
        if(!(o instanceof BaseQualityBin))
        {
            return false;
        }

        final BaseQualityBin that = (BaseQualityBin) o;

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
        if(trinucleotideContext0 != that.trinucleotideContext0)
        {
            return false;
        }
        if(trinucleotideContext1 != that.trinucleotideContext1)
        {
            return false;
        }
        if(trinucleotideContext2 != that.trinucleotideContext2)
        {
            return false;
        }
        return rawBaseQuality == that.rawBaseQuality;
    }

    @Override
    public int hashCode()
    {
        final int PRIME = 127;
        int result = (firstOfPair ? 1 : 0);
        result = PRIME * result + readPosition;
        result = PRIME * result + (int) ref;
        result = PRIME * result + (int) alt;
        result = PRIME * result + (int) trinucleotideContext0;
        result = PRIME * result + (int) trinucleotideContext1;
        result = PRIME * result + (int) trinucleotideContext2;
        result = PRIME * result + rawBaseQuality;
        return result;
    }
}
