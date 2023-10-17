package com.hartwig.hmftools.errorprofile;

import java.util.Arrays;

import com.hartwig.hmftools.errorprofile.utils.IRefGenomeCache;

import org.apache.commons.lang3.Validate;

import htsjdk.samtools.util.StringUtil;

// a simplified version of RefGenomeCache that does not require an indexed ref genome fasta
// for testing purpose
public class FakeRefGenomeCache implements IRefGenomeCache
{
    private long mChunkStartPosition = 0;
    private byte[] mRefChunk;

    public FakeRefGenomeCache(long refStartPosition, String refGenomeBases)
    {
        mChunkStartPosition = refStartPosition;
        mRefChunk = StringUtil.stringToBytes(refGenomeBases);
    }

    @Override public byte getBase(long refGenomePosition)
    {
        Validate.isTrue(refGenomePosition > 0);
        return mRefChunk[(int)(refGenomePosition - mChunkStartPosition)];
    }

    // positions are 1 base and end inclusive
    @Override public byte[] getBases(long refStart, long refEndInclusive)
    {
        Validate.isTrue(refStart > 0);
        return Arrays.copyOfRange(mRefChunk, (int) (refStart - mChunkStartPosition), (int) (refEndInclusive - mChunkStartPosition + 1));
    }

    @Override public String getBasesAsString(long refStart, long refEndInclusive)
    {
        return Arrays.toString(getBases(refStart, refEndInclusive));
    }
}
