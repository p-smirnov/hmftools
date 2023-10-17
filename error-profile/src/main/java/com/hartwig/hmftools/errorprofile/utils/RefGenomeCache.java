package com.hartwig.hmftools.errorprofile.utils;

import java.util.Arrays;

import org.apache.commons.lang3.Validate;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

// cache chunks of the reference genome
// a class to avoid querying reference genome too frequently
// relies on the fact that we usually request the same bases
public class RefGenomeCache implements IRefGenomeCache
{
    // by default cache 100k
    private static final int DEFAULT_CHUNK_SIZE = 100_000;
    private static final int LEFT_FLANK_BASES = 2_000;

    private IndexedFastaSequenceFile mRefGenome;

    private String mContig;
    private int mChunkSize;

    private long mChunkStartPosition = 0;
    private long mContigLength;
    private byte[] mRefChunk = null;

    public RefGenomeCache(IndexedFastaSequenceFile refGenome, String contig, int chunkSize)
    {
        mRefGenome = refGenome;
        mChunkSize = chunkSize;
        mContig = contig;
        mContigLength = refGenome.getSequence(contig).length();
    }

    public RefGenomeCache(IndexedFastaSequenceFile refGenome, String contig)
    {
        this(refGenome, contig, DEFAULT_CHUNK_SIZE);
    }

    // position is 1 based
    @Override public byte getBase(long refGenomePosition)
    {
        Validate.isTrue(refGenomePosition > 0);
        if(mRefChunk == null || refGenomePosition < mChunkStartPosition || refGenomePosition >= mChunkStartPosition + mRefChunk.length)
        {
            getNewChunk(refGenomePosition);
        }
        return mRefChunk[(int)(refGenomePosition - mChunkStartPosition)];
    }

    // positions are 1 base and end inclusive
    @Override public byte[] getBases(long refStart, long refEndInclusive)
    {
        Validate.isTrue(refStart > 0);
        if(refStart < mChunkStartPosition || refEndInclusive >= mChunkStartPosition + mRefChunk.length)
        {
            getNewChunk(refStart);
        }
        return Arrays.copyOfRange(mRefChunk, (int) (refStart - mChunkStartPosition), (int) (refEndInclusive - mChunkStartPosition + 1));
    }

    @Override public String getBasesAsString(long refStart, long refEndInclusive)
    {
        return Arrays.toString(getBases(refStart, refEndInclusive));
    }

    private void getNewChunk(long requestedPosition)
    {
        // get a new chunk, we assume that we are processing the reads in sequential manner.
        // so we want to do it such that there is more right side than left side
        mChunkStartPosition = Math.max(requestedPosition - LEFT_FLANK_BASES, 1);
        long end = Math.min(mChunkStartPosition + mChunkSize - 1, mContigLength);
        mRefChunk = mRefGenome.getSubsequenceAt(mContig, mChunkStartPosition, end).getBases();
    }
}
