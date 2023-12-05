package com.hartwig.hmftools.errorprofile;

import static htsjdk.samtools.util.SequenceUtil.A;
import static htsjdk.samtools.util.SequenceUtil.C;
import static htsjdk.samtools.util.SequenceUtil.G;
import static htsjdk.samtools.util.SequenceUtil.N;
import static htsjdk.samtools.util.SequenceUtil.T;

import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicIntegerFieldUpdater;

import com.google.common.collect.Interner;
import com.google.common.collect.Interners;

import org.apache.commons.lang3.Validate;
import org.apache.commons.lang3.tuple.ImmutableTriple;

// class to look at all the bam error stats
// comprehensive
//
public class BaseQualityBinCounter
{
    // We need to make sure the count can be incremented atomically.
    // use AtomicIntegerFieldUpdater instead of AtomicInteger, it is much more efficient
    // https://stackoverflow.com/questions/58863128/unexpected-varhandle-performance-4x-slower-than-alternatives
    public static class Count
    {
        private volatile int errorCount = 0;
        private volatile int totalCount = 0;

        public int getErrorCount()
        {
            return errorCount;
        }

        public int getTotalCount()
        {
            return totalCount;
        }

        void incrementErrorCount()
        {
            ERROR_COUNT_UPDATER.getAndAdd(this, 1);
        }

        void incrementTotalCount()
        {
            TOTAL_COUNT_UPDATER.getAndAdd(this, 1);
        }

        void setErrorCount(int c)
        {
            ERROR_COUNT_UPDATER.set(this, c);
        }

        void setTotalCount(int c)
        {
            TOTAL_COUNT_UPDATER.set(this, c);
        }

        private static final AtomicIntegerFieldUpdater<Count> ERROR_COUNT_UPDATER =
                AtomicIntegerFieldUpdater.newUpdater(Count.class, "errorCount");
        private static final AtomicIntegerFieldUpdater<Count> TOTAL_COUNT_UPDATER =
                AtomicIntegerFieldUpdater.newUpdater(Count.class, "totalCount");
    }

    public static class BaseCount
    {
        int aCount = 0;
        int tCount = 0;
        int cCount = 0;
        int gCount = 0;

        void increment(byte base)
        {
            switch(base)
            {
                case A: ++aCount; break;
                case T: ++tCount; break;
                case C: ++cCount; break;
                case G: ++gCount; break;
                default: throw new IllegalArgumentException("illegal base: " + base);
            }
        }

        int getCount(byte base)
        {
            switch(base)
            {
                case A: return aCount;
                case T: return tCount;
                case C: return cCount;
                case G: return gCount;
                default: throw new IllegalArgumentException("illegal base: " + base);
            }
        }
    }

    Map<BaseQualityBin, Count> mBaseQualityCountMap = new ConcurrentHashMap<>(100_000, 0.5f);
    //Map<TileBaseQualityBin, Count> mTileBaseQualityCountMap = new ConcurrentHashMap<>(200_000, 0.5f);
    TileBaseQualityBinMap mFastTileBaseQualityCountMap = new TileBaseQualityBinMap();

    Interner<String> mStringInterner = Interners.newStrongInterner();

    public Map<BaseQualityBin, Count> getBaseQualityCountMap() { return mBaseQualityCountMap; }
    public Map<TileBaseQualityBin, Count> getTileBaseQualityCountMap() { return mFastTileBaseQualityCountMap.toTileBaseQualityMap(); }
    //public Map<TileBaseQualityBin, Count> getTileBaseQualityCountMap() { return mTileBaseQualityCountMap; }

    void onReadProfile(ReadProfile readProfile)
    {
        // add the read profile to the counts

    }

    /*
    final boolean firstOfPair, final int readPosition,
                final char ref, final char alt,
                final char trinucleotideContext0, final char trinucleotideContext1, final char trinucleotideContext2,
                final int rawBaseQuality
     */
    void onReadBaseSupport(ReadBaseSupport readBaseSupport)
    {
        // get the read id
        ImmutableTriple<String, Integer, Integer> flowcellLaneTile = ErrorProfileUtils.parseFlowcellLaneTile(readBaseSupport.read.getReadName());
        Validate.notNull(flowcellLaneTile);
        String flowcell = mStringInterner.intern(flowcellLaneTile.getLeft());
        int lane = flowcellLaneTile.getMiddle();
        int tile = flowcellLaneTile.getRight();

        // add the readBaseSupport to the counts
        for(ReadBaseSupport.PositionSupport posSupport : readBaseSupport.positionSupports)
        {
            if(!posSupport.isAlignment())
            {
                continue;
            }

            Validate.isTrue(posSupport.ref.length == 1);

            byte ref = posSupport.ref[0];

            if(ref == N || posSupport.alt == N)
            {
                continue;
            }

            if(posSupport.trinucleotideContext[0] == N ||
                posSupport.trinucleotideContext[1] == N ||
                posSupport.trinucleotideContext[2] == N)
            {
                continue;
            }

            BaseQualityBin bqBin = new BaseQualityBin(readBaseSupport.read.getFirstOfPairFlag(),
                    posSupport.readPosition5To3,
                    ref,
                    posSupport.alt,
                    posSupport.trinucleotideContext[0],
                    posSupport.trinucleotideContext[1],
                    posSupport.trinucleotideContext[2],
                    posSupport.baseQuality);

            TileBaseQualityBin tileBaseQualityBin = new TileBaseQualityBin(
                    flowcell,
                    lane,
                    tile,
                    readBaseSupport.read.getFirstOfPairFlag(),
                    posSupport.readPosition5To3,
                    //ref,
                    //posSupport.alt,
                    posSupport.baseQuality);

            // TODO: remove this
            //if (posSupport.readPosition5To3 != 0 || posSupport.baseQuality != 37)
            //    continue;

            Count c = mBaseQualityCountMap.computeIfAbsent(bqBin, (k) -> new Count());
            //Count tileCount = mTileBaseQualityCountMap.computeIfAbsent(tileBaseQualityBin, (k) -> new Count());

            mFastTileBaseQualityCountMap.incrementCount(tileBaseQualityBin, false);
            c.incrementTotalCount();
            //tileCount.incrementTotalCount();

            if(ref != posSupport.alt && !BaseQualCalcs.likelyRealVariant(posSupport))
            {
                c.incrementErrorCount();
                //tileCount.incrementErrorCount();
                mFastTileBaseQualityCountMap.incrementCount(tileBaseQualityBin, true);
            }
        }
    }
}
