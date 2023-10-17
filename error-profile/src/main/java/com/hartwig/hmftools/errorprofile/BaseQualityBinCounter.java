package com.hartwig.hmftools.errorprofile;

import static htsjdk.samtools.util.SequenceUtil.A;
import static htsjdk.samtools.util.SequenceUtil.N;
import static htsjdk.samtools.util.SequenceUtil.T;
import static htsjdk.samtools.util.SequenceUtil.C;
import static htsjdk.samtools.util.SequenceUtil.G;

import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

import org.apache.commons.lang3.Validate;
import org.apache.commons.lang3.tuple.ImmutableTriple;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.lang3.tuple.Triple;

// class to look at all the bam error stats
// comprehensive
//
public class BaseQualityBinCounter
{
    public static class Count
    {
        AtomicInteger nonVariantCount = new AtomicInteger(0);
        AtomicInteger realVariantCount = new AtomicInteger(0);
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

    Map<BaseQualityBin, Count> mBaseQualityCountMap = new ConcurrentHashMap<>();
    Map<TileBaseQualityBin, Count> mTileBaseQualityCountMap = new ConcurrentHashMap<>();

    Map<BaseQualityBin, Count> getBaseQualityCountMap() { return mBaseQualityCountMap; }
    Map<TileBaseQualityBin, Count> getTileBaseQualityCountMap() { return mTileBaseQualityCountMap; }

    synchronized void onReadProfile(ReadProfile readProfile)
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
        String flowcell = flowcellLaneTile.getLeft();
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

            if(posSupport.ref[0] == N || posSupport.alt == N)
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
                    posSupport.ref[0],
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
                    posSupport.ref[0],
                    posSupport.alt,
                    posSupport.baseQuality);

            Count c = mBaseQualityCountMap.computeIfAbsent(bqBin, (k) -> new Count());
            Count tileCount = mTileBaseQualityCountMap.computeIfAbsent(tileBaseQualityBin, (k) -> new Count());

            if(ErrorProfileCalcs.likelyRealVariant(posSupport))
            {
                c.realVariantCount.incrementAndGet();
                tileCount.realVariantCount.incrementAndGet();
            }
            else
            {
                c.nonVariantCount.incrementAndGet();
                tileCount.nonVariantCount.incrementAndGet();
            }
        }
    }

    void getOrCreateBaseQualityBin()
    {

    }
}
