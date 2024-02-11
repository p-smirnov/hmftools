package com.hartwig.hmftools.errorprofile.repeat;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.commons.lang3.Validate;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.SAMRecord;
import it.unimi.dsi.fastutil.ints.Int2IntArrayMap;
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.IntBinaryOperator;

// works with one repeat
// find the number of reads of each number of repeat across the section
public class RepeatAnalyser
{
    public static final Logger sLogger = LogManager.getLogger(RepeatAnalyser.class);

    final RefGenomeMicrosatellite refGenomeMicrosatellite;

    private final List<ReadRepeatMatch> mReadRepeatMatches = new ArrayList<>();

    public List<ReadRepeatMatch> getReadRepeatMatches() { return mReadRepeatMatches; }

    public List<ReadRepeatMatch> getPassingReadRepeatMatches()
    {
        return mReadRepeatMatches.stream().filter(o -> !o.shouldDropRead).collect(Collectors.toList());
    }

    public int numReadRejected()
    {
        return (int)mReadRepeatMatches.stream().filter(o -> o.shouldDropRead).count();
    }

    public RepeatAnalyser(final RefGenomeMicrosatellite refGenomeMicrosatellite)
    {
        this.refGenomeMicrosatellite = refGenomeMicrosatellite;
    }

    public void addReadToStats(final SAMRecord read)
    {
        if(read.getReadUnmappedFlag() || read.getDuplicateReadFlag())
            return;

        mReadRepeatMatches.add(ReadRepeatMatch.from(refGenomeMicrosatellite, read));
    }

    public int getCountWithRepeatUnits(int numRepeatUnits)
    {
        return (int)getPassingReadRepeatMatches().stream().filter(o -> o.numRepeatUnits() == numRepeatUnits).count();
    }

    public boolean isRealVariant(final double altCountFractionInit,
            final double altCountFractionCutoffStep,
            final double rejectedReadFractionCutoff)
    {
        // debug: use old logic
        IntBinaryOperator sumFunc = Integer::sum;

        /*
        # filter out cases of real variants, this is defined by more than 30% of the reads are in other counts
        filter_df = df[[ x for x in df.columns if x.startswith("count")]]
        filter_df = filter_df.divide(filter_df.sum(axis=1), axis="index")
        filter_df = filter_df.drop(columns="count+0") > 0.3
        df["realVariant"] = filter_df.any(axis=1)
         */
        Int2IntArrayMap repeatReadCounts = new Int2IntArrayMap();

        for(ReadRepeatMatch readRepeatMatch : getPassingReadRepeatMatches())
        {
            if(readRepeatMatch.numRepeatUnits() != refGenomeMicrosatellite.numRepeat)
            {
                repeatReadCounts.mergeInt(readRepeatMatch.numRepeatUnits(), 1, sumFunc);
            }
        }

        double cutoff = getPassingReadRepeatMatches().size() * 0.3;

        return repeatReadCounts.values().stream().anyMatch(i -> i > cutoff);


        // new logic below
        /*

        Validate.isTrue(altCountFractionCutoffStep <= 0.0);

        double fractionRejected = 1.0 - getPassingReadRepeatMatches().size() / (double)getReadRepeatMatches().size();

        if(Doubles.greaterOrEqual(fractionRejected, rejectedReadFractionCutoff))
        {
            return true;
        }

        Int2IntArrayMap repeatReadCounts = new Int2IntArrayMap();

        for(ReadRepeatMatch readRepeatMatch : getPassingReadRepeatMatches())
        {
            int repeatDiff = refGenomeMicrosatellite.numRepeat - readRepeatMatch.numRepeatUnits();

            if(repeatDiff != 0)
            {
                repeatReadCounts.mergeInt(repeatDiff, 1, Integer::sum);
            }
        }

        for(Int2IntMap.Entry entry : repeatReadCounts.int2IntEntrySet())
        {
            int repeatDiff = entry.getIntKey();
            int readCount = entry.getIntValue();

            double countCutoff = (altCountFractionInit + (Math.abs(repeatDiff) - 1) * altCountFractionCutoffStep) * getPassingReadRepeatMatches().size();
            if(Doubles.greaterThan(readCount, countCutoff))
            {
                return true;
            }
        }

        return false;
         */
    }

    /*
    public void completeRegion(
            final AtomicLong readTagGenerator,
            final Consumer<ReadProfile> readProfileConsumer,
            final Consumer<ReadBaseSupport> readBaseSupportConsumer)
    {
        sLogger.info("region: {} read count: {}", , mReads.size());

        List<ReadProfile> readProfiles = new ArrayList<>();
        List<ReadBaseSupport> readBaseSupports = new ArrayList<>();

        // process all reads in this region
        for(SAMRecord read : mReads)
        {
            // we skip reads that are not entirely in this region
            if(read.getAlignmentStart() < mGenomeRegion.start() || read.getAlignmentStart() + read.getReadLength() > mGenomeRegion.end() + 1)
            {
                continue;
            }

            // also skip reads that matches ref genome

            long readTag = readTagGenerator.incrementAndGet();
            ReadBaseSupport readBaseSupport = calcReadBaseSupport(read, readTag);
            readBaseSupportConsumer.accept(readBaseSupport);
            readProfileConsumer.accept(ReadProfiler.profileRead(read, readTag));
        }

        sLogger.info("finished processing region: {}", mGenomeRegion);
    }

    // check if this read matches ref genome
    // TODO: write a ref genome cache to make it more generic
    public boolean isAllMatch(final SAMRecord record)
    {
        if (record.getCigar().numCigarElements() > 1 || record.getCigar().getFirstCigarElement().getOperator() != CigarOperator.M)
            return false;

        // this record has 151M, check it against the ref
        // we need to
        byte[] readBases = record.getReadBases();
        int startOffset = record.getAlignmentStart() - mGenomeRegion.start();
        int readOffset = Math.max(-startOffset, 0);
        int length = record.getReadLength();

        boolean allMatch = true;

        for(int i = 0; i < length; ++i)
        {
            byte refBase = mRefGenomeCache.getBase(record.getAlignmentStart() + i);
            if(refBase == N)
                continue;
            byte readBase = readBases[readOffset + i];
            if(readBase != N && readBase != refBase)
            {
                allMatch = false;
                break;
            }
        }

        if(allMatch)
        {
            // sLogger.debug("read all match ref, cigar={}", record.getCigarString());
        }
        return allMatch;
    }*/
}
