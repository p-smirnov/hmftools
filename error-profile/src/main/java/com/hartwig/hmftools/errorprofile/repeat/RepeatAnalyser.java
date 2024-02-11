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

    // have threshold for ALT site differ depending on INDEL length (e.g. 30% for INDEL=1, 25% for INDEL=2, ..., 10% for INDEL=5)
    public boolean isRealVariant(final double altCountFractionInit,
            final double altCountFractionCutoffStep,
            final double rejectedReadFractionCutoff)
    {
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
    }

    public boolean isRealVariantOld(final double altCountFractionInit,
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
    }
}
