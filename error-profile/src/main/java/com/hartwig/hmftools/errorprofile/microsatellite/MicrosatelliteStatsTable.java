package com.hartwig.hmftools.errorprofile.microsatellite;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import it.unimi.dsi.fastutil.ints.Int2IntArrayMap;

// store the data for each MS table
public class MicrosatelliteStatsTable
{
    public class Row
    {
        public final int refNumUnits;
        public int totalReadCount = 0;
        public final Int2IntArrayMap jitterCounts = new Int2IntArrayMap();

        public Row(final int refNumUnits)
        {
            this.refNumUnits = refNumUnits;
        }

        void addRead(int jitter)
        {
            totalReadCount++;
            jitterCounts.mergeInt(jitter, 1, Integer::sum);
        }

        int getJitterReadCount(int jitter)
        {
            return jitterCounts.get(jitter);
        }

        String getRepeatUnit() { return repeatUnit; }
    }

    public final String repeatUnit;

    // ref num unit to rows
    private final List<Row> mRows = new ArrayList<>();

    public MicrosatelliteStatsTable(final String repeatUnit)
    {
        this.repeatUnit = repeatUnit;
    }

    void summarise(@NotNull final Collection<MicrosatelliteSiteAnalyser> microsatelliteSiteAnalysers)
    {
        for(MicrosatelliteSiteAnalyser microsatelliteSiteAnalyser : microsatelliteSiteAnalysers)
        {
            if(microsatelliteSiteAnalyser.isRealVariant(MicrosatelliteAnalyserConstants.ALT_COUNT_FRACTION_INIT, MicrosatelliteAnalyserConstants.ALT_COUNT_FRACTION_STEP,
                    MicrosatelliteAnalyserConstants.MAX_REJECTED_READ_FRACTION))
            {
                continue;
            }

            // add all the counts
            for(MicrosatelliteRead microsatelliteRead : microsatelliteSiteAnalyser.getPassingReadRepeatMatches())
            {
                int refNumUnits = microsatelliteSiteAnalyser.refGenomeMicrosatellite.numRepeat;
                int numRepeatUnits = microsatelliteRead.numRepeatUnits();
                int jitter = numRepeatUnits - refNumUnits;
                getOrCreateRow(refNumUnits).addRead(jitter);
            }
        }
    }

    public List<Row> getRows()
    {
        return mRows;
    }

    public int getReadCount(int numRepeats)
    {
        Row row = getRow(numRepeats);
        return row == null ? 0 : row.totalReadCount;
    }

    @NotNull
    public Row getOrCreateRow(int refNumUnits)
    {
        for(int i = 0; i < mRows.size(); ++i)
        {
            Row row = mRows.get(i);
            if(row.refNumUnits == refNumUnits)
            {
                return row;
            }
            if(row.refNumUnits > refNumUnits)
            {
                // make a row object, insert but keep list sorted
                Row newRow = new Row(refNumUnits);
                mRows.add(i, newRow);
                return newRow;
            }
        }
        // make new object and add to end
        Row row = new Row(refNumUnits);
        mRows.add(row);
        return row;
    }

    @Nullable
    public Row getRow(int refNumUnits)
    {
        for(Row row : mRows)
        {
            if(row.refNumUnits == refNumUnits)
            {
                return row;
            }
        }
        return null;
    }
}
