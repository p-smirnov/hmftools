package com.hartwig.hmftools.errorprofile.microsatellite;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import org.jetbrains.annotations.NotNull;

import it.unimi.dsi.fastutil.ints.Int2IntArrayMap;

// store the data for each MS table
public class MicrosatelliteStatsTable
{
    public class Row
    {
        public final int refNumUnits;
        public int totalReadCount = 0;
        public Int2IntArrayMap repeatDiffCounts = new Int2IntArrayMap();

        public Row(final int refNumUnits)
        {
            this.refNumUnits = refNumUnits;
        }

        void addRead(int repeatDiff)
        {
            totalReadCount++;
            repeatDiffCounts.mergeInt(repeatDiff, 1, Integer::sum);
        }

        int getRepeatDiffReadCount(int repeatDiff)
        {
            return repeatDiffCounts.get(repeatDiff);
        }

        String getRepeatUnit() { return repeatUnit; }
    }

    public String repeatUnit;

    // ref num unit to rows
    public Map<Integer, Row> rows = new HashMap<>();

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
                int repeatDiff = numRepeatUnits - refNumUnits;
                rows.computeIfAbsent(refNumUnits, k -> new Row(refNumUnits)).addRead(repeatDiff);
            }
        }
    }
}
