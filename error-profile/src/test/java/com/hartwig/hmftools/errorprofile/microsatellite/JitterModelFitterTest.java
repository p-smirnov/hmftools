package com.hartwig.hmftools.errorprofile.microsatellite;

import java.util.List;

import org.junit.Assert;
import org.junit.Test;

/*
unit numUnits readCount count-1   count+0   count+1
A/T     4       10000     1000     8000     1000
A/T     5       10000     1000     8000     1000
A/T     6       10000     1000     8000     1000
A/T     7       10000     1000     8000     1000
A/T     8       10000     1000     8000     1000
A/T     9       10000     1000     8000     1000
A/T     10      10000     1000     8000     1000
A/T     11      10000     1000     8000     1000
A/T     12      10000     1000     8000     1000
A/T     13      10000     1000     8000     1000
A/T     14      10000     1000     8000     1000
A/T     15      10000     1000     8000     1000
 */

public class JitterModelFitterTest
{
    private static final double EPS = 1e-5;

    @Test
    public void testJitterModelFitterSimple()
    {
        // create a stats table
        MicrosatelliteStatsTable statsTable = new MicrosatelliteStatsTable("A/T");

        // populate some stats
        for(int refNumRepeat = 4; refNumRepeat <= 15; ++refNumRepeat)
        {
            MicrosatelliteStatsTable.Row row = statsTable.getOrCreateRow(refNumRepeat);
            row.totalReadCount = 10000;
            row.jitterCounts.put(0, 8000);
            row.jitterCounts.put(-1, 1000);
            row.jitterCounts.put(1, 1000);
        }

        // do fitter
        JitterModelFitter fitter = new JitterModelFitter(statsTable);
        fitter.performFit();

        JitterModelParams jitterModelParams = fitter.getJitterModelParams();
        Assert.assertEquals(0.5, jitterModelParams.OptimalScaleRepeat4, EPS);
        Assert.assertEquals(0.5, jitterModelParams.OptimalScaleRepeat5, EPS);
        Assert.assertEquals(0.5, jitterModelParams.OptimalScaleRepeat6, EPS);
        Assert.assertEquals(1.0, jitterModelParams.MicrosatelliteSkew, EPS);
        Assert.assertEquals(0.04, jitterModelParams.ScaleFitGradient, EPS);
    }

    @Test
    public void testJitterModelFitterSkewRegression()
    {
        // create a stats table
        MicrosatelliteStatsTable statsTable = new MicrosatelliteStatsTable("A/T");

        // populate some stats, make sure there is some skew
        for(int refNumRepeat = 4; refNumRepeat <= 6; ++refNumRepeat)
        {
            MicrosatelliteStatsTable.Row row = statsTable.getOrCreateRow(refNumRepeat);
            row.totalReadCount = 10000;
            row.jitterCounts.put(0, 8500);
            row.jitterCounts.put(-1, 500);
            row.jitterCounts.put(1, 1000);
        }

        for(int refNumRepeat = 7; refNumRepeat <= 15; ++refNumRepeat)
        {
            MicrosatelliteStatsTable.Row row = statsTable.getOrCreateRow(refNumRepeat);
            row.totalReadCount = 10000;
            int jitter1Count = refNumRepeat * 100;
            row.jitterCounts.put(0, row.totalReadCount - 3 * jitter1Count);
            row.jitterCounts.put(-1, jitter1Count);
            row.jitterCounts.put(1, 2 * jitter1Count);
        }

        // do fitter
        JitterModelFitter fitter = new JitterModelFitter(statsTable);
        fitter.performFit();

        JitterModelParams jitterModelParams = fitter.getJitterModelParams();
        Assert.assertEquals(0.4, jitterModelParams.OptimalScaleRepeat4, EPS);
        Assert.assertEquals(0.4, jitterModelParams.OptimalScaleRepeat5, EPS);
        Assert.assertEquals(0.4, jitterModelParams.OptimalScaleRepeat6, EPS);
        Assert.assertEquals(0.7444444, jitterModelParams.MicrosatelliteSkew, EPS);
        Assert.assertEquals(0.0808333, jitterModelParams.ScaleFitGradient, EPS);
        Assert.assertEquals(-0.055833, jitterModelParams.ScaleFitIntercept, EPS);

        // MicrosatelliteStatsTableFile.write("/Users/hongwingl/hmftools/error-profile/src/test/java/com/hartwig/hmftools/errorprofile/microsatellite/test.tsv",
           //     List.of(statsTable));
    }

    @Test
    public void testJitterModelFitterSkewRegression1()
    {
        // create a stats table
        MicrosatelliteStatsTable statsTable = new MicrosatelliteStatsTable("A/T");

        // populate some stats, make sure there is some skew
        for(int refNumRepeat = 4; refNumRepeat <= 6; ++refNumRepeat)
        {
            MicrosatelliteStatsTable.Row row = statsTable.getOrCreateRow(refNumRepeat);
            row.totalReadCount = 10000;
            row.jitterCounts.put(0, 7750);
            row.jitterCounts.put(-2, 250);
            row.jitterCounts.put(-1, 500);
            row.jitterCounts.put(1, 1000);
            row.jitterCounts.put(2, 500);
        }

        for(int refNumRepeat = 7; refNumRepeat <= 15; ++refNumRepeat)
        {
            MicrosatelliteStatsTable.Row row = statsTable.getOrCreateRow(refNumRepeat);
            row.totalReadCount = 10000;
            int jitter1Count = refNumRepeat * 100;
            int jitter2Count = refNumRepeat * 50;
            row.jitterCounts.put(0, row.totalReadCount - 2 * jitter1Count - 2 * jitter2Count);
            row.jitterCounts.put(-1, jitter1Count);
            row.jitterCounts.put(1, jitter1Count);
            row.jitterCounts.put(-2, jitter2Count);
            row.jitterCounts.put(2, jitter2Count);

        }

        // do fitter
        JitterModelFitter fitter = new JitterModelFitter(statsTable);
        fitter.performFit();

        JitterModelParams jitterModelParams = fitter.getJitterModelParams();
        Assert.assertEquals(0.7, jitterModelParams.OptimalScaleRepeat4, EPS);
        Assert.assertEquals(0.7, jitterModelParams.OptimalScaleRepeat5, EPS);
        Assert.assertEquals(0.7, jitterModelParams.OptimalScaleRepeat6, EPS);
        Assert.assertEquals(1.0, jitterModelParams.MicrosatelliteSkew, EPS);
        Assert.assertEquals(0.05916667, jitterModelParams.ScaleFitGradient, EPS);
        Assert.assertEquals(0.26583333, jitterModelParams.ScaleFitIntercept, EPS);

        MicrosatelliteStatsTableFile.write("/Users/hongwingl/hmftools/error-profile/src/test/java/com/hartwig/hmftools/errorprofile/microsatellite/test.tsv",
                List.of(statsTable));
    }

    @Test
    public void testJitterModelFitterFallback()
    {
        // this tests the fallback case where there are not enough data to do the linear regression
        // create a stats table
        MicrosatelliteStatsTable statsTable = new MicrosatelliteStatsTable("A/T");

        // populate some stats, make sure there is some skew
        for(int refNumRepeat = 4; refNumRepeat <= 6; ++refNumRepeat)
        {
            MicrosatelliteStatsTable.Row row = statsTable.getOrCreateRow(refNumRepeat);
            row.totalReadCount = 10000;
            row.jitterCounts.put(0, 8000);
            row.jitterCounts.put(-1, 1000);
            row.jitterCounts.put(1, 1000);
        }

        // not enough read count in 7-15
        for(int refNumRepeat = 7; refNumRepeat <= 15; ++refNumRepeat)
        {
            MicrosatelliteStatsTable.Row row = statsTable.getOrCreateRow(refNumRepeat);
            row.totalReadCount = 1000;
            row.jitterCounts.put(0, 700);
            row.jitterCounts.put(-1, 100);
            row.jitterCounts.put(1, 200);
        }

        // do fitter
        JitterModelFitter fitter = new JitterModelFitter(statsTable);
        fitter.performFit();

        JitterModelParams jitterModelParams = fitter.getJitterModelParams();
        Assert.assertEquals(0.5, jitterModelParams.OptimalScaleRepeat4, EPS);
        Assert.assertEquals(0.5, jitterModelParams.OptimalScaleRepeat5, EPS);
        Assert.assertEquals(0.5, jitterModelParams.OptimalScaleRepeat6, EPS);
        Assert.assertEquals(1.0, jitterModelParams.MicrosatelliteSkew, EPS);
        Assert.assertEquals(0.04, jitterModelParams.ScaleFitGradient, EPS);

        // scale intercept = OptimalScaleRepeat6 - 0.04 * 6 = 0.5 - 0.24 = 0.26
        Assert.assertEquals(0.26, jitterModelParams.ScaleFitIntercept, EPS);

        // MicrosatelliteStatsTableFile.write("test.tsv", List.of(statsTable));
    }
}

/*
unit numUnits readCount count-2 count-1 count+0 count+1 count+2
A/T     4       10000     250     500     7750    1000    500
A/T     5       10000     250     500     7750    1000    500
A/T     6       10000     250     500     7750    1000    500
A/T     7       10000     350     700     7900    700     350
A/T     8       10000     400     800     7600    800     400
A/T     9       10000     450     900     7300    900     450
A/T     10      10000     500     1000    7000    1000    500
A/T     11      10000     550     1100    6700    1100    550
A/T     12      10000     600     1200    6400    1200    600
A/T     13      10000     650     1300    6100    1300    650
A/T     14      10000     700     1400    5800    1400    700
A/T     15      10000     750     1500    5500    1500    750
 */