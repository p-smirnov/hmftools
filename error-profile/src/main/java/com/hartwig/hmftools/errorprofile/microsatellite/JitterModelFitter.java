package com.hartwig.hmftools.errorprofile.microsatellite;

import java.util.HashMap;
import java.util.Map;

import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.commons.lang3.Validate;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;
import org.jetbrains.annotations.Nullable;

import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.doubles.DoubleList;

public class JitterModelFitter
{
    public static class ScaleSkew
    {
        public final double scale;
        public final double skew;
        public final double loss;

        private ScaleSkew(final double scale, final double skew, final double loss)
        {
            this.scale = scale;
            this.skew = skew;
            this.loss = loss;
        }

        public static ScaleSkew of(final double scale, final double skew, final double loss)
        {
            return new ScaleSkew(scale, skew, loss);
        }
    }

    private final MicrosatelliteStatsTable mStatsTable;
    private final Map<Integer, ScaleSkew> msRepeatFittedScaleSkew = new HashMap<>();

    private JitterModelParams mJitterModelParams = null;

    public JitterModelParams getJitterModelParams()
    {
        return mJitterModelParams;
    }

    public JitterModelFitter(MicrosatelliteStatsTable statsTable)
    {
        this.mStatsTable = statsTable;
    }

    public void performFit()
    {
        msRepeatFittedScaleSkew.clear();
        mJitterModelParams = null;

        for(int numRepeats = 4; numRepeats <= 15; ++numRepeats)
        {
            ScaleSkew bestScaleSkew = gridSearch(numRepeats);
            if(bestScaleSkew != null)
            {
                msRepeatFittedScaleSkew.put(numRepeats, bestScaleSkew);
            }
        }

        double readCountSum7_15 = 0.0;
        for(int numRepeats = 7; numRepeats <= 15; ++numRepeats)
        {
            readCountSum7_15 += mStatsTable.getReadCount(numRepeats);
        }

        double scaleFitGradient = Double.NaN;
        double scaleFitIntercept = Double.NaN;
        double microsatelliteSkew = 1.0;

        if(readCountSum7_15 > 10_000)
        {
            // next we do a weighted least square fitting linear regression through the repeat lengths
            // 7-15 to find the scale
            // https://stackoverflow.com/questions/5684282/weighted-linear-regression-in-java

            DoubleList xUnweighted = new DoubleArrayList();
            DoubleList yUnweighted = new DoubleArrayList();
            DoubleList weights = new DoubleArrayList();
            DoubleList skews = new DoubleArrayList();

            for(int numRepeats = 7; numRepeats <= 15; ++numRepeats)
            {
                double x = numRepeats;
                ScaleSkew scaleSkew = msRepeatFittedScaleSkew.get(numRepeats);
                if(scaleSkew == null)
                {
                    continue;
                }
                double y = scaleSkew.scale;
                double w = Math.min(mStatsTable.getReadCount(numRepeats) / readCountSum7_15, 0.2);

                xUnweighted.add(x);
                yUnweighted.add(y);
                weights.add(w);
                skews.add(scaleSkew.skew);
            }

            if(xUnweighted.size() >= 2)
            {
                double[] y = new double[yUnweighted.size()];
                double[][] x = new double[xUnweighted.size()][2];

                for (int i = 0; i < y.length; i++) {
                    y[i] = Math.sqrt(weights.getDouble(i)) * yUnweighted.getDouble(i);
                    x[i][0] = Math.sqrt(weights.getDouble(i)) * xUnweighted.getDouble(i);
                    x[i][1] = Math.sqrt(weights.getDouble(i));
                }

                OLSMultipleLinearRegression regression = new OLSMultipleLinearRegression();
                regression.setNoIntercept(true);
                regression.newSampleData(y, x);

                double[] regressionParameters = regression.estimateRegressionParameters();

                // gradient must be positive
                if(Doubles.greaterThan(regressionParameters[0], 0.0))
                {
                    scaleFitGradient = regressionParameters[0];
                    scaleFitIntercept = regressionParameters[1];

                    // also calculate the microsatellite skew as the weighted average of the skews
                    microsatelliteSkew = 0.0;
                    for(int i = 0; i < skews.size(); ++i)
                    {
                        microsatelliteSkew += weights.getDouble(i) * skews.getDouble(i);
                    }

                    microsatelliteSkew /= weights.doubleStream().sum();
                }
            }
        }

        // fall back
        if(Double.isNaN(scaleFitGradient))
        {
            ScaleSkew scaleSkewAt6 = msRepeatFittedScaleSkew.get(6);
            if(scaleSkewAt6 != null)
            {
                scaleFitGradient = 0.04;
                scaleFitIntercept = scaleSkewAt6.scale - 0.04 * 6;
            }
        }

        double optimalScaleRepeat4 = Double.NaN;
        double optimalScaleRepeat5 = Double.NaN;
        double optimalScaleRepeat6 = Double.NaN;

        ScaleSkew scaleSkewAt4 = msRepeatFittedScaleSkew.get(4);
        if(scaleSkewAt4 != null)
        {
            optimalScaleRepeat4 = scaleSkewAt4.scale;
        }

        ScaleSkew scaleSkewAt5 = msRepeatFittedScaleSkew.get(5);
        if(scaleSkewAt5 != null)
        {
            optimalScaleRepeat5 = scaleSkewAt5.scale;
        }

        ScaleSkew scaleSkewAt6 = msRepeatFittedScaleSkew.get(6);
        if(scaleSkewAt6 != null)
        {
            optimalScaleRepeat6 = scaleSkewAt6.scale;
        }

        // assign the params
        mJitterModelParams = new JitterModelParams(
                mStatsTable.repeatUnit,
                optimalScaleRepeat4, optimalScaleRepeat5, optimalScaleRepeat6,
                scaleFitGradient, scaleFitIntercept, microsatelliteSkew);
    }

    @Nullable
    public ScaleSkew gridSearch(int numRepeats)
    {
        double minLoss = Double.MAX_VALUE;
        double bestScale = Double.NaN;
        double bestSkew = Double.NaN;

        MicrosatelliteStatsTable.Row row = mStatsTable.getRow(numRepeats);

        if(row == null)
        {
            return null;
        }

        // we need to find the fitted scale of the previous number of repeated
        double lengthMinusOneScale = 0;
        if (msRepeatFittedScaleSkew.containsKey(numRepeats - 1))
        {
            lengthMinusOneScale = msRepeatFittedScaleSkew.get(numRepeats - 1).scale;
        }

        JitterModelLoss lossCalc = new JitterModelLoss(row, numRepeats, lengthMinusOneScale);

        for(double scale = 0.05; scale <= 5.001; scale = searchIncrement(scale))
        {
            for(double skew = 0.05; skew <= 5.001; skew = searchIncrement(skew))
            {
                double loss = lossCalc.totalLoss(scale, skew);

                if(loss < minLoss)
                {
                    minLoss = loss;
                    bestScale = scale;
                    bestSkew = skew;
                }
            }
        }

        return ScaleSkew.of(bestScale, bestSkew, minLoss);
    }

    private static double searchIncrement(double last)
    {
        if(Doubles.lessOrEqual(last, 0.25))
            return last + 0.01;
        if(Doubles.lessThan(last, 0.5))
            return last + 0.02;
        if(Doubles.lessOrEqual(last, 2.0))
            return last + 0.05;
        return last + 0.2;
    }
}
