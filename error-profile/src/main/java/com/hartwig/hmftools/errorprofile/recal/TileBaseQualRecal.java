package com.hartwig.hmftools.errorprofile.recal;

import static com.hartwig.hmftools.errorprofile.ErrorProfileConstants.BASEQ_OUTLIER_CUTOFF;
import static com.hartwig.hmftools.errorprofile.ErrorProfileConstants.BASEQ_OUTLIER_MEDIAN_COUNT;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.errorprofile.BaseQualityBinCounter;
import com.hartwig.hmftools.errorprofile.TileAdjustmentUnit;
import com.hartwig.hmftools.errorprofile.TileBaseQualityBin;

import org.apache.commons.math3.stat.regression.SimpleRegression;

import it.unimi.dsi.fastutil.doubles.DoubleDoublePair;

public class TileBaseQualRecal
{
    BaseQualityBinCounter mBaseQualityBinCounter;

    // first find all the tiles
    Map<TileAdjustmentUnit, TileBaseQualStats> computeTileAdjustmentStats()
    {
        Map<TileAdjustmentUnit, TileBaseQualStats> tileAdjustmentStats = new HashMap<>();

        for(Map.Entry<TileBaseQualityBin, BaseQualityBinCounter.Count> entry : mBaseQualityBinCounter.getTileBaseQualityCountMap().entrySet())
        {
            TileBaseQualityBin tileBaseQualityBin = entry.getKey();
            BaseQualityBinCounter.Count count = entry.getValue();

            TileAdjustmentUnit tileAdjustmentUnit = new TileAdjustmentUnit(tileBaseQualityBin.flowcell, tileBaseQualityBin.lane, tileBaseQualityBin.tile,
                    tileBaseQualityBin.rawBaseQuality,
                    tileBaseQualityBin.firstOfPair);
            TileBaseQualStats stats = tileAdjustmentStats.computeIfAbsent(tileAdjustmentUnit, k -> new TileBaseQualStats(tileAdjustmentUnit));

            if(tileBaseQualityBin.ref != tileBaseQualityBin.alt)
            {
                stats.addToCount(tileBaseQualityBin.readPosition, count.getNonVariantCount(), count.getNonVariantCount() + count.getRealVariantCount());
            }
            else
            {
                stats.addToCount(tileBaseQualityBin.readPosition, 0, count.getNonVariantCount() + count.getRealVariantCount());
            }
        }

        return tileAdjustmentStats;
    }

    // for each tile, we find the functions, ignoring ref / alt
    // i.e. the empirical error rate would be the sum (non variant | ref != alt) / sum
    void computeAdjustment(TileBaseQualStats stats)
    {
        // find the lowest point in the stats
        double rawBaseQuality = stats.tileAdjustmentUnit.rawBaseQuality;
        double highestBaseQual = 0;
        double lowestBaseQual = Double.MAX_VALUE;
        int posWithHighestBaseQual = -1;
        int posWithLowestBaseQual = -1;

        List<Double> empiricalBaseQualities = stats.positionCounts.stream().map(TileBaseQualStats.Count::empiricalBaseQuality).collect(Collectors.toList());
        List<Integer> outlierPositions = new ArrayList<>();
        SimpleRegression midRegression = new SimpleRegression();

        for(int readPosition = 0; readPosition < stats.positionCounts.size(); ++readPosition)
        {
            // calculate the median of points around it, to find outliers
            ArrayList<Double> medianOf = new ArrayList<>();

            for(int i = Math.max(0, readPosition - BASEQ_OUTLIER_MEDIAN_COUNT);
                    i < Math.min(stats.positionCounts.size(), readPosition + BASEQ_OUTLIER_MEDIAN_COUNT + 1);
                    ++i)
            {
                if(i != readPosition)
                {
                    medianOf.add(empiricalBaseQualities.get(i));
                }
            }

            double median = Doubles.median(medianOf);

            double empiricalBaseQuality = empiricalBaseQualities.get(readPosition);

            // we only exclude for low quality values, as outliers are basically always bad ones with
            // problems
            if(empiricalBaseQuality > median - BASEQ_OUTLIER_CUTOFF)
            {
                outlierPositions.add(readPosition);
                continue;
            }

            // update lowest and highest
            if(empiricalBaseQuality > highestBaseQual)
            {
                posWithHighestBaseQual = readPosition;
                highestBaseQual = empiricalBaseQuality;
            }

            if(empiricalBaseQuality < lowestBaseQual)
            {
                posWithLowestBaseQual = readPosition;
                lowestBaseQual = empiricalBaseQuality;
            }

            midRegression.addData(readPosition, empiricalBaseQuality);
        }

        // use the regression line to find how we should model that start and the end part of the curve
        // usually the first few bases has poor quality
        // what we want to do is to get all the points that are below the regression line
        // and use those to form another regression line at the start
        // Note here we do not remove outliers
        SimpleRegression leftRegression = new SimpleRegression();
        int midStart = 0;
        for(int readPosition = 0; readPosition < stats.positionCounts.size(); ++readPosition)
        {
            double empiricalBaseQuality = empiricalBaseQualities.get(readPosition);
            if(midRegression.predict(readPosition) > empiricalBaseQuality)
            {
                leftRegression.addData(readPosition, empiricalBaseQuality);
            }
            else
            {
                midStart = readPosition;
                break;
            }
        }

        // do same regression for points at the end. Note here we do not remove outliers
        SimpleRegression rightRegression = new SimpleRegression();
        int rightStart = 0;
        for(int readPosition = stats.positionCounts.size() - 1; readPosition >= 0; --readPosition)
        {
            double empiricalBaseQuality = empiricalBaseQualities.get(readPosition);
            if(midRegression.predict(readPosition) > empiricalBaseQuality)
            {
                rightRegression.addData(readPosition, empiricalBaseQuality);
                rightStart = readPosition;
            }
            else
            {
                break;
            }
        }

        // now put everything together, we have left, mid and right regression curves
        // write it out to CSV for plotting
        TileBaseQualAdjustment tileBaseQualAdjustment = new TileBaseQualAdjustment(leftRegression.getSlope(), leftRegression.getIntercept(),
                midStart, midRegression.getSlope(), midRegression.getIntercept(),
                rightStart, rightRegression.getSlope(), rightRegression.getIntercept());

        stats.adjustmentFunction = tileBaseQualAdjustment;
    }

    // perform first linear regression to get the point points
    // then we fix up the ends by doing another two linear regressions
    SimpleRegression performRegression(List<DoubleDoublePair> points)
    {
        SimpleRegression regression = new SimpleRegression();
        points.forEach(p -> regression.addData(p.leftDouble(), p.rightDouble()));
        //regression.regress()
        return regression;
    }
}
