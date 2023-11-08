package com.hartwig.hmftools.errorprofile.recal;

import static com.hartwig.hmftools.errorprofile.ErrorProfileConstants.BASEQ_OUTLIER_CUTOFF;
import static com.hartwig.hmftools.errorprofile.ErrorProfileConstants.BASEQ_OUTLIER_MEDIAN_COUNT;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.errorprofile.BaseQualCalcs;
import com.hartwig.hmftools.errorprofile.BaseQualityBinCounter;
import com.hartwig.hmftools.errorprofile.ErrorProfileConstants;
import com.hartwig.hmftools.errorprofile.TileAdjustmentKey;
import com.hartwig.hmftools.errorprofile.TileBaseQualityBin;

import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import it.unimi.dsi.fastutil.doubles.DoubleDoublePair;

public class TileBaseQualRecal
{
    public static final Logger sLogger = LogManager.getLogger(TileBaseQualRecal.class);

    // first find all the tiles
    Map<TileAdjustmentKey, TileBaseQualStats> computeTileAdjustmentStats(Map<TileBaseQualityBin, BaseQualityBinCounter.Count> tileBaseQualityCountMap)
    {
        sLogger.info("start collating tile adjustment stats");

        Map<TileAdjustmentKey, TileBaseQualStats> tileAdjustmentStats = new HashMap<>();

        for(Map.Entry<TileBaseQualityBin, BaseQualityBinCounter.Count> entry : tileBaseQualityCountMap.entrySet())
        {
            TileBaseQualityBin tileBaseQualityBin = entry.getKey();
            BaseQualityBinCounter.Count count = entry.getValue();

            TileAdjustmentKey tileAdjustmentKey = new TileAdjustmentKey(tileBaseQualityBin.flowcell, tileBaseQualityBin.lane, tileBaseQualityBin.tile,
                    tileBaseQualityBin.firstOfPair,
                    tileBaseQualityBin.rawBaseQuality);
            TileBaseQualStats stats = tileAdjustmentStats.computeIfAbsent(tileAdjustmentKey, k -> new TileBaseQualStats(tileAdjustmentKey));

            stats.addToCount(tileBaseQualityBin.readPosition, count.getErrorCount(), count.getTotalCount());
        }

        sLogger.info("start calculating adjustments");

        // calculate the adjustments
        tileAdjustmentStats.values().forEach(this::computeAdjustment);

        return tileAdjustmentStats;
    }

    // empirical base quality but if error count is 0 we assume it to be 0.1
    private static double empiricalBaseQualityAdj(int errorCount, int totalCount)
    {
        double denom = errorCount;

        if(denom == 0)
            denom = ErrorProfileConstants.SPARSE_BIN_ERROR_COUNT;

        return BaseQualCalcs.empiricalBaseQuality(denom, totalCount);
    }

    // for each tile, we find the functions, ignoring ref / alt
    // i.e. the empirical error rate would be the sum (non variant | ref != alt) / sum
    void computeAdjustment(TileBaseQualStats stats)
    {
        // first calculate the sum error rate of this tile
        int tileErrorCount = stats.positionCounts.stream().mapToInt(c -> c.errorCount).sum();
        int tileTotalCount = stats.positionCounts.stream().mapToInt(c -> c.totalCount).sum();
        double refBaseQual = BaseQualCalcs.empiricalBaseQuality(tileErrorCount, tileTotalCount);

        List<Double> positionBaseQualDiffs = stats.positionCounts.stream()
                .map(c -> empiricalBaseQualityAdj(c.errorCount, c.totalCount) - refBaseQual)
                .collect(Collectors.toList());

        SimpleRegression midRegression = new SimpleRegression();

        for(int readPosition = 0; readPosition < stats.positionCounts.size(); ++readPosition)
        {
            double positionBqAdj = positionBaseQualDiffs.get(readPosition);

            if(!Double.isFinite(positionBqAdj))
            {
                // filter out NaNs, 0 values and infinity
                continue;
            }

            // calculate the median of points around it, to find outliers
            ArrayList<Double> medianOf = new ArrayList<>();

            for(int i = Math.max(0, readPosition - BASEQ_OUTLIER_MEDIAN_COUNT);
                    i < Math.min(stats.positionCounts.size(), readPosition + BASEQ_OUTLIER_MEDIAN_COUNT + 1);
                    ++i)
            {
                if(i != readPosition)
                {
                    medianOf.add(positionBaseQualDiffs.get(i));
                }
            }

            double median = Doubles.median(medianOf);

            // we only exclude for low quality values, as outliers are basically always bad ones with
            // problems
            if(positionBqAdj < median - BASEQ_OUTLIER_CUTOFF)
            {
                TileBaseQualStats.Count count = stats.positionCounts.get(readPosition);
                assert(empiricalBaseQualityAdj(count.errorCount, count.totalCount) == positionBqAdj);
                //outlierPositions.add(new TileBaseQualOutlier(stats.tileAdjustmentKey, readPosition, count.errorCount, count.totalCount));
                continue;
            }

            midRegression.addData(readPosition, positionBqAdj);
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
            double positionBqAdj = positionBaseQualDiffs.get(readPosition);

            if(!Double.isFinite(positionBqAdj))
            {
                // filter out NaNs and infinity
                continue;
            }
            if(midRegression.predict(readPosition) > positionBqAdj || leftRegression.getN() < 2)
            {
                leftRegression.addData(readPosition, positionBqAdj);
            }
            else
            {
                midStart = readPosition;
                break;
            }
        }

        if(Double.isNaN(leftRegression.getSlope()))
        {
            sLogger.error("left regression slope is NaN");

            leftRegression = new SimpleRegression();
            midStart = 0;
            for(int readPosition = 0; readPosition < stats.positionCounts.size(); ++readPosition)
            {
                double positionBqAdj = positionBaseQualDiffs.get(readPosition);

                if(midRegression.predict(readPosition) > positionBqAdj || readPosition <= 2)
                {
                    leftRegression.addData(readPosition, positionBqAdj);
                }
                else
                {
                    midStart = readPosition;
                    break;
                }
            }
        }

        // do same regression for points at the end. Note here we do not remove outliers
        SimpleRegression rightRegression = new SimpleRegression();
        int rightStart = 0;
        for(int readPosition = stats.positionCounts.size() - 1; readPosition >= 0; --readPosition)
        {
            double positionBqAdj = positionBaseQualDiffs.get(readPosition);

            if(!Double.isFinite(positionBqAdj))
            {
                // filter out NaNs and infinity
                continue;
            }
            if(midRegression.predict(readPosition) > positionBqAdj || rightRegression.getN() < 2)
            {
                rightRegression.addData(readPosition, positionBqAdj);
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
                rightStart, rightRegression.getSlope(), rightRegression.getIntercept(),
                refBaseQual);

        List<TileBaseQualOutlier> outlierPositions = new ArrayList<>();

        for(int readPosition = 0; readPosition < stats.positionCounts.size(); ++readPosition)
        {
            double positionBqAdj = positionBaseQualDiffs.get(readPosition);

            if(!Double.isFinite(positionBqAdj))
            {
                // filter out NaNs and infinity
                continue;
            }

            // we only exclude for low quality values, as outliers are basically always bad ones with
            // problems
            if(positionBqAdj < tileBaseQualAdjustment.calcBaseQualAdjustment(readPosition) - BASEQ_OUTLIER_CUTOFF)
            {
                TileBaseQualStats.Count count = stats.positionCounts.get(readPosition);
                assert(empiricalBaseQualityAdj(count.errorCount, count.totalCount) == positionBqAdj);
                outlierPositions.add(new TileBaseQualOutlier(stats.tileAdjustmentKey, readPosition, count.errorCount, count.totalCount));
            }
        }

        stats.adjustmentFunction = tileBaseQualAdjustment;
        stats.outliers = outlierPositions;
        stats.empiricalBaseQuality = refBaseQual;
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
