package com.hartwig.hmftools.errorprofile.utils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

import org.apache.commons.math3.distribution.NormalDistribution;

public class KernelDensitySmoother
{
    private static final int SIGMA_CUTOFF = 5;

    private static class Point
    {
        public double x;
        public double y;

        public Point(final double x, final double y)
        {
            this.x = x;
            this.y = y;
        }
    }

    public KernelDensitySmoother(double sigma)
    {
        mSigma = sigma;
        mNormDist = new NormalDistribution(0.0, mSigma);
    }

    public void addPoint(double x, double y)
    {
        Point point = new Point(x, y);

        // must be sorted
        int index = Collections.binarySearch(mPointList, point, mXComparator);

        // this value already exist
        if (index >= 0)
        {
            throw new IllegalArgumentException("value already exist");
        }

        // insert it
        mPointList.add(-index - 1, point);
    }

    public double getSmoothedY(double x)
    {
        int index = Collections.binarySearch(mPointList, new Point(x, 0.0), mXComparator);
        double kernelSum = 0.0;

        // we go backwards and forwards
        if (index < 0)
        {
            index = -index - 1;
        }

        // go forward
        for (int i = index; i < mPointList.size(); ++i)
        {
            double d = mPointList.get(i).x - x;

            // we cut off at 10 sigma
            if (d >= SIGMA_CUTOFF * mSigma)
            {
                break;
            }

            kernelSum += mNormDist.density(d) * mPointList.get(i).y;
        }

        // go backwards
        for (int i = index - 1; i >= 0; --i)
        {
            double d = x - mPointList.get(i).x;

            // we cut off at 10 sigma
            if (d >= SIGMA_CUTOFF * mSigma)
            {
                break;
            }

            kernelSum += mNormDist.density(d) * mPointList.get(i).y;
        }

        return kernelSum;
    }

    private final double mSigma;

    private final NormalDistribution mNormDist;

    private final ArrayList<Point> mPointList = new ArrayList<>();

    private static final Comparator<Point> mXComparator = Comparator.comparingDouble(p -> p.x);
}
