package com.hartwig.hmftools.errorprofile.recal;

// essentially regression function for base qual adjustment
public class TileBaseQualAdjustment
{
    public double leftSlope;
    public double leftIntersect;
    public int midStart;
    public double midSlope;
    public double midIntersect;
    public int rightStart;
    public double rightSlope;
    public double rightIntersect;

    public TileBaseQualAdjustment(final double leftSlope, final double leftIntersect, final int midStart, final double midSlope,
            final double midIntersect, final int rightStart, final double rightSlope, final double rightIntersect)
    {
        this.leftSlope = leftSlope;
        this.leftIntersect = leftIntersect;
        this.midStart = midStart;
        this.midSlope = midSlope;
        this.midIntersect = midIntersect;
        this.rightStart = rightStart;
        this.rightSlope = rightSlope;
        this.rightIntersect = rightIntersect;
    }

    double calcBaseQualAdjustment(int position)
    {
        double slope;
        double intersect;
        if(position >= rightStart)
        {
            slope = leftSlope;
            intersect = leftIntersect;
        }
        else if(position >= midStart)
        {
            slope = midSlope;
            intersect = midIntersect;
        }
        else
        {
            slope = rightSlope;
            intersect = rightIntersect;
        }

        return intersect + slope * position;
    }
}
