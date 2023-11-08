package com.hartwig.hmftools.errorprofile;

public class BaseQualCalcs
{
    public static double empiricalBaseQuality(double errorCount, double totalCount)
    {
        if(totalCount == 0)
            return Double.NaN;

        if(errorCount == 0)
            return Double.POSITIVE_INFINITY;

        double p = errorCount / totalCount;
        return -10 * Math.log10(p);
    }

    public static boolean likelyRealVariant(ReadBaseSupport.PositionSupport posSupport)
    {
        if(posSupport.ref.length == 1 && posSupport.alt == posSupport.ref[0])
        {
            // no variant
            return false;
        }
        if(posSupport.posStrandSupport == 0 && posSupport.posStrandDepth >= 10)
        {
            return false;
        }
        if(posSupport.negStrandSupport == 0 && posSupport.negStrandDepth >= 10)
        {
            return false;
        }
        if(posSupport.totalSupport() < 0.05 * posSupport.totalDepth())
        {
            return false;
        }

        // real variant
        return true;
    }
}
