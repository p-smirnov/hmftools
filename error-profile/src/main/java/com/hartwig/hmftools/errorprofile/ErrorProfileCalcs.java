package com.hartwig.hmftools.errorprofile;

public class ErrorProfileCalcs
{
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
