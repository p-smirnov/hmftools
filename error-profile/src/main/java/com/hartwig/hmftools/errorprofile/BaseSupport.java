package com.hartwig.hmftools.errorprofile;

// all supports are divided by strand
public class BaseSupport
{
    public int posStrandDepth;
    public int negStrandDepth;
    public int posStrandSupport;
    public int negStrandSupport;

    public int totalDepth() { return posStrandDepth + negStrandDepth; }
    public int totalSupport() { return posStrandSupport + negStrandSupport; }

    public BaseSupport(int posStrandDepth, int negStrandDepth, int posStrandSupport, int negStrandSupport)
    {
        this.posStrandDepth = posStrandDepth;
        this.negStrandDepth = negStrandDepth;
        this.posStrandSupport = posStrandSupport;
        this.negStrandSupport = negStrandSupport;
    }
}
