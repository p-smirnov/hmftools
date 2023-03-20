package com.hartwig.hmftools.bee;

public class BeeInput
{
    public String readString;
    public String baseQualityString;

    // we use N if it is not mapped
    public String referenceGenomeBases;

    public double[] bqrBaseQualities;

    public boolean isRead1;
    public boolean isMapped;
    public boolean isMateMapped;
    public double gcContent;
    public double mappability;
    public double mapQuality;

    public String cigar;

    // target vector to predict
    public boolean[] baseError;
}
