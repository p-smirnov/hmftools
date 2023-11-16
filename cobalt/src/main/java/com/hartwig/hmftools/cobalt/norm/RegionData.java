package com.hartwig.hmftools.cobalt.norm;

import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;

public class RegionData
{
    public final int Position;

    // reference data
    private double mMappability;

    // calculate values
    private double mRelativeEnrichment;

    private final List<SampleRegionData> mSampleRegionData;

    public RegionData(int position)
    {
        Position = position;
        mMappability = 0;
        mSampleRegionData = Lists.newArrayList();
    }

    public double mappability() { return mMappability; }

    public void setGcProfile(double mappability)
    {
        mMappability = mappability;
    }

    public double relativeEnrichment() { return mRelativeEnrichment; }
    public void setRelativeEnrichment(double relativeEnrichment) { mRelativeEnrichment = relativeEnrichment; }

    public void addSampleRegionData(final SampleRegionData sampleRegionData)
    {
        mSampleRegionData.add(sampleRegionData);
    }

    public int sampleCount() { return mSampleRegionData.size(); }
    public SampleRegionData getSampleData(final int sampleIndex) { return mSampleRegionData.get(sampleIndex); }
    public List<SampleRegionData> getSamples() { return mSampleRegionData; }

    public String toString()
    {
        return format("%d: relEnrichment(%.4f)", Position, mRelativeEnrichment);
    }
}
