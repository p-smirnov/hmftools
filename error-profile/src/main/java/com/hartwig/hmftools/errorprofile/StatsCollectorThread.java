package com.hartwig.hmftools.errorprofile;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import htsjdk.samtools.SAMRecord;

/*
// we use a thread to collect the stats we need, to keep things fast
// we combine them later
public class StatsCollectorThread implements AsyncBamReader.IReadHandler
{
    BaseQualityBinCounter mBaseQualityBinCounter;

    @Override
    public void processRead(SAMRecord read, ChrBaseRegion baseRegion)
    {
        // mReadClassifier.classifyRead(read);
        // mReadBaseSupportAnalyser.processRead(read, baseRegion);
    }

    @Override
    public void regionComplete(ChrBaseRegion baseRegion)
    {
        // mReadBaseSupportAnalyser.regionComplete(baseRegion);
        // mReadClassifier.classifyRead(read);
    }

    public void processReadProfile(ReadProfile readProfile)
    {
        // write the read profiles
        //mReadProfileFileWriter.write(readProfile);
        mBaseQualityBinCounter.onReadProfile(readProfile);
    }

    public void processReadBaseSupport(ReadBaseSupport readBaseSupport)
    {
        // write the read bases supports
        //mReadBaseSupportFileWriter.write(readBaseSupport);

        mBaseQualityBinCounter.onReadBaseSupport(readBaseSupport);
    }
}
*/