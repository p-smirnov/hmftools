package com.hartwig.hmftools.bamtools.compare;

import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.samtools.BamSlicer;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

public class SAMRecordProvider implements ISAMRecordProvider
{
    private final SamReader mSamReader;

    public SAMRecordProvider(final SamReader samReader)
    {
        mSamReader = samReader;
    }

    @Override
    public void processSlice(final BamSlicer bamSlicer, final ChrBaseRegion region, final Consumer<SAMRecord> consumer)
    {
        bamSlicer.slice(mSamReader, Lists.newArrayList(region), consumer);
    }
}
