package com.hartwig.hmftools.bamtools.compare;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.samtools.BamSlicer;

import htsjdk.samtools.SAMRecord;

public interface ISAMRecordProvider
{
    void processSlice(final BamSlicer bamSlicer, final ChrBaseRegion region, final Consumer<SAMRecord> consumer);
}
