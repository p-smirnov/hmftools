package com.hartwig.hmftools.bee.train;

import java.util.ArrayList;
import java.util.Collection;

import htsjdk.samtools.SAMRecord;

public class UmiReadGroup
{
    public final String umi;
    public final String chromosome;
    public final int alignmentStart;
    public final boolean firstOfPair;

    SAMRecord consensusRead = null;
    Collection<SAMRecord> duplicateReads = new ArrayList<>();

    public UmiReadGroup(final String umi, String chromosome, int alignmentStart, boolean firstOfPair)
    {
        this.umi = umi;
        this.chromosome = chromosome;
        this.alignmentStart = alignmentStart;
        this.firstOfPair = firstOfPair;
    }
}
