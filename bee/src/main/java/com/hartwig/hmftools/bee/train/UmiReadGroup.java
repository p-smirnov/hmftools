package com.hartwig.hmftools.bee.train;

import java.util.ArrayList;
import java.util.Collection;

import htsjdk.samtools.SAMRecord;

public class UmiReadGroup
{
    public final String umi;
    public final String chromosome;
    public final int read5PrimeReferencePos;
    public final boolean negativeStrand;

    SAMRecord consensusRead = null;
    Collection<SAMRecord> duplicateReads = new ArrayList<>();

    public UmiReadGroup(final String umi, String chromosome, int read5PrimeReferencePos, boolean negativeStrand)
    {
        this.umi = umi;
        this.chromosome = chromosome;
        this.read5PrimeReferencePos = read5PrimeReferencePos;
        this.negativeStrand = negativeStrand;
    }
}
