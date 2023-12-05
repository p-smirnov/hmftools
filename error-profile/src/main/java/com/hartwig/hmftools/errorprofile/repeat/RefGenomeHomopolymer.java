package com.hartwig.hmftools.errorprofile.repeat;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class RefGenomeHomopolymer
{
    public final ChrBaseRegion genomeRegion;
    public final byte base;
    public final int numRepeat;

    public RefGenomeHomopolymer(final ChrBaseRegion genomeRegion, final byte base, final int numRepeat)
    {
        this.genomeRegion = genomeRegion;
        this.base = base;
        this.numRepeat = numRepeat;
    }

    public RefGenomeHomopolymer(final String chromosome, int start, int end, final byte base, final int numRepeat)
    {
        this(new ChrBaseRegion(chromosome, start, end), base, numRepeat);
    }

    public String chromosome()
    {
        return genomeRegion.chromosome();
    }

    public int referenceStart()
    {
        return genomeRegion.start();
    }

    public int referenceEnd()
    {
        return genomeRegion.end();
    }
}
