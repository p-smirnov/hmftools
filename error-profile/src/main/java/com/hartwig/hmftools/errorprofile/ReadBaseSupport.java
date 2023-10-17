package com.hartwig.hmftools.errorprofile;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.genome.gc.GcCalcs;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

// for a read, detail the base support for each position
// Each position can be
// 1.
public class ReadBaseSupport
{
    public static class PositionSupport extends BaseSupport
    {
        public CigarOperator cigarOperator;

        public byte[] ref;
        public byte alt;

        public int readPosition5To3;
        public int refPosition;

        public int baseQuality;

        public byte[] trinucleotideContext;

        /*
        public PositionSupport(
                Type type,
                String ref, char alt, int readPosition,
                int posStrandDepth, int negStrandDepth, int posStrandSupport, int negStrandSupport)
        {
            super(posStrandDepth, negStrandDepth, posStrandSupport, negStrandSupport);
            this.type = type;
            this.ref = ref;
            this.alt = alt;
            this.readPosition = readPosition;
        }*/

        public PositionSupport(CigarOperator cigarOp, byte[] ref, byte alt, int readPosition5To3, int refPosition, byte baseQuality,
                byte[] trinucleotideContext,
                BaseSupport baseSupport)
        {
            super(baseSupport.posStrandDepth, baseSupport.negStrandDepth, baseSupport.posStrandSupport, baseSupport.negStrandSupport);
            this.cigarOperator = cigarOp;
            this.ref = ref;
            this.alt = alt;
            this.readPosition5To3 = readPosition5To3;
            this.refPosition = refPosition;
            this.baseQuality = baseQuality;
            this.trinucleotideContext = trinucleotideContext;
        }

        public boolean isAlignment()
        {
            return cigarOperator.isAlignment();
        }
    }

    public SAMRecord read;

    public long readTag;

    public double gc;

    // base support for the whole read
    public List<PositionSupport> positionSupports = new ArrayList<>();

    public ReadBaseSupport(SAMRecord read, long readTag)
    {
        this.read = read;
        this.readTag = readTag;
        this.gc = GcCalcs.calcGcPercent(read.getReadString());
    }

    public void addPositionSupport(PositionSupport positionSupport)
    {
        positionSupports.add(positionSupport);
    }
}
