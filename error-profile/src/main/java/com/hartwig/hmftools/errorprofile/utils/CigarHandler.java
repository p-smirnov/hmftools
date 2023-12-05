package com.hartwig.hmftools.errorprofile.utils;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public interface CigarHandler
{
    default void handleLeftSoftClip(final SAMRecord record, final CigarElement element) {}

    default void handleRightSoftClip(final SAMRecord record, final CigarElement element, int readIndex, int refPosition) {}

    default void handleAlignment(final SAMRecord record, final CigarElement element, int readIndex, int refPosition) {}

    default void handleInsert(final SAMRecord record, final CigarElement element, int readIndex, int refPosition) {}

    default void handleDelete(final SAMRecord record, final CigarElement element, int readIndex, int refPosition) {}

    default void handleSkippedReference(final SAMRecord record, final CigarElement element, int readIndex, int refPosition) {}

    default void traverseCigar(final SAMRecord record)
    {
        final Cigar cigar = record.getCigar();

        int readIndex = 0;
        int refBase = record.getAlignmentStart();

        for(int i = 0; i < cigar.numCigarElements(); i++)
        {
            final CigarElement e = cigar.getCigarElement(i);
            switch(e.getOperator())
            {
                case H:
                    break; // ignore hard clips - no need to skip either bases or positions
                case P:
                    break; // ignore pads
                case S:
                    if(i == 0)
                    {
                        handleLeftSoftClip(record, e);
                    }
                    else if(i == cigar.numCigarElements() - 1)
                    {
                        handleRightSoftClip(record, e, readIndex, refBase);
                    }
                    readIndex += e.getLength();
                    break; // soft clip read bases
                case N:
                    //handleSkippedReference(record, e, readIndex, refBase);
                    refBase += e.getLength();
                    break;  // reference skip
                case D:
                    handleDelete(record, e, readIndex, refBase);
                    refBase += e.getLength();
                    break;
                case I:
                    handleInsert(record, e, readIndex, refBase);
                    readIndex += e.getLength();
                    break;
                case M:
                case EQ:
                case X:
                    handleAlignment(record, e, readIndex, refBase);
                    readIndex += e.getLength();
                    refBase += e.getLength();
                    break;
                default:
                    throw new IllegalStateException("Case statement didn't deal with op: " + e.getOperator() + "in CIGAR: " + cigar);
            }
        }
    }
}
