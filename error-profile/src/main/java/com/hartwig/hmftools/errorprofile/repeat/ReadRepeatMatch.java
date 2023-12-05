package com.hartwig.hmftools.errorprofile.repeat;

import static com.hartwig.hmftools.errorprofile.repeat.RepeatProfileConstant.MAX_DISTANCE;

import com.hartwig.hmftools.errorprofile.utils.CigarHandler;

import org.apache.commons.lang3.Validate;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

// we need to look through the cigar to decide whether this read has matches
// we can either
// 1. check the M elements and see if they add up
// 2. check the D elements to see if any was deleted
// 3. we must check the I elements in case the polymer was lengthened.
class ReadRepeatMatch implements CigarHandler
{
    public static final Logger sLogger = LogManager.getLogger(ReadRepeatMatch.class);

    final RefGenomeHomopolymer refGenomeHomopolymer;
    boolean shouldDropRead = false;
    int numAligned = 0;
    int numInserted = 0;
    int numDeleted = 0;

    private ReadRepeatMatch(final RefGenomeHomopolymer refGenomeHomopolymer, final SAMRecord record)
    {
        this.refGenomeHomopolymer = refGenomeHomopolymer;

        // this read needs to wholly contain the homopolymer to be counted
        if(record.getAlignmentStart() > refGenomeHomopolymer.referenceStart() ||
        record.getAlignmentEnd() < refGenomeHomopolymer.referenceEnd())
        {
            shouldDropRead = true;
        }
        else
        {
            traverseCigar(record);
        }

        // do some validation and logging
        if(!shouldDropRead)
        {
            int readRepeatLength = numAligned + numInserted;
            if(readRepeatLength != refGenomeHomopolymer.numRepeat - numDeleted + numInserted)
            {
                sLogger.error("read({}) {}, incorrect read repeat length({}) numAligned({}) numInserted({}) numDeleted({}) " +
                        "homopolymer({})",
                        record, record.getCigarString(), readRepeatLength, numAligned, numInserted, numDeleted, refGenomeHomopolymer.genomeRegion);
            }
        }
    }

    public int readRepeatLength()
    {
        Validate.isTrue(!shouldDropRead);
        // calculate the repeat length
        int readRepeatLength = numAligned + numInserted;

        // check that deleted makes sense
        // Validate.isTrue(readRepeatLength == refGenomeHomopolymer.numRepeat - numDeleted + numInserted);

        return readRepeatLength;
    }

    public static ReadRepeatMatch from(final RefGenomeHomopolymer refGenomeHomopolymer, final SAMRecord record)
    {
        return new ReadRepeatMatch(refGenomeHomopolymer, record);
    }

    public void handleAlignment(final SAMRecord record, final CigarElement e, final int startReadIndex, final int startRefPos)
    {
        // check if this alignment spans the repeat
        int repeatAlignStart = Math.max(startRefPos, refGenomeHomopolymer.referenceStart());
        int repeatAlignEnd = Math.min(startRefPos + e.getLength() - 1, refGenomeHomopolymer.referenceEnd());

        if(repeatAlignStart <= repeatAlignEnd)
        {
            numAligned += repeatAlignEnd - repeatAlignStart + 1;
        }
    }

    public void handleInsert(final SAMRecord record, final CigarElement e, final int readIndex, final int refPos)
    {
        if(refPos >= refGenomeHomopolymer.referenceStart() && refPos <= refGenomeHomopolymer.referenceEnd() + 1)
        {
            sLogger.trace("read: {} inserted {} bases", record, e.getLength());

            // the insert is in the same base, we should check if the base is the repeat base
            for(int i = readIndex; i < readIndex + e.getLength(); ++i)
            {
                if(record.getReadBases()[i] != refGenomeHomopolymer.base)
                {
                    shouldDropRead = true;
                    // should drop this base
                    sLogger.trace("read: {} inserted bases: {} vs ref: {} mismatch, dropping read",
                            record, (char)record.getReadBases()[i], (char)refGenomeHomopolymer.base);
                }
            }

            numInserted += e.getLength();
        }
        else if(refPos + MAX_DISTANCE >= refGenomeHomopolymer.referenceStart() && refPos - MAX_DISTANCE <= refGenomeHomopolymer.referenceEnd() + 1)
        {
            // there is an insert very close to the homopolymer, drop this read
            shouldDropRead = true;
        }
    }

    public void handleDelete(final SAMRecord record, final CigarElement e, final int readIndex, final int startRefPos)
    {
        int endRefPos = startRefPos + e.getLength() - 1;
        if(startRefPos >= refGenomeHomopolymer.referenceStart() && endRefPos <= refGenomeHomopolymer.referenceEnd())
        {
            // the whole delete is inside the repeat, this is nice simple case
            numDeleted += e.getLength();
        }
        else if((startRefPos < refGenomeHomopolymer.referenceStart() && endRefPos >= refGenomeHomopolymer.referenceStart() - 1) ||
                (startRefPos <= refGenomeHomopolymer.referenceEnd() + 1 && endRefPos > refGenomeHomopolymer.referenceEnd()))
        {
            // if the delete cross over the start or the end, drop the read
            // also drop if the delete is just before the start of the polymer
            shouldDropRead = true;
        }
    }

    public void handleLeftSoftClip(final SAMRecord record, final CigarElement element)
    {
        // drop this read completely if the soft clip is near the repeat
        if(record.getAlignmentStart() >= refGenomeHomopolymer.referenceStart())
        {
            shouldDropRead = true;
        }
    }

    public void handleRightSoftClip(final SAMRecord record, final CigarElement element, int startReadIndex, int startRefPosition)
    {
        // drop this read completely if the soft clip is inside the repeat
        if(startRefPosition <= refGenomeHomopolymer.referenceEnd())
        {
            shouldDropRead = true;
        }
    }
}
