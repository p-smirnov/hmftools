package com.hartwig.hmftools.bee.train;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.aligner.AlignmentOperator;
import com.hartwig.hmftools.common.aligner.GlobalSequenceAligner;
import com.hartwig.hmftools.common.utils.Strings;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.util.PeekIterator;

// using cigar to align reads
// this works in following steps
// 1. using the intersecting M regions to map out the parts which both reads match the ref genome.
// 2. expand the M regions where only 1 read matches ref genome.
// 3. finally deal with parts where neither read matches ref genome.
// It could make sense to use local aligner to align the regions outside intersecting M.
public class ReadCigarAligner
{
    private static final Logger sLogger = LogManager.getLogger(ReadCigarAligner.class);

    // walk the cigar by the number of bases
    public static class CigarBaseWalker
    {
        private final Cigar mCigar;
        private int mCurrentElemIndex = 0;
        private int mConsumedBases = 0;

        public CigarBaseWalker(Cigar cigar)
        {
            mCigar = cigar;
        }

        public CigarElement current()
        {
            return mCigar.getCigarElement(mCurrentElemIndex);
        }

        public CigarOperator currentOperator()
        {
            return current().getOperator();
        }

        public int currentRemaining()
        {
            Preconditions.checkArgument(current().getLength() > mConsumedBases);
            return current().getLength() - mConsumedBases;
        }

        public boolean reachedEnd()
        {
            return mCurrentElemIndex == mCigar.numCigarElements();
        }

        public void consumeBases(int numBases)
        {
            while (numBases > 0)
            {
                Preconditions.checkArgument(mCurrentElemIndex <= mCigar.numCigarElements());

                CigarElement cigarElement = current();

                if (mConsumedBases + numBases < cigarElement.getLength())
                {
                    mConsumedBases += numBases;
                }
                else
                {
                    numBases -= cigarElement.getLength() - mConsumedBases;
                    mCurrentElemIndex++;
                    mConsumedBases = 0;
                }
            }
        }
    }

    // indices are 0 based, easier for Java code
    public static class ReadRegion
    {
        public final int subjectReadStart;
        public final int subjectReadEnd;
        public final int refReadStart;
        public final int refReadEnd;

        // could be -1
        public final int referenceGenomeStart;

        ReadRegion(int subjectReadStart, int subjectReadEnd,
                   int refReadStart, int refReadEnd,
                   int referenceGenomeStart)
        {
            this.subjectReadStart = subjectReadStart;
            this.subjectReadEnd = subjectReadEnd;
            this.refReadStart = refReadStart;
            this.refReadEnd = refReadEnd;
            this.referenceGenomeStart = referenceGenomeStart;
        }

        int getSubjectReadLength() { return subjectReadEnd - subjectReadStart; }
        int getRefReadLength() { return refReadEnd - refReadStart; }
    }

    // Go through the align blocks of the two reads, and find the intersecting match regions to the ref genomes
    static List<ReadRegion> findAlignedRegions(Cigar subjectCigar, int subjectAlignStart,
            Cigar refCigar, int refAlignStart)
    {
        List<AlignmentBlock> subjectAlignBs = SAMUtils.getAlignmentBlocks(subjectCigar, subjectAlignStart, "read");
        List<AlignmentBlock> refAlignBs = SAMUtils.getAlignmentBlocks(refCigar, refAlignStart, "read");

        PeekIterator<AlignmentBlock> subjectItr = new PeekIterator<>(subjectAlignBs.listIterator());
        PeekIterator<AlignmentBlock> refItr = new PeekIterator<>(refAlignBs.listIterator());

        List<ReadRegion> result = new ArrayList<>();

        while (subjectItr.hasNext() && refItr.hasNext())
        {
            AlignmentBlock subjectBlock = subjectItr.peek();
            AlignmentBlock refBlock = refItr.peek();

            int subjectBlockEnd = subjectBlock.getReferenceStart() + subjectBlock.getLength();
            int refBlockEnd = refBlock.getReferenceStart() + refBlock.getLength();

            // test if they overlap
            if (subjectBlock.getReferenceStart() < refBlockEnd &&
                    refBlock.getReferenceStart() < subjectBlockEnd)
            {
                //
                int overlapStart = Math.max(subjectBlock.getReferenceStart(), refBlock.getReferenceStart());
                int overlapLength = Math.min(subjectBlockEnd, refBlockEnd) - overlapStart;
                int subjectReadStart = subjectBlock.getReadStart() - 1 + overlapStart - subjectBlock.getReferenceStart();
                int refReadStart = refBlock.getReadStart() - 1 + overlapStart - refBlock.getReferenceStart();

                result.add(new ReadRegion(subjectReadStart, subjectReadStart + overlapLength,
                        refReadStart, refReadStart + overlapLength,
                        overlapStart));
            }

            // choose which one to advance
            if (subjectBlockEnd > refBlockEnd)
            {
                refItr.next();
            }
            else if (subjectBlockEnd < refBlockEnd)
            {
                subjectItr.next();
            }
            else
            {
                refItr.next();
                subjectItr.next();
            }
        }

        return result;
    }

    public static List<ReadRegion> buildReadRegions(Cigar subjectCigar, int subjectAlignStart,
            Cigar refCigar, int refAlignStart)
    {
        List<ReadRegion> alignedRegions = findAlignedRegions(subjectCigar, subjectAlignStart, refCigar, refAlignStart);

        List<ReadRegion> result = new ArrayList<>();

        // find regions that are not covered
        if (alignedRegions.isEmpty())
        {
            // whole read is not aligned.
            result.add(new ReadRegion(0, subjectCigar.getReadLength(),
                    0, refCigar.getReadLength(),
                    -1));
            return result;
        }

        int subjectNonAlignedStart = 0;
        int refNonAlignedStart = 0;

        for (ReadRegion alignedRegion : alignedRegions)
        {
            if (alignedRegion.subjectReadStart > subjectNonAlignedStart ||
                alignedRegion.refReadStart > refNonAlignedStart)
            {
                result.add(new ReadRegion(subjectNonAlignedStart, alignedRegion.subjectReadStart,
                        refNonAlignedStart, alignedRegion.refReadStart,
                        -1));
            }

            result.add(alignedRegion);

            subjectNonAlignedStart = alignedRegion.subjectReadEnd;
            refNonAlignedStart = alignedRegion.refReadEnd;
        }

        if (subjectNonAlignedStart < subjectCigar.getReadLength() ||
            refNonAlignedStart < refCigar.getReadLength())
        {
            result.add(new ReadRegion(subjectNonAlignedStart, subjectCigar.getReadLength(),
                    refNonAlignedStart, refCigar.getReadLength(),
                    -1));
        }

        return result;
    }

    public static String longDebugString(SAMRecord read)
    {
        return String.format("%s cigar(%s) reversed?(%s) ummapped?(%s) proper paired?(%s) mate align(%s:%d) isize(%d)",
                read, read.getCigar(), read.getReadNegativeStrandFlag(),
                read.getReadUnmappedFlag(), read.getProperPairFlag(),
                read.getMateReferenceName(),
                read.getMateAlignmentStart(),
                read.getInferredInsertSize());
    }

    public static List<AlignmentOperator> alignReads(SAMRecord subjectRead, SAMRecord refRead)
    {
        sLogger.trace("subject read({}), ref read({})", subjectRead, refRead);
        List<AlignmentOperator> alignOps = alignReads(subjectRead.getReadString(), subjectRead.getCigar(), subjectRead.getAlignmentStart(),
                refRead.getReadString(), refRead.getCigar(), refRead.getAlignmentStart());

        // count the number of mismatches, if it is more than 30, log it as error
        long numMismatches = alignOps.stream().filter(op -> op == AlignmentOperator.MISMATCH).count();

        if (numMismatches > 40)
        {
            sLogger.error("too many mismatches({}) subject({}) ref({})", numMismatches,
                    longDebugString(subjectRead), longDebugString(refRead));
            var alignment = new GlobalSequenceAligner.Alignment(alignOps, 0);
            GlobalSequenceAligner.Alignment.logAlignment(sLogger, Level.ERROR, subjectRead.getReadString(),
                    refRead.getReadString(), alignment);
        }
        else if (sLogger.isTraceEnabled())
        {
            var alignment = new GlobalSequenceAligner.Alignment(alignOps, 0);
            GlobalSequenceAligner.Alignment.logAlignment(sLogger, Level.TRACE, subjectRead.getReadString(),
                    refRead.getReadString(), alignment);
        }

        return alignOps;
    }

    //
    public static List<AlignmentOperator> alignReads(String subjectReadBases, Cigar subjectCigar, int subjectAlignStart,
            String refReadBases, Cigar refCigar, int refAlignStart)
    {
        var alignOps = new ArrayList<AlignmentOperator>();

        sLogger.trace("subject cigar({}), ref cigar({})", subjectCigar, refCigar);

        // if the number of bases different are less than 5, we assume there is no indel.
        if (subjectReadBases.length() == refReadBases.length() &&
            Strings.countCharDiff(subjectReadBases, refReadBases) < 5)
        {
            // they both align each other, assume both to be 151M
            // walk through both and work out the alignment operators
            appendAlignOps(alignOps, subjectReadBases, 0, refReadBases, 0, subjectReadBases.length());
            return alignOps;
        }

        List<ReadRegion> readRegions = buildReadRegions(subjectCigar, subjectAlignStart, refCigar, refAlignStart);

        for (ReadRegion region : readRegions)
        {
            sLogger.trace("region: sub({}:{}) ref({}:{}) match?({})",
                    region.subjectReadStart, region.subjectReadEnd,
                    region.refReadStart, region.refReadEnd,
                    region.referenceGenomeStart >= 1);

            if (region.referenceGenomeStart >= 1)
            {
                Preconditions.checkArgument(region.getRefReadLength() == region.getSubjectReadLength());
                // this part is aligned
                appendAlignOps(alignOps, subjectReadBases, region.subjectReadStart,
                        refReadBases, region.refReadStart, region.getSubjectReadLength());
                continue;
            }
            // if this part is not aligned, we do some simple checks
            if (region.getRefReadLength() == region.getSubjectReadLength())
            {
                // same length, they might already match each other. But we cannot be certain, test it out
                // if we compare every base and there are not too many mismatch, there is no indel
                if (Strings.countCharDiff(subjectReadBases, region.subjectReadStart,
                        refReadBases, region.refReadStart, region.getSubjectReadLength()) < 3)
                {
                    // the bases all match
                    appendAlignOps(alignOps, subjectReadBases, region.subjectReadStart,
                            refReadBases, region.refReadStart, region.getSubjectReadLength());
                    continue;
                }
            }

            // if one of the read region is empty then we set it as either insert or delete
            if (region.getRefReadLength() == 0)
            {
                for (int i = 0; i < region.getSubjectReadLength(); ++i)
                {
                    alignOps.add(AlignmentOperator.INSERTION);
                }
                continue;
            }
            if (region.getSubjectReadLength() == 0)
            {
                for (int i = 0; i < region.getRefReadLength(); ++i)
                {
                    alignOps.add(AlignmentOperator.DELETION);
                }
                continue;
            }

            // otherwise use aligner to build it up
            String subject = subjectReadBases.substring(region.subjectReadStart, region.subjectReadEnd);
            String ref = refReadBases.substring(region.refReadStart, region.refReadEnd);

            // we heavily favour mismatch over gap. We are dealing with sequencing errors, not
            // biological insert, therefore gap extension should have same cost as gap opening
            // note this could be different and need configuration if we use Utima
            var aligner = new GlobalSequenceAligner(1, -1, -6, -6);
            //aligner.setLogWorkMatrix(true);
            GlobalSequenceAligner.Alignment alignment = aligner.alignSequence(subject, ref);
            alignOps.addAll(alignment.getOperators());
        }

        return alignOps;
    }

    static void appendAlignOps(Collection<AlignmentOperator> ops, String leftSeq, int leftSeqOffset,
            String rightSeq, int rightSeqOffset, int length)
    {
        for (int i = 0; i < length; ++i)
        {
            ops.add(leftSeq.charAt(leftSeqOffset + i) == rightSeq.charAt(rightSeqOffset + i) ?
                    AlignmentOperator.MATCH:
                    AlignmentOperator.MISMATCH);
        }
    }
}
