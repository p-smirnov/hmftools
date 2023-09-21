package com.hartwig.hmftools.errorprofile;

import static com.hartwig.hmftools.errorprofile.ErrorProfileUtils.getReadStrand;

import static htsjdk.samtools.util.SequenceUtil.N;

import java.util.Comparator;

import com.hartwig.hmftools.common.samtools.CigarHandler;
import com.hartwig.hmftools.common.samtools.CigarTraversal;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.apache.commons.lang3.Validate;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.SequenceUtil;

// go through all the reads in a certain region
// and for each position we collect statistics based on
// 1. number of each mutation
// 2. insertion and types of each
// 3. deletion and types of each
//
// classifying deletion is a little bit more tricky.
// Also deciding whether it is real or not
// Also need to record the strandedness of the types of variant
//
// How do I handle soft clips?
//
// for each read, output
//
public class GenomeRegionStats
{
    public static final Logger sLogger = LogManager.getLogger(GenomeRegionStats.class);

    private final ChrBaseRegion mGenomeRegion;

    private final IndexedFastaSequenceFile mRefGenome;

    private final String mRefBases;

    private final GenomePositionStats[] mGenomePositionStats;

    private final StatsCollector statsCollector = new StatsCollector();

    public GenomeRegionStats(ChrBaseRegion genomeRegion, IndexedFastaSequenceFile refGenome)
    {
        mGenomeRegion = genomeRegion;
        mRefGenome = refGenome;

        // TODO: make sure this is not too large
        mRefBases = refGenome.getSubsequenceAt(genomeRegion.chromosome(), genomeRegion.start(), genomeRegion.end()).getBaseString();

        // this doesn't seem to work too well
        mGenomePositionStats = new GenomePositionStats[genomeRegion.end() - genomeRegion.start() + 1];
    }

    public void addReadToStats(final SAMRecord read)
    {
        CigarTraversal.traverseCigar(read, statsCollector);

        // if(mReadCounter > 0 && (mReadCounter % 1000) == 0)
           // purgeBaseDataList(record.getAlignmentStart());
    }

    // check if this read matches ref genome
    // TODO: write a ref genome cache to make it more generic
    public boolean isAllMatch(final SAMRecord record)
    {
        if (record.getCigar().numCigarElements() > 1 || record.getCigar().getFirstCigarElement().getOperator() != CigarOperator.M)
            return false;

        // this record has 151M, check it against the ref
        // we need to
        String readString = record.getReadString();
        int startOffset = record.getAlignmentStart() - mGenomeRegion.start();
        int readOffset = Math.max(-startOffset, 0);
        int length = record.getReadLength() - readOffset - Math.max(record.getAlignmentEnd() - mGenomeRegion.end(), 0);
        int refOffset = Math.max(startOffset, 0);

        boolean allMatch = true;

        for(int i = 0; i < length; ++i)
        {
            char refBase = mRefBases.charAt(refOffset + i);
            if(refBase == N)
                continue;
            char readBase = readString.charAt(readOffset + i);
            if(readBase != N && readBase != refBase)
            {
                allMatch = false;
                break;
            }
        }

        if(allMatch)
        {
            // sLogger.debug("read all match ref, cigar={}", record.getCigarString());
        }
        return allMatch;
    }

    public ReadBaseSupport calcReadBaseSupport(SAMRecord read, long readTag)
    {
        ReadBaseSupportBuilder readBaseSupportBuilder = new ReadBaseSupportBuilder(read, readTag);
        // which base look correct?
        CigarTraversal.traverseCigar(read, readBaseSupportBuilder);
        readBaseSupportBuilder.convertToReadStrand();

        return readBaseSupportBuilder.mBaseSupports;
    }

    protected @Nullable GenomePositionStats getPositionStats(int position)
    {
        if(position > mGenomeRegion.end() || position < mGenomeRegion.start())
            return null;

        int posIndex = position - mGenomeRegion.start();
        return mGenomePositionStats[posIndex];
    }

    protected GenomePositionStats getOrCreatePositionStats(int position)
    {
        int posIndex = position - mGenomeRegion.start();
        GenomePositionStats genomePositionStats = mGenomePositionStats[posIndex];

        if(genomePositionStats == null)
        {
            char ref = mRefBases.charAt(posIndex);

            genomePositionStats = new GenomePositionStats();
            genomePositionStats.refBase = ref;
            mGenomePositionStats[posIndex] = genomePositionStats;
        }

        return genomePositionStats;
    }

    class StatsCollector implements CigarHandler
    {
        @Override
        public void handleInsert(final SAMRecord record, final CigarElement e, final int readIndex, final int refPos)
        {
            if(refPos < mGenomeRegion.start() || refPos > mGenomeRegion.end())
                return;

            // find the inserted sequence
            // readIndex is 0 based
            String insertedSeq = record.getReadString().substring(readIndex, readIndex + e.getLength());

            getOrCreatePositionStats(refPos).addInsert(getReadStrand(record), insertedSeq);
        }

        @Override
        public void handleDelete(final SAMRecord record, final CigarElement e, final int readIndex, final int startRefPos)
        {
            for(int i = 0; i < e.getLength(); i++)
            {
                int refPos = startRefPos + i;

                if(refPos > mGenomeRegion.end())
                    return;

                if(refPos < mGenomeRegion.start())
                    continue;

                getOrCreatePositionStats(refPos).addDelete(getReadStrand(record));
            }
        }

        @Override
        public void handleAlignment(final SAMRecord record, final CigarElement cigarElement, boolean beforeIndel, final int startReadIndex, final int startRefPos)
        {
            for(int i = 0; i < cigarElement.getLength(); i++)
            {
                int refPos = startRefPos + i;

                if(refPos > mGenomeRegion.end())
                    return;

                if(refPos < mGenomeRegion.start())
                    continue;

                int readIndex = startReadIndex + i;

                char alt = record.getReadString().charAt(readIndex);

                if(alt == N)
                    continue;

                getOrCreatePositionStats(refPos).addAlignedBase(getReadStrand(record), alt);
            }
        }

        @Override
        public void handleLeftSoftClip(final SAMRecord record, final CigarElement element)
        {
            // we get the extrapolated ref position
            int startRefPos = record.getAlignmentStart() - element.getLength();

            for(int readIndex = 0; readIndex < element.getLength(); readIndex++)
            {
                int refPos = startRefPos + readIndex;

                if(refPos > mGenomeRegion.end())
                    return;

                if(refPos < mGenomeRegion.start())
                    continue;

                char alt = record.getReadString().charAt(readIndex);

                if(alt == N)
                    continue;

                getOrCreatePositionStats(refPos).addSoftClippedBase(getReadStrand(record));
            }
        }

        @Override
        public void handleRightSoftClip(final SAMRecord record, final CigarElement element, int startReadIndex, int startRefPosition)
        {
            for(int i = 0; i < element.getLength(); i++)
            {
                int refPos = startRefPosition + i;

                if(refPos > mGenomeRegion.end())
                    return;

                if(refPos < mGenomeRegion.start())
                    continue;

                int readIndex = startReadIndex + i;

                char alt = record.getReadString().charAt(readIndex);

                if(alt == N)
                    continue;

                getOrCreatePositionStats(refPos).addSoftClippedBase(getReadStrand(record));
            }
        }
    }

    class ReadBaseSupportBuilder implements CigarHandler
    {
        private SAMRecord mRead;
        private ReadBaseSupport mBaseSupports;

        ReadBaseSupport getBaseSupport() { return mBaseSupports; }

        ReadBaseSupportBuilder(SAMRecord read, long readTag)
        {
            mRead = read;
            mBaseSupports = new ReadBaseSupport(read, readTag);
        }

        // if it is read negative strand then fix it
        private void convertToReadStrand()
        {
            if(mRead.getReadNegativeStrandFlag())
            {
                for(ReadBaseSupport.PositionSupport posSupport : mBaseSupports.positionSupports)
                {
                    posSupport.readPosition5To3 = mRead.getReadLength() - posSupport.readPosition5To3 - 1;
                    posSupport.ref = SequenceUtil.reverseComplement(posSupport.ref);
                    posSupport.alt = (char) SequenceUtil.complement((byte) posSupport.alt);
                }

                // resort the position support by their position
                // we rely on this being a stable sort
                mBaseSupports.positionSupports.sort(Comparator.comparingInt(o -> o.readPosition5To3));
            }
        }

        @Override
        public void handleInsert(final SAMRecord record, final CigarElement e, final int startReadIndex, final int refPos)
        {
            if(refPos < mGenomeRegion.start() || refPos > mGenomeRegion.end())
                return;

            // find the inserted sequence
            // readIndex is 0 based
            String insertedSeq = record.getReadString().substring(startReadIndex, startReadIndex + e.getLength());

            GenomePositionStats genomePositionStats = getPositionStats(refPos);
            if(genomePositionStats != null)
            {
                BaseSupport support = genomePositionStats.getInsertSupport(insertedSeq);

                for(int i = 0; i < e.getLength(); i++)
                {
                    int readIndex = startReadIndex + i;
                    char alt = record.getReadString().charAt(readIndex);

                    if(alt == N)
                        continue;

                    mBaseSupports.addPositionSupport(new ReadBaseSupport.PositionSupport(CigarOperator.I,
                            "", alt, readIndex, refPos, record.getBaseQualities()[readIndex],
                            support));
                }
            }
        }

        @Override
        public void handleDelete(final SAMRecord record, final CigarElement e, final int readIndex, final int startRefPos)
        {
            // a delete is represented as a deleted sequence at a read position.
            // We do not actually keep track of how many times a multi base delete has happened at a particular point in genome
            // we proxy by finding the position with the lowest supporting count
            BaseSupport leastSupport = null;
            String refBases = mRefBases.substring(startRefPos - mGenomeRegion.start(), startRefPos - mGenomeRegion.start() + e.getLength());

            for(int i = 0; i < e.getLength(); i++)
            {
                int refPos = startRefPos + i;

                if(refPos > mGenomeRegion.end())
                    return;

                if(refPos < mGenomeRegion.start())
                    continue;

                GenomePositionStats genomePositionStats = getPositionStats(refPos);
                if(genomePositionStats != null)
                {
                    Validate.isTrue(refBases.charAt(i) == genomePositionStats.refBase);
                    BaseSupport baseSupport = genomePositionStats.getDeleteBaseSupport();
                    if (leastSupport == null || baseSupport.totalSupport() < leastSupport.totalSupport())
                    {
                        leastSupport = baseSupport;
                    }
                }
            }

            if (leastSupport != null)
            {
                mBaseSupports.addPositionSupport(new ReadBaseSupport.PositionSupport(CigarOperator.D,
                    refBases, '-', readIndex, startRefPos, record.getBaseQualities()[readIndex],
                        leastSupport));
            }
        }

        @Override
        public void handleAlignment(final SAMRecord record, final CigarElement e, boolean beforeIndel, final int startReadIndex, final int startRefPos)
        {
            for(int i = 0; i < e.getLength(); i++)
            {
                int refPos = startRefPos + i;

                if(refPos > mGenomeRegion.end())
                    return;

                if(refPos < mGenomeRegion.start())
                    continue;

                int readIndex = startReadIndex + i;

                char alt = record.getReadString().charAt(readIndex);

                if(alt == N)
                    continue;

                GenomePositionStats genomePositionStats = getPositionStats(refPos);
                if(genomePositionStats != null)
                {
                    Validate.isTrue(mRefBases.charAt(refPos - mGenomeRegion.start()) == genomePositionStats.refBase);

                    CigarOperator cigarOp = genomePositionStats.refBase == alt ? CigarOperator.EQ : CigarOperator.X;

                    mBaseSupports.addPositionSupport(new ReadBaseSupport.PositionSupport(cigarOp,
                            String.valueOf(genomePositionStats.refBase), alt, readIndex,
                            refPos, record.getBaseQualities()[readIndex],
                            genomePositionStats.getAlignedBaseSupport(alt)));
                }
            }
        }

        @Override
        public void handleLeftSoftClip(final SAMRecord record, final CigarElement element)
        {
            // we get the extrapolated ref position
            int startRefPos = record.getAlignmentStart() - element.getLength();

            for(int readIndex = 0; readIndex < element.getLength(); readIndex++)
            {
                int refPos = startRefPos + readIndex;
                char alt = record.getReadString().charAt(readIndex);

                if(alt == N)
                    continue;

                GenomePositionStats genomePositionStats = getPositionStats(refPos);
                if(genomePositionStats != null)
                {
                    Validate.isTrue(mRefBases.charAt(refPos - mGenomeRegion.start()) == genomePositionStats.refBase);
                    mBaseSupports.addPositionSupport(new ReadBaseSupport.PositionSupport(CigarOperator.S,
                            String.valueOf(genomePositionStats.refBase), alt, readIndex, refPos,
                            record.getBaseQualities()[readIndex],
                            genomePositionStats.getSoftClippedBaseSupport(alt)));
                }
            }
        }

        @Override
        public void handleRightSoftClip(final SAMRecord record, final CigarElement element, int startReadIndex, int startRefPosition)
        {
            for(int i = 0; i < element.getLength(); i++)
            {
                int refPos = startRefPosition + i;

                if(refPos > mGenomeRegion.end())
                    return;

                if(refPos < mGenomeRegion.start())
                    continue;

                int readIndex = startReadIndex + i;

                char alt = record.getReadString().charAt(readIndex);

                if(alt == N)
                    continue;

                GenomePositionStats genomePositionStats = getPositionStats(refPos);
                if(genomePositionStats != null)
                {
                    Validate.isTrue(mRefBases.charAt(refPos - mGenomeRegion.start()) == genomePositionStats.refBase);
                    mBaseSupports.addPositionSupport(new ReadBaseSupport.PositionSupport(CigarOperator.S,
                            String.valueOf(genomePositionStats.refBase), alt, readIndex, refPos,
                            record.getBaseQualities()[readIndex],
                            genomePositionStats.getSoftClippedBaseSupport(alt)));
                }
            }
        }
    }
}
