package com.hartwig.hmftools.errorprofile;

import static com.hartwig.hmftools.errorprofile.ErrorProfileUtils.getReadStrand;

import static htsjdk.samtools.util.SequenceUtil.N;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.atomic.AtomicLong;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.samtools.CigarHandler;
import com.hartwig.hmftools.common.samtools.CigarTraversal;
import com.hartwig.hmftools.errorprofile.utils.IRefGenomeCache;

import org.apache.commons.lang3.Validate;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
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
public class GenomeRegionReadQualAnalyser
{
    public static final Logger sLogger = LogManager.getLogger(GenomeRegionReadQualAnalyser.class);

    private final ChrBaseRegion mGenomeRegion;

    private final IRefGenomeCache mRefGenomeCache;

    private final GenomePositionStats[] mGenomePositionStats;

    private final GenomePositionStatsAdder mGenomePositionStatsAdder = new GenomePositionStatsAdder();

    private final List<SAMRecord> mReads = new ArrayList<>();

    public GenomeRegionReadQualAnalyser(ChrBaseRegion genomeRegion, IRefGenomeCache refGenomeCache)
    {
        mGenomeRegion = genomeRegion;
        mRefGenomeCache = refGenomeCache;

        // using one array for now
        mGenomePositionStats = new GenomePositionStats[genomeRegion.end() - genomeRegion.start() + 1];
    }

    public void addReadToStats(final SAMRecord read)
    {
        if(read.getReadUnmappedFlag() || read.getDuplicateReadFlag())
            return;

        mReads.add(read);
        CigarTraversal.traverseCigar(read, mGenomePositionStatsAdder);

        // if(mReadCounter > 0 && (mReadCounter % 1000) == 0)
           // purgeBaseDataList(record.getAlignmentStart());
    }

    public void completeRegion(
            final AtomicLong readTagGenerator,
            final Consumer<ReadProfile> readProfileConsumer,
            final Consumer<ReadBaseSupport> readBaseSupportConsumer)
    {
        sLogger.info("region: {} read count: {}", mGenomeRegion, mReads.size());

        List<ReadProfile> readProfiles = new ArrayList<>();
        List<ReadBaseSupport> readBaseSupports = new ArrayList<>();

        // process all reads in this region
        for(SAMRecord read : mReads)
        {
            // we skip reads that are not entirely in this region
            if(read.getAlignmentStart() < mGenomeRegion.start() || read.getAlignmentStart() + read.getReadLength() > mGenomeRegion.end() + 1)
            {
                continue;
            }

            // also skip reads that matches ref genome
                /* if(regionData.stats.isAllMatch(read))
                {
                    continue;
                }*/

            long readTag = readTagGenerator.incrementAndGet();
            ReadBaseSupport readBaseSupport = calcReadBaseSupport(read, readTag);
            readBaseSupportConsumer.accept(readBaseSupport);
            readProfileConsumer.accept(ReadProfiler.profileRead(read, readTag));
        }

        sLogger.info("region: {}", mGenomeRegion);
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
        int length = record.getReadLength();

        boolean allMatch = true;

        for(int i = 0; i < length; ++i)
        {
            byte refBase = mRefGenomeCache.getBase(record.getAlignmentStart() + i);
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

    protected GenomePositionStats getOrCreatePositionStats(int refPosition)
    {
        int posIndex = refPosition - mGenomeRegion.start();
        GenomePositionStats genomePositionStats = mGenomePositionStats[posIndex];

        if(genomePositionStats == null)
        {
            byte ref = mRefGenomeCache.getBase(refPosition);
            genomePositionStats = new GenomePositionStats();
            genomePositionStats.refBase = ref;
            mGenomePositionStats[posIndex] = genomePositionStats;
        }

        return genomePositionStats;
    }

    protected byte[] getTrinucleotideContext(int refPosition)
    {
        Validate.isTrue(refPosition > 1);
        return mRefGenomeCache.getBases(refPosition - 1, refPosition + 1);
    }

    class GenomePositionStatsAdder implements CigarHandler
    {
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

                byte alt = record.getReadBases()[readIndex];

                if(alt == N)
                    continue;

                getOrCreatePositionStats(refPos).addAlignedBase(getReadStrand(record), alt);
            }
        }

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

    // this class builds the base support plus the context
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
                    SequenceUtil.reverseComplement(posSupport.ref);
                    posSupport.alt = SequenceUtil.complement(posSupport.alt);
                    SequenceUtil.reverseComplement(posSupport.trinucleotideContext);
                }

                // resort the position support by their position
                // we rely on this being a stable sort
                mBaseSupports.positionSupports.sort(Comparator.comparingInt(o -> o.readPosition5To3));
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

                if(refPos < mGenomeRegion.start() || refPos <= 1)
                    continue;

                int readIndex = startReadIndex + i;

                byte alt = record.getReadBases()[readIndex];

                if(alt == N)
                    continue;

                GenomePositionStats genomePositionStats = getPositionStats(refPos);
                if(genomePositionStats != null)
                {
                    Validate.isTrue(mRefGenomeCache.getBase(refPos) == genomePositionStats.refBase);

                    CigarOperator cigarOp = genomePositionStats.refBase == alt ? CigarOperator.EQ : CigarOperator.X;

                    mBaseSupports.addPositionSupport(new ReadBaseSupport.PositionSupport(cigarOp,
                            new byte[] { genomePositionStats.refBase }, alt, readIndex,
                            refPos, record.getBaseQualities()[readIndex],
                            getTrinucleotideContext(refPos),
                            genomePositionStats.getAlignedBaseSupport(alt)));
                }
            }
        }

        @Override
        public void handleInsert(final SAMRecord record, final CigarElement e, final int startReadIndex, final int refPos)
        {
            if(refPos < mGenomeRegion.start() || refPos > mGenomeRegion.end() || refPos <= 1)
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
                    byte alt = record.getReadBases()[readIndex];

                    if(alt == N)
                        continue;

                    mBaseSupports.addPositionSupport(new ReadBaseSupport.PositionSupport(CigarOperator.I,
                            new byte[]{}, alt, readIndex, refPos, record.getBaseQualities()[readIndex],
                            getTrinucleotideContext(refPos),
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
            byte[] refBases = mRefGenomeCache.getBases(startRefPos, startRefPos + e.getLength() - 1);

            for(int i = 0; i < e.getLength(); i++)
            {
                int refPos = startRefPos + i;

                if(refPos > mGenomeRegion.end())
                    return;

                if(refPos < mGenomeRegion.start() || refPos <= 1)
                    continue;

                GenomePositionStats genomePositionStats = getPositionStats(refPos);
                if(genomePositionStats != null)
                {
                    Validate.isTrue(refBases[i] == genomePositionStats.refBase);
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
                    refBases, (byte) '-', readIndex, startRefPos, record.getBaseQualities()[readIndex],
                        getTrinucleotideContext(startRefPos),
                        leastSupport));
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

                if(refPos < mGenomeRegion.start() || refPos <= 1)
                    continue;

                byte alt = record.getReadBases()[readIndex];

                if(alt == N)
                    continue;

                GenomePositionStats genomePositionStats = getPositionStats(refPos);
                if(genomePositionStats != null)
                {
                    Validate.isTrue(mRefGenomeCache.getBase(refPos) == genomePositionStats.refBase);
                    mBaseSupports.addPositionSupport(new ReadBaseSupport.PositionSupport(CigarOperator.S,
                            new byte[] { genomePositionStats.refBase }, alt, readIndex, refPos,
                            record.getBaseQualities()[readIndex],
                            getTrinucleotideContext(refPos),
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

                if(refPos < mGenomeRegion.start() || refPos <= 1)
                    continue;

                int readIndex = startReadIndex + i;

                byte alt = record.getReadBases()[readIndex];

                if(alt == N)
                    continue;

                GenomePositionStats genomePositionStats = getPositionStats(refPos);
                if(genomePositionStats != null)
                {
                    Validate.isTrue(mRefGenomeCache.getBase(refPos) == genomePositionStats.refBase);
                    mBaseSupports.addPositionSupport(new ReadBaseSupport.PositionSupport(CigarOperator.S,
                            new byte[] { genomePositionStats.refBase }, alt, readIndex, refPos,
                            record.getBaseQualities()[readIndex],
                            getTrinucleotideContext(refPos),
                            genomePositionStats.getSoftClippedBaseSupport(alt)));
                }
            }
        }
    }
}
