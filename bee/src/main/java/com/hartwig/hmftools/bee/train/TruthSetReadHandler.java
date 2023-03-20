package com.hartwig.hmftools.bee.train;

import java.io.Closeable;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Random;

import com.hartwig.hmftools.bee.BeeInput;
import com.hartwig.hmftools.bee.FragmentGCProfile;
import com.hartwig.hmftools.bee.RefGenomeCompare;
import com.hartwig.hmftools.common.aligner.AlignmentOperator;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.samtools.CigarUtils;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;

// from a bam with consensus reads, get the read id
public class TruthSetReadHandler implements Closeable
{
    private static final Logger sLogger = LogManager.getLogger(TruthSetReadHandler.class);

    private final BeeInputWriter mBeeInputWriter;

    private final int mReadLength;

    private final int mMinUmiGroupSize;

    private final FragmentGCProfile mFragmentGCProfile;

    private final RefGenomeCompare mRefGenomeCompare;

    private final BqrLogic mBqrLogic;

    private final UmiReadGroupDatastore mUmiReadGroupDatastore = new UmiReadGroupDatastore();

    private final ReadKeyTracker mReadKeyTracker = new ReadKeyTracker();

    // track which chromosome position we are up to
    private final ChromosomeProgress mChromosomeProgress = new ChromosomeProgress();

    private final Random mRandom = new Random();

    public TruthSetReadHandler(BeeInputWriter beeInputWriter,
            final Collection<GCProfile> gcProfileList,
            RefGenomeCompare refGenomeCompare,
            BqrLogic bqrLogic,
            int readLength,
            int minUmiGroupSize)
    {
        mBeeInputWriter = beeInputWriter;
        mReadLength = readLength;
        mMinUmiGroupSize = minUmiGroupSize;
        mFragmentGCProfile = new FragmentGCProfile(gcProfileList, readLength);
        mRefGenomeCompare = refGenomeCompare;
        mBqrLogic = bqrLogic;
    }

    @Override
    public void close() throws IOException
    {
        mBeeInputWriter.close();
    }

    public void handleRead1(SAMRecord record) {}

    // we create a data point for each read
    public void handleRead(SAMRecord record)
    {
        if (record.getReferenceIndex().equals(-1))
        {
            // don't process read that are unmapped
            return;
        }

        if (record.isSecondaryOrSupplementary())
        {
            // don't process secondary or supplementary alignments
            return;
        }

        if (!mReadKeyTracker.tryAdd(record.getReadName(), record.getFirstOfPairFlag()))
        {
            return;
        }

        // see if this is a consensus read
        String[] tok = StringUtils.split(record.getReadName(), ':');
        String umi = tok[tok.length - 1];
        boolean isConsensusRead = false;

        if (umi.startsWith("CNS_"))
        {
            // this is consensus read
            isConsensusRead = true;
            umi = umi.substring(4);
        }

        // we only process either consensus reads or duplicate reads
        if (!isConsensusRead && !record.getDuplicateReadFlag())
        {
            return;
        }

        // TODO: what about hard clip?
        // start position is the aligned position of the read start
        int read5PrimeReferencePos;

        if (record.getReadNegativeStrandFlag())
        {
            read5PrimeReferencePos = record.getAlignmentEnd() + CigarUtils.rightSoftClipLength(record);
        }
        else
        {
            read5PrimeReferencePos = record.getAlignmentStart() - CigarUtils.leftSoftClipLength(record);
        }

        UmiReadGroup umiGroup = mUmiReadGroupDatastore.getOrCreate(umi, record.getReferenceName(),
                read5PrimeReferencePos, record.getReadNegativeStrandFlag());

        if (isConsensusRead)
        {
            umiGroup.consensusRead = record;
        }
        else
        {
            umiGroup.duplicateReads.add(record);
        }
    }

    synchronized public void onRegionComplete(ChrBaseRegion chrBaseRegion)
    {
        // we have to avoid processing reads that are not in this designated region, not sure how to do that?
        // a region is completed, we can process those groups
        if (mChromosomeProgress.onRegionCompleted(chrBaseRegion))
        {
            int upToPos = mChromosomeProgress.getChromosomePosition(chrBaseRegion.chromosome()).completedPosition;

            sLogger.info("chromosome region completed: {}, progress head: {}", chrBaseRegion, upToPos);

            Collection<UmiReadGroup> completedUmiGroups = mUmiReadGroupDatastore.popCompleted(
                    chrBaseRegion.chromosome(), upToPos - 10 * mReadLength);

            for (UmiReadGroup g : completedUmiGroups)
            {
                sLogger.info("completed group umi({}) num reads({})", g.umi, g.duplicateReads.size());

                // check that it has consensus read
                if (g.consensusRead == null)
                {
                    if (g.duplicateReads.size() > 1)
                    {
                        sLogger.error("consensus read for umi group(umi({}), aligned({}:{})) not found: ", g.umi, g.chromosome, g.read5PrimeReferencePos);

                        for (SAMRecord r : g.duplicateReads)
                        {
                            sLogger.error("read: {}", r);
                        }
                    }

                    continue;
                }

                if (g.duplicateReads.size() < mMinUmiGroupSize)
                {
                    // not enough reads to be 100% sure of the consensus
                    continue;
                }

                if (mRandom.nextInt(100) != 1)
                {
                    // only 1% of groups needed for now
                    continue;
                }

                var cigarErrors = g.consensusRead.getCigar().isValid(g.consensusRead.getReadName(), 0);

                if (cigarErrors != null && !cigarErrors.isEmpty())
                {
                    sLogger.error("consensus read({}): invalid cigar({})", g.consensusRead, g.consensusRead.getCigarString());
                    continue;
                }

                for (SAMRecord dupRead: g.duplicateReads)
                {
                    // each read we need to store them somehow
                    BeeInput beeInput = new BeeInput();

                    beeInput.readString = dupRead.getReadString();
                    beeInput.baseQualityString = dupRead.getBaseQualityString();
                    beeInput.isRead1 = dupRead.getFirstOfPairFlag();
                    beeInput.isMapped = !dupRead.getReadUnmappedFlag();
                    beeInput.isMateMapped = !dupRead.getMateUnmappedFlag();
                    beeInput.mapQuality = dupRead.getMappingQuality();

                    beeInput.referenceGenomeBases = mRefGenomeCompare.getReferenceGenomeBases(dupRead);
                    beeInput.bqrBaseQualities = mBqrLogic.calcBqrQualities(dupRead);

                    // for ML it is easier to learn the patterns by the direction a read is sequenced
                    if (dupRead.getReadNegativeStrandFlag())
                    {
                        beeInput.readString = SequenceUtil.reverseComplement(beeInput.readString);
                        beeInput.baseQualityString = StringUtils.reverse(beeInput.baseQualityString);
                        beeInput.referenceGenomeBases = StringUtils.reverse(beeInput.referenceGenomeBases);
                        ArrayUtils.reverse(beeInput.bqrBaseQualities);
                    }

                    @Nullable
                    GCProfile gcProfile = mFragmentGCProfile.selectFragmentGCProfile(dupRead);

                    if (gcProfile != null)
                    {
                        beeInput.gcContent = gcProfile.gcContent();
                        beeInput.mappability = gcProfile.mappablePercentage();
                    }
                    else
                    {
                        beeInput.gcContent = Double.NaN;
                        beeInput.mappability = Double.NaN;
                    }

                    beeInput.cigar = dupRead.getCigarString();

                    beeInput.baseError = new boolean[dupRead.getReadLength()];

                    List<AlignmentOperator> alignOps = ReadCigarAligner.alignReads(dupRead, g.consensusRead);

                    int i = 0;

                    // go through the align ops
                    for (AlignmentOperator op : alignOps)
                    {
                        switch (op)
                        {
                            case MATCH:
                                beeInput.baseError[i++] = false;
                                break;
                            case MISMATCH:
                            case INSERTION:
                                beeInput.baseError[i++] = true;
                                break;
                            case DELETION:
                                // do nothing for now
                                break;
                        }
                    }

                    // deep debugging: log out the alignment to make sure we got it correct
                    if (false)
                    {
                        AlignmentOperator.logAlignment(sLogger, Level.TRACE, dupRead.getReadString(),
                                g.consensusRead.getReadString(), alignOps);

                        var sb = new StringBuilder();

                        for (var b : beeInput.baseError)
                        {
                            sb.append(b ? 1 : 0);
                        }

                        sLogger.trace("{}", sb);
                    }

                    try
                    {
                        mBeeInputWriter.write(beeInput);
                    }
                    catch (IOException e)
                    {
                        // changed to unchecked exception to get around Consumer signature
                        throw new RuntimeException(String.format("IO Exception: %s", e));
                    }
                }
            }

            sLogger.info("chromosome region: {} completed umi groups: {}", chrBaseRegion, completedUmiGroups.size());
        }
    }
}
