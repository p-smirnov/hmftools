package com.hartwig.hmftools.bee.train;

import java.io.IOException;

import com.hartwig.hmftools.bee.BeeInput;
import com.hartwig.hmftools.bee.BeeUtils;
import com.hartwig.hmftools.common.genome.gc.GCProfile;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;

public class TrainDataWriter
{
    private static final Logger sLogger = LogManager.getLogger(TrainDataWriter.class);

    public void processUmiReadGroup(UmiReadGroup umiGroup)
    {
        // TODO: check that consensus read is there

        sLogger.debug("handle umi group: umi({}), mapped({}:{})",
                umiGroup.umi, umiGroup.consensusRead.getReferenceName(),
                umiGroup.consensusRead.getAlignmentStart());

        // we got a umi read group, compare each base to consensus
        // get the consensus read
        SAMRecord consensus = umiGroup.consensusRead;
        String consensusReadString = consensus.getReadString();

        for (SAMRecord dupRead : umiGroup.duplicateReads)
        {
            // we need to measure the read bases
            String readString = dupRead.getReadString();

            // they should align perfectly, but if they do not, say with
            // insertions then I am not sure what I can do?

        }
    }

    /*
    public void something(SAMRecord record)
    {
        sLogger.debug("handle read");

        // each read we need to store them somehow
        BeeInput beeInput = new BeeInput();

        beeInput.readString = record.getReadString();
        beeInput.baseQualityString = record.getBaseQualityString();
        beeInput.isRead1 = record.getFirstOfPairFlag();
        beeInput.isMapped = !record.getReadUnmappedFlag();
        beeInput.isMateMapped = !record.getMateUnmappedFlag();
        beeInput.mapQuality = record.getMappingQuality();

        beeInput.matchesRef = mRefGenomeCompare.getMatchRef(record);

        // for ML it is easier to learn the patterns by the direction a read is sequenced
        if (record.getReadNegativeStrandFlag())
        {
            beeInput.readString = SequenceUtil.reverseComplement(beeInput.readString);
            beeInput.baseQualityString = StringUtils.reverse(beeInput.baseQualityString);
            ArrayUtils.reverse(beeInput.matchesRef);
        }

        @Nullable
        GCProfile gcProfile = mFragmentGCProfile.selectFragmentGCProfile(record);

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

        beeInput.cigar = record.getCigarString();

        beeInput.target = new int[record.getReadLength()];

        // create the target vector
        int polyGcCount = BeeUtils.polyGCTail(beeInput.readString, false);

        if (polyGcCount > 0)
        {
            for (int i = 0; i < polyGcCount; ++i)
            {
                beeInput.target[beeInput.target.length - 1 - i] = 1;
            }
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
     */
}
