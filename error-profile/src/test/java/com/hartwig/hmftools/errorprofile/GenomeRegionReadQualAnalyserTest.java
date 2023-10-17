package com.hartwig.hmftools.errorprofile;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicLong;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.junit.Assert;
import org.junit.Test;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;

public class GenomeRegionReadQualAnalyserTest
{
    @Test
    public void tesReadQualAnalyser()
    {
        String refBases =
                "CAAACACAGCAAGCAGCTTCCCTCCCTGCTTTGGGGCCTGGAAGGGATAGCAGGAAGTTGACTGGACCAGGGAGATGACCACAGCTGCTGACCTCTCACTCCAGGGAGATGACCACAGCTGCTGACCTCTCACTC";
        FakeRefGenomeCache refGenomeCache = new FakeRefGenomeCache(10000, refBases);

        GenomeRegionReadQualAnalyser analyser = new GenomeRegionReadQualAnalyser(new ChrBaseRegion("chr1", 10000, 10150), refGenomeCache);
        final AtomicLong readTagGenerator = new AtomicLong();

        // Create a SAM record with soft clip on both sides
        SAMRecord record = new SAMRecord(null);
        record.setAlignmentStart(10001);

        // 30M
        record.setCigarString("30M");
        // 30 bases here
        record.setReadString("AAACACAGCAAGCAGCTTCCCTCCCTGCTT");
        record.setReadNegativeStrandFlag(false);
        record.setBaseQualityString("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF");
        record.setMappingQuality(100);
        record.setDuplicateReadFlag(false);
        record.setReadUnmappedFlag(false);
        record.setProperPairFlag(true);
        record.setReadPairedFlag(true);

        List<ReadProfile> readProfileList = new ArrayList<>();
        List<ReadBaseSupport> readBaseSupportList = new ArrayList<>();
        analyser.addReadToStats(record);
        analyser.completeRegion(readTagGenerator, readProfileList::add, readBaseSupportList::add);

        Assert.assertEquals(1, readProfileList.size());

        int refBasesPos = 1;

        // check that read base support contains 30 EQ
        for(ReadBaseSupport.PositionSupport positionSupport : readBaseSupportList.get(0).positionSupports)
        {
            Assert.assertEquals(CigarOperator.EQ, positionSupport.cigarOperator);

            Assert.assertEquals(refBases.charAt(refBasesPos - 1), positionSupport.trinucleotideContext[0]);
            Assert.assertEquals(refBases.charAt(refBasesPos), positionSupport.trinucleotideContext[1]);
            Assert.assertEquals(refBases.charAt(refBasesPos + 1), positionSupport.trinucleotideContext[2]);

            refBasesPos++;
        }
    }

    @Test
    public void tesReadQualAnalyserNegStrand()
    {
        String refBases =
                "CAAACACAGCAAGCAGCTTCCCTCCCTGCTTTGGGGCCTGGAAGGGATAGCAGGAAGTTGACTGGACCAGGGAGATGACCACAGCTGCTGACCTCTCACTCCAGGGAGATGACCACAGCTGCTGACCTCTCACTC";
        FakeRefGenomeCache refGenomeCache = new FakeRefGenomeCache(10000, refBases);

        GenomeRegionReadQualAnalyser analyser = new GenomeRegionReadQualAnalyser(new ChrBaseRegion("chr1", 10000, 10150), refGenomeCache);
        final AtomicLong readTagGenerator = new AtomicLong();

        // Create a SAM record with soft clip on both sides
        SAMRecord record = new SAMRecord(null);
        record.setAlignmentStart(10001);

        // 30M
        record.setCigarString("30M");
        // 30 bases here
        record.setReadString("AAACACAGCAAGCAGCTTCCCTCCCTGCTT");
        record.setReadNegativeStrandFlag(true);
        record.setBaseQualityString("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF");
        record.setMappingQuality(100);
        record.setDuplicateReadFlag(false);
        record.setReadUnmappedFlag(false);
        record.setProperPairFlag(true);
        record.setReadPairedFlag(true);

        List<ReadProfile> readProfileList = new ArrayList<>();
        List<ReadBaseSupport> readBaseSupportList = new ArrayList<>();
        analyser.addReadToStats(record);
        analyser.completeRegion(readTagGenerator, readProfileList::add, readBaseSupportList::add);

        Assert.assertEquals(1, readProfileList.size());
        String refBasesReverseComp = SequenceUtil.reverseComplement(refBases);

        // since we reversed the read, the first is last
        int refBasesPos = refBases.length() - 31;

        // check that read base support contains 30 EQ
        for(ReadBaseSupport.PositionSupport positionSupport : readBaseSupportList.get(0).positionSupports)
        {
            Assert.assertEquals(CigarOperator.EQ, positionSupport.cigarOperator);

            // difference is the trinucleotide context
            Assert.assertEquals(refBasesReverseComp.charAt(refBasesPos - 1), positionSupport.trinucleotideContext[0]);
            Assert.assertEquals(refBasesReverseComp.charAt(refBasesPos), positionSupport.trinucleotideContext[1]);
            Assert.assertEquals(refBasesReverseComp.charAt(refBasesPos + 1), positionSupport.trinucleotideContext[2]);

            refBasesPos++;
        }
    }
}
