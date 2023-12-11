// TODO: Not in final
//package com.hartwig.hmftools.sage.sagevis;
//
//import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
//import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
//
//import static org.junit.Assert.assertEquals;
//import static org.junit.Assert.assertNull;
//import static org.junit.Assert.assertTrue;
//
//import java.io.File;
//
//import com.hartwig.hmftools.common.test.SamRecordTestUtils;
//
//import org.junit.Test;
//
//import htsjdk.samtools.DefaultSAMRecordFactory;
//import htsjdk.samtools.SAMFileHeader;
//import htsjdk.samtools.SAMLineParser;
//import htsjdk.samtools.SAMRecord;
//import htsjdk.samtools.SamReader;
//import htsjdk.samtools.ValidationStringency;
//
//public class SageVisTest
//{
//    @Test
//    public void testSamStringDeserialization()
//    {
//        final String readId = "READ_001";
//        String readBases = "A".repeat(100);
//        String readCigar = "100M";
//        int numMutations = 5;
//        SAMRecord expectedRead = SamRecordTestUtils.createSamRecord(
//                readId, CHR_1, 100, readBases, readCigar, CHR_1, 600, false,
//                false, null, true, readCigar);
//        expectedRead.setAttribute(NUM_MUTATONS_ATTRIBUTE, numMutations);
//
//        String expectedSamString = expectedRead.getSAMString();
//
//        // samFileHeader only used for setting ref index and mate ref index
//        SAMLineParser parser = new SAMLineParser(new DefaultSAMRecordFactory(), ValidationStringency.SILENT, new SAMFileHeader(), (SamReader)null, (File)null);
//        SAMRecord actualRead = parser.parseLine(expectedSamString.trim());
//        String actualSamString = actualRead.getSAMString();
//
//        assertEquals(expectedSamString, actualSamString);
//    }
//}