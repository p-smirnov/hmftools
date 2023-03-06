package com.hartwig.hmftools.bee.train;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.aligner.AlignmentOperator;
import com.hartwig.hmftools.common.samtools.CigarUtils;

import org.junit.Before;
import org.junit.Test;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;

public class ReadCigarAlignerTest
{
    @Before
    public void setUp()
    {
        org.apache.logging.log4j.core.config.Configurator.setRootLevel(org.apache.logging.log4j.Level.TRACE);
    }

    @Test
    public void testFindIntersectingMRegion()
    {
        Cigar subjectCigar = TextCigarCodec.decode("151M");
        int subjectAlignStart = 10000;
        Cigar refCigar = TextCigarCodec.decode("100M");
        int refAlignStart = 10025;

        List<ReadCigarAligner.ReadRegion> readRegions = ReadCigarAligner.findAlignedRegions(
                subjectCigar, subjectAlignStart, refCigar, refAlignStart);

        assertEquals(1, readRegions.size());
        assertEquals(10025, readRegions.get(0).referenceGenomeStart);
        assertEquals(25, readRegions.get(0).subjectReadStart);
        assertEquals(0, readRegions.get(0).refReadStart);
        assertEquals(100, readRegions.get(0).getRefReadLength());

        subjectCigar = TextCigarCodec.decode("151M");
        subjectAlignStart = 10000;
        refCigar = TextCigarCodec.decode("10S100M");
        refAlignStart = 10025;

        readRegions = ReadCigarAligner.findAlignedRegions(
                subjectCigar, subjectAlignStart, refCigar, refAlignStart);

        assertEquals(1, readRegions.size());
        assertEquals(10025, readRegions.get(0).referenceGenomeStart);
        assertEquals(25, readRegions.get(0).subjectReadStart);
        assertEquals(10, readRegions.get(0).refReadStart);
        assertEquals(100, readRegions.get(0).getRefReadLength());

        // split into two sections, first one is 25, second one is 40
        //
        // sub read           |-----50M-----|10I|---40M----|
        // ref read                  |--------------100M-----------|
        //                    |      |      |              |       |
        // ref genome pos   10000  10025  10050          10090   10125
        // intersect                 |--25--|   |----40----|
        subjectCigar = TextCigarCodec.decode("50M10I40M");
        subjectAlignStart = 10000;
        refCigar = TextCigarCodec.decode("10S100M");
        refAlignStart = 10025;

        readRegions = ReadCigarAligner.findAlignedRegions(
                subjectCigar, subjectAlignStart, refCigar, refAlignStart);

        assertEquals(2, readRegions.size());
        assertEquals(10025, readRegions.get(0).referenceGenomeStart);
        assertEquals(25, readRegions.get(0).subjectReadStart);
        assertEquals(10, readRegions.get(0).refReadStart);
        assertEquals(25, readRegions.get(0).getRefReadLength());

        assertEquals(10050, readRegions.get(1).referenceGenomeStart);
        assertEquals(60, readRegions.get(1).subjectReadStart);
        assertEquals(35, readRegions.get(1).refReadStart);
        assertEquals(40, readRegions.get(1).getRefReadLength());

        // test that reversing the subject and ref produces equivalent results
        readRegions = ReadCigarAligner.findAlignedRegions(
                refCigar, refAlignStart, subjectCigar, subjectAlignStart);

        assertEquals(2, readRegions.size());
        assertEquals(10025, readRegions.get(0).referenceGenomeStart);
        assertEquals(25, readRegions.get(0).refReadStart);
        assertEquals(10, readRegions.get(0).subjectReadStart);
        assertEquals(25, readRegions.get(0).getRefReadLength());

        assertEquals(10050, readRegions.get(1).referenceGenomeStart);
        assertEquals(60, readRegions.get(1).refReadStart);
        assertEquals(35, readRegions.get(1).subjectReadStart);
        assertEquals(40, readRegions.get(1).getRefReadLength());


        // split into two sections, first one is 25, second one is 40
        //
        // sub read           |-----50M-----|10D|---40M----|
        // ref read                  |--------------100M-----------|
        //                    |      |      |   |          |       |
        // ref genome pos   10000  10025  10050 60       10090   10125
        // intersect                 |--25--|   |----40----|

        subjectCigar = TextCigarCodec.decode("50M10D40M");
        subjectAlignStart = 10000;
        refCigar = TextCigarCodec.decode("10S100M");
        refAlignStart = 10025;

        readRegions = ReadCigarAligner.findAlignedRegions(
                subjectCigar, subjectAlignStart, refCigar, refAlignStart);

        assertEquals(2, readRegions.size());
        assertEquals(10025, readRegions.get(0).referenceGenomeStart);
        assertEquals(25, readRegions.get(0).subjectReadStart);
        assertEquals(10, readRegions.get(0).refReadStart);
        assertEquals(25, readRegions.get(0).getRefReadLength());

        assertEquals(10060, readRegions.get(1).referenceGenomeStart);
        assertEquals(50, readRegions.get(1).subjectReadStart);
        assertEquals(45, readRegions.get(1).refReadStart);
        assertEquals(40, readRegions.get(1).getRefReadLength());

        // test that reversing the subject and ref produces equivalent results
        readRegions = ReadCigarAligner.findAlignedRegions(
                refCigar, refAlignStart, subjectCigar, subjectAlignStart);

        assertEquals(2, readRegions.size());
        assertEquals(10025, readRegions.get(0).referenceGenomeStart);
        assertEquals(25, readRegions.get(0).refReadStart);
        assertEquals(10, readRegions.get(0).subjectReadStart);
        assertEquals(25, readRegions.get(0).getRefReadLength());

        assertEquals(10060, readRegions.get(1).referenceGenomeStart);
        assertEquals(50, readRegions.get(1).refReadStart);
        assertEquals(45, readRegions.get(1).subjectReadStart);
        assertEquals(40, readRegions.get(1).getRefReadLength());
    }

    @Test
    public void testSimpleReadAlign()
    {
        String subjectRead = "ACTTCTAGGCTTAGA";
        Cigar subjectCigar = TextCigarCodec.decode("15M");
        int subjectAlignStart = 10000;
        String refRead = "ACTTCTAGGCTTAGA";
        Cigar refCigar = TextCigarCodec.decode("15M");
        int refAlignStart = 10000;

        List<AlignmentOperator> ops = ReadCigarAligner.alignReads(subjectRead, subjectCigar, subjectAlignStart, refRead, refCigar, refAlignStart);

        assertEquals("MMMMMMMMMMMMMMM", AlignmentOperator.toString(ops));
    }

    @Test
    public void testReadAlignMismatch()
    {
        String subjectRead = "ACTTGTAGGCTTAGA";
        String refRead =     "ACTTCTAGGATTAGA";
        Cigar subjectCigar = TextCigarCodec.decode("15M");
        int subjectAlignStart = 10000;
        Cigar refCigar = TextCigarCodec.decode("15M");
        int refAlignStart = 10000;

        List<AlignmentOperator> ops = ReadCigarAligner.alignReads(subjectRead, subjectCigar, subjectAlignStart, refRead, refCigar, refAlignStart);

        // 2 mismatch
        assertEquals("MMMMSMMMMSMMMMM", AlignmentOperator.toString(ops));
    }

    @Test
    public void testReadAlignComplexAlignment()
    {
        // there are left and right soft clip and also middle alignment part is split into two
        String subjectRead = "ACTGGTATGCACTTGTAGGCTTAGATAGGC";
        String refRead =     "ACCGGTATGCACTTGTAGGCTTAGATAGGGC";
        Cigar subjectCigar = TextCigarCodec.decode("10S-7M-5D-8M-5S".replace("-", ""));
        int subjectAlignStart = 10010;
        Cigar refCigar = TextCigarCodec.decode("17M-1D-10M-4S".replace("-", ""));
        int refAlignStart = 10000;

        assertEquals(subjectRead.length(), subjectCigar.getReadLength());
        assertEquals(refRead.length(), refCigar.getReadLength());

        List<AlignmentOperator> ops = ReadCigarAligner.alignReads(subjectRead, subjectCigar, subjectAlignStart, refRead, refCigar, refAlignStart);

        // 2 mismatch
        assertEquals("MMSMMMMMMMMMMMMMM----SSSSMMM+++MMM", AlignmentOperator.toString(ops));
    }

}
