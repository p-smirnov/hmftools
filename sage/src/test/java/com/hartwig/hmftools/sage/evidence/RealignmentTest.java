package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildDefaultBaseQuals;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_FLANK_LENGTH;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.FULL;
import static com.hartwig.hmftools.sage.common.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.sage.common.TestUtils.REF_SEQUENCE_200;
import static com.hartwig.hmftools.sage.common.TestUtils.buildCigarString;
import static com.hartwig.hmftools.sage.common.TestUtils.buildSamRecord;
import static com.hartwig.hmftools.sage.common.TestUtils.createSamRecord;
import static com.hartwig.hmftools.sage.common.VariantUtils.createSimpleVariant;
import static com.hartwig.hmftools.sage.evidence.Realignment.realigned;
import static com.hartwig.hmftools.sage.evidence.RealignedType.EXACT;
import static com.hartwig.hmftools.sage.evidence.RealignedType.LENGTHENED;
import static com.hartwig.hmftools.sage.evidence.RealignedType.NONE;
import static com.hartwig.hmftools.sage.evidence.RealignedType.SHORTENED;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.sage.common.ReadContextMatch;
import com.hartwig.hmftools.sage.common.ReadContextMatcher;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.RegionTaskTester;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.common.VariantReadContextBuilder;
import com.hartwig.hmftools.sage.pipeline.RegionTask;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class RealignmentTest
{
    @Test
    public void testIndelRealignment()
    {
        String refBases = "X" + REF_BASES_200.substring(0, 100)
                // 0123456789012345678901234567890123456789
                //             ->
                + "GTTGTTGTTGTCGTTGTTGTTGTTGGTTTTTCTGAGACAGAGTC" + REF_BASES_200.substring(0, 100);

        RefSequence refSequence = new RefSequence(0, refBases.getBytes());

        String insertedBases = "TTGTTGTTGTTGGTTTTTC";

        String varBuildReadBases = refBases.substring(91, 101) + refBases.substring(101, 114) + insertedBases + refBases.substring(114, 150);

        String readCigar = "33M19I36M";

        SAMRecord varBuildRead = buildSamRecord(91, readCigar, varBuildReadBases);
        String ref = refBases.substring(113, 114);

        VariantReadContextBuilder builder = new VariantReadContextBuilder(DEFAULT_FLANK_LENGTH);

        SimpleVariant var = createSimpleVariant(113, ref, ref + insertedBases);
        VariantReadContext readContext = builder.createContext(var, varBuildRead, 22, refSequence);

        assertEquals(113, readContext.variant().position());
        assertEquals(11, readContext.VarReadIndex);
        assertEquals("CGTTGTTGTTGTTGGTTTTTCTTGTTGTTGTTGGTTTTTCTGA", readContext.coreStr());

        // variant is at read index 41 when should be 22
        int realignedStartPos = 72;
        readCigar = buildCigarString(varBuildReadBases.length());
        SAMRecord realignedRead = buildSamRecord(realignedStartPos, readCigar, varBuildReadBases);

        ReadContextMatcher matcher = new ReadContextMatcher(readContext);

        int realignedReadIndex = Realignment.realignedReadIndexPosition(readContext, realignedRead);
        ReadContextMatch matchType = matcher.determineReadMatch(realignedRead, realignedReadIndex);

        assertEquals(FULL, matchType);
    }

    /* REALIGN

    @Test
    public void testIndelRealignedLeft()
    {
        ChrBaseRegion region = new ChrBaseRegion(CHR_1, 1, 450);

        RegionTaskTester tester = new RegionTaskTester();

        RegionTask task = tester.createRegionTask(region);

        // insert at position 288
        String refBases = "X" + generateRandomBases(187)
                + "CTACCCAAATCTTATAAAATGGCCCCACCCATATCTCCCTTCGCTGACTCTCTATTCGGACTCAGCCCGCCTGCACCCAGGTGAAATAAACAGCCTTGTT"
                + "GCACACACACACACACACACACACAAAGAATTCATTGCTCTTAGCTTTTGCACTGTTTCCTTTTTTAGACCTTCCTACATAAGTATTTGTGTATATGTCTGTATTT";

        tester.RefGenome.RefGenomeMap.put(CHR_1, refBases + generateRandomBases(1500));

        SAMRecord read0 = createSamRecord("READ_00", CHR_1, 257, // 24220, actual pos = 257
                "CCTGCACCCAGTTGAAATAAACAGCCTTGTTGCACACACACACACACACACACACACACACACACACACAAAGAATTCATTGCTCTTAGCTTTTGCACTGTTTCC"
                        + "TTTTTTAGACCTTCCTACATAAGTATTTGTGTATATGTCTGTATTT",
                "32M14I105M");

        tester.TumorSamSlicer.ReadRecords.add(read0);
        tester.TumorSamSlicer.ReadRecords.add(read0);

        SAMRecord read1 = createSamRecord("READ_01", CHR_1, 188, // 24220, actual pos =
                "CTACCCAAATCTTATAAAATGGCCCCACCCATATCTCCCTTCGCTGACTCTCTATTCGGACTCAGCCCGCCTGCACCCAGGTGAAATAAACAGCCTTGTTGC"
                        + "ACACACACACACACACACACACACACACACACACACAAAGAATTCATTG",
                "125M26S");

        tester.TumorSamSlicer.ReadRecords.add(read1);
        tester.TumorSamSlicer.ReadRecords.add(read1);

        task.run();

        SageVariant var1 = task.getVariants().stream().filter(x -> x.position() == 288).findFirst().orElse(null);
        SageVariant var2 = task.getVariants().stream().filter(x -> x.position() == 312).findFirst().orElse(null);
        assertNotNull(var1);
        assertNotNull(var2);
        assertEquals(4, var1.tumorReadCounters().get(0).readSupportCounts().Full);
        assertEquals(2, var2.tumorReadCounters().get(0).readSupportCounts().Full);
        assertEquals(2, var2.tumorReadCounters().get(0).readSupportCounts().Realigned);
    }

    @Test
    public void testIndelRealignedRight()
    {
        ChrBaseRegion region = new ChrBaseRegion(CHR_1, 1, 450);

        RegionTaskTester tester = new RegionTaskTester();

        RegionTask task = tester.createRegionTask(region);

        String refBases = "X"
                + "AATTTGAACCTAATTTTTTTTTTTCTTCTGCACTAACATGCCTGTTGAACCATTTGGACTTAACTTTTGTGCATGGTGTGAAATAGGTGCCCAGCCTCATTCTTTTGCATGTAGATAT"
                + "CCTTTTCCCAGCACCATTCGTTGAATGGAGACTATTCTTTCCCCACTGAATAGTCTTGGTACCCTCTTTGAAAATCAATTGATGATAAATAGATGTGTTTATTTCTGAACTCTCCATTT"
                + "TATTCCATTGACCTATATCTCTCCTTATGCCAGTTTTTATTACTGTGCAGTTTTGATTACTAC";

        tester.RefGenome.RefGenomeMap.put(CHR_1, refBases + generateRandomBases(1500));

        SAMRecord read0 = createSamRecord("READ_00", CHR_1, 28, // 11490
                "CTGCACTAACATGCCTGTTGAACCATTTGGACTTAACTTTTGTGCATGGTGTGAAATAGGTGCCCAGCCTCTTTTGTGCATGGTGTGAAATAGGTGCCCAGCCTCATT"
                        + "CTTTTGCATGTAGATATCCTTTTCCCAGCACCATTCGTTGAAT",
                "36M34I81M");

        tester.TumorSamSlicer.ReadRecords.add(read0);
        tester.TumorSamSlicer.ReadRecords.add(read0);

        SAMRecord read1 = createSamRecord("READ_01", CHR_1, 14, // 15718
                "ATTTTTTTTTTCTTCTGCACTAACATGCCTGTTGAACCATTTGGACTTAACTTTTGTGCATGGTGTGAAATAGGTGCCCAGCCTCTTTTGTGCATGGTGTGAAATAGGTG"
                        + "CCCAGCCTCATTCTTTTGCATGTAGATATCCTTTTCCCAGC",
                "85M66S");

        tester.TumorSamSlicer.ReadRecords.add(read1);
        tester.TumorSamSlicer.ReadRecords.add(read1);

        SAMRecord read2 = createSamRecord("READ_02", CHR_1, 64, // 22780
                "TTGGACTTAACTTTTGTGCATGGTGTGAAATAGGTGCCCAGCCTCTTTTGTGCATGGTGTGAAATAGGTGCCCAGCCTCATTCTTTTGCATGTAGATATCCTTTTCCCAG"
                        + "CACCATTCGTTGAATGGAGACTATTCTTTCCCCACTGAATA",
                "44S107M");

        tester.TumorSamSlicer.ReadRecords.add(read2);
        tester.TumorSamSlicer.ReadRecords.add(read2);

        task.run();

        SageVariant var1 = task.getVariants().stream().filter(x -> x.position() == 63).findFirst().orElse(null);
        SageVariant var2 = task.getVariants().stream().filter(x -> x.position() == 98).findFirst().orElse(null);
        assertNotNull(var1);
        assertNotNull(var2);
        assertEquals(4, var1.tumorReadCounters().get(0).readSupportCounts().Full);
        assertEquals(2, var1.tumorReadCounters().get(0).readSupportCounts().Partial);

        assertEquals(2, var2.tumorReadCounters().get(0).readSupportCounts().Full);
        assertEquals(4, var2.tumorReadCounters().get(0).readSupportCounts().Realigned);
    }

    @Test
    public void testRealignedCoreInDelete()
    {
        ChrBaseRegion region = new ChrBaseRegion(CHR_1, 1, 450);

        RegionTaskTester tester = new RegionTaskTester();

        RegionTask task = tester.createRegionTask(region);

        String refBases = "TGGGGCCTCGCTTCTTCGAGCTGTCCCCCGGTGAGCTGGGCTTGTCATCCACTCTGCTGGTACCCCGCTCTCGTTCCCTCTCACGTTCCCGCTCCCTCTCTCGCTCCCTCTCCCTTTCT"
                + "CGATCCCGCTCCCGGTCCCTATCCCGGTCTCGGTCCCGGTCCCGAGGGGCCTCTGCGCGGCAGGTGCTGCTGGTGCTGGTGCTGCTCGCTGGCGGGGACAAACCCATGTTTCGAAGGCC"
                + "CTCCCGGGCCCGGGAAGTGGAGGCCAGGTCCTGGGGCGAGCGGCCAGGGTAGGGGATCCGGAT";

        tester.RefGenome.RefGenomeMap.put(CHR_1, refBases + generateRandomBases(1500));

        SAMRecord read0 = createSamRecord("READ_00", CHR_1, 4, // 22498
                "GCCTCGCTTCTTCGAGCTGTCCCCCGGTGAGCTGGGCTTGTCATCCACTCTGCTGGTACCCCGCTCTCGTTCCCTCTCACGTTCCCGCTCCCTCTCTCGCTCCCTCTCCCT"
                        + "TTCTCGATCCCGGTCTCGGTCCCGGTCCCGAGGGGCCTCT",
                "117M18D34M");

        tester.TumorSamSlicer.ReadRecords.add(read0);
        tester.TumorSamSlicer.ReadRecords.add(read0);

        SAMRecord read1 = createSamRecord("READ_01", CHR_1, 1, // 10473
                "GGGGCCTCGCTTCTTCGAGCTGTCCCCCGGTGAGCTGGGCTTGTCATCCACTCTGCTGGTACCCCGCTCTCGTTCCCTCTCACGTTCCCGCTCCCTCTCTCGCTCCCTCTC"
                        + "CCTTTCTCGATCCCGGTCTCGGTCCCGGTCCCG",
                "144M");

        tester.TumorSamSlicer.ReadRecords.add(read1);
        tester.TumorSamSlicer.ReadRecords.add(read1);

        // readBases(GGGGCTTGGGGCCTCGCTTCTTCGAGCTGTCCCCCGGTGAGCTGGGCTTGTCATCCACTCTGCTGGTACCCCGCTCTCGTTCCCTCTCACGTTCCCGCTCCCTCTCTCGCTCCCTCTCCCTTTCTCGATCCCGGTCTCGGTCCCGGTCCCG)
        // 11:29:20.883 [TRACE] var(17:77769130 C>T) readContext(132-134-136) support(FULL) read(idx=136 posStart=77768994 cigar=151M
        // id=A00260:46:HGWWNDSXX:4:2324:12310:10473)
        // readBases(GGGGCTTGGGGCCTCGCTTCTTCGAGCTGTCCCCCGGTGAGCTGGGCTTGTCATCCACTCTGCTGGTACCCCGCTCTCGTTCCCTCTCACGTTCCCGCTCCCTCTCTCGCTCCCTCTCCCTTTCTCGATCCCGGTCTCGGTCCCGGTCCCG)

        task.run();

        SageVariant var1 = task.getVariants().stream().filter(x -> x.position() == 120).findFirst().orElse(null);
        SageVariant var2 = task.getVariants().stream().filter(x -> x.position() == 127).findFirst().orElse(null);
        SageVariant var3 = task.getVariants().stream().filter(x -> x.position() == 130).findFirst().orElse(null);
        assertNotNull(var1);
        assertNotNull(var2);
        assertNotNull(var3);
        assertEquals(4, var1.tumorReadCounters().get(0).readSupportCounts().Full);

        assertEquals(2, var2.tumorReadCounters().get(0).readSupportCounts().Full);
        assertEquals(2, var2.tumorReadCounters().get(0).readSupportCounts().Realigned);
        assertEquals(2, var3.tumorReadCounters().get(0).readSupportCounts().Full);
        assertEquals(2, var3.tumorReadCounters().get(0).readSupportCounts().Realigned);
    }

    @Test
    public void testRealignedTooShort()
    {
        String sequence = "GAGAGTGTGTGTGTGTGTCTGTGTGTATGTATATATATATATATATATATCACATTTTT";
        String truncatedAtEnd = "GAGAGTGTGTGTGTGTGTCTGTGTGTATGTATATATATATATATATATATCACATTTT";
        String truncatedAtStart = "GAGTGTGTGTGTGTGTCTGTGTGTATGTATATATATATATATATATATCACATTTTT";
        int startIndex = 0;
        int endIndex = sequence.length() - 1;

        assertRealigned(EXACT, realigned(startIndex, endIndex, sequence.getBytes(), startIndex, sequence.getBytes(), 0));
        assertRealigned(NONE, realigned(startIndex, endIndex, sequence.getBytes(), startIndex, truncatedAtEnd.getBytes(), 0));
        assertRealigned(NONE, realigned(startIndex, endIndex, sequence.getBytes(), -2, truncatedAtStart.getBytes(), 0));
    }

    @Test
    public void testRealigned()
    {
        String sequence = "GAGAGTGTGTGTGTGTGTCTGTGTGTATGTATATATATATATATATATATCACATTTTTATTATTG";
        int startIndex = 3;
        int endIndex = startIndex + 55;

        assertRealigned(EXACT, realigned(startIndex, endIndex, sequence.getBytes(), startIndex, sequence.getBytes(), 0));

        assertRealigned(NONE, realigned(startIndex, endIndex, sequence.getBytes(), startIndex - 1, sequence.getBytes(), 0));
        assertRealigned(NONE, realigned(startIndex, endIndex, sequence.getBytes(), startIndex + 1, sequence.getBytes(), 0));

        assertRealigned(NONE, realigned(startIndex, endIndex, sequence.getBytes(), startIndex - 2, sequence.getBytes(), 1));
        assertRealigned(EXACT, realigned(startIndex, endIndex, sequence.getBytes(), startIndex - 1, sequence.getBytes(), 1));
        assertRealigned(EXACT, realigned(startIndex, endIndex, sequence.getBytes(), startIndex + 1, sequence.getBytes(), 1));
        assertRealigned(NONE, realigned(startIndex, endIndex, sequence.getBytes(), startIndex + 2, sequence.getBytes(), 1));

        assertRealigned(EXACT, realigned(startIndex, endIndex, sequence.getBytes(), startIndex - 2, sequence.getBytes(), 2));
        assertRealigned(EXACT, realigned(startIndex, endIndex, sequence.getBytes(), startIndex + 2, sequence.getBytes(), 2));
    }

    @Test
    public void testPolyA()
    {
        String shorter = "GATCAAAAAAAAAGATC";
        String ref = "GATCAAAAAAAAAAGATC";
        String longer = "GATCAAAAAAAAAAAGATC";

        int startIndex = 0;
        int endIndex = ref.length() - 1;

        assertRealigned(EXACT, realigned(startIndex, endIndex, ref.getBytes(), startIndex, ref.getBytes(), 10));
        assertRealigned(SHORTENED, 9, realigned(startIndex, endIndex, ref.getBytes(), startIndex, shorter.getBytes(), 10));
        assertRealigned(LENGTHENED, 11, realigned(startIndex, endIndex, ref.getBytes(), startIndex, longer.getBytes(), 10));
    }

    @Test
    public void testDiNucleotideRepeat()
    {
        String shorter = "GATCATATATATGATC";
        String ref = "GATCATATATATATGATC";
        String longer = "GATCATATATATATATGATC";

        int startIndex = 0;
        int endIndex = ref.length() - 1;

        assertRealigned(EXACT, realigned(startIndex, endIndex, ref.getBytes(), startIndex, ref.getBytes(), 10));
        assertRealigned(SHORTENED, 4, realigned(startIndex, endIndex, ref.getBytes(), startIndex, shorter.getBytes(), 10));
        assertRealigned(LENGTHENED, 6, realigned(startIndex, endIndex, ref.getBytes(), startIndex, longer.getBytes(), 10));
    }

    private static void assertRealigned(RealignedType expectedType, int expectedCount, RealignedContext context)
    {
        assertEquals(expectedCount, context.RepeatCount);
        assertEquals(expectedType, context.Type);
    }

    private static void assertRealigned(RealignedType expectedType, RealignedContext context)
    {
        assertEquals(expectedType, context.Type);
    }
    */
}
