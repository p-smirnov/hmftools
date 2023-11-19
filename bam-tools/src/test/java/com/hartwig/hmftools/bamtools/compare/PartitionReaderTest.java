package com.hartwig.hmftools.bamtools.compare;

import static com.hartwig.hmftools.bamtools.compare.MismatchType.NEW_ONLY;
import static com.hartwig.hmftools.bamtools.compare.MismatchType.REF_ONLY;
import static com.hartwig.hmftools.bamtools.compare.MismatchType.VALUE;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NO_CHROMOSOME_NAME;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SupplementaryReadData.ALIGNMENTS_DELIM;
import static com.hartwig.hmftools.common.samtools.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_3;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.ExcludedRegions;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NO_CHROMOSOME_NAME;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SupplementaryReadData.ALIGNMENTS_DELIM;
import static com.hartwig.hmftools.common.samtools.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_3;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;


public class PartitionReaderTest
{
    private static final ChrBaseRegion REGION = new ChrBaseRegion(CHR_1, 1, 1_000_000);

    private class MockDiffConsumer implements IDiffConsumer
    {
        public final List<DiffRecord> Diffs;

        public MockDiffConsumer()
        {
            Diffs = Lists.newArrayList();
        }

        @Override
        public void addDiff(final DiffRecord diff)
        {
            Diffs.add(diff);
        }
    }

    private class MockSAMRecordProvider implements ISAMRecordProvider
    {
        private final List<SAMRecord> mRecords;

        public MockSAMRecordProvider(final List<SAMRecord> records)
        {
            mRecords = records;
        }

        @Override
        public void processSlice(final BamSlicer bamSlicer, final ChrBaseRegion region, final Consumer<SAMRecord> consumer)
        {
            for (SAMRecord record : mRecords)
            {
                if (bamSlicer.getConsumerHalt())
                    return;

                if (!region.overlapsAlignment(record))
                    continue;

                if (!bamSlicer.passesFilters(record))
                    continue;

                consumer.accept(record);
            }
        }
    }

    @Test
    public void testSingleRefRead()
    {
        // TODO: extract common
        SAMRecord refRecord1 = SamRecordTestUtils.createSamRecord("READ_001", CHR_1, 100, "", "151M", CHR_2, 100, false, false, null, true, "151M");
        List<SAMRecord> refRecords = Lists.newArrayList(refRecord1);

        MockSAMRecordProvider refRecordProver = new MockSAMRecordProvider(refRecords);
        MockSAMRecordProvider newRecordProvider = new MockSAMRecordProvider(Lists.newArrayList());

        CompareConfig config = new CompareConfig(0,false);
        MockDiffConsumer diffs = new MockDiffConsumer();
        PartitionReader reader = new PartitionReader(REGION, config, refRecordProver, newRecordProvider, diffs);
        reader.run();
        Statistics stats = reader.stats();

        assertEquals(1, diffs.Diffs.size());
        DiffRecord diff = diffs.Diffs.get(0);

        assertEquals(refRecord1, diff.Read);
        assertEquals(REF_ONLY, diff.Type);
        assertNull(diff.DiffList);

        assertEquals(1, stats.RefReadCount);
        assertEquals(0, stats.NewReadCount);
        assertEquals(1, stats.DiffCount);
    }

    @Test
    public void testSingleNewRead()
    {
        // TODO: extract common
        SAMRecord newRecord1 = SamRecordTestUtils.createSamRecord("READ_001", CHR_1, 100, "", "151M", CHR_2, 100, false, false, null, true, "151M");
        List<SAMRecord> newRecords = Lists.newArrayList(newRecord1);

        MockSAMRecordProvider refRecordProver = new MockSAMRecordProvider(Lists.newArrayList());
        MockSAMRecordProvider newRecordProvider = new MockSAMRecordProvider(newRecords);

        CompareConfig config = new CompareConfig(0,false);
        MockDiffConsumer diffs = new MockDiffConsumer();
        PartitionReader reader = new PartitionReader(REGION, config, refRecordProver, newRecordProvider, diffs);
        reader.run();
        Statistics stats = reader.stats();

        assertEquals(1, diffs.Diffs.size());
        DiffRecord diff = diffs.Diffs.get(0);

        assertEquals(newRecord1, diff.Read);
        assertEquals(NEW_ONLY, diff.Type);
        assertNull(diff.DiffList);

        assertEquals(0, stats.RefReadCount);
        assertEquals(1, stats.NewReadCount);
        assertEquals(1, stats.DiffCount);
    }

    @Test
    public void testExcludeRefConsensusRead()
    {
        // TODO: extract common
        SAMRecord refRecord1 = SamRecordTestUtils.createSamRecord("READ_001", CHR_1, 100, "", "151M", CHR_2, 100, false, false, null, true, "151M");
        refRecord1.setAttribute(CONSENSUS_READ_ATTRIBUTE, "");
        List<SAMRecord> refRecords = Lists.newArrayList(refRecord1);

        MockSAMRecordProvider refRecordProver = new MockSAMRecordProvider(refRecords);
        MockSAMRecordProvider newRecordProvider = new MockSAMRecordProvider(Lists.newArrayList());

        CompareConfig config = new CompareConfig(0,true);
        MockDiffConsumer diffs = new MockDiffConsumer();
        PartitionReader reader = new PartitionReader(REGION, config, refRecordProver, newRecordProvider, diffs);
        reader.run();
        Statistics stats = reader.stats();

        assertTrue(diffs.Diffs.isEmpty());
        assertEquals(0, stats.RefReadCount);
        assertEquals(0, stats.NewReadCount);
        assertEquals(0, stats.DiffCount);
    }

    @Test
    public void testExcludeNewConsensusRead()
    {
        // TODO: extract common
        SAMRecord newRecord1 = SamRecordTestUtils.createSamRecord("READ_001", CHR_1, 100, "", "151M", CHR_2, 100, false, false, null, true, "151M");
        newRecord1.setAttribute(CONSENSUS_READ_ATTRIBUTE, "");
        List<SAMRecord> newRecords = Lists.newArrayList(newRecord1);

        MockSAMRecordProvider refRecordProver = new MockSAMRecordProvider(Lists.newArrayList());
        MockSAMRecordProvider newRecordProvider = new MockSAMRecordProvider(newRecords);

        CompareConfig config = new CompareConfig(0,true);
        MockDiffConsumer diffs = new MockDiffConsumer();
        PartitionReader reader = new PartitionReader(REGION, config, refRecordProver, newRecordProvider, diffs);
        reader.run();
        Statistics stats = reader.stats();

        assertTrue(diffs.Diffs.isEmpty());
        assertEquals(0, stats.RefReadCount);
        assertEquals(0, stats.NewReadCount);
        assertEquals(0, stats.DiffCount);
    }

    @Test
    public void testSingleMatchingRead()
    {
        // TODO: extract common
        SAMRecord refRead1 = SamRecordTestUtils.createSamRecord("READ_001", CHR_1, 100, "", "151M", CHR_2, 100, false, false, null, true, "151M");
        SAMRecord newRead1 = SamRecordTestUtils.createSamRecord("READ_001", CHR_1, 100, "", "151M", CHR_2, 100, false, false, null, true, "151M");

        MockSAMRecordProvider refRecordProver = new MockSAMRecordProvider(Lists.newArrayList(refRead1));
        MockSAMRecordProvider newRecordProvider = new MockSAMRecordProvider(Lists.newArrayList(newRead1));

        CompareConfig config = new CompareConfig(0,false);
        MockDiffConsumer diffs = new MockDiffConsumer();
        PartitionReader reader = new PartitionReader(REGION, config, refRecordProver, newRecordProvider, diffs);
        reader.run();
        Statistics stats = reader.stats();

        assertTrue(diffs.Diffs.isEmpty());
        assertEquals(1, stats.RefReadCount);
        assertEquals(1, stats.NewReadCount);
        assertEquals(0, stats.DiffCount);
    }

    @Test
    public void testRefReadBeforeNewRead()
    {
        // TODO: extract common
        SAMRecord refRead1 = SamRecordTestUtils.createSamRecord("READ_001", CHR_1, 100, "", "151M", CHR_2, 100, false, false, null, true, "151M");
        SAMRecord newRead1 = SamRecordTestUtils.createSamRecord("READ_001", CHR_1, 101, "", "151M", CHR_2, 100, false, false, null, true, "151M");

        MockSAMRecordProvider refRecordProver = new MockSAMRecordProvider(Lists.newArrayList(refRead1));
        MockSAMRecordProvider newRecordProvider = new MockSAMRecordProvider(Lists.newArrayList(newRead1));

        CompareConfig config = new CompareConfig(0,false);
        MockDiffConsumer diffs = new MockDiffConsumer();
        PartitionReader reader = new PartitionReader(REGION, config, refRecordProver, newRecordProvider, diffs);
        reader.run();
        Statistics stats = reader.stats();

        assertEquals(2, diffs.Diffs.size());

        DiffRecord refDiff = diffs.Diffs.get(0);
        assertEquals(refRead1, refDiff.Read);
        assertEquals(REF_ONLY, refDiff.Type);
        assertNull(refDiff.DiffList);

        DiffRecord newDiff = diffs.Diffs.get(1);
        assertEquals(newRead1, newDiff.Read);
        assertEquals(NEW_ONLY, newDiff.Type);
        assertNull(newDiff.DiffList);

        assertEquals(1, stats.RefReadCount);
        assertEquals(1, stats.NewReadCount);
        assertEquals(2, stats.DiffCount);
    }

    @Test
    public void testNewReadBeforeRefRead()
    {
        // TODO: extract common
        SAMRecord refRead1 = SamRecordTestUtils.createSamRecord("READ_001", CHR_1, 101, "", "151M", CHR_2, 100, false, false, null, true, "151M");
        SAMRecord newRead1 = SamRecordTestUtils.createSamRecord("READ_001", CHR_1, 100, "", "151M", CHR_2, 100, false, false, null, true, "151M");

        MockSAMRecordProvider refRecordProver = new MockSAMRecordProvider(Lists.newArrayList(refRead1));
        MockSAMRecordProvider newRecordProvider = new MockSAMRecordProvider(Lists.newArrayList(newRead1));

        CompareConfig config = new CompareConfig(0,false);
        MockDiffConsumer diffs = new MockDiffConsumer();
        PartitionReader reader = new PartitionReader(REGION, config, refRecordProver, newRecordProvider, diffs);
        reader.run();
        Statistics stats = reader.stats();

        assertEquals(2, diffs.Diffs.size());

        DiffRecord newDiff = diffs.Diffs.get(0);
        assertEquals(newRead1, newDiff.Read);
        assertEquals(NEW_ONLY, newDiff.Type);
        assertNull(newDiff.DiffList);

        DiffRecord refDiff = diffs.Diffs.get(1);
        assertEquals(refRead1, refDiff.Read);
        assertEquals(REF_ONLY, refDiff.Type);
        assertNull(refDiff.DiffList);

        assertEquals(1, stats.RefReadCount);
        assertEquals(1, stats.NewReadCount);
        assertEquals(2, stats.DiffCount);
    }

    @Test
    public void testAlignedMismatch()
    {
        // TODO: extract common
        SAMRecord refRead1 = SamRecordTestUtils.createSamRecord("READ_001", CHR_1, 100, "", "151M", CHR_2, 100, false, false, null, true, "151M");
        SAMRecord newRead1 = SamRecordTestUtils.createSamRecord("READ_002", CHR_1, 100, "", "151M", CHR_2, 100, false, false, null, true, "151M");

        MockSAMRecordProvider refRecordProver = new MockSAMRecordProvider(Lists.newArrayList(refRead1));
        MockSAMRecordProvider newRecordProvider = new MockSAMRecordProvider(Lists.newArrayList(newRead1));

        CompareConfig config = new CompareConfig(0,false);
        MockDiffConsumer diffs = new MockDiffConsumer();
        PartitionReader reader = new PartitionReader(REGION, config, refRecordProver, newRecordProvider, diffs);
        reader.run();
        Statistics stats = reader.stats();

        assertEquals(2, diffs.Diffs.size());

        DiffRecord newDiff = diffs.Diffs.get(0);
        assertEquals(newRead1, newDiff.Read);
        assertEquals(NEW_ONLY, newDiff.Type);
        assertNull(newDiff.DiffList);

        DiffRecord refDiff = diffs.Diffs.get(1);
        assertEquals(refRead1, refDiff.Read);
        assertEquals(REF_ONLY, refDiff.Type);
        assertNull(refDiff.DiffList);

        assertEquals(1, stats.RefReadCount);
        assertEquals(1, stats.NewReadCount);
        assertEquals(2, stats.DiffCount);
    }

    @Test
    public void testFlagMismatch()
    {
        // TODO: extract common
        SAMRecord refRead1 = SamRecordTestUtils.createSamRecord("READ_001", CHR_1, 100, "", "151M", CHR_2, 100, false, false, null, true, "151M");
        SAMRecord newRead1 = SamRecordTestUtils.createSamRecord("READ_001", CHR_1, 100, "", "151M", CHR_2, 100, false, false, null, true, "151M");
        newRead1.setDuplicateReadFlag(true);

        MockSAMRecordProvider refRecordProver = new MockSAMRecordProvider(Lists.newArrayList(refRead1));
        MockSAMRecordProvider newRecordProvider = new MockSAMRecordProvider(Lists.newArrayList(newRead1));

        CompareConfig config = new CompareConfig(0,false);
        MockDiffConsumer diffs = new MockDiffConsumer();
        PartitionReader reader = new PartitionReader(REGION, config, refRecordProver, newRecordProvider, diffs);
        reader.run();
        Statistics stats = reader.stats();

        assertEquals(1, diffs.Diffs.size());
        DiffRecord newDiff = diffs.Diffs.get(0);
        assertEquals(refRead1, newDiff.Read);
        assertEquals(VALUE, newDiff.Type);
        assertEquals(1, newDiff.DiffList.size());
        assertEquals("duplicate(false/true)", newDiff.DiffList.get(0));

        assertEquals(1, stats.RefReadCount);
        assertEquals(1, stats.NewReadCount);
        assertEquals(1, stats.DiffCount);
    }

    @Test
    public void testExcludedRegion()
    {
        ChrBaseRegion excludeRegion = ExcludedRegions.getPolyGRegion(V37);

        // TODO: extract common
        SAMRecord refRead1 = SamRecordTestUtils.createSamRecord("READ_001", excludeRegion.Chromosome, excludeRegion.start(), "", "151M", CHR_2, 100, false, false, null, true, "151M");
        SAMRecord newRead1 = SamRecordTestUtils.createSamRecord("READ_001", excludeRegion.Chromosome, excludeRegion.start(), "", "151M", CHR_2, 100, false, false, null, true, "151M");

        MockSAMRecordProvider refRecordProver = new MockSAMRecordProvider(Lists.newArrayList(refRead1));
        MockSAMRecordProvider newRecordProvider = new MockSAMRecordProvider(Lists.newArrayList(newRead1));

        CompareConfig config = new CompareConfig(0,false);
        MockDiffConsumer diffs = new MockDiffConsumer();
        PartitionReader reader = new PartitionReader(excludeRegion, config, refRecordProver, newRecordProvider, diffs);
        reader.run();
        Statistics stats = reader.stats();

        assertEquals(0, diffs.Diffs.size());

        assertEquals(0, stats.RefReadCount);
        assertEquals(0, stats.NewReadCount);
        assertEquals(0, stats.DiffCount);
    }

    @Test
    public void testHaltRefProcessing()
    {
        // TODO: extract common
        SAMRecord refRead1 = SamRecordTestUtils.createSamRecord("READ_001", CHR_1, 100, "", "151M", CHR_2, 100, false, false, null, true, "151M");
        SAMRecord refRead2 = SamRecordTestUtils.createSamRecord("READ_001", CHR_1, 101, "", "151M", CHR_2, 100, false, false, null, true, "151M");
        SAMRecord newRead1 = SamRecordTestUtils.createSamRecord("READ_001", CHR_1, 100, "", "151M", CHR_2, 100, false, false, null, true, "151M");

        MockSAMRecordProvider refRecordProver = new MockSAMRecordProvider(Lists.newArrayList(refRead1, refRead2));
        MockSAMRecordProvider newRecordProvider = new MockSAMRecordProvider(Lists.newArrayList(newRead1));

        CompareConfig config = new CompareConfig(1,false);
        MockDiffConsumer diffs = new MockDiffConsumer();
        PartitionReader reader = new PartitionReader(REGION, config, refRecordProver, newRecordProvider, diffs);
        reader.run();
        Statistics stats = reader.stats();

        assertEquals(1, diffs.Diffs.size());
        DiffRecord refDiff = diffs.Diffs.get(0);
        assertEquals(refRead1, refDiff.Read);
        assertEquals(REF_ONLY, refDiff.Type);
        assertNull(refDiff.DiffList);

        assertEquals(1, stats.RefReadCount);
        assertEquals(0, stats.NewReadCount);
        assertEquals(1, stats.DiffCount);
    }

    @Test
    public void testHaltNewProcessing()
    {
        // TODO: extract common
        SAMRecord newRead1 = SamRecordTestUtils.createSamRecord("READ_001", CHR_1, 100, "", "151M", CHR_2, 100, false, false, null, true, "151M");
        SAMRecord newRead2 = SamRecordTestUtils.createSamRecord("READ_001", CHR_1, 101, "", "151M", CHR_2, 100, false, false, null, true, "151M");

        MockSAMRecordProvider refRecordProver = new MockSAMRecordProvider(Lists.newArrayList());
        MockSAMRecordProvider newRecordProvider = new MockSAMRecordProvider(Lists.newArrayList(newRead1, newRead2));

        CompareConfig config = new CompareConfig(1,false);
        MockDiffConsumer diffs = new MockDiffConsumer();
        PartitionReader reader = new PartitionReader(REGION, config, refRecordProver, newRecordProvider, diffs);
        reader.run();
        Statistics stats = reader.stats();

        assertEquals(1, diffs.Diffs.size());
        DiffRecord newDiff = diffs.Diffs.get(0);
        assertEquals(newRead1, newDiff.Read);
        assertEquals(NEW_ONLY, newDiff.Type);
        assertNull(newDiff.DiffList);

        assertEquals(0, stats.RefReadCount);
        assertEquals(1, stats.NewReadCount);
        assertEquals(1, stats.DiffCount);
    }

//    private static final String READ_BASES = "A".repeat(100);
//    private static final String READ_CIGAR = "100M";
//    private static final String CHR_4 = "4";
//    private static final String CHR_5 = "5";
//    private static final String SOFT_CLIPPED_READ_CIGAR = "21S79M";
//    private static final String READ_ID = "READ_001";
//    private static final Map<String, List<HighDepthRegion>> CHR_LOCATION_MAP;
//    private static final ReadUnmapper READ_UNMAPPER;
//
//    static
//    {
//        CHR_LOCATION_MAP = Maps.newHashMap();
//        CHR_LOCATION_MAP.put(CHR_1, Lists.newArrayList(new HighDepthRegion(500, 700, 0)));
//        CHR_LOCATION_MAP.put(CHR_2, Lists.newArrayList());
//        CHR_LOCATION_MAP.put(CHR_3, Lists.newArrayList(new HighDepthRegion(500, 700, UNMAP_MIN_HIGH_DEPTH)));
//
//        CHR_LOCATION_MAP.put(CHR_4, Lists.newArrayList(
//                new HighDepthRegion(1000, 2000, UNMAP_MIN_HIGH_DEPTH),
//                new HighDepthRegion(3000, 4000, UNMAP_MIN_HIGH_DEPTH),
//                new HighDepthRegion(5000, 6000, UNMAP_MIN_HIGH_DEPTH),
//                new HighDepthRegion(7000, 8000, UNMAP_MIN_HIGH_DEPTH)));
//        READ_UNMAPPER = new ReadUnmapper(CHR_LOCATION_MAP);
//    }
//
//    @Test
//    public void testGetSoftClipCountFromCigarStr()
//    {
//        assertEquals(0, ReadUnmapper.getSoftClipCountFromCigarStr("151M"));
//        assertEquals(10, ReadUnmapper.getSoftClipCountFromCigarStr("10S141M"));
//        assertEquals(10, ReadUnmapper.getSoftClipCountFromCigarStr("141M10S"));
//        assertEquals(12, ReadUnmapper.getSoftClipCountFromCigarStr("2S138M10S"));
//    }
//
//    @Test
//    public void testNoUnmap()
//    {
//        // neither read in unmapp regions
//        SAMRecord read = SamRecordTestUtils.createSamRecord(
//                READ_ID, CHR_1, 100, READ_BASES, READ_CIGAR, CHR_2, 600, false,
//                false, null, true, READ_CIGAR);
//
//        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
//        assertFalse(checkTransformRead(read, CHR_1));
//        assertFalse(read.getReadUnmappedFlag());
//        assertFalse(read.getMateUnmappedFlag());
//        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
//    }
//
//    @Test
//    public void testNoUnmapSoftClipped()
//    {
//        // soft-clipped but not in an unmap region
//        SAMRecord read = SamRecordTestUtils.createSamRecord(
//                READ_ID, CHR_1, 100, READ_BASES, SOFT_CLIPPED_READ_CIGAR, CHR_2, 600, false,
//                false, null, true, SOFT_CLIPPED_READ_CIGAR);
//
//        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
//        assertFalse(checkTransformRead(read, CHR_1));
//        assertFalse(read.getReadUnmappedFlag());
//        assertFalse(read.getMateUnmappedFlag());
//        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
//    }
//
//    @Test
//    public void testUnmapReadSoftClipped()
//    {
//        // read unmapped since chimeric, sets read coords to mate
//        SAMRecord read = SamRecordTestUtils.createSamRecord(
//                READ_ID, CHR_1, 550, READ_BASES, SOFT_CLIPPED_READ_CIGAR, CHR_2, 600, false,
//                false, null, true, READ_CIGAR);
//
//        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
//        assertTrue(checkTransformRead(read, CHR_1));
//        assertTrue(read.getReadUnmappedFlag());
//        assertFalse(read.getMateUnmappedFlag());
//        assertEquals(600, read.getAlignmentStart());
//        assertEquals(CHR_2, read.getReferenceName());
//    }
//
//    @Test
//    public void testUnmapMateSoftClipped()
//    {
//        // read unmapped since soft-cliped
//        SAMRecord read = SamRecordTestUtils.createSamRecord(
//                READ_ID, CHR_1, 100, READ_BASES, READ_CIGAR, CHR_1, 600, false,
//                false, null, true, SOFT_CLIPPED_READ_CIGAR);
//
//        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
//        assertTrue(checkTransformRead(read, CHR_1));
//        assertFalse(read.getReadUnmappedFlag());
//        assertTrue(read.getMateUnmappedFlag());
//        assertEquals(100, read.getMateAlignmentStart());
//        assertEquals(CHR_1, read.getMateReferenceName());
//        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
//    }
//
//    @Test
//    public void testUnmapMateAndSupplementarySoftClipped()
//    {
//        // mate unmapped and supp removed since both are in unmap regions and both have long SCs
//        SAMRecord read = SamRecordTestUtils.createSamRecord(
//                READ_ID, CHR_1, 100, READ_BASES, READ_CIGAR, CHR_1, 600, false,
//                false,
//                new SupplementaryReadData(CHR_1, 550, SUPP_POS_STRAND, SOFT_CLIPPED_READ_CIGAR, 60),
//                true, SOFT_CLIPPED_READ_CIGAR);
//
//        assertTrue(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
//        assertTrue(checkTransformRead(read, CHR_1));
//        assertFalse(read.getReadUnmappedFlag());
//        assertTrue(read.getMateUnmappedFlag());
//        assertEquals(100, read.getMateAlignmentStart());
//        assertEquals(CHR_1, read.getMateReferenceName());
//        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
//    }
//
//    @Test
//    public void testUnmapReadSoftClippedMateUnmapped()
//    {
//        SAMRecord read = SamRecordTestUtils.createSamRecord(
//                READ_ID, CHR_1, 550, READ_BASES, SOFT_CLIPPED_READ_CIGAR, CHR_1, 550, false,
//                false, null, true, null);
//        read.setMateUnmappedFlag(true);
//
//        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
//        assertTrue(checkTransformRead(read, CHR_1));
//        assertTrue(read.getReadUnmappedFlag());
//        assertTrue(read.getMateUnmappedFlag());
//        assertEquals(NO_POSITION, read.getAlignmentStart());
//        assertEquals(NO_CHROMOSOME_NAME, read.getReferenceName());
//        assertEquals(NO_POSITION, read.getMateAlignmentStart());
//        assertEquals(NO_CHROMOSOME_NAME, read.getMateReferenceName());
//        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
//    }
//
//    @Test
//    public void testUnmapMateSoftClippedReadUnmapped()
//    {
//        SAMRecord read = SamRecordTestUtils.createSamRecord(
//                READ_ID, CHR_1, 550, READ_BASES, NO_CIGAR, CHR_1, 550, false,
//                false, null, true, SOFT_CLIPPED_READ_CIGAR);
//        read.setReadUnmappedFlag(true);
//
//        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
//        assertTrue(checkTransformRead(read, CHR_1));
//        assertTrue(read.getReadUnmappedFlag());
//        assertTrue(read.getMateUnmappedFlag());
//        assertEquals(NO_POSITION, read.getAlignmentStart());
//        assertEquals(NO_CHROMOSOME_NAME, read.getReferenceName());
//        assertEquals(NO_POSITION, read.getMateAlignmentStart());
//        assertEquals(NO_CHROMOSOME_NAME, read.getMateReferenceName());
//        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
//    }
//
//    @Test
//    public void testUnmapSupplementaryWhenPrimaryIsUnmappedSoftClipped()
//    {
//        // Note that the suppReadData for a supplementary read refers to the associated primary read
//        SupplementaryReadData suppReadData = new SupplementaryReadData(
//                CHR_1, 550, SUPP_POS_STRAND, SOFT_CLIPPED_READ_CIGAR, 60);
//
//        SAMRecord read = SamRecordTestUtils.createSamRecord(
//                READ_ID, CHR_2, 100, READ_BASES, READ_CIGAR, CHR_2, 1000, false,
//                true, suppReadData, true, READ_CIGAR);
//
//        // when supplementaries are unmapped not all properties are unset since it will be dropped immediately
//        assertTrue(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
//        assertTrue(checkTransformRead(read, CHR_2));
//        assertTrue(read.getReadUnmappedFlag());
//    }
//
//    @Test
//    public void testUnmapReadChimeric()
//    {
//        SAMRecord read = SamRecordTestUtils.createSamRecord(
//                READ_ID, CHR_1, 550, READ_BASES, READ_CIGAR, CHR_2, 600, false,
//                false, null, true, READ_CIGAR);
//
//        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
//        assertTrue(checkTransformRead(read, CHR_1));
//        assertTrue(read.getReadUnmappedFlag());
//        assertFalse(read.getMateUnmappedFlag());
//        assertEquals(600, read.getAlignmentStart());
//        assertEquals(CHR_2, read.getReferenceName());
//    }
//
//    @Test
//    public void testUnmapMateChimeric()
//    {
//        SAMRecord read = SamRecordTestUtils.createSamRecord(
//                READ_ID, CHR_2, 100, READ_BASES, READ_CIGAR, CHR_1, 600, false,
//                false, null, true, READ_BASES);
//
//        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
//        assertTrue(checkTransformRead(read, CHR_1));
//        assertFalse(read.getReadUnmappedFlag());
//        assertTrue(read.getMateUnmappedFlag());
//        assertEquals(100, read.getMateAlignmentStart());
//        assertEquals(CHR_2, read.getMateReferenceName());
//    }
//
//    @Test
//    public void testSupplementaryNotUnmappedChimeric()
//    {
//        // Note that the suppReadData for a supplementary read refers to the associated primary read.
//        SupplementaryReadData suppReadData = new SupplementaryReadData(CHR_2, 100, SUPP_POS_STRAND, READ_CIGAR, 60);
//
//        SAMRecord read = SamRecordTestUtils.createSamRecord(
//                READ_ID, CHR_1, 100, READ_BASES, READ_CIGAR, CHR_2, 600, false, true,
//                suppReadData, true, READ_BASES);
//
//        assertTrue(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
//        assertFalse(checkTransformRead(read, CHR_1));
//        assertFalse(read.getReadUnmappedFlag());
//        assertFalse(read.getMateUnmappedFlag());
//        assertTrue(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
//    }
//
//    @Test
//    public void testUnmapSupplementaryWhenPrimaryIsUnmappedChimeric()
//    {
//        // Note that the suppReadData for a supplementary read refers to the associated primary read.
//        SupplementaryReadData suppReadData = new SupplementaryReadData(CHR_1, 550, SUPP_POS_STRAND, READ_CIGAR, 60);
//
//        SAMRecord read = SamRecordTestUtils.createSamRecord(
//                READ_ID, CHR_2, 100, READ_BASES, READ_CIGAR, CHR_2, 1000, false,
//                true, suppReadData, true, READ_CIGAR);
//
//        // Note that when we are unmapping a supplementary we do not bother unsetting all of its properties, because we will drop an unmapped
//        // supplementary read immediately.
//        assertTrue(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
//        assertTrue(checkTransformRead(read, CHR_2));
//        assertTrue(read.getReadUnmappedFlag());
//    }
//
//    @Test
//    public void testUnmapReadHighDepth()
//    {
//        SAMRecord read = SamRecordTestUtils.createSamRecord(
//                READ_ID, CHR_3, 600, READ_BASES, READ_CIGAR, CHR_3, 701, false,
//                false, null, true, READ_CIGAR);
//
//        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
//        assertTrue(checkTransformRead(read, CHR_3));
//        assertTrue(read.getReadUnmappedFlag());
//        assertFalse(read.getMateUnmappedFlag());
//        assertEquals(701, read.getAlignmentStart());
//        assertEquals(CHR_3, read.getReferenceName());
//    }
//
//    @Test
//    public void testUnmapMateHighDepth()
//    {
//        SAMRecord read = SamRecordTestUtils.createSamRecord(
//                READ_ID, CHR_3, 400, READ_BASES, READ_CIGAR, CHR_3, 500, false,
//                false, null, true, READ_CIGAR);
//
//        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
//        assertTrue(checkTransformRead(read, CHR_3));
//        assertFalse(read.getReadUnmappedFlag());
//        assertTrue(read.getMateUnmappedFlag());
//        assertEquals(400, read.getMateAlignmentStart());
//        assertEquals(CHR_3, read.getMateReferenceName());
//    }
//
//    @Test
//    public void testUnmapMateAndSupplementaryHighDepth()
//    {
//        SAMRecord read = SamRecordTestUtils.createSamRecord(
//                READ_ID, CHR_3, 400, READ_BASES, READ_CIGAR, CHR_3, 500, false,
//                false,
//                new SupplementaryReadData(CHR_3, 600, SUPP_POS_STRAND, READ_CIGAR, 60),
//                true, READ_CIGAR);
//
//        assertTrue(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
//        assertTrue(checkTransformRead(read, CHR_3));
//        assertFalse(read.getReadUnmappedFlag());
//        assertTrue(read.getMateUnmappedFlag());
//        assertEquals(400, read.getMateAlignmentStart());
//        assertEquals(CHR_3, read.getMateReferenceName());
//        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
//    }
//
//    @Test
//    public void testUnmapMultipleSupplementaryAlignments()
//    {
//        SAMRecord read = SamRecordTestUtils.createSamRecord(
//                READ_ID, CHR_1, 100, READ_BASES, READ_CIGAR, CHR_3, 200, false,
//                false, null, true, READ_CIGAR);
//
//
//        SupplementaryReadData suppData1 = new SupplementaryReadData(CHR_1, 1000, SUPP_POS_STRAND, READ_CIGAR, 60);
//        SupplementaryReadData suppData2 = new SupplementaryReadData(CHR_2, 6000, SUPP_POS_STRAND, READ_CIGAR, 60);
//        SupplementaryReadData suppData3 = new SupplementaryReadData(CHR_3, 600, SUPP_POS_STRAND, READ_CIGAR, 60);
//
//        StringJoiner saJoiner = new StringJoiner(ALIGNMENTS_DELIM);
//        saJoiner.add(suppData1.asSamTag());
//        saJoiner.add(suppData2.asSamTag());
//        saJoiner.add(suppData3.asSamTag());
//
//        read.setAttribute(SUPPLEMENTARY_ATTRIBUTE, saJoiner.toString());
//        assertTrue(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
//
//        assertTrue(checkTransformRead(read, CHR_1));
//        assertFalse(read.getReadUnmappedFlag());
//        assertFalse(read.getMateUnmappedFlag());
//        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
//    }
//
//    @Test
//    public void testUnmapSupplementaryWhenUnmappingAssociatedSupplementary()
//    {
//        SAMRecord read = SamRecordTestUtils.createSamRecord(
//                READ_ID, CHR_1, 100, READ_BASES, READ_CIGAR, CHR_3, 300, false,
//                true, null, true, READ_CIGAR);
//
//        SupplementaryReadData suppDataPrimary = new SupplementaryReadData(CHR_3, 100, SUPP_POS_STRAND, READ_CIGAR, 60);
//        SupplementaryReadData suppDataOther = new SupplementaryReadData(CHR_3, 500, SUPP_POS_STRAND, READ_CIGAR, 60);
//
//        StringJoiner saJoiner = new StringJoiner(ALIGNMENTS_DELIM);
//        saJoiner.add(suppDataPrimary.asSamTag());
//        saJoiner.add(suppDataOther.asSamTag());
//
//        read.setAttribute(SUPPLEMENTARY_ATTRIBUTE, saJoiner.toString());
//        assertTrue(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
//
//        // Note that when we are unmapping a supplementary we do not bother unsetting all of its properties, because we will drop an
//        // unmapped supplementary read immediately.
//        assertTrue(checkTransformRead(read, CHR_3));
//        assertTrue(read.getReadUnmappedFlag());
//    }
//
//    @Test
//    public void testUnmapRegionState()
//    {
//        UnmapRegionState regionState = new UnmapRegionState(new ChrBaseRegion(
//                CHR_4, 1, 1000000), CHR_LOCATION_MAP.get(CHR_4));
//
//        // primary read
//        SAMRecord read = SamRecordTestUtils.createSamRecord(
//                READ_ID, CHR_4, 1000, READ_BASES, READ_CIGAR, CHR_5, 100, false,
//                false, null, true, READ_CIGAR);
//
//        assertTrue(READ_UNMAPPER.checkTransformRead(read, regionState));
//        assertTrue(regionState.LastMatchedRegionIndex != null);
//        assertEquals(0, regionState.LastMatchedRegionIndex.intValue());
//
//        read = SamRecordTestUtils.createSamRecord(
//                READ_ID, CHR_4, 2100, READ_BASES, READ_CIGAR, CHR_5, 100, false,
//                false, null, true, READ_CIGAR);
//
//        assertFalse(READ_UNMAPPER.checkTransformRead(read, regionState));
//        assertEquals(0, regionState.LastMatchedRegionIndex.intValue());
//
//        read = SamRecordTestUtils.createSamRecord(
//                READ_ID, CHR_4, 5100, READ_BASES, READ_CIGAR, CHR_5, 100, false,
//                false, null, true, READ_CIGAR);
//
//        assertTrue(READ_UNMAPPER.checkTransformRead(read, regionState));
//        assertEquals(2, regionState.LastMatchedRegionIndex.intValue());
//
//        read = SamRecordTestUtils.createSamRecord(
//                READ_ID, CHR_4, 10000, READ_BASES, READ_CIGAR, CHR_5, 100, false,
//                false, null, true, READ_CIGAR);
//
//        assertFalse(READ_UNMAPPER.checkTransformRead(read, regionState));
//        assertEquals(3, regionState.LastMatchedRegionIndex.intValue());
//    }
//
//    @Test
//    public void testNoUnmapUnpaired()
//    {
//        // read not in unmap regions
//        SAMRecord read = SamRecordTestUtils.createSamRecordUnpaired(READ_ID, CHR_1, 100, READ_BASES, READ_CIGAR, false, false, null);
//
//        assertFalse(checkTransformRead(read, CHR_1));
//        assertFalse(read.getReadUnmappedFlag());
//
//    }
//
//    @Test
//    public void testNoUnmapSoftClippedUnpaired()
//    {
//        // soft-clipped but not in an unmap region
//        SAMRecord read = SamRecordTestUtils.createSamRecordUnpaired(
//                READ_ID, CHR_1, 100, READ_BASES, SOFT_CLIPPED_READ_CIGAR, false, false, null);
//
//        assertFalse(checkTransformRead(read, CHR_1));
//        assertFalse(read.getReadUnmappedFlag());
//    }
//
//    @Test
//    public void testUnmapReadSoftClippedUnpaired()
//    {
//        SAMRecord read = SamRecordTestUtils.createSamRecordUnpaired(
//                READ_ID, CHR_1, 550, READ_BASES, SOFT_CLIPPED_READ_CIGAR, false, false, null);
//
//        assertTrue(checkTransformRead(read, CHR_1));
//        assertTrue(read.getReadUnmappedFlag());
//        assertEquals(NO_POSITION, read.getAlignmentStart());
//        assertEquals(NO_CHROMOSOME_NAME, read.getReferenceName());
//    }
//
//    @Test
//    public void testSupplementarySoftClippedUnpaired()
//    {
//        // supp removed since it is in an unmap region and has a long SC
//        SAMRecord read = SamRecordTestUtils.createSamRecordUnpaired(
//                READ_ID, CHR_1, 100, READ_BASES, READ_CIGAR, false,
//                false,
//                new SupplementaryReadData(CHR_1, 550, SUPP_POS_STRAND, SOFT_CLIPPED_READ_CIGAR, 60));
//
//        assertTrue(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
//        assertTrue(checkTransformRead(read, CHR_1));
//        assertFalse(read.getReadUnmappedFlag());
//        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
//    }
//
//    @Test
//    public void testUnmapSupplementaryWhenPrimaryIsUnmappedSoftClippedUnpaired()
//    {
//        // Note that the suppReadData for a supplementary read refers to the associated primary read
//        SupplementaryReadData suppReadData = new SupplementaryReadData(
//                CHR_1, 550, SUPP_POS_STRAND, SOFT_CLIPPED_READ_CIGAR, 60);
//
//        SAMRecord read = SamRecordTestUtils.createSamRecordUnpaired(
//                READ_ID, CHR_2, 100, READ_BASES, READ_CIGAR, false,
//                true, suppReadData);
//
//        // when supplementaries are unmapped not all properties are unset since it will be dropped immediately
//        assertTrue(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
//        assertTrue(checkTransformRead(read, CHR_2));
//        assertTrue(read.getReadUnmappedFlag());
//    }
//
//    @Test
//    public void testUnpairedChimericNotUnmapped()
//    {
//        SAMRecord read = SamRecordTestUtils.createSamRecord(
//                READ_ID, CHR_1, 550, READ_BASES, READ_CIGAR, CHR_2, 600, false,
//                false, null, true, READ_CIGAR);
//        read.setReadPairedFlag(false);
//
//        assertFalse(checkTransformRead(read, CHR_1));
//        assertFalse(read.getReadUnmappedFlag());
//    }
//
//    @Test
//    public void testUnmapReadHighDepthUnpaired()
//    {
//        SAMRecord read = SamRecordTestUtils.createSamRecordUnpaired(READ_ID, CHR_3, 600, READ_BASES, READ_CIGAR, false, false, null);
//
//        assertTrue(checkTransformRead(read, CHR_3));
//        assertTrue(read.getReadUnmappedFlag());
//        assertEquals(NO_POSITION, read.getAlignmentStart());
//        assertEquals(NO_CHROMOSOME_NAME, read.getReferenceName());
//    }
//
//    @Test
//    public void testUnmapSupplementaryHighDepthUnpaired()
//    {
//        SAMRecord read = SamRecordTestUtils.createSamRecordUnpaired(
//                READ_ID, CHR_3, 400, READ_BASES, READ_CIGAR, false,
//                false,
//                new SupplementaryReadData(CHR_3, 600, SUPP_POS_STRAND, READ_CIGAR, 60));
//
//        assertTrue(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
//        assertTrue(checkTransformRead(read, CHR_3));
//        assertFalse(read.getReadUnmappedFlag());
//        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
//    }
//
//    @Test
//    public void testUnmapMultipleSupplementaryAlignmentsUnpaired()
//    {
//        SAMRecord read = SamRecordTestUtils.createSamRecordUnpaired(
//                READ_ID, CHR_1, 100, READ_BASES, READ_CIGAR, false,
//                false, null);
//
//        SupplementaryReadData suppData1 = new SupplementaryReadData(CHR_1, 1000, SUPP_POS_STRAND, READ_CIGAR, 60);
//        SupplementaryReadData suppData2 = new SupplementaryReadData(CHR_2, 6000, SUPP_POS_STRAND, READ_CIGAR, 60);
//        SupplementaryReadData suppData3 = new SupplementaryReadData(CHR_3, 600, SUPP_POS_STRAND, READ_CIGAR, 60);
//
//        StringJoiner saJoiner = new StringJoiner(ALIGNMENTS_DELIM);
//        saJoiner.add(suppData1.asSamTag());
//        saJoiner.add(suppData2.asSamTag());
//        saJoiner.add(suppData3.asSamTag());
//
//        read.setAttribute(SUPPLEMENTARY_ATTRIBUTE, saJoiner.toString());
//
//        assertTrue(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
//        assertTrue(checkTransformRead(read, CHR_1));
//        assertFalse(read.getReadUnmappedFlag());
//        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
//    }
//
//    @Test
//    public void testUnmapSupplementaryWhenUnmappingAssociatedSupplementaryUnpaired()
//    {
//        SAMRecord read = SamRecordTestUtils.createSamRecordUnpaired(
//                READ_ID, CHR_1, 100, READ_BASES, READ_CIGAR, false,
//                true, null);
//
//        SupplementaryReadData suppDataPrimary = new SupplementaryReadData(CHR_3, 100, SUPP_POS_STRAND, READ_CIGAR, 60);
//        SupplementaryReadData suppDataOther = new SupplementaryReadData(CHR_3, 500, SUPP_POS_STRAND, READ_CIGAR, 60);
//
//        StringJoiner saJoiner = new StringJoiner(ALIGNMENTS_DELIM);
//        saJoiner.add(suppDataPrimary.asSamTag());
//        saJoiner.add(suppDataOther.asSamTag());
//
//        read.setAttribute(SUPPLEMENTARY_ATTRIBUTE, saJoiner.toString());
//        assertTrue(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
//
//        // Note that when we are unmapping a supplementary we do not bother unsetting all of its properties, because we will drop an
//        // unmapped supplementary read immediately.
//        assertTrue(checkTransformRead(read, CHR_3));
//        assertTrue(read.getReadUnmappedFlag());
//    }
//
//    @Test
//    public void testUnmapSupplementaryWhenPrimaryIsUnmappedUnpaired()
//    {
//        // Note that the suppReadData for a supplementary read refers to the associated primary read.
//        SupplementaryReadData suppReadData = new SupplementaryReadData(CHR_3, 550, SUPP_POS_STRAND, READ_CIGAR, 60);
//
//        SAMRecord read = SamRecordTestUtils.createSamRecordUnpaired(READ_ID, CHR_2, 100, READ_BASES, READ_CIGAR, false, true, suppReadData);
//
//        // Note that when we are unmapping a supplementary we do not bother unsetting all of its properties, because we will drop an unmapped
//        // supplementary read immediately.
//        assertTrue(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
//        assertTrue(checkTransformRead(read, CHR_2));
//        assertTrue(read.getReadUnmappedFlag());
//    }
//
//    private boolean checkTransformRead(final SAMRecord read, final String chromosome)
//    {
//        UnmapRegionState regionState = new UnmapRegionState(new ChrBaseRegion(
//                chromosome, 1, 1000000), CHR_LOCATION_MAP.get(chromosome));
//
//        return READ_UNMAPPER.checkTransformRead(read, regionState);
//    }
}
