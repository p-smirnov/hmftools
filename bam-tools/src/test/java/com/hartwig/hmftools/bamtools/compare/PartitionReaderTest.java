package com.hartwig.hmftools.bamtools.compare;

import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.compare.MismatchType.NEW_ONLY;
import static com.hartwig.hmftools.bamtools.compare.MismatchType.REF_ONLY;
import static com.hartwig.hmftools.bamtools.compare.MismatchType.VALUE;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
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
import static org.junit.Assert.assertNotNull;
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
import com.hartwig.hmftools.common.samtools.SamRecordUtils;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;

import org.jetbrains.annotations.Nullable;
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
        SAMRecord refRecord1 = createSamRecord("READ_001", CHR_1, 100);
        List<SAMRecord> refRecords = Lists.newArrayList(refRecord1);

        MockSAMRecordProvider refRecordProver = new MockSAMRecordProvider(refRecords);
        MockSAMRecordProvider newRecordProvider = new MockSAMRecordProvider(Lists.newArrayList());

        CompareConfig config = new CompareConfig(0,false, false);
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
        SAMRecord newRecord1 = createSamRecord("READ_001", CHR_1, 100);
        List<SAMRecord> newRecords = Lists.newArrayList(newRecord1);

        MockSAMRecordProvider refRecordProver = new MockSAMRecordProvider(Lists.newArrayList());
        MockSAMRecordProvider newRecordProvider = new MockSAMRecordProvider(newRecords);

        CompareConfig config = new CompareConfig(0, false, false);
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
        SAMRecord refRecord1 = createSamRecord("READ_001", CHR_1, 100, false, true);
        List<SAMRecord> refRecords = Lists.newArrayList(refRecord1);

        MockSAMRecordProvider refRecordProver = new MockSAMRecordProvider(refRecords);
        MockSAMRecordProvider newRecordProvider = new MockSAMRecordProvider(Lists.newArrayList());

        CompareConfig config = new CompareConfig(0, false, true);
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
        SAMRecord newRecord1 = createSamRecord("READ_001", CHR_1, 100, false, true);
        List<SAMRecord> newRecords = Lists.newArrayList(newRecord1);

        MockSAMRecordProvider refRecordProver = new MockSAMRecordProvider(Lists.newArrayList());
        MockSAMRecordProvider newRecordProvider = new MockSAMRecordProvider(newRecords);

        CompareConfig config = new CompareConfig(0, false, true);
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
        SAMRecord refRead1 = createSamRecord("READ_001", CHR_1, 100);
        SAMRecord newRead1 = createSamRecord("READ_001", CHR_1, 100);

        MockSAMRecordProvider refRecordProver = new MockSAMRecordProvider(Lists.newArrayList(refRead1));
        MockSAMRecordProvider newRecordProvider = new MockSAMRecordProvider(Lists.newArrayList(newRead1));

        CompareConfig config = new CompareConfig(0, false, false);
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
        SAMRecord refRead1 = createSamRecord("READ_001", CHR_1, 100);
        SAMRecord newRead1 = createSamRecord("READ_001", CHR_1, 101);

        MockSAMRecordProvider refRecordProver = new MockSAMRecordProvider(Lists.newArrayList(refRead1));
        MockSAMRecordProvider newRecordProvider = new MockSAMRecordProvider(Lists.newArrayList(newRead1));

        CompareConfig config = new CompareConfig(0, false, false);
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
        SAMRecord refRead1 = createSamRecord("READ_001", CHR_1, 101);
        SAMRecord newRead1 = createSamRecord("READ_001", CHR_1, 100);

        MockSAMRecordProvider refRecordProver = new MockSAMRecordProvider(Lists.newArrayList(refRead1));
        MockSAMRecordProvider newRecordProvider = new MockSAMRecordProvider(Lists.newArrayList(newRead1));

        CompareConfig config = new CompareConfig(0, false, false);
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
        SAMRecord refRead1 = createSamRecord("READ_001", CHR_1, 100);
        SAMRecord newRead1 = createSamRecord("READ_002", CHR_1, 100);

        MockSAMRecordProvider refRecordProver = new MockSAMRecordProvider(Lists.newArrayList(refRead1));
        MockSAMRecordProvider newRecordProvider = new MockSAMRecordProvider(Lists.newArrayList(newRead1));

        CompareConfig config = new CompareConfig(0, false, false);
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
        SAMRecord refRead1 = createSamRecord("READ_001", CHR_1, 100);
        SAMRecord newRead1 = createSamRecord("READ_001", CHR_1, 100, true, false);

        MockSAMRecordProvider refRecordProver = new MockSAMRecordProvider(Lists.newArrayList(refRead1));
        MockSAMRecordProvider newRecordProvider = new MockSAMRecordProvider(Lists.newArrayList(newRead1));

        CompareConfig config = new CompareConfig(0, false, false);
        MockDiffConsumer diffs = new MockDiffConsumer();
        PartitionReader reader = new PartitionReader(REGION, config, refRecordProver, newRecordProvider, diffs);
        reader.run();
        Statistics stats = reader.stats();

        assertEquals(1, diffs.Diffs.size());
        DiffRecord diff = diffs.Diffs.get(0);
        assertEquals(refRead1, diff.Read);
        assertEquals(VALUE, diff.Type);
        assertEquals(1, diff.DiffList.size());
        assertEquals("duplicate(false/true)", diff.DiffList.get(0));

        assertEquals(1, stats.RefReadCount);
        assertEquals(1, stats.NewReadCount);
        assertEquals(1, stats.DiffCount);
    }

    @Test
    public void testIgnoreDuplicateFlagMismatch()
    {
        // TODO: extract common
        SAMRecord refRead1 = createSamRecord("READ_001", CHR_1, 100);
        SAMRecord newRead1 = createSamRecord("READ_001", CHR_1, 100, true, false);

        MockSAMRecordProvider refRecordProver = new MockSAMRecordProvider(Lists.newArrayList(refRead1));
        MockSAMRecordProvider newRecordProvider = new MockSAMRecordProvider(Lists.newArrayList(newRead1));

        CompareConfig config = new CompareConfig(0, true, false);
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
    public void testExcludedRegion()
    {
        ChrBaseRegion excludeRegion = ExcludedRegions.getPolyGRegion(V37);

        // TODO: extract common
        SAMRecord refRead1 = createSamRecord("READ_001", excludeRegion.Chromosome, excludeRegion.start());
        SAMRecord newRead1 = createSamRecord("READ_001", excludeRegion.Chromosome, excludeRegion.start());

        MockSAMRecordProvider refRecordProver = new MockSAMRecordProvider(Lists.newArrayList(refRead1));
        MockSAMRecordProvider newRecordProvider = new MockSAMRecordProvider(Lists.newArrayList(newRead1));

        CompareConfig config = new CompareConfig(0, false, false);
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
        SAMRecord refRead1 = createSamRecord("READ_001", CHR_1, 100);
        SAMRecord refRead2 = createSamRecord("READ_001", CHR_1, 101);
        SAMRecord newRead1 = createSamRecord("READ_001", CHR_1, 100);

        MockSAMRecordProvider refRecordProver = new MockSAMRecordProvider(Lists.newArrayList(refRead1, refRead2));
        MockSAMRecordProvider newRecordProvider = new MockSAMRecordProvider(Lists.newArrayList(newRead1));

        CompareConfig config = new CompareConfig(1, false, false);
        MockDiffConsumer diffs = new MockDiffConsumer();
        PartitionReader reader = new PartitionReader(REGION, config, refRecordProver, newRecordProvider, diffs);
        reader.run();
        Statistics stats = reader.stats();

        assertEquals(1, diffs.Diffs.size());
        DiffRecord diff = diffs.Diffs.get(0);
        assertEquals(refRead1, diff.Read);
        assertEquals(REF_ONLY, diff.Type);
        assertNull(diff.DiffList);

        assertEquals(1, stats.RefReadCount);
        assertEquals(0, stats.NewReadCount);
        assertEquals(1, stats.DiffCount);
    }

    @Test
    public void testHaltNewProcessing()
    {
        // TODO: extract common
        SAMRecord newRead1 = createSamRecord("READ_001", CHR_1, 100);
        SAMRecord newRead2 = createSamRecord("READ_001", CHR_1, 101);

        MockSAMRecordProvider refRecordProver = new MockSAMRecordProvider(Lists.newArrayList());
        MockSAMRecordProvider newRecordProvider = new MockSAMRecordProvider(Lists.newArrayList(newRead1, newRead2));

        CompareConfig config = new CompareConfig(1, false, false);
        MockDiffConsumer diffs = new MockDiffConsumer();
        PartitionReader reader = new PartitionReader(REGION, config, refRecordProver, newRecordProvider, diffs);
        reader.run();
        Statistics stats = reader.stats();

        assertEquals(1, diffs.Diffs.size());
        DiffRecord diff = diffs.Diffs.get(0);
        assertEquals(newRead1, diff.Read);
        assertEquals(NEW_ONLY, diff.Type);
        assertNull(diff.DiffList);

        assertEquals(0, stats.RefReadCount);
        assertEquals(1, stats.NewReadCount);
        assertEquals(1, stats.DiffCount);
    }

    @Test
    public void testCigarMismatch()
    {
        // TODO: extract common
        SAMRecord refRead1 = createSamRecord("READ_001", CHR_1, 100, false, false, "151M", "151M");
        SAMRecord newRead1 = createSamRecord("READ_001", CHR_1, 100, false, false, "1S150M", "151M");

        MockSAMRecordProvider refRecordProver = new MockSAMRecordProvider(Lists.newArrayList(refRead1));
        MockSAMRecordProvider newRecordProvider = new MockSAMRecordProvider(Lists.newArrayList(newRead1));

        CompareConfig config = new CompareConfig(0, false, false);
        MockDiffConsumer diffs = new MockDiffConsumer();
        PartitionReader reader = new PartitionReader(REGION, config, refRecordProver, newRecordProvider, diffs);
        reader.run();
        Statistics stats = reader.stats();

        assertEquals(1, diffs.Diffs.size());
        DiffRecord diff = diffs.Diffs.get(0);
        assertEquals(refRead1, diff.Read);
        assertEquals(VALUE, diff.Type);
        assertEquals(1, diff.DiffList.size());
        assertEquals("cigar(151M/1S150M)", diff.DiffList.get(0));

        assertEquals(1, stats.RefReadCount);
        assertEquals(1, stats.NewReadCount);
        assertEquals(1, stats.DiffCount);
    }

    @Test
    public void testMateCigarAttributeMismatch()
    {
        // TODO: extract common
        SAMRecord refRead1 = createSamRecord("READ_001", CHR_1, 100, false, false, "151M", "1S150M");
        SAMRecord newRead1 = createSamRecord("READ_001", CHR_1, 100, false, false, "151M", "151M");

        assertNotNull(refRead1.getStringAttribute(MATE_CIGAR_ATTRIBUTE));
        assertEquals(refRead1.getStringAttribute(MATE_CIGAR_ATTRIBUTE), "1S150M");

        assertNotNull(newRead1.getStringAttribute(MATE_CIGAR_ATTRIBUTE));
        assertEquals(newRead1.getStringAttribute(MATE_CIGAR_ATTRIBUTE), "151M");

        MockSAMRecordProvider refRecordProver = new MockSAMRecordProvider(Lists.newArrayList(refRead1));
        MockSAMRecordProvider newRecordProvider = new MockSAMRecordProvider(Lists.newArrayList(newRead1));

        CompareConfig config = new CompareConfig(0, false, false);
        MockDiffConsumer diffs = new MockDiffConsumer();
        PartitionReader reader = new PartitionReader(REGION, config, refRecordProver, newRecordProvider, diffs);
        reader.run();
        Statistics stats = reader.stats();

        assertEquals(1, diffs.Diffs.size());
        DiffRecord diff = diffs.Diffs.get(0);
        assertEquals(refRead1, diff.Read);
        assertEquals(VALUE, diff.Type);
        assertEquals(1, diff.DiffList.size());
        assertEquals(format("attrib_%s(1S150M/151M)", MATE_CIGAR_ATTRIBUTE), diff.DiffList.get(0));

        assertEquals(1, stats.RefReadCount);
        assertEquals(1, stats.NewReadCount);
        assertEquals(1, stats.DiffCount);
    }

    @Test
    public void testMateCigarAttributeBothNull()
    {
        // TODO: extract common
        SAMRecord refRead1 = createSamRecord("READ_001", CHR_1, 100, false, false, "151M", null);
        SAMRecord newRead1 = createSamRecord("READ_001", CHR_1, 100, false, false, "151M", null);

        assertNull(refRead1.getStringAttribute(MATE_CIGAR_ATTRIBUTE));
        assertNull(newRead1.getStringAttribute(MATE_CIGAR_ATTRIBUTE));

        MockSAMRecordProvider refRecordProver = new MockSAMRecordProvider(Lists.newArrayList(refRead1));
        MockSAMRecordProvider newRecordProvider = new MockSAMRecordProvider(Lists.newArrayList(newRead1));

        CompareConfig config = new CompareConfig(0, false, false);
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
    public void testMateCigarAttributeRefIsNull()
    {
        // TODO: extract common
        SAMRecord refRead1 = createSamRecord_("READ_001", CHR_1, 100, "", "151M", CHR_2, 100, false, false, null, true, null);
        SAMRecord newRead1 = createSamRecord_("READ_001", CHR_1, 100, "", "151M", CHR_2, 100, false, false, null, true, "151M");

        assertNull(refRead1.getStringAttribute(MATE_CIGAR_ATTRIBUTE));

        assertNotNull(newRead1.getStringAttribute(MATE_CIGAR_ATTRIBUTE));
        assertEquals(newRead1.getStringAttribute(MATE_CIGAR_ATTRIBUTE), "151M");

        MockSAMRecordProvider refRecordProver = new MockSAMRecordProvider(Lists.newArrayList(refRead1));
        MockSAMRecordProvider newRecordProvider = new MockSAMRecordProvider(Lists.newArrayList(newRead1));

        CompareConfig config = new CompareConfig(0, false, false);
        MockDiffConsumer diffs = new MockDiffConsumer();
        PartitionReader reader = new PartitionReader(REGION, config, refRecordProver, newRecordProvider, diffs);
        reader.run();
        Statistics stats = reader.stats();

        assertEquals(1, diffs.Diffs.size());
        DiffRecord diff = diffs.Diffs.get(0);
        assertEquals(refRead1, diff.Read);
        assertEquals(VALUE, diff.Type);
        assertEquals(1, diff.DiffList.size());
        assertEquals(format("attrib_%s(missing/151M)", MATE_CIGAR_ATTRIBUTE), diff.DiffList.get(0));

        assertEquals(1, stats.RefReadCount);
        assertEquals(1, stats.NewReadCount);
        assertEquals(1, stats.DiffCount);
    }

    @Test
    public void testMateCigarAttributeNewIsNull()
    {
        // TODO: extract common
        SAMRecord refRead1 = createSamRecord_("READ_001", CHR_1, 100, "", "151M", CHR_2, 100, false, false, null, true, "151M");
        SAMRecord newRead1 = createSamRecord_("READ_001", CHR_1, 100, "", "151M", CHR_2, 100, false, false, null, true, null);

        assertNotNull(refRead1.getStringAttribute(MATE_CIGAR_ATTRIBUTE));
        assertEquals(refRead1.getStringAttribute(MATE_CIGAR_ATTRIBUTE), "151M");

        assertNull(newRead1.getStringAttribute(MATE_CIGAR_ATTRIBUTE));

        MockSAMRecordProvider refRecordProver = new MockSAMRecordProvider(Lists.newArrayList(refRead1));
        MockSAMRecordProvider newRecordProvider = new MockSAMRecordProvider(Lists.newArrayList(newRead1));

        CompareConfig config = new CompareConfig(0, false, false);
        MockDiffConsumer diffs = new MockDiffConsumer();
        PartitionReader reader = new PartitionReader(REGION, config, refRecordProver, newRecordProvider, diffs);
        reader.run();
        Statistics stats = reader.stats();

        assertEquals(1, diffs.Diffs.size());
        DiffRecord diff = diffs.Diffs.get(0);
        assertEquals(refRead1, diff.Read);
        assertEquals(VALUE, diff.Type);
        assertEquals(1, diff.DiffList.size());
        assertEquals(format("attrib_%s(151M/missing)", MATE_CIGAR_ATTRIBUTE), diff.DiffList.get(0));

        assertEquals(1, stats.RefReadCount);
        assertEquals(1, stats.NewReadCount);
        assertEquals(1, stats.DiffCount);
    }

    private static SAMRecord createSamRecord(final String readName, final String chromosome, int posStart)
    {
        return createSamRecord(readName, chromosome, posStart, false, false);
    }

    private static SAMRecord createSamRecord(final String readName, final String chromosome, int posStart, boolean isDuplicate, boolean isConsensus)
    {
        return createSamRecord(readName, chromosome, posStart, false, false, "151M", "151M");
    }

    private static SAMRecord createSamRecord(final String readName, final String chromosome, int posStart, boolean isDuplicate, boolean isConsensus, final String cigar, @Nullable final String mateCigar)
    {
        SAMRecord read = SamRecordTestUtils.createSamRecord(readName, chromosome, posStart, "", cigar, CHR_2, 100, false, false, null, true, mateCigar);

        if (isDuplicate)
            read.setDuplicateReadFlag(true);

        if (isConsensus)
            read.setAttribute(CONSENSUS_READ_ATTRIBUTE, "");

        return read;
    }
}
