package com.hartwig.hmftools.sage.common;

import static java.util.Map.entry;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_READ_CONTEXT_FLANK_SIZE;
import static com.hartwig.hmftools.sage.SageConstants.MIN_CORE_DISTANCE;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.bed.NamedBed;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.SageCallConfig;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.bqr.BqrRecord;
import com.hartwig.hmftools.sage.bqr.BqrRecordMap;
import com.hartwig.hmftools.sage.candidate.AltContext;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.candidate.Candidates;
import com.hartwig.hmftools.sage.candidate.RefContextCache;
import com.hartwig.hmftools.sage.candidate.RefContextConsumer;
import com.hartwig.hmftools.sage.coverage.Coverage;
import com.hartwig.hmftools.sage.evidence.FragmentLengths;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;
import com.hartwig.hmftools.sage.old.IndexedBases;
import com.hartwig.hmftools.sage.old.ReadContext;
import com.hartwig.hmftools.sage.phase.PhaseSetCounter;
import com.hartwig.hmftools.sage.pipeline.RegionResults;
import com.hartwig.hmftools.sage.pipeline.RegionTask;
import com.hartwig.hmftools.sage.quality.QualityCalculator;
import com.hartwig.hmftools.sage.vcf.VcfWriter;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;

public class TestUtils
{
    public static final SageConfig TEST_CONFIG = createSageConfig();

    public static final BqrRecordMap RECALIBRATION = new BqrRecordMap(Collections.emptyList());

    public static final MockRefGenome MOCK_REF_GENOME = new MockRefGenome();

    private static final RefSequence REF_BASES = new RefSequence(550, "TGTTTCTGTTTC".getBytes());

    public static final QualityCalculator QUALITY_CALCULATOR = new QualityCalculator(TEST_CONFIG, RECALIBRATION, REF_BASES, MOCK_REF_GENOME );

    public static SageConfig createSageConfig()
    {
        // add input arguments as necessary or take the defaults
        return new SageConfig(false);
    }

    public static SimpleVariant createSimpleVariant(int position)
    {
        return new SimpleVariant(CHR_1, position, "A", "C");
    }

    public static SimpleVariant createSimpleVariant(int position, final String ref, final String alt)
    {
        return new SimpleVariant(CHR_1, position, ref, alt);
    }

    public static SageVariant createVariant(int position, final String ref, final String alt)
    {
        String readBases = buildReadContextBases(alt);

        // LF          L   I   RC          RF
        // 0123456789  01  2  34  0123456789
        int leftCoreIndex = DEFAULT_READ_CONTEXT_FLANK_SIZE;
        int index = leftCoreIndex + MIN_CORE_DISTANCE;
        int rightCoreIndex = index + alt.length() - 1 + MIN_CORE_DISTANCE;

        IndexedBases indexBases = new IndexedBases(
                position, index, leftCoreIndex, rightCoreIndex, DEFAULT_READ_CONTEXT_FLANK_SIZE, readBases.getBytes());

        return createVariant(CHR_1, position, ref, alt, indexBases);
    }

    public static SageVariant createVariant(final ReadContextCounter readContextCounter)
    {
        Candidate candidate = new Candidate(
                VariantTier.HIGH_CONFIDENCE, readContextCounter.variant(), readContextCounter.readContext(), 1, 1);

        return new SageVariant(candidate, Collections.emptyList(), Lists.newArrayList(readContextCounter));
    }

    public static SageVariant createVariant(final String chromosome, int position, final String ref, final String alt, final IndexedBases indexBases)
    {
        SimpleVariant variant = new SimpleVariant(chromosome, position, ref, alt);

        ReadContext readContext = new ReadContext(position, "", 0, "", indexBases, false);

        ReadContextCounter readCounter = new ReadContextCounter(
                0, variant, readContext, VariantTier.LOW_CONFIDENCE, 100, 1,
                TEST_CONFIG, QUALITY_CALCULATOR, null);

        List<ReadContextCounter> tumorCounters = Lists.newArrayList(readCounter);

        Candidate candidate = new Candidate(
                VariantTier.HIGH_CONFIDENCE, variant, tumorCounters.get(0).readContext(), 1, 1);

        List<ReadContextCounter> normalCounters = Lists.newArrayList();

        return new SageVariant(candidate, normalCounters, tumorCounters);
    }

    public static void setBaseQualities(final SAMRecord record, int baseQual)
    {
        for(int i = 0; i < record.getBaseQualities().length; ++i)
        {
            record.getBaseQualities()[i] = (byte)baseQual;
        }
    }

    public static void setTumorQuality(final SageVariant variant, int count, int quality)
    {
        variant.tumorReadCounters().get(0).readSupportCounts().Full = count;
        variant.tumorReadCounters().get(0).readSupportQualityCounts().Full = quality;
    }

    public static void addLocalPhaseSet(final SageVariant variant, int lps, int readCount)
    {
        variant.tumorReadCounters().get(0).addLocalPhaseSet(lps, readCount, 0);
    }

    public static void clearFilters(final List<SageVariant> variants)
    {
        variants.forEach(x -> x.filters().clear());
    }

    public static void clearFilters(final SageVariant variant)
    {
        variant.filters().clear();
    }

    public static boolean hasFilter(final SageVariant variant, final String filter)
    {
        return variant.filters().contains(filter);
    }

    public static String buildReadContextBases(final String alt)
    {
        String flank = generateRandomBases(DEFAULT_READ_CONTEXT_FLANK_SIZE);
        String core = generateRandomBases(MIN_CORE_DISTANCE);
        return flank + core + alt + core + flank;
    }

    public static SimpleVariant createSimpleVariant(final String chromosome, int position, final String ref, final String alt)
    {
        return new SimpleVariant(chromosome, position, ref, alt);
    }

    public static ReadContext createReadContext(
            int refPosition, int readIndex, int leftCentreIndex, int rightCentreIndex, String readBases, String microhomology)
    {
        int adjLeftCentreIndex = Math.max(leftCentreIndex, 0);
        int adjRightCentreIndex = Math.min(rightCentreIndex, readBases.length() - 1);
        boolean incompleteCore = adjLeftCentreIndex != leftCentreIndex || adjRightCentreIndex != rightCentreIndex;

        IndexedBases readBasesIndexed = new IndexedBases(refPosition, readIndex, adjLeftCentreIndex, adjRightCentreIndex, 0, readBases.getBytes());

        return new ReadContext(refPosition, "", 0, microhomology, readBasesIndexed, incompleteCore);
    }

    public static SAMRecord createSamRecord(
            final String readId, final String chrStr, int readStart, final String readBases, final String cigar)
    {
        SAMRecordSetBuilder recordBuilder = new SAMRecordSetBuilder();
        recordBuilder.setUnmappedHasBasesAndQualities(false);

        HumanChromosome chromosome = HumanChromosome.fromString(chrStr);

        SAMRecord record = recordBuilder.addFrag(
                readId, chromosome.ordinal(), readStart, false, false,
                cigar, readBases, 37, false);

        record.setReadBases(readBases.getBytes());

        final byte[] qualities = new byte[readBases.length()];

        for(int i = 0; i < readBases.length(); ++i)
            qualities[i] = 37;

        record.setBaseQualities(qualities);
        record.setReferenceName(chrStr);
        record.setReferenceIndex(chromosome.ordinal()); // need to override since no header is present

        record.setMateReferenceName(chrStr);
        record.setMateReferenceIndex(chromosome.ordinal());
        record.setMateAlignmentStart(readStart + 300);
        record.setMateNegativeStrandFlag(true);
        record.setAttribute(MATE_CIGAR_ATTRIBUTE, cigar);

        // to be correct this should match the cigar element count
        record.setAttribute(NUM_MUTATONS_ATTRIBUTE, 1);
        record.setFirstOfPairFlag(true);

        record.setReadPairedFlag(true);
        record.setProperPairFlag(true);

        return record;
    }

    public static SAMRecord buildSamRecord(final int alignmentStart, final String cigar, final String readString, final String qualities)
    {
        return buildSamRecord(alignmentStart, cigar, readString, qualities.getBytes());
    }

    public static SAMRecord buildSamRecord(final int alignmentStart, final String cigar, final String readString, final byte[] qualities)
    {
        final SAMRecord record = new SAMRecord(null);
        record.setReferenceName(CHR_1);
        record.setAlignmentStart(alignmentStart);
        record.setCigarString(cigar);
        record.setReadString(readString);
        record.setReadNegativeStrandFlag(false);
        record.setBaseQualities(qualities);
        record.setMappingQuality(20);
        record.setDuplicateReadFlag(false);
        record.setReadUnmappedFlag(false);
        record.setProperPairFlag(true);
        record.setReadPairedFlag(true);
        record.setInferredInsertSize(600);
        return record;
    }

    public static List<Candidate> findInitialCandidates(final ChrBaseRegion region, final RefGenomeInterface refGenome,
            final List<SAMRecord> reads)
    {
        List<VariantHotspot> hotspots = Lists.newArrayList();
        List<BaseRegion> panelRegions = Lists.newArrayList();
        List<BaseRegion> highConfidenceRegions = Lists.newArrayList();
        Set<Integer> specificPositions = Sets.newHashSet();
        SageConfig config = new SageConfig(false);
        RefSequence refSequence = new RefSequence(region, refGenome);

        RefContextCache refContextCache = new RefContextCache(config, hotspots, panelRegions);
        RefContextConsumer refContextConsumer = new RefContextConsumer(config, region, refSequence, refContextCache, hotspots);

        // generate read contexts
        List<AltContext> altContexts = Lists.newArrayList();
        reads.forEach(refContextConsumer::processRead);
        altContexts.addAll(refContextCache.altContexts());

        // choose candidates from the read contexts
        Candidates initialCandidates = new Candidates(hotspots, panelRegions, highConfidenceRegions);
        initialCandidates.addSingleSample(altContexts);
        return initialCandidates.candidates(specificPositions);
    }

    public static void dumpVariantVis(final String sampleId, final List<SAMRecord> reads,
            final ChrBaseRegion region, final MockRefGenome refGenome, final String fullVisOutputDir)
    {
        int taskId = 0;

        SamSlicerInterface samSlicer = (final Consumer<SAMRecord> consumer) -> reads.forEach(consumer);
        SamSlicerFactory samSlicerFactory = new SamSlicerFactory();
        samSlicerFactory.addSamSlicer(sampleId, samSlicer);

        VcfWriter vcfWriter = null;
        RegionResults results = new RegionResults(vcfWriter);

        SageCallConfig config = new SageCallConfig(fullVisOutputDir);
        config.TumorIds.add(sampleId);

        List<BqrRecord> bqrRecords = Lists.newArrayList();
        Map<String, BqrRecordMap> qualityRecalibrationMap = Map.ofEntries(
                entry(sampleId, new BqrRecordMap(bqrRecords))
        );

        Collection<NamedBed> panel = Lists.newArrayList();
        Coverage coverage = new Coverage(config.TumorIds, panel, config.Common);

        List<VariantHotspot> hotspots = Lists.newArrayList();
        List<BaseRegion> panelRegions = Lists.newArrayList();
        List<BaseRegion> highConfidenceRegions = Lists.newArrayList();
        List<TranscriptData> transcripts = null;
        PhaseSetCounter phaseSetCounter = new PhaseSetCounter();
        FragmentLengths fragmentLengths = null;

        RegionTask regionTask = new RegionTask(taskId, region, results, config, refGenome, hotspots, panelRegions, transcripts,
                highConfidenceRegions, qualityRecalibrationMap, phaseSetCounter, coverage, samSlicerFactory, fragmentLengths);
        regionTask.run();
    }
}
