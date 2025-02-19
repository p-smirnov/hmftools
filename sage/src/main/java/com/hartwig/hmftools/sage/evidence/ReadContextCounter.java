package com.hartwig.hmftools.sage.evidence;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.extractUmiType;
import static com.hartwig.hmftools.common.variant.SageVcfTags.UMI_TYPE_COUNT;
import static com.hartwig.hmftools.common.variant.VariantReadSupport.CORE;
import static com.hartwig.hmftools.common.variant.VariantReadSupport.FULL;
import static com.hartwig.hmftools.common.variant.VariantReadSupport.REALIGNED;
import static com.hartwig.hmftools.common.variant.VariantReadSupport.REF;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.SageConstants.EVIDENCE_MIN_MAP_QUAL;
import static com.hartwig.hmftools.sage.SageConstants.MQ_RATIO_SMOOTHING;
import static com.hartwig.hmftools.sage.SageConstants.REALIGN_READ_CONTEXT_MIN_SEARCH_BUFFER;
import static com.hartwig.hmftools.sage.SageConstants.REALIGN_READ_CONTEXT_MIN_SEARCH_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.REALIGN_READ_MIN_INDEL_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.REQUIRED_UNIQUE_FRAG_COORDS;
import static com.hartwig.hmftools.sage.SageConstants.SC_READ_EVENTS_FACTOR;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.NONE;
import static com.hartwig.hmftools.sage.evidence.ReadEdgeDistance.calcAdjustedVariantPosition;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.ALT_SUPPORT;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.CHIMERIC;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.IN_SPLIT;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.MAP_QUAL;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.MAX_COVERAGE;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.NON_CORE;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.REF_SUPPORT;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.SOFT_CLIP;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.UNRELATED;
import static com.hartwig.hmftools.sage.evidence.RealignedType.EXACT;
import static com.hartwig.hmftools.sage.evidence.RealignedType.LENGTHENED;
import static com.hartwig.hmftools.sage.evidence.RealignedType.SHORTENED;
import static com.hartwig.hmftools.sage.filter.ReadFilters.isChimericRead;
import static com.hartwig.hmftools.sage.quality.QualityCalculator.isImproperPair;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;

import java.util.List;

import javax.annotation.Nullable;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.UmiReadType;
import com.hartwig.hmftools.common.sage.FragmentLengthCounts;
import com.hartwig.hmftools.common.variant.VariantReadSupport;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.common.ReadContextMatcher;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.common.ReadContextMatch;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantTier;
import com.hartwig.hmftools.sage.filter.FragmentCoords;
import com.hartwig.hmftools.sage.filter.StrandBiasData;
import com.hartwig.hmftools.sage.quality.QualityCalculator;
import com.hartwig.hmftools.sage.quality.UltimaQualModel;
import com.hartwig.hmftools.sage.common.NumberEvents;
import com.hartwig.hmftools.sage.vis.VariantVis;
import com.hartwig.hmftools.sage.sync.FragmentData;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class ReadContextCounter
{
    private final int mId;
    private final VariantTier mTier;
    private final SimpleVariant mVariant;
    private final VariantReadContext mReadContext;
    private final ReadContextMatcher mReadContextMatcher;
    private final SageConfig mConfig;
    private final QualityCalculator mQualityCalculator;
    private final String mSample;
    private final int mMaxCoverage;
    private final VariantVis mVariantVis;

    // local variant-related state
    private final int mMinNumberOfEvents;
    private final boolean mIsMnv;
    private final boolean mIsIndel;
    private final BqrQualCache mBqrQualCache;

    // counts and quals by support type
    private final ReadSupportCounts mQualities;
    private final ReadSupportCounts mCounts;

    private final StrandBiasData mAltFragmentStrandBias;
    private final StrandBiasData mRefFragmentStrandBias;
    private final StrandBiasData mAltReadStrandBias;
    private final StrandBiasData mRefReadStrandBias;

    private final JitterData mJitterData;
    private int mImproperPairCount;

    private final QualCounters mQualCounters;

    private int mMaxCandidateDeleteLength;
    private final ReadEdgeDistance mReadEdgeDistance;

    private List<Integer> mLocalPhaseSets;
    private List<Integer> mLpsCounts;
    private int[] mUmiTypeCounts;
    private FragmentLengthCounts mFragmentLengthData;
    private FragmentCoords mFragmentCoords;

    public ReadContextCounter(
            final int id, final VariantReadContext readContext, final VariantTier tier, final int maxCoverage, final int minNumberOfEvents,
            final SageConfig config, final QualityCalculator qualityCalculator, final String sampleId)
    {
        mId = id;

        mTier = tier;
        mMaxCoverage = maxCoverage;
        mMinNumberOfEvents = minNumberOfEvents;
        mSample = sampleId;
        mQualityCalculator = qualityCalculator;
        mConfig = config;

        mReadContext = readContext;
        mReadContextMatcher = new ReadContextMatcher(mReadContext);
        mVariant = readContext.variant();

        mVariantVis = config.Visualiser.Enabled && config.Visualiser.processVariant(mVariant)
                ? new VariantVis(mConfig, mSample, mVariant, mReadContext, mTier) : null;

        // set local state to avoid testing on each read
        mIsMnv = mVariant.isMNV();
        mIsIndel = mVariant.isIndel();
        mBqrQualCache = !mIsIndel ? new BqrQualCache(mVariant.position(), mVariant.alt()) : null;

        mQualities = new ReadSupportCounts();
        mCounts = new ReadSupportCounts();

        mJitterData = new JitterData();

        mAltFragmentStrandBias = new StrandBiasData(true);
        mRefFragmentStrandBias = new StrandBiasData(false);
        mAltReadStrandBias = new StrandBiasData(true);
        mRefReadStrandBias = new StrandBiasData(false);

        mImproperPairCount = 0;

        mQualCounters = new QualCounters();
        mMaxCandidateDeleteLength = 0;

        mReadEdgeDistance = new ReadEdgeDistance(calcAdjustedVariantPosition(mVariant.position(), indelLength()));

        mLocalPhaseSets = null;
        mLpsCounts = null;
        mUmiTypeCounts = null;
        mFragmentLengthData = mConfig.WriteFragmentLengths ? new FragmentLengthCounts() : null;
        mFragmentCoords = mConfig.Quality.HighDepthMode ? new FragmentCoords(REQUIRED_UNIQUE_FRAG_COORDS) : null;
    }

    public int id() { return mId; }
    public SimpleVariant variant() { return mVariant; }
    public VariantReadContext readContext() { return mReadContext; }
    public ReadContextMatcher readContextMatcher() { return mReadContextMatcher; }
    public VariantTier tier() { return mTier; }
    public int indelLength() { return mVariant.isIndel() ? max(mVariant.alt().length(), mVariant.ref().length()) : 0; }
    public boolean isSnv() { return mVariant.isSNV(); }
    public boolean isIndel() { return mIsIndel; }
    public final BqrQualCache bqrQualCache() { return mBqrQualCache; }
    public String chromosome() { return mVariant.chromosome(); }
    public int position() { return mVariant.position(); }
    public String ref() { return mVariant.ref(); }
    public String alt() { return mVariant.alt(); }

    public int altSupport() { return mCounts.altSupport(); }
    public int strongAltSupport() { return mCounts.strongSupport(); }
    public int refSupport() { return mCounts.Ref; }

    public int depth() { return mCounts.Total; }

    public QualCounters qualCounters() { return mQualCounters; }

    public int baseQualityTotal() { return mQualCounters.BaseQualityTotal; }
    public int altBaseQualityTotal() { return mQualCounters.AltBaseQualityTotal; }

    public long mapQualityTotal() { return mQualCounters.MapQualityTotal; }
    public long altMapQualityTotal() { return mQualCounters.AltMapQualityTotal; }

    public double vaf()
    {
        return mCounts.Total == 0 ? 0d : mCounts.altSupport() / (double)mCounts.Total;
    }

    public int tumorQuality()
    {
        int tumorQuality = mQualities.Full + mQualities.PartialCore + mQualities.Realigned;
        return Math.max(0, tumorQuality - mJitterData.penalty());
    }

    public int[] counts() { return mCounts.toArray(); }
    public int[] quality() { return mQualities.toArray(); }

    public int[] jitter()
    {
        return mJitterData.summary();
    }

    public ArtefactContext artefactContext() { return mReadContext.artefactContext(); }
    public UltimaQualModel ultimaQualModel() { return mReadContext.ultimaQualModel(); }

    public StrandBiasData fragmentStrandBiasAlt() { return mAltFragmentStrandBias; }
    public StrandBiasData fragmentStrandBiasRef() { return mRefFragmentStrandBias; }
    public StrandBiasData readStrandBiasAlt() { return mAltReadStrandBias; }
    public StrandBiasData readStrandBiasRef() { return mRefReadStrandBias; }
    public int improperPairCount() { return mImproperPairCount; }

    public ReadEdgeDistance readEdgeDistance() { return mReadEdgeDistance; }
    public int minNumberOfEvents() { return mMinNumberOfEvents; }

    @Nullable
    public VariantVis variantVis() { return mVariantVis; }

    public double averageAltBaseQuality()
    {
        // excludes realigned
        int supportCount = mCounts.Full + mCounts.PartialCore + mCounts.Core;
        return supportCount > 0 ? mQualCounters.AltBaseQualityTotal / (double)supportCount : 0;
    }

    public void setMaxCandidateDeleteLength(int length) { mMaxCandidateDeleteLength = length; }
    public int maxCandidateDeleteLength() { return mMaxCandidateDeleteLength; }

    public List<Integer> localPhaseSets() { return mLocalPhaseSets; }
    public List<Integer> lpsCounts() { return mLpsCounts; }

    public int[] umiTypeCounts() { return mUmiTypeCounts; }
    public FragmentLengthCounts fragmentLengths() { return mFragmentLengthData; }

    public boolean exceedsMaxCoverage() { return mCounts.Total >= mMaxCoverage; }

    public boolean belowMinFragmentCoords() { return mFragmentCoords != null && !mFragmentCoords.atCapacity(); }

    public String toString()
    {
        return format("id(%d) var(%s) core(%s) counts(f=%d p=%d c=%d)",
                mId, varString(), mReadContext.toString(), mCounts.Full, mCounts.PartialCore, mCounts.Core);
    }

    public String varString()
    {
        return format("%s:%d %s>%s", mVariant.chromosome(), mVariant.position(), mVariant.ref(), mVariant.alt());
    }

    public ReadMatchType processRead(final SAMRecord record, int numberOfEvents, @Nullable final FragmentData fragmentData)
    {
        if(exceedsMaxCoverage())
            return MAX_COVERAGE;

        if(mTier != VariantTier.HOTSPOT && record.getMappingQuality() < EVIDENCE_MIN_MAP_QUAL)
        {
            addVariantVisRecord(record, ReadContextMatch.NONE, null, fragmentData);
            return MAP_QUAL;
        }

        if(mConfig.Quality.HighDepthMode && isChimericRead(record))
        {
            addVariantVisRecord(record, ReadContextMatch.NONE, null, fragmentData);
            return CHIMERIC;
        }

        RawContext rawContext = RawContext.create(mVariant, record);

        if(mConfig.Quality.HighDepthMode && rawContext.PositionType == VariantReadPositionType.SOFT_CLIP)
        {
            addVariantVisRecord(record, ReadContextMatch.NONE, null, fragmentData);
            return SOFT_CLIP;
        }

        // TODO: need to check if the soft-clip bases cover the variant and allow for realignment
        if(rawContext.ReadIndex < 0)
        {
            addVariantVisRecord(record, ReadContextMatch.NONE, null, fragmentData);
            return UNRELATED;
        }

        if(rawContext.PositionType == VariantReadPositionType.SKIPPED)
        {
            addVariantVisRecord(record, ReadContextMatch.NONE, null, fragmentData);
            return IN_SPLIT;
        }

        int readIndex = rawContext.ReadIndex;

        boolean covered = mReadContextMatcher.coversVariant(record, readIndex);

        if(!covered)
        {
            addVariantVisRecord(record, ReadContextMatch.NONE, null, fragmentData);
            return NON_CORE;
        }

        double adjustedNumOfEvents = numberOfEvents;

        if(mIsMnv)
            adjustedNumOfEvents = NumberEvents.calcWithMnvRaw(numberOfEvents, mVariant.ref(), mVariant.alt());

        if(max(mVariant.ref().length(), mVariant.alt().length()) <= SC_READ_EVENTS_FACTOR)
        {
            // penalise variants except long INDELs for their soft-clipped bases
            adjustedNumOfEvents += NumberEvents.calcSoftClipAdjustment(record);
        }

        adjustedNumOfEvents = max(mMinNumberOfEvents, adjustedNumOfEvents);

        double rawBaseQuality = mQualityCalculator.rawBaseQuality(this, readIndex, record);

        if(mConfig.Quality.HighDepthMode && rawBaseQuality < mConfig.Quality.HighBaseQualLimit)
        {
            // TODO: check if strand bias and frag coords need to be updated in this scenario - cannot now rely on raw AltSupport
            // if(rawContext.AltSupport)
            //     countAltSupportMetrics(record, fragmentData);

            addVariantVisRecord(record, ReadContextMatch.NONE, null, fragmentData);
            return UNRELATED;
        }

        QualityCalculator.QualityScores qualityScores = mQualityCalculator.calculateQualityScores(
                this, readIndex, record, adjustedNumOfEvents, rawBaseQuality);

        double quality = qualityScores.ModifiedQuality;

        // Check if FULL, PARTIAL, OR CORE
        ReadContextMatch matchType = rawContext.PositionType != VariantReadPositionType.DELETED ?
                determineReadContextMatch(record, readIndex, true) : NONE;

        if(matchType.SupportsAlt)
        {
            VariantReadSupport readSupport = matchType.toReadSupport();

            registerReadSupport(record, readSupport, quality);

            mQualCounters.MapQualityTotal += record.getMappingQuality();
            mQualCounters.AltMapQualityTotal += record.getMappingQuality();
            mQualCounters.AltBaseQualityTotal += qualityScores.RecalibratedBaseQuality;
            mQualCounters.BaseQualityTotal += qualityScores.RecalibratedBaseQuality;

            mReadEdgeDistance.update(record, fragmentData, true);

            addVariantVisRecord(record, matchType, qualityScores, fragmentData);
            logReadEvidence(record, matchType, readIndex, quality);
            countAltSupportMetrics(record, fragmentData);

            checkImproperCount(record);
            return ALT_SUPPORT;
        }

        boolean canRealign = matchType != ReadContextMatch.REF
                && (abs(mVariant.indelLength()) >= REALIGN_READ_MIN_INDEL_LENGTH || readHasIndelInCore(record));

        RealignedContext realignment = canRealign ? checkRealignment(record, readIndex) : RealignedContext.NONE;

        if(realignment.Type == EXACT)
        {
            matchType = ReadContextMatch.REALIGNED;
            registerReadSupport(record, REALIGNED, quality);

            mQualCounters.MapQualityTotal += record.getMappingQuality();
            mQualCounters.AltMapQualityTotal += record.getMappingQuality();

            addVariantVisRecord(record, matchType, qualityScores, fragmentData);
            logReadEvidence(record, matchType, readIndex,quality);

            return ALT_SUPPORT;
        }

        mQualCounters.BaseQualityTotal += qualityScores.RecalibratedBaseQuality;

        RealignedContext jitterRealign = Realignment.realignedAroundIndex(mReadContext, readIndex, record);

        if(rawContext.PositionType == VariantReadPositionType.SOFT_CLIP)
        {
            if(jitterRealign.Type != LENGTHENED && jitterRealign.Type != SHORTENED)
            {
                addVariantVisRecord(record, ReadContextMatch.NONE, qualityScores, fragmentData);
                return SOFT_CLIP;
            }
        }

        mQualCounters.MapQualityTotal += record.getMappingQuality();

        VariantReadSupport readSupport = null;

        if(matchType == ReadContextMatch.REF)
        {
            readSupport = REF;

            mRefFragmentStrandBias.registerFragment(record);
            mRefReadStrandBias.registerRead(record, fragmentData, this);

            mReadEdgeDistance.update(record, fragmentData, false);
        }

        registerReadSupport(record, readSupport, quality);

        // add to jitter penalty as a function of the number of repeats found
        mJitterData.update(jitterRealign, mConfig.Quality);

        addVariantVisRecord(record, matchType, qualityScores, fragmentData);
        logReadEvidence(record, matchType, readIndex, quality);

        return matchType == ReadContextMatch.REF ? REF_SUPPORT : UNRELATED;
    }

    private ReadContextMatch determineReadContextMatch(final SAMRecord record, int readIndex, boolean allowCoreVariation)
    {
        // CLEAN-UP: fix this for RNA
        // better approaches would be to have the read matcher stop checking if it is in a N-section,
        // or to avoid affecting all DNA samples, make a new SAMRecord with the Ns filled out sufficiently

        /*
        ReadIndexBases readIndexBases;
        if(record.getCigar().containsOperator(CigarOperator.N))
        {
            readIndexBases = SplitReadUtils.expandSplitRead(readIndex, record);
        }

        final ReadContextMatch match = mReadContext.indexedBases().matchAtPosition(
                readIndexBases, record.getBaseQualities(),
                allowCoreVariation ? mAllowWildcardMatchInCore : false,
                allowCoreVariation ? mMaxCoreMismatches : 0);
        */

        // REALIGN: realignment did not allow low-qual mismatches or wildcards - is that still a necessary condition?

        return mReadContextMatcher.determineReadMatch(record, readIndex);
    }

    private void registerReadSupport(final SAMRecord record, @Nullable final VariantReadSupport support, final double quality)
    {
        mCounts.addSupport(support, 1);
        mQualities.addSupport(support, (int)quality);

        boolean supportsVariant = support != null
                && (support == FULL || support == VariantReadSupport.PARTIAL_CORE || support == CORE || support == REALIGNED);

        if(mConfig.Sequencing.HasUMIs)
        {
            countUmiType(record, supportsVariant);
        }

        if(mFragmentLengthData != null && (support == REF || supportsVariant))
        {
            mFragmentLengthData.addLength(abs(record.getInferredInsertSize()), supportsVariant);
        }
    }

    private void countUmiType(final SAMRecord record, final boolean supportsVariant)
    {
        if(mUmiTypeCounts == null)
        {
            // 3 total depth values followed by the 3 variant support values
            mUmiTypeCounts = new int[UMI_TYPE_COUNT];
        }

        UmiReadType umiReadType = extractUmiType(record);

        // add to total and variant support if applicable
        ++mUmiTypeCounts[umiReadType.ordinal()];

        if(supportsVariant)
            ++mUmiTypeCounts[umiReadType.ordinal() + 3];
    }

    private void addVariantVisRecord(
            final SAMRecord record, final ReadContextMatch matchType,
            @Nullable QualityCalculator.QualityScores modifiedQualities, @Nullable final FragmentData fragmentData)
    {
        if(mVariantVis != null)
            mVariantVis.addEvidence(record, fragmentData, matchType, modifiedQualities);
    }

    private void logReadEvidence(final SAMRecord record, final ReadContextMatch matchType, int readIndex, double quality)
    {
        if(!mConfig.LogEvidenceReads || !SG_LOGGER.isTraceEnabled())
            return;

        // mQualityCalculator.logReadQualCalcs(this, readIndex, record, adjustedNumOfEvents);

        // Variant,MatchType,ReadId,ReadStart,Cigar,LeftCore,Index,RightCore,ReadIndex,Quality
        SG_LOGGER.trace("READ_EV,{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                mSample, chromosome(), position(), ref(), alt(),
                matchType, record.getReadName(), record.getAlignmentStart(), record.getCigarString(),
                mReadContext.CoreIndexStart, mReadContext.VarReadIndex,
                mReadContext.CoreIndexEnd, readIndex, format("%.1f", quality));
    }

    private boolean readHasIndelInCore(final SAMRecord record)
    {
        if(!record.getCigar().containsOperator(D) && !record.getCigar().containsOperator(I))
            return false;

        int variantLeftCorePos = mReadContext.CorePositionStart;
        int variantRightCorePos = mReadContext.CorePositionEnd;

        int currentPos = record.getAlignmentStart() - 1;

        // eg 2S10M2D10M starting at 100: first non-SC element, in this case a delete, starts at 109
        for(CigarElement element : record.getCigar())
        {
            if(element.getOperator() == S)
                continue;

            if(element.getOperator() == I || element.getOperator() == D)
            {
                int indelLowerPos = currentPos;
                int indelUpperPos = indelLowerPos + (element.getOperator() == D ? element.getLength() : 1);

                if(positionsOverlap(variantLeftCorePos, variantRightCorePos, indelLowerPos, indelUpperPos))
                    return true;
            }
            else if(element.getOperator() == M)
            {
                currentPos += element.getLength();
            }
        }

        return false;
    }

    private RealignedContext checkRealignment(final SAMRecord record, int readIndex)
    {
        // do a comparison of the bases going from this calculated RI - coreLength - flankLength vs the variant's full read base sequence

        // the read index corresponding to the ref position at the end of the core

        int realignedReadIndex = Realignment.realignedReadIndexPosition(mReadContext, record);

        if(realignedReadIndex < 0 || realignedReadIndex >= record.getReadBases().length)
            return RealignedContext.NONE;

        ReadContextMatch match = determineReadContextMatch(record, realignedReadIndex, true);

        if(match == ReadContextMatch.FULL || match == ReadContextMatch.PARTIAL_CORE)
            return new RealignedContext(EXACT, mReadContext.totalLength(), realignedReadIndex);

        return RealignedContext.NONE;
    }

    public void addLocalPhaseSet(int lps, int readCount, double allocCount)
    {
        if(mLocalPhaseSets == null)
        {
            mLocalPhaseSets = Lists.newArrayList();
            mLpsCounts = Lists.newArrayList();
        }

        // add in order of highest counts
        int index = 0;
        while(index < mLpsCounts.size())
        {
            int existingCount = mLpsCounts.get(index);
            if(readCount + allocCount > existingCount)
                break;

            ++index;
        }

        mLocalPhaseSets.add(index, lps);
        int lpsTotalCount = readCount + (int)allocCount;
        mLpsCounts.add(index, lpsTotalCount);
    }

    private void countAltSupportMetrics(final SAMRecord record, final FragmentData fragmentData)
    {
        mAltReadStrandBias.registerRead(record, fragmentData, this);
        mAltFragmentStrandBias.registerFragment(record);

        if(mFragmentCoords != null)
        {
            if(fragmentData != null)
                mFragmentCoords.addRead(fragmentData.First, fragmentData.Second);
            else
                mFragmentCoords.addRead(record, null);
        }
    }

    private void checkImproperCount(final SAMRecord record)
    {
        if(isImproperPair(record) || record.getSupplementaryAlignmentFlag())
        {
            mImproperPairCount++;
        }
    }

    public void applyMapQualityRatio()
    {
        int depth = depth();
        int avgTotalMapQuality = depth > 0 ? (int)round(mapQualityTotal() / (double)depth) : 0;

        if(avgTotalMapQuality == 0)
            return;

        int altSupport = altSupport();
        int avgAltMapQuality = altSupport > 0 ? (int)round(altMapQualityTotal() / (double)altSupport) : 0;

        double ratioRaw = (avgAltMapQuality + MQ_RATIO_SMOOTHING) / (avgTotalMapQuality + MQ_RATIO_SMOOTHING);
        double calcRatio = pow(min(1, ratioRaw), mConfig.Quality.MapQualityRatioFactor);

        mQualities.applyRatio(calcRatio);
    }

    public boolean logEvidence() { return mConfig.LogEvidenceReads; }

    @VisibleForTesting
    public ReadSupportCounts readSupportQualityCounts() { return mQualities; };
    public ReadSupportCounts readSupportCounts() { return mCounts; }
    public FragmentCoords fragmentCoords() { return mFragmentCoords; }


    // CLEAN-UP: old realignment methods
    private RealignedContext checkRealignmentOld(final SAMRecord record)
    {
        // try left and right alignment in turn
        int realignLeftReadIndex = calcLeftAlignmentIndex(record);

        if(realignLeftReadIndex >= 0) //  && realignLeftReadIndex != readIndex
        {
            ReadContextMatch match = determineReadContextMatch(record, realignLeftReadIndex, false);

            if(match == ReadContextMatch.FULL || match == ReadContextMatch.PARTIAL_CORE)
                return new RealignedContext(EXACT, mReadContext.totalLength(), realignLeftReadIndex);
        }

        int realignRightReadIndex = calcRightAlignmentIndex(record);

        if(realignRightReadIndex >= 0)
        {
            // still need to test even if this index matches the original readIndex since if the readIndex was in a delete
            // it will have been skipped above
            ReadContextMatch match = determineReadContextMatch(record, realignRightReadIndex, false);

            if(match == ReadContextMatch.FULL || match == ReadContextMatch.PARTIAL_CORE)
                return new RealignedContext(RealignedType.EXACT, mReadContext.totalLength(), realignRightReadIndex);
            // else if(match == ReadContextMatch.CORE_PARTIAL) // REALIGN
            //    return new RealignedContext(RealignedType.CORE_PARTIAL, mReadContext.indexedBases().length(), realignRightReadIndex);
        }

        // try a simple string search and take it as exact if the matched index is within the expected range
        String readContext = mReadContext.readBases();

        if(readContext.length() >= REALIGN_READ_CONTEXT_MIN_SEARCH_LENGTH)
        {
            int matchedReadIndex = record.getReadString().indexOf(readContext);

            if(matchedReadIndex >= 0)
            {
                int matchedIndex = matchedReadIndex + mReadContext.leftLength();
                if(abs(matchedIndex - realignLeftReadIndex) <= REALIGN_READ_CONTEXT_MIN_SEARCH_BUFFER
                        || abs(matchedIndex - realignRightReadIndex) <= REALIGN_READ_CONTEXT_MIN_SEARCH_BUFFER)
                {
                    return new RealignedContext(RealignedType.EXACT, mReadContext.totalLength(), matchedIndex);
                }
            }
        }

        return RealignedContext.NONE;
    }

    private int calcLeftAlignmentIndex(final SAMRecord record)
    {
        // Left alignment: Match full read context starting at base = pos - rc_index
        int leftCoreOffset = mReadContext.leftCoreLength();
        int realignLeftCorePos = position() - leftCoreOffset;
        int realignLeftCoreIndex = record.getReadPositionAtReferencePosition(realignLeftCorePos);

        if(realignLeftCoreIndex > 0)
        {
            int realignLeftReadIndex = realignLeftCoreIndex - 1 + leftCoreOffset;
            return realignLeftReadIndex;
        }

        int deleteCount = (int)record.getCigar().getCigarElements().stream().filter(x -> x.getOperator() == D).count();
        if(deleteCount == 0)
            return -1;

        int deleteStartPos = record.getAlignmentStart();
        int deleteStartIndex = 0;
        for(CigarElement element : record.getCigar())
        {
            if(element.getOperator() == S || element.getOperator() == I)
                continue;

            if(element.getOperator() == M)
            {
                deleteStartPos += element.getLength();
                deleteStartIndex += element.getLength();
            }
            else if(element.getOperator() == D)
            {
                --deleteCount;

                if(deleteCount == 0)
                    break;

                deleteStartPos += element.getLength();
            }
        }

        --deleteStartPos;

        int posDiff = realignLeftCorePos - deleteStartPos;
        int realignLeftReadIndex = deleteStartIndex + posDiff - 1 + leftCoreOffset;
        return realignLeftReadIndex;
    }

    private int calcRightAlignmentIndex(final SAMRecord record)
    {
        // Right alignment: Match full read context ending at base = pos + length[RC} - rc_index - 1 - length(alt) + length(ref)
        int rightCoreOffset = mReadContext.rightCoreLength();
        int realignRightPos = position() + rightCoreOffset - mVariant.alt().length() + mVariant.ref().length();
        int realignRightCoreIndex = record.getReadPositionAtReferencePosition(realignRightPos);

        if(realignRightCoreIndex > 0)
        {
            int realignRightReadIndex = realignRightCoreIndex - 1 - rightCoreOffset;
            return realignRightReadIndex;
        }

        return -1;
    }
}
