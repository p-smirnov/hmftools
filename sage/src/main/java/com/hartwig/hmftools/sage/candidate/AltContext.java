package com.hartwig.hmftools.sage.candidate;

import static java.lang.Math.min;

import static com.hartwig.hmftools.sage.SageConstants.MIN_SECOND_CANDIDATE_FULL_READS;
import static com.hartwig.hmftools.sage.SageConstants.MIN_SECOND_CANDIDATE_FULL_READS_PERC;
import static com.hartwig.hmftools.sage.common.ReadContextMatcher.compareReadContexts;

import java.util.Collections;
import java.util.List;

import javax.annotation.Nullable;

import com.google.common.collect.Lists;
import com.google.common.primitives.Ints;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.common.ReadContextMatch;

import org.jetbrains.annotations.NotNull;

public class AltContext extends SimpleVariant
{
    public final RefContext RefContext;
    
    private final List<ReadContextCandidate> mReadContextCandidates;

    private int mRawSupportAlt;
    private int mRawBaseQualityAlt;
    private ReadContextCandidate mCandidate;
    private AltContext mSecondCandidate; // relevant if has a different read context and sufficient support

    public AltContext(final RefContext refContext, final String ref, final String alt)
    {
        super(refContext.Chromosome, refContext.Position, ref, alt);
        RefContext = refContext;

        mReadContextCandidates = Lists.newArrayList();
        mCandidate = null;
        mSecondCandidate = null;
    }

    public AltContext(
            final RefContext refContext, final String ref, final String alt, final ReadContextCandidate candidate,
            int rawSupportAlt, int rawBaseQualAlt)
    {
        super(refContext.Chromosome, refContext.Position, ref, alt);
        RefContext = refContext;

        mCandidate = candidate;
        mRawSupportAlt = rawSupportAlt;
        mRawBaseQualityAlt = rawBaseQualAlt;
        mReadContextCandidates = null;
    }

    public void incrementAltRead(int baseQuality)
    {
        mRawSupportAlt++;
        mRawBaseQualityAlt += baseQuality;
    }

    public void addReadContext(int numberOfEvents, final VariantReadContext newReadContext)
    {
        if(!newReadContext.isValid())
            return;

        int coreMatch = 0;
        ReadContextCandidate fullMatchCandidate = null;

        for(ReadContextCandidate candidate : mReadContextCandidates)
        {
            ReadContextMatch match = compareReadContexts(candidate.readContext(), newReadContext);

            switch(match)
            {
                case FULL:
                    candidate.incrementFull(1, numberOfEvents);
                    fullMatchCandidate = candidate;
                    break;

                case CORE:
                    candidate.CoreMatch++;
                    coreMatch++;
                    break;
            }
        }

        if(fullMatchCandidate == null)
        {
            final ReadContextCandidate candidate = new ReadContextCandidate(numberOfEvents, newReadContext);
            candidate.CoreMatch += coreMatch;
            mReadContextCandidates.add(candidate);
        }
    }

    public int readContextSupport() { return mCandidate.FullMatch; }
    public int minNumberOfEvents()
    {
        return mCandidate.minNumberOfEvents();
    }

    public boolean hasValidCandidate() { return mCandidate != null; }
    public boolean hasSecondCandidate() { return mSecondCandidate != null; }
    public AltContext secondCandidate() { return mSecondCandidate; }

    public void selectCandidates()
    {
        if(mReadContextCandidates.isEmpty())
            return;

        // sort by full, then partial then core read counts
        Collections.sort(mReadContextCandidates);

        if(mReadContextCandidates.isEmpty())
            return;

        mCandidate = mReadContextCandidates.get(0);

        if(mReadContextCandidates.size() > 1)
        {
            final String topCore = mReadContextCandidates.get(0).readContext().coreStr();
            double topCandidateRcThreshold = mReadContextCandidates.get(0).FullMatch * MIN_SECOND_CANDIDATE_FULL_READS_PERC;

            // add a second if its core is different and it has sufficient support
            for(int i = 0; i < mReadContextCandidates.size(); ++i)
            {
                ReadContextCandidate candidate = mReadContextCandidates.get(i);

                if(candidate.FullMatch < MIN_SECOND_CANDIDATE_FULL_READS || candidate.FullMatch < topCandidateRcThreshold)
                    break;

                String coreStr = candidate.readContext().coreStr();
                if(coreStr.contains(topCore) || topCore.contains(coreStr))
                    continue;

                mSecondCandidate = new AltContext(RefContext, Ref, Alt, candidate, mRawSupportAlt, mRawBaseQualityAlt);
                break;
            }
        }

        mReadContextCandidates.clear();
    }

    public VariantReadContext readContext() { return mCandidate.readContext(); }

    @Override
    public String ref() { return Ref; }

    @Override
    public String alt() { return Alt; }

    @Override
    public String chromosome() { return RefContext.chromosome(); }

    @Override
    public int position() { return RefContext.position(); }

    public int rawAltSupport() { return mRawSupportAlt; }

    public int rawAltBaseQuality() { return mRawBaseQualityAlt; }

    @Override
    public boolean equals(@Nullable Object another)
    {
        if(this == another)
            return true;

        return another instanceof SimpleVariant && equalTo((SimpleVariant) another);
    }

    private boolean equalTo(final SimpleVariant other)
    {
        return Ref.equals(other.Ref) && Alt.equals(other.alt()) && Chromosome.equals(other.Chromosome) && Position == other.Position;
    }

    @Override
    public int hashCode()
    {
        int h = 5381;
        h += (h << 5) + Ref.hashCode();
        h += (h << 5) + Alt.hashCode();
        h += (h << 5) + chromosome().hashCode();
        h += (h << 5) + Ints.hashCode(position());
        return h;
    }

    public String toString()
    {
        return String.format("var(%s:%d %s->%s) readCandidates(%d)",
                chromosome(), position(), Ref, Alt, mReadContextCandidates != null ? mReadContextCandidates.size() : 0);
    }

    protected class ReadContextCandidate implements Comparable<ReadContextCandidate>
    {
        private final VariantReadContext mReadContext;

        public int FullMatch;
        public int CoreMatch;
        public int MinNumberOfEvents;

        ReadContextCandidate(int numberOfEvents, final VariantReadContext readContext)
        {
            mReadContext = readContext;
            MinNumberOfEvents = numberOfEvents;
        }

        public void incrementFull(int count, int numberOfEvents)
        {
            FullMatch += count;
            MinNumberOfEvents = min(MinNumberOfEvents, numberOfEvents);
        }

        public int minNumberOfEvents() { return MinNumberOfEvents; }

        public VariantReadContext readContext() { return mReadContext; }

        @Override
        public int compareTo(@NotNull final ReadContextCandidate other)
        {
            int fullCompare = -Integer.compare(FullMatch, other.FullMatch);

            if(fullCompare != 0)
                return fullCompare;

            return -Integer.compare(CoreMatch, other.CoreMatch);
        }

        public String toString()
        {
            return String.format("matches(full=%d core=%d)", FullMatch, CoreMatch);
        }
    }
}
