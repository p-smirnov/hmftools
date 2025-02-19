package com.hartwig.hmftools.sage.candidate;

import com.hartwig.hmftools.sage.common.VariantReadContext;

import org.jetbrains.annotations.Nullable;

public class AltRead
{
    public final String Ref;
    public final String Alt;
    public final int BaseQuality;
    public final int NumberOfEvents;
    public final boolean SufficientMapQuality;

    private final RefContext mRefContext;

    @Nullable
    private VariantReadContext mReadContext;

    public AltRead(
            final RefContext refContext, final String ref, final String alt, final int baseQuality, final int numberOfEvents,
            final boolean sufficientMapQuality, final VariantReadContext readContext)
    {
        mRefContext = refContext;
        Ref = ref;
        Alt = alt;
        BaseQuality = baseQuality;
        NumberOfEvents = numberOfEvents;
        SufficientMapQuality = sufficientMapQuality;

        mReadContext = readContext;
    }

    public int position()
    {
        return mRefContext.position();
    }

    public boolean isIndel()
    {
        return Ref.length() != Alt.length();
    }

    public int length()
    {
        return Math.abs(Ref.length() - Alt.length());
    }

    public void updateRefContext()
    {
        mRefContext.processAltRead(Ref, Alt, BaseQuality, NumberOfEvents, mReadContext);
    }

    public String toString() { return String.format("%d: %s>%s", position(), Ref, Alt); }
}
