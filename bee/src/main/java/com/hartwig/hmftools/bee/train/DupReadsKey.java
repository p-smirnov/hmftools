package com.hartwig.hmftools.bee.train;

import java.util.Objects;

import org.jetbrains.annotations.NotNull;

// currently a duplicate group is identified by the
// UMI + extrapolated map start
public class DupReadsKey
{
    @NotNull
    public final String umi;
    public final int referenceIndex;
    public final int startPosition;

    public DupReadsKey(String umi, int referenceIndex, int startPosition)
    {
        this.umi = umi;
        this.referenceIndex = referenceIndex;
        this.startPosition = startPosition;
    }

    @Override
    public boolean equals(final Object o)
    {
        if (this == o)
        {
            return true;
        }
        if (o == null || getClass() != o.getClass())
        {
            return false;
        }
        final DupReadsKey that = (DupReadsKey) o;
        return referenceIndex == that.referenceIndex && startPosition == that.startPosition && umi.equals(that.umi);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(umi, referenceIndex, startPosition);
    }
}
