package com.hartwig.hmftools.esvee.caller.annotation;

import static java.lang.String.format;

import com.hartwig.hmftools.common.gripss.RepeatMaskData;

public class RepeatMaskAnnotation
{
    public final RepeatMaskData RmData;
    public final double Coverage;
    public final AlignmentData Alignment;

    public RepeatMaskAnnotation(final RepeatMaskData rmData, final double coverage, final AlignmentData alignment)
    {
        RmData = rmData;
        Coverage = coverage;
        Alignment = alignment;
    }

    public String toString() { return format("id(%d) type(%s, %s) coverage(%.3f) align(%s)",
            RmData.Id, RmData.ClassType, RmData.Repeat, Coverage, Alignment.toString()); }
}
