package com.hartwig.hmftools.bamtools.compare;

import java.util.List;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;

public class DiffRecord
{
    public final SAMRecord Read;
    public final MismatchType Type;
    public final List<String> DiffList;

    public DiffRecord(final SAMRecord read, final MismatchType type, @Nullable final List<String> diffList)
    {
        Read = read;
        Type = type;
        DiffList = diffList;
    }
}
