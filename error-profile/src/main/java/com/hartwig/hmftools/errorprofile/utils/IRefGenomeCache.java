package com.hartwig.hmftools.errorprofile.utils;

public interface IRefGenomeCache
{
    byte getBase(long refGenomePosition);

    // positions are 1 base and end inclusive
    byte[] getBases(long refStart, long refEndInclusive);

    String getBasesAsString(long refStart, long refEndInclusive);
}
