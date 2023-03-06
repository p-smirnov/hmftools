package com.hartwig.hmftools.bee.train;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;

// given the consensus read and consensus read cigar
// walk through the read and return the error state
// of each position
public class ReadPositionWalker
{
    enum TargetState
    {
        INSERT,
        DELETE,
        MATCH
    }

    Cigar mReadCigar;
    int mAlignStart;
    Cigar mConsensusReadCigar;
    int mConsensusReadAlignStart;

    int mCurrentPos = -1;

    CigarElement mCurrentReadCigarElement;
    CigarElement mCurrentConsensusCigarElement;

    public ReadPositionWalker(Cigar readCigar, int alignStart, Cigar consensusReadCigar, int consensusReadAlignStart)
    {
        mReadCigar = readCigar;
        mAlignStart = alignStart;
        mConsensusReadCigar = consensusReadCigar;
        mConsensusReadAlignStart = consensusReadAlignStart;
    }

}