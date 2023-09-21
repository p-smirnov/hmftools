package com.hartwig.hmftools.errorprofile;

/*
Read_Id
Chr
Pos
GC content (0-151)
Avg base qual
# of low qual bases ( base qual < 30)
# of 3’ soft clipped bases
# of 5’ soft clipped bases
Total bases in homopolymer / dinucleotide repeats wth repeat length > 4
Longest homopolymer length
Longest homopolymer base {A,G,C,T}
# of bases 5’ side of longest homopolymer
# of bases 3’ side of longest homopolymer
# of low qual bases on 5’ side of longest homopolymer
# of low qual bases on 3’ side of longest homopolymer
Longest dinucleotide repeat length
Longest dinucleotide repeat bases {AC,AG,AT,CG,CT,GT}
# of bases 5’ side of longest  dinucleotide repeat
# of bases 3’ side of longest  dinucleotide repeat
# of low qual bases on 5’ side of longest  dinucleotide repeat
# of low qual bases on 3’ side of longest  dinucleotide repeat
 */

import com.hartwig.hmftools.common.genome.region.Strand;

// error profile of each read
public class ReadProfile
{
    public long readTag;
    public String readId;
    public boolean firstOfPair;
    public String chromosome;
    public int position;
    public Strand strand;
    public int readLength;
    public String cigar;
    public int mapQ;
    public double gcContent;
    public double averageBaseQual;
    public int numLowQualBases;
    public int insertSize = 0;
    public int num5PrimeSoftClipped = 0;
    public int num3PrimeSoftClipped = 0;
    public int polyGLength = 0;
    public char homopolymerBase = 'N';
    public int homopolymerStart = -1;
    public int homopolymerEnd = -1;
    public int numLowQualBeforeHomopolymer = 0;
    public int numLowQualAfterHomopolymer = 0;
    public TandemRepeat tandemRepeat = null;
    public int numLowQualBeforeTandemRepeat = 0;
    public int numLowQualAfterTandemRepeat = 0;

    public boolean hasHomopolymer() { return homopolymerStart != -1; }
    public int getHomopolymerLength() { return homopolymerEnd - homopolymerStart; }
}
