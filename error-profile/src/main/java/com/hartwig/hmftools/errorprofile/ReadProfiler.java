package com.hartwig.hmftools.errorprofile;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.genome.gc.GcCalcs;
import com.hartwig.hmftools.common.genome.region.Strand;
import com.hartwig.hmftools.common.samtools.CigarUtils;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.tuple.Triple;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;

///
// Starting point:
//
//What % of fragments have soft clips and are not adapters
//
//Need a definition of:
//
//Poor alignment
//
//PolyG
//
//Likely phasing error (has to be 3' AND some low qual)
//
//After we mark all those, what proportion is left, and start to manually curate
//
//
public class ReadProfiler
{
    public static final Logger sLogger = LogManager.getLogger(ReadProfiler.class);

    // see https://teichlab.github.io/scg_lib_structs/methods_html/Illumina.html
    // we use Truseq?
    public static final String read1Adaptor = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT";
    public static final String read2Adaptor = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";

    ErrorProfileStats mStats = new ErrorProfileStats();

    final List<ReadProfile> readProfiles = new ArrayList<>();

    // count of how many G / Cs at the end
    final Map<Integer, Integer> mPolyGCounts = new HashMap<>();

    // homopolymer lengths
    // key is the base and the length, value is the count
    final Map<Map.Entry<Character, Integer>, Integer> mHomopolymerLengthCounts = new HashMap<>();

    // tandem repeats
    final Map<Map.Entry<String, Integer>, Integer> mTandemRepeatCounts = new HashMap<>();

    public Map<Integer, Integer> getPolyGCounts()
    {
        return mPolyGCounts;
    }

    public Map<Map.Entry<Character, Integer>, Integer> getHomopolymerLengthCounts()
    {
        return mHomopolymerLengthCounts;
    }

    public Map<Map.Entry<String, Integer>, Integer> getTandemRepeatCounts()
    {
        return mTandemRepeatCounts;
    }

    // we only keep reads that have homopolymer, tandem repeat, short fragment, or soft clips
    // the rest are not kept to keep the data smaller
    public static ReadProfile profileRead(SAMRecord read, long readTag)
    {
        // mStats.numReads++;

        String readString = read.getReadString();
        byte[] baseQualities = read.getBaseQualities();
        if(read.getReadNegativeStrandFlag())
        {
            readString = SequenceUtil.reverseComplement(readString);
            ArrayUtils.reverse(baseQualities);
        }

        ReadProfile readProfile = new ReadProfile();

        readProfile.readTag = readTag;
        readProfile.readId = read.getReadName();
        readProfile.firstOfPair = read.getFirstOfPairFlag();
        readProfile.chromosome = read.getReferenceName();
        readProfile.position = read.getAlignmentStart();
        readProfile.strand = read.getReadNegativeStrandFlag() ? Strand.REVERSE : Strand.FORWARD;
        readProfile.readLength = read.getReadLength();
        readProfile.cigar = read.getCigarString();
        readProfile.mapQ = read.getMappingQuality();
        readProfile.gcContent = GcCalcs.calcGcPercent(readString);
        readProfile.averageBaseQual = averageBaseQual(baseQualities);
        readProfile.numLowQualBases = numLowQualBases(baseQualities);

        readProfile.polyGLength = polyGTailLength(read);

        // add to the poly g stats
        // mPolyGCounts.compute(polyGCount, (k, v) -> (v == null) ? 1 : v + 1);

        Triple<Character, Integer, Integer> homopolymer = ErrorProfileUtils.findHomopolymerLength(readString);

        if(homopolymer != null)
        {
            readProfile.homopolymerBase = homopolymer.getLeft();
            readProfile.homopolymerStart = homopolymer.getMiddle();
            readProfile.homopolymerEnd = homopolymer.getRight();

            readProfile.numLowQualBeforeHomopolymer = numLowQualBases(baseQualities, 0, readProfile.homopolymerStart);
            readProfile.numLowQualAfterHomopolymer = numLowQualBases(baseQualities, readProfile.homopolymerEnd, baseQualities.length);
        }

        @Nullable TandemRepeat tandemRepeat = TandemRepeatFinder.findTandemRepeat(readString);

        if(tandemRepeat != null)
        {
            readProfile.tandemRepeat = tandemRepeat;
            readProfile.numLowQualBeforeTandemRepeat = numLowQualBases(baseQualities, 0, tandemRepeat.startIndex);
            readProfile.numLowQualAfterTandemRepeat = numLowQualBases(baseQualities, tandemRepeat.endIndex, baseQualities.length);
        }

        readProfile.insertSize = read.getInferredInsertSize();

        // check against the adaptor
        // check insert size
        if(read.getInferredInsertSize() != 0 && Math.abs(read.getInferredInsertSize()) < read.getReadLength())
        {
            sLogger.trace("short fragment read: {} insert size: {} contains adaptor: {}",
                    read, read.getInferredInsertSize(), containsAdaptor(read));
        }

        // work out which direction the split is at
        if(read.getReadNegativeStrandFlag())
        {
            readProfile.num5PrimeSoftClipped = CigarUtils.rightSoftClipLength(read);
            readProfile.num3PrimeSoftClipped = CigarUtils.leftSoftClipLength(read);
        }
        else
        {
            readProfile.num5PrimeSoftClipped = CigarUtils.leftSoftClipLength(read);
            readProfile.num3PrimeSoftClipped = CigarUtils.rightSoftClipLength(read);
        }

        return readProfile;
    }

    public static boolean shouldKeepErrorProfile(ReadProfile errorProfile)
    {
        return true;
        // return errorProfile.isShortFragment || errorProfile.is3PrimeSoftClipped || errorProfile.is5PrimeSoftClipped ||
           //     errorProfile.polyGLength > 0 || errorProfile.getHomopolymerLength() > 0 || errorProfile.tandemRepeat != null;
    }

    public static boolean detectPhasingError(SAMRecord read)
    {
        // what is a phasing error?
        // a phasing error is when some fragments in the flowcell become out of sync
        // causing downstream bases to have much lower base qual
        // this most frequently is caused by homopolymer or tandem repeats
        return false;
    }

    public static boolean containsAdaptor(SAMRecord read)
    {
        int adaptorIndex = -1;

        String seq = read.getReadString();

        if(read.getReadNegativeStrandFlag())
        {
            // do a reverse complement
            seq = SequenceUtil.reverseComplement(seq);
        }

        if(read.getFirstOfPairFlag())
        {
            // try the R1 adatpor
            adaptorIndex = seq.indexOf(read1Adaptor);
        }
        else
        {
            // try the R2 adatpor
            adaptorIndex = seq.indexOf(read2Adaptor);
        }

        return adaptorIndex != -1;
    }

    public static int polyGTailLength(SAMRecord read)
    {
        if(read.getReadNegativeStrandFlag())
        {
            return numLeadingPolyC(read.getReadString());
        }
        else
        {
            return numTrailingPolyG(read.getReadString());
        }
    }

    public static int numTrailingPolyG(String seq)
    {
        for(int i = 0; i < seq.length(); i++)
        {
            if(seq.charAt(seq.length() - i - 1) != 'G')
            {
                return i;
            }
        }
        return seq.length();
    }

    public static int numLeadingPolyC(String seq)
    {
        for(int i = 0; i < seq.length(); i++)
        {
            if(seq.charAt(i) != 'C')
            {
                return i;
            }
        }
        return seq.length();
    }

    public static double averageBaseQual(byte[] baseQualities)
    {
        double sumBaseQual = 0;
        for(byte bq : baseQualities)
        {
            sumBaseQual += bq;
        }
        return sumBaseQual / baseQualities.length;
    }

    public static int numLowQualBases(byte[] baseQualities)
    {
        return numLowQualBases(baseQualities, 0, baseQualities.length);
    }

    public static int numLowQualBases(byte[] baseQualities, int startIndex, int endIndex)
    {
        int count = 0;
        for(int i = startIndex; i < endIndex; ++i)
        {
            if(baseQualities[i] < ErrorProfileConstants.BASE_QUAL_CUTOFF)
            {
                ++count;
            }
        }
        return count;
    }
}
