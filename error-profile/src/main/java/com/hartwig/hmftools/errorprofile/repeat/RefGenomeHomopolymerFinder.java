package com.hartwig.hmftools.errorprofile.repeat;

import static htsjdk.samtools.util.SequenceUtil.N;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;

// simple class to find homopolymer regions in reference genome
public class RefGenomeHomopolymerFinder
{
    public static final Logger sLogger = LogManager.getLogger(RefGenomeHomopolymerFinder.class);

    public static List<RefGenomeHomopolymer> findHomopolymer(ReferenceSequenceFile referenceSequenceFile, int minRepeatLength)
    {
        int chunkSize = 100_000;
        return findHomopolymer(referenceSequenceFile, minRepeatLength, chunkSize);
    }

    static List<RefGenomeHomopolymer> findHomopolymer(ReferenceSequenceFile referenceSequenceFile, int minRepeatLength, int chunkSize)
    {
        List<RefGenomeHomopolymer> refGenomeHomopolymers = new ArrayList<>();

        for(SAMSequenceRecord sequenceRecord : referenceSequenceFile.getSequenceDictionary().getSequences())
        {
            if(!HumanChromosome.contains(sequenceRecord.getContig()))
            {
                continue;
            }

            int length = sequenceRecord.getSequenceLength();
            byte currentBase = N;
            int homopolymerLength = 1;

            for(int start = 1; start <= length; start = start + chunkSize)
            {
                int endInclusive = Math.min(start + chunkSize, length) - 1;
                byte[] seq = referenceSequenceFile.getSubsequenceAt(sequenceRecord.getContig(), start, endInclusive).getBases();

                for(int i = 0; i < seq.length; ++i)
                {
                    byte base = seq[i];
                    if(currentBase == base)
                    {
                        homopolymerLength++;
                    }
                    else
                    {
                        if(homopolymerLength >= minRepeatLength && currentBase != N)
                        {
                            int homopolymerStart = start + i - homopolymerLength;
                            int homopolymerEnd = start + i - 1;
                            RefGenomeHomopolymer refGenomeHomopolymer = new RefGenomeHomopolymer(sequenceRecord.getContig(), homopolymerStart,
                                    homopolymerEnd, currentBase, homopolymerLength);
                            refGenomeHomopolymers.add(refGenomeHomopolymer);
                        }

                        currentBase = base;
                        homopolymerLength = 1;
                    }
                }
            }
        }

        sLogger.info("found {} homopolymer regions in ref genome", refGenomeHomopolymers.size());

        return refGenomeHomopolymers;
    }
}
