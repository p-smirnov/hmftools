package com.hartwig.hmftools.errorprofile;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.genome.region.Strand;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.apache.commons.lang3.tuple.ImmutableTriple;
import org.checkerframework.checker.units.qual.C;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ErrorProfileUtils
{
    public static int DEFAULT_PARTITION_SIZE = 1000000;

    public static Strand getReadStrand(SAMRecord read)
    {
        return read.getReadNegativeStrandFlag() ? Strand.REVERSE : Strand.FORWARD;
    }

    public static List<ChrBaseRegion> createPartitions(ErrorProfileConfig config)
    {
        List<SAMSequenceRecord> samSequences;
        try(SamReader samReader = openSamReader(config.BamPath, config.RefGenomeFile))
        {
            samSequences = samReader.getFileHeader().getSequenceDictionary().getSequences();
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }
        List<ChrBaseRegion> partitions = new ArrayList<>();
        int partitionSize = DEFAULT_PARTITION_SIZE;

        List<ChrBaseRegion> sequences;

        if(config.SpecificRegions.isEmpty())
        {
            sequences = new ArrayList<>();
            for(SAMSequenceRecord seq : samSequences)
            {
                String chrStr = seq.getSequenceName();
                int chromosomeLength = seq.getSequenceLength();
                sequences.add(new ChrBaseRegion(chrStr, 1, chromosomeLength));
            }
        }
        else
        {
            sequences = config.SpecificRegions;
        }

        for(ChrBaseRegion seq : sequences)
        {
            int chromosomeLength = seq.end();
            int startPos = seq.start();

            while(startPos < chromosomeLength)
            {
                int endPos = startPos + partitionSize - 1;
                if(endPos + partitionSize * 0.2 > chromosomeLength)
                {
                    endPos = chromosomeLength;
                }
                partitions.add(new ChrBaseRegion(seq.chromosome(), startPos, endPos));
                startPos = endPos + 1;
            }
        }

        return partitions;
    }

    public static SamReader openSamReader(@NotNull String bamFile, @Nullable String refGenomeFile)
    {
        SamReaderFactory factory = SamReaderFactory.makeDefault();

        if(refGenomeFile != null && !refGenomeFile.isEmpty())
        {
            factory = factory.referenceSequence(new File(refGenomeFile));
        }

        return factory.open(new File(bamFile));
    }

    public static @Nullable ImmutableTriple<Character, Integer, Integer> findHomopolymerLength(String seq)
    {
        char currentBase = '.';
        int currentStart = 0;

        char homopolymerBase = '.';
        int homopolymerStart = 0;
        int homopolymerLength = 0;

        for(int i = 0; i < seq.length(); i++)
        {
            if(seq.charAt(i) == currentBase)
            {
                int currentLen = i - currentStart + 1;
                if(currentLen > homopolymerLength)
                {
                    homopolymerBase = currentBase;
                    homopolymerStart = currentStart;
                    homopolymerLength = currentLen;
                }
            }
            else
            {
                currentBase = seq.charAt(i);
                currentStart = i;
            }
        }

        if(homopolymerLength >= ErrorProfileConstants.MIN_HOMOPOLYMER_LENGTH)
        {
            return ImmutableTriple.of(homopolymerBase, homopolymerStart, homopolymerStart + homopolymerLength);
        }
        else
        {
            return null;
        }
    }
}
