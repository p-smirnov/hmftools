package com.hartwig.hmftools.errorprofile;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.genome.region.Strand;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.apache.commons.lang3.tuple.ImmutableTriple;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.lang3.tuple.Triple;
import org.apache.logging.log4j.util.BiConsumer;
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

    public static void processBamAsync(ErrorProfileConfig config, BiConsumer<SAMRecord, ChrBaseRegion> asyncRecordHandler,
            Consumer<ChrBaseRegion> regionCompleteHandler) throws InterruptedException
    {
        int numBamReaders = Math.max(config.Threads - 1, 1);
        List<ChrBaseRegion> partitions = createPartitions(config);

        SamReaderFactory readerFactory = SamReaderFactory.make().validationStringency(config.BamStringency);
        if(config.RefGenomeFile != null)
        {
            readerFactory = readerFactory.referenceSource(new ReferenceSource(new File(config.RefGenomeFile)));
        }

        AsyncBamReader.processBam(config.BamPath, readerFactory, partitions, asyncRecordHandler, regionCompleteHandler,
                numBamReaders, config.MinMappingQuality);
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

    // parse the flowcell, lane and tile for Illumina read id
    // @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>
    // see https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm
    @Nullable public static ImmutableTriple<String, Integer, Integer> parseFlowcellLaneTile(String readName)
    {
        // we want to find the 4th : and the 5th :
        int colonCount = 0;
        int thirdColonIndex = -1;
        int fourthColonIndex = -1;
        int fifthColonIndex = -1;
        for(int i = 0; i < readName.length(); ++i)
        {
            if(readName.charAt(i) == ':')
            {
                colonCount++;
                if(colonCount == 3)
                {
                    thirdColonIndex = i;
                }
                if(colonCount == 4)
                {
                    fourthColonIndex = i;
                }
                else if(colonCount == 5)
                {
                    fifthColonIndex = i;
                    break;
                }
            }
        }

        if(fifthColonIndex == -1)
        {
            return null;
        }

        String flowcell = readName.substring(0, fourthColonIndex);
        Integer lane = Integer.parseInt(readName.substring(thirdColonIndex + 1, fourthColonIndex));
        Integer tile = Integer.parseInt(readName.substring(fourthColonIndex + 1, fifthColonIndex));
        return ImmutableTriple.of(flowcell, lane, tile);
    }
}
