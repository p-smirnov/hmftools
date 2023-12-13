package com.hartwig.hmftools.errorprofile.repeat;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.SAMRecord;

// works with one repeat
// find the number of reads of each number of repeat across the section
public class RepeatAnalyser
{
    public static final Logger sLogger = LogManager.getLogger(RepeatAnalyser.class);

    final RefGenomeMicrosatellite refGenomeMicrosatellite;

    private final List<ReadRepeatMatch> mReadRepeatMatches = new ArrayList<>();

    public List<ReadRepeatMatch> getReadRepeatMatches() { return mReadRepeatMatches; }

    public List<ReadRepeatMatch> getPassingReadRepeatMatches()
    {
        return mReadRepeatMatches.stream().filter(o -> !o.shouldDropRead).collect(Collectors.toList());
    }

    public int numReadRejected()
    {
        return (int)mReadRepeatMatches.stream().filter(o -> o.shouldDropRead).count();
    }

    public RepeatAnalyser(final RefGenomeMicrosatellite refGenomeMicrosatellite)
    {
        this.refGenomeMicrosatellite = refGenomeMicrosatellite;
    }

    public void addReadToStats(final SAMRecord read)
    {
        if(read.getReadUnmappedFlag() || read.getDuplicateReadFlag())
            return;

        mReadRepeatMatches.add(ReadRepeatMatch.from(refGenomeMicrosatellite, read));
    }

    /*
    public void completeRegion(
            final AtomicLong readTagGenerator,
            final Consumer<ReadProfile> readProfileConsumer,
            final Consumer<ReadBaseSupport> readBaseSupportConsumer)
    {
        sLogger.info("region: {} read count: {}", , mReads.size());

        List<ReadProfile> readProfiles = new ArrayList<>();
        List<ReadBaseSupport> readBaseSupports = new ArrayList<>();

        // process all reads in this region
        for(SAMRecord read : mReads)
        {
            // we skip reads that are not entirely in this region
            if(read.getAlignmentStart() < mGenomeRegion.start() || read.getAlignmentStart() + read.getReadLength() > mGenomeRegion.end() + 1)
            {
                continue;
            }

            // also skip reads that matches ref genome

            long readTag = readTagGenerator.incrementAndGet();
            ReadBaseSupport readBaseSupport = calcReadBaseSupport(read, readTag);
            readBaseSupportConsumer.accept(readBaseSupport);
            readProfileConsumer.accept(ReadProfiler.profileRead(read, readTag));
        }

        sLogger.info("finished processing region: {}", mGenomeRegion);
    }

    // check if this read matches ref genome
    // TODO: write a ref genome cache to make it more generic
    public boolean isAllMatch(final SAMRecord record)
    {
        if (record.getCigar().numCigarElements() > 1 || record.getCigar().getFirstCigarElement().getOperator() != CigarOperator.M)
            return false;

        // this record has 151M, check it against the ref
        // we need to
        byte[] readBases = record.getReadBases();
        int startOffset = record.getAlignmentStart() - mGenomeRegion.start();
        int readOffset = Math.max(-startOffset, 0);
        int length = record.getReadLength();

        boolean allMatch = true;

        for(int i = 0; i < length; ++i)
        {
            byte refBase = mRefGenomeCache.getBase(record.getAlignmentStart() + i);
            if(refBase == N)
                continue;
            byte readBase = readBases[readOffset + i];
            if(readBase != N && readBase != refBase)
            {
                allMatch = false;
                break;
            }
        }

        if(allMatch)
        {
            // sLogger.debug("read all match ref, cigar={}", record.getCigarString());
        }
        return allMatch;
    }*/
}
