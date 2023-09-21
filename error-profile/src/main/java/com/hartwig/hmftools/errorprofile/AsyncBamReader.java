package com.hartwig.hmftools.errorprofile;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.BiConsumer;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class AsyncBamReader
{
    private static final Logger logger = LogManager.getLogger(AsyncBamReader.class);

    static class BamReaderThread extends Thread
    {
        final Queue<ChrBaseRegion> mTaskQ;
        final SamReader mSamReader;

        private final BamSlicer mBamSlicer;

        private final int mMinMappingQuality;
        BiConsumer<SAMRecord, ChrBaseRegion> mConsumer;
        Consumer<ChrBaseRegion> mRegionComplete;

        BamReaderThread(final String bamFile, final SamReaderFactory samReaderFactory, final Queue<ChrBaseRegion> inTaskQ,
                BiConsumer<SAMRecord, ChrBaseRegion> consumer, Consumer<ChrBaseRegion> regionComplete,
                int minMappingQuality)
        {
            mTaskQ = inTaskQ;
            mSamReader = samReaderFactory.open(new File(bamFile));
            // mMinMappingQuality = minMappingQuality;
            mConsumer = consumer;
            mRegionComplete = regionComplete;
            mMinMappingQuality = minMappingQuality;
            mBamSlicer = new BamSlicer(0, false, false, false);
        }

        @Override
        public void run()
        {
            logger.debug("bam reader thread start");

            while (true)
            {
                ChrBaseRegion region;
                try
                {
                    region = mTaskQ.remove();
                }
                catch (NoSuchElementException e)
                {
                    // finished processing
                    break;
                }

                try (final SAMRecordIterator iterator = mSamReader.queryOverlapping(region.chromosome(), region.start(),
                        region.end()))
                {
                    while (iterator.hasNext())
                    {
                        final SAMRecord record = iterator.next();

                        if (!passesFilters(record))
                        {
                            continue;
                        }

                        int alignmentEnd = record.getAlignmentEnd();

                        // check overlap
                        if (region.end() >= record.getAlignmentStart() && region.start() <= alignmentEnd)
                        {
                            mConsumer.accept(record, region);
                        }
                    }

                    mRegionComplete.accept(region);
                }
            }

            try {
                mSamReader.close();
            }
            catch (IOException e)
            {
                logger.error("IO exception in SamReader::close: {}", e.getMessage());
            }

            logger.debug("bam reader thread finish");
        }

        private boolean passesFilters(final SAMRecord record)
        {
            if(record.getMappingQuality() < mMinMappingQuality || record.getReadUnmappedFlag())
                return false;

            if(record.isSecondaryAlignment())
                return false;

            if(record.getSupplementaryAlignmentFlag())
                return false;

            return !record.getDuplicateReadFlag();
        }
    }

    public static void processBam(final String bamFile,
            final SamReaderFactory samReaderFactory,
            List<ChrBaseRegion> chrBaseRegions,
            BiConsumer<SAMRecord, ChrBaseRegion> asyncRecordHandler,
            Consumer<ChrBaseRegion> regionCompleteHandler,
            int threadCount,
            int minMappingQuality)
            throws InterruptedException
    {
        logger.debug("Processing {} regions in bam {}", chrBaseRegions.size(), bamFile);

        final Queue<ChrBaseRegion> taskQ = new ConcurrentLinkedQueue<>(chrBaseRegions);

        // we create the consumer and producer
        var bamReaders = new ArrayList<BamReaderThread>();

        for (int i = 0; i < Math.max(threadCount, 1); ++i)
        {
            var t = new BamReaderThread(bamFile, samReaderFactory, taskQ, asyncRecordHandler, regionCompleteHandler, minMappingQuality);
            t.setName(String.format("worker-%d", i));
            t.start();
            bamReaders.add(t);
        }

        logger.info("{} bam reader threads started", bamReaders.size());

        //AmberTaskCompletion taskCompletion = new AmberTaskCompletion(taskQ.size());
        for (BamReaderThread t : bamReaders)
        {
            while (t.isAlive())
            {
                // check status every 30 seconds
                t.join(30_000);

                // check status
                //taskCompletion.progress(taskQ.size());
            }
        }

        logger.info("{} bam reader threads finished", bamReaders.size());
    }
}
