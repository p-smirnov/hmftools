package com.hartwig.hmftools.errorprofile;

import org.apache.logging.log4j.LogManager;

import java.util.*;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.LinkedBlockingDeque;

import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

public class BamProcessor
{
    /*
    private static final org.apache.logging.log4j.Logger logger = LogManager.getLogger(BamProcessor.class);

    public static void processBam(ErrorProfileConfig config) throws InterruptedException
    {
        int numBamReaders = Math.max(config.getThreadCount() - 1, 1);
        logger.info("processing bam: {} with {} bam reader threads", config.getBamFile(), numBamReaders);
        Queue<BamReader.Task> bamReaderTaskQ = new ConcurrentLinkedQueue<>();
        BlockingQueue<TelBamRecord> telBamRecordQ = new LinkedBlockingDeque<>();
        Set<String> incompleteReadNames = Collections.newSetFromMap(new ConcurrentHashMap<>());
        BamRecordWriter writer = new BamRecordWriter(config, telBamRecordQ, incompleteReadNames);
        writer.setProcessingMissingReadRegions(false);
        List<BamReader> bamReaders = new ArrayList<>();

        // create all the bam readers
        for(int i = 0; i < numBamReaders; i++)
        {
            BamReader bamReader = new BamReader(config, bamReaderTaskQ, telBamRecordQ, Collections.unmodifiableSet(incompleteReadNames));
            bamReaders.add(bamReader);
        }

        // add unmapped read first cause it is slower to process
        bamReaderTaskQ.add(BamReader.Task.fromQueryUnmapped());
        List<ChrBaseRegion> partitions = createPartitions(config);

        // put all into the queue
        for(ChrBaseRegion baseRegion : partitions)
        {
            bamReaderTaskQ.add(BamReader.Task.fromBaseRegion(baseRegion));
        }

        // start processing threads and run til completion
        runThreadsTillCompletion(telBamRecordQ, bamReaders, writer);
        logger.info("initial BAM file processing complete");

        // we want to process until no new reads have been accepted
        while(!writer.incompleteReadGroups.isEmpty())
        {
            int numAcceptedReads = writer.numAcceptedReads;

            // add unmapped read first cause it is slower to process
            bamReaderTaskQ.add(BamReader.Task.fromQueryUnmapped());

            // now we process the incomplete groups, since the threads have already finished there is no
            // concurrency problem
            List<ChrBaseRegion> missingReadRegions = getMissingReadRegions(writer.incompleteReadGroups, config.getSpecificChromosomes());
            if(!missingReadRegions.isEmpty())
            {
                logger.info("processing {} missing read regions", missingReadRegions.size());
                for(ChrBaseRegion baseRegion : missingReadRegions)
                {
                    bamReaderTaskQ.add(BamReader.Task.fromBaseRegion(baseRegion));
                }
            }
            writer.setProcessingMissingReadRegions(true);

            // start processing threads and run til completion
            runThreadsTillCompletion(telBamRecordQ, bamReaders, writer);
            if(writer.numAcceptedReads == numAcceptedReads)
            {
                // no change, so we did not find any more reads
                break;
            }
        }
        writer.finish();
    }

    // run the bam reader and record writer threads till completion
    private static void runThreadsTillCompletion(Queue<TelBamRecord> telBamRecordQ, List<BamReader> bamReaders, BamRecordWriter writer)
            throws InterruptedException
    {
        List<Thread> bamReaderThreads = new ArrayList<>();
        for(BamReader bamReader : bamReaders)
        {
            Thread thread = new Thread(bamReader);
            bamReaderThreads.add(thread);
            thread.start();
        }

        Thread writerThread = new Thread(writer);
        writerThread.start();

        // wait for reader to finish
        for(Thread thread : bamReaderThreads)
        {
            thread.join();
        }
        logger.info("{} bam reader threads finished", bamReaderThreads.size());

        // now tell writer to finish also
        TelBamRecord telBamRecord = new TelBamRecord();
        telBamRecord.setPoison(true);
        telBamRecordQ.add(telBamRecord);
        writerThread.join();
        logger.info("writer thread finished");
    }

    private static List<ChrBaseRegion> getMissingReadRegions(Map<String, ReadGroup> incompleteReadGroups, List<String> specificChromosomes)
    {
        List<ChrBaseRegion> missingReadRegions = new ArrayList<>();
        for(ReadGroup readGroup : incompleteReadGroups.values())
        {
            assert readGroup.invariant();
            missingReadRegions.addAll(readGroup.findMissingReadBaseRegions());
            if(!specificChromosomes.isEmpty())
            {
                List<ChrBaseRegion> list = new ArrayList<>();
                for(ChrBaseRegion baseRegion : missingReadRegions)
                {
                    if(specificChromosomes.contains(baseRegion.chromosome()))
                    {
                        list.add(baseRegion);
                    }
                }
                missingReadRegions = list;
            }
        }

        // sort the mate regions
        missingReadRegions.sort(Comparator.naturalOrder());
        if(!missingReadRegions.isEmpty())
        {
            List<ChrBaseRegion> compactBaseRegions = new ArrayList<>();

            // now we go through the mate regions and create a new condense list
            ChrBaseRegion overlapRegion = null;
            for(ChrBaseRegion br : missingReadRegions)
            {
                // see if still overlap
                if(overlapRegion != null && br.chromosome().equals(overlapRegion.chromosome()))
                {
                    assert br.start() >= overlapRegion.start();
                    if(overlapRegion.end() >= br.start())
                    {
                        // still overlapping, we can set it
                        overlapRegion.setStart(br.start());
                        continue;
                    }
                }
                // no longer can reuse old one, make a new one
                overlapRegion = (ChrBaseRegion) br.clone();
                compactBaseRegions.add(overlapRegion);
            }
            missingReadRegions = compactBaseRegions;
        }
        return missingReadRegions;
    }
     */
}
