package com.hartwig.hmftools.errorprofile;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SAM_LOGGER;
import static com.hartwig.hmftools.errorprofile.ErrorProfileUtils.createPartitions;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.atomic.AtomicLong;
import java.util.function.Supplier;
import java.util.stream.Collectors;

import com.google.common.util.concurrent.UncheckedExecutionException;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.errorprofile.utils.RefGenomeCache;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

// For performance reasons, this class is not decoupled from the threading
// We want to avoid locks, queues as much as possible
public class ReadQualAnalyser
{
    private static final Logger sLogger = LogManager.getLogger(ReadQualAnalyser.class);

    /*
    static class StatsCollectorThread implements AsyncBamReader.IReadHandler
    {
        // open it once per thread
        final IndexedFastaSequenceFile refGenome;
        final AtomicLong readTagGenerator;

        // each region should only be processed by one thread
        // we combine them all back at the end
        final Map<ChrBaseRegion, GenomeRegionReadQualAnalyser> genomeRegionDataMap = new HashMap<>();

        BaseQualityBinCounter mBaseQualityBinCounter;

        public StatsCollectorThread(final String indexedRefGenomePath, final AtomicLong readTagGenerator, BaseQualityBinCounter baseQualityBinCounter)
        {
            try
            {
                this.refGenome = new IndexedFastaSequenceFile(new File(indexedRefGenomePath));
            }
            catch(FileNotFoundException e)
            {
                throw new RuntimeException(e);
            }
            this.readTagGenerator = readTagGenerator;
            this.mBaseQualityBinCounter = baseQualityBinCounter;
        }

        @Override
        public void processRead(SAMRecord read, ChrBaseRegion baseRegion)
        {
            if(read.getReadUnmappedFlag())
            {
                return;
            }

            GenomeRegionReadQualAnalyser analyser = getOrCreateGenomRegionData(baseRegion);
            analyser.addReadToStats(read);
        }

        @Override
        public void regionComplete(ChrBaseRegion baseRegion)
        {
            GenomeRegionReadQualAnalyser regionData = genomeRegionDataMap.get(baseRegion);
            if(regionData != null)
            {
                regionData.completeRegion(readTagGenerator, this::processReadProfile, this::processReadBaseSupport);

                // remove it to free up memory
                genomeRegionDataMap.remove(baseRegion);
            }
        }

        public GenomeRegionReadQualAnalyser getOrCreateGenomRegionData(ChrBaseRegion baseRegion)
        {
            return genomeRegionDataMap.computeIfAbsent(baseRegion, k -> new GenomeRegionReadQualAnalyser(k, new RefGenomeCache(refGenome, baseRegion.chromosome())));
        }

        public void processReadProfile(ReadProfile readProfile)
        {
            // write the read profiles
            //mReadProfileFileWriter.write(readProfile);
            mBaseQualityBinCounter.onReadProfile(readProfile);
        }

        public void processReadBaseSupport(ReadBaseSupport readBaseSupport)
        {
            // write the read bases supports
            //mReadBaseSupportFileWriter.write(readBaseSupport);

            mBaseQualityBinCounter.onReadBaseSupport(readBaseSupport);
        }
    }*/

    final AtomicLong mReadTagGenerator = new AtomicLong();

    BaseQualityBinCounter mBaseQualityBinCounter = new BaseQualityBinCounter();

    BaseQualityBinCounter getBaseQualityBinCounter() { return mBaseQualityBinCounter; }

    final Map<ChrBaseRegion, GenomeRegionReadQualAnalyser> genomeRegionDataMap = new ConcurrentHashMap<>();

    public ReadQualAnalyser()
    {
    }

    public void processBam(final ErrorProfileConfig config, ExecutorService executorService) throws InterruptedException
    {
        SamReaderFactory readerFactory = SamReaderFactory.make().validationStringency(config.BamStringency);
        if(config.RefGenomeFile != null)
        {
            readerFactory = readerFactory.referenceSource(new ReferenceSource(new File(config.RefGenomeFile)));
        }

        List<ChrBaseRegion> partitions = createPartitions(config);

        Random rand = new Random(0);

        // use sampling frac to filter the regions
        partitions = partitions.stream()
                .filter(o -> rand.nextDouble() <= config.SamplingFraction)
                .collect(Collectors.toList());

        List<IndexedFastaSequenceFile> refGenomeList = Collections.synchronizedList(new ArrayList<>());

        // one ref genome fasta file per thread
        final ThreadLocal<IndexedFastaSequenceFile> threadRefGenome = ThreadLocal.withInitial(() -> {
            IndexedFastaSequenceFile f = null;
            try
            {
                 f = new IndexedFastaSequenceFile(new File(config.RefGenomeFile));
            }
            catch(FileNotFoundException e)
            {
                throw new RuntimeException(e);
            }
            refGenomeList.add(f);
            return f;
        });

        BamSlicer bamSlicer = new BamSlicer(config.MinMappingQuality, false, false, false);
        CompletableFuture<Void> bamSliceTasks = bamSlicer.queryAsync(new File(config.BamPath), readerFactory, partitions,
                false, executorService,
                (read, chrBaseRegion) -> processRead(read, chrBaseRegion, threadRefGenome),
                this::regionComplete);

        // add another task to close all fasta files
        bamSliceTasks = bamSliceTasks.thenRun(() ->
        {
            for(IndexedFastaSequenceFile fastaSequenceFile : refGenomeList)
            {
                SAM_LOGGER.info("cleaning up");
                try
                {
                    fastaSequenceFile.close();
                }
                catch(IOException e)
                {
                    throw new UncheckedIOException(e);
                }
            }
        });

        try
        {
            // wait for all to complete
            bamSliceTasks.get();
        }
        catch(ExecutionException e)
        {
            throw new UncheckedExecutionException(e);
        }
    }

    public void processRead(SAMRecord read, ChrBaseRegion baseRegion, ThreadLocal<IndexedFastaSequenceFile> threadRefGenome)
    {
        if(read.getReadUnmappedFlag())
        {
            return;
        }

        GenomeRegionReadQualAnalyser analyser = getOrCreateGenomRegionData(baseRegion, threadRefGenome);
        analyser.addReadToStats(read);
    }

    public void regionComplete(ChrBaseRegion baseRegion)
    {
        GenomeRegionReadQualAnalyser regionData = genomeRegionDataMap.get(baseRegion);
        if(regionData != null)
        {
            regionData.completeRegion(mReadTagGenerator, this::processReadProfile, this::processReadBaseSupport);

            // remove it to free up memory
            genomeRegionDataMap.remove(baseRegion);
        }
    }

    public GenomeRegionReadQualAnalyser getOrCreateGenomRegionData(ChrBaseRegion baseRegion, ThreadLocal<IndexedFastaSequenceFile> threadRefGenome)
    {
        return genomeRegionDataMap.computeIfAbsent(baseRegion,
                k -> new GenomeRegionReadQualAnalyser(k, new RefGenomeCache(threadRefGenome.get(), baseRegion.chromosome())));
    }

    public void processReadProfile(ReadProfile readProfile)
    {
        // write the read profiles
        //mReadProfileFileWriter.write(readProfile);
        mBaseQualityBinCounter.onReadProfile(readProfile);
    }

    public void processReadBaseSupport(ReadBaseSupport readBaseSupport)
    {
        // write the read bases supports
        //mReadBaseSupportFileWriter.write(readBaseSupport);

        mBaseQualityBinCounter.onReadBaseSupport(readBaseSupport);
    }
}
