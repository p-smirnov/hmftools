package com.hartwig.hmftools.errorprofile;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicLong;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.errorprofile.utils.RefGenomeCache;

import org.apache.commons.lang3.Validate;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

// For performance reasons, this class is not decoupled from the threading
// We want to avoid locks, queues as much as possible
public class ReadQualAnalyser
{
    private static final Logger sLogger = LogManager.getLogger(ReadQualAnalyser.class);

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
    }

    final List<StatsCollectorThread> mThreadHandlers = new ArrayList<>();

    final AtomicLong mReadTagGenerator = new AtomicLong();

    BaseQualityBinCounter mBaseQualityBinCounter = new BaseQualityBinCounter();

    BaseQualityBinCounter getBaseQualityBinCounter() { return mBaseQualityBinCounter; }

    public ReadQualAnalyser()
    {
    }

    public void processBam(final ErrorProfileConfig config) throws InterruptedException
    {
        ErrorProfileUtils.processBamAsync(config, () -> createThreadReadhandler(config.RefGenomeFile));
    }

    private StatsCollectorThread createThreadReadhandler(final String refGenomeFilePath)
    {
        StatsCollectorThread threadReadHandler = new StatsCollectorThread(refGenomeFilePath, mReadTagGenerator, mBaseQualityBinCounter);
        mThreadHandlers.add(threadReadHandler);
        return threadReadHandler;
    }
}
