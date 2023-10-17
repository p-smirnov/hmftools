package com.hartwig.hmftools.errorprofile;

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

public class ReadQualAnalyser
{
    private static final Logger sLogger = LogManager.getLogger(ReadQualAnalyser.class);

    static class GenomeRegionData
    {
        GenomeRegionReadQualAnalyser analyser;
        private final Thread expectedThread = Thread.currentThread();

        GenomeRegionData(ChrBaseRegion genomeRegion, IndexedFastaSequenceFile refGenome)
        {
            analyser = new GenomeRegionReadQualAnalyser(genomeRegion, new RefGenomeCache(refGenome, genomeRegion.chromosome()));
        }

        void checkThread()
        {
            Validate.isTrue(Thread.currentThread() == expectedThread);
        }
    }

    final IndexedFastaSequenceFile refGenome;
    final Consumer<ReadBaseSupport> readBaseSupportConsumer;

    final Consumer<ReadProfile> readProfileConsumer;
    final Map<ChrBaseRegion, GenomeRegionData> genomeRegionDataMap = new ConcurrentHashMap<>();

    final AtomicLong readTagGenerator = new AtomicLong();

    public ReadQualAnalyser(IndexedFastaSequenceFile refGenome,
            Consumer<ReadProfile> readProfileConsumer,
            Consumer<ReadBaseSupport> readBaseSupportConsumer)
    {
        this.refGenome = refGenome;
        this.readBaseSupportConsumer = readBaseSupportConsumer;
        this.readProfileConsumer = readProfileConsumer;
    }

    public void processBam(final ErrorProfileConfig config) throws InterruptedException
    {
        ErrorProfileUtils.processBamAsync(config, this::processRead, this::regionComplete);
    }

    private void processRead(SAMRecord read, ChrBaseRegion baseRegion)
    {
        if(read.getReadUnmappedFlag())
        {
            return;
        }

        GenomeRegionData regionData = getOrCreateGenomRegionData(baseRegion);
        regionData.checkThread();
        regionData.analyser.addReadToStats(read);
    }

    private void regionComplete(ChrBaseRegion baseRegion)
    {
        GenomeRegionData regionData = genomeRegionDataMap.get(baseRegion);
        if(regionData != null)
        {
            regionData.checkThread();
            regionData.analyser.completeRegion(readTagGenerator, readProfileConsumer, readBaseSupportConsumer);

            // remove it to free up memory
            genomeRegionDataMap.remove(baseRegion);
        }
    }

    public GenomeRegionData getOrCreateGenomRegionData(ChrBaseRegion baseRegion)
    {
        return genomeRegionDataMap.computeIfAbsent(baseRegion, k -> new GenomeRegionData(k, refGenome));
    }
}
