package com.hartwig.hmftools.errorprofile;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicLong;

import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.apache.commons.lang3.Validate;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class GenomeRegionStatsCollector
{
    private static final Logger sLogger = LogManager.getLogger(GenomeRegionStatsCollector.class);

    static class GenomeRegionData
    {
        GenomeRegionStats stats;
        List<SAMRecord> reads = new ArrayList<>();
        private final Thread expectedThread = Thread.currentThread();

        GenomeRegionData(ChrBaseRegion genomeRegion, IndexedFastaSequenceFile refGenome)
        {
            stats = new GenomeRegionStats(genomeRegion, refGenome);
        }

        void checkThread()
        {
            Validate.isTrue(Thread.currentThread() == expectedThread);
        }
    }

    final IndexedFastaSequenceFile refGenome;
    final ReadBaseSupportFileWriter readBaseSupportFileWriter;

    final ReadProfileFileWriter readProfileFileWriter;
    final Map<ChrBaseRegion, GenomeRegionData> genomeRegionDataMap = new ConcurrentHashMap<>();

    final AtomicLong readTagGenerator = new AtomicLong();

    public GenomeRegionStatsCollector(IndexedFastaSequenceFile refGenome,
            ReadBaseSupportFileWriter readBaseSupportFileWriter,
            ReadProfileFileWriter readProfileFileWriter)
    {
        this.refGenome = refGenome;
        this.readBaseSupportFileWriter = readBaseSupportFileWriter;
        this.readProfileFileWriter = readProfileFileWriter;
    }

    public void processRead(SAMRecord read, ChrBaseRegion baseRegion)
    {
        GenomeRegionData regionData = getOrCreateGenomRegionData(baseRegion);
        regionData.checkThread();
        regionData.reads.add(read);
        regionData.stats.addReadToStats(read);
    }

    public void regionComplete(ChrBaseRegion baseRegion)
    {
        GenomeRegionData regionData = genomeRegionDataMap.get(baseRegion);
        if(regionData != null)
        {
            regionData.checkThread();

            sLogger.info("region: {} read count: {}", baseRegion, regionData.reads.size());

            List<ReadProfile> readProfiles = new ArrayList<>();
            List<ReadBaseSupport> readBaseSupports = new ArrayList<>();

            // process all reads in this region
            for(SAMRecord read : regionData.reads)
            {
                // we skip reads that are not entirely in this region
                if(read.getAlignmentStart() < baseRegion.start() || read.getAlignmentStart() + read.getReadLength() > baseRegion.end() + 1)
                {
                    continue;
                }

                // also skip reads that matches ref genome
                /* if(regionData.stats.isAllMatch(read))
                {
                    continue;
                } */

                long readTag = readTagGenerator.incrementAndGet();
                ReadBaseSupport readBaseSupport = regionData.stats.calcReadBaseSupport(read, readTag);
                readBaseSupports.add(readBaseSupport);
                readProfiles.add(ReadProfiler.profileRead(read, readTag));
            }

            // write the read profiles
            readProfileFileWriter.write(readProfiles);

            // write the read bases supports
            readBaseSupportFileWriter.write(readBaseSupports);

            sLogger.info("region: {} finished writing {} reads out of {}", baseRegion, readBaseSupports.size(), regionData.reads.size());

            // remove it to free up memory
            genomeRegionDataMap.remove(baseRegion);
        }
    }

    public GenomeRegionData getOrCreateGenomRegionData(ChrBaseRegion baseRegion)
    {
        return genomeRegionDataMap.computeIfAbsent(baseRegion, k -> new GenomeRegionData(k, refGenome));
    }
}
