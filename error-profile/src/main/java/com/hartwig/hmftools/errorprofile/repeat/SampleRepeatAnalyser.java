package com.hartwig.hmftools.errorprofile.repeat;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.Random;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.stream.Collectors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ImmutableListMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multimaps;
import com.google.common.util.concurrent.UncheckedExecutionException;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.samtools.BamSlicer;

import org.apache.commons.lang3.Validate;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;

// For performance reasons, this class is not decoupled from the threading
// We want to avoid locks, queues as much as possible
public class SampleRepeatAnalyser
{
    private static final Logger sLogger = LogManager.getLogger(SampleRepeatAnalyser.class);

    private final List<RefGenomeMicrosatellite> mRefGenomeMicrosatellites;

    // for speed reasons we need to consolidate the chr base regions into bigger chunks
    private final Multimap<ChrBaseRegion, RepeatAnalyser> mRepeatAnalysers = ArrayListMultimap.create();

    public Collection<RepeatAnalyser> getRepeatAnalysers() { return mRepeatAnalysers.values(); }

    public SampleRepeatAnalyser(List<RefGenomeMicrosatellite> refGenomeMicrosatellites, double samplingFraction)
    {
        Random random = new Random(0);
        mRefGenomeMicrosatellites = refGenomeMicrosatellites.stream().filter(o -> random.nextDouble() <= samplingFraction).collect(Collectors.toList());
        partitionGenomeRepeats();
    }

    private void partitionGenomeRepeats()
    {
        final int PARTITION_SIZE = 1_000_000;

        mRepeatAnalysers.clear();

        ImmutableListMultimap<String, RefGenomeMicrosatellite> chromosomeRepeats = Multimaps.index(mRefGenomeMicrosatellites, RefGenomeMicrosatellite::chromosome);

        List<RepeatAnalyser> regionAnalysers = new ArrayList<>();

        for(String chromosome : chromosomeRepeats.keySet())
        {
            ChrBaseRegion currentRegion = null;

            List<RefGenomeMicrosatellite> refGenomeMicrosatellites = new ArrayList<>(chromosomeRepeats.get(chromosome));
            refGenomeMicrosatellites.sort(Comparator.comparing(RefGenomeMicrosatellite::referenceStart));

            for(RefGenomeMicrosatellite refGenomeMicrosatellite : refGenomeMicrosatellites)
            {
                if(currentRegion == null || currentRegion.baseLength() >= PARTITION_SIZE)
                {
                    if(currentRegion != null)
                    {
                        mRepeatAnalysers.putAll(currentRegion, regionAnalysers);
                    }

                    currentRegion = refGenomeMicrosatellite.genomeRegion.clone();
                    regionAnalysers.clear();
                }
                currentRegion.setEnd(refGenomeMicrosatellite.genomeRegion.end());
                regionAnalysers.add(new RepeatAnalyser(refGenomeMicrosatellite));
            }

            // final one
            if(currentRegion != null)
            {
                mRepeatAnalysers.putAll(currentRegion, regionAnalysers);
            }
        }
    }

    public void queryBam(final RepeatProfileConfig config, ExecutorService executorService) throws InterruptedException
    {
        SamReaderFactory readerFactory = SamReaderFactory.make().validationStringency(config.BamStringency);
        if(config.RefGenomeFile != null)
        {
            readerFactory = readerFactory.referenceSource(new ReferenceSource(new File(config.RefGenomeFile)));
        }

        Collection<ChrBaseRegion> partitions = mRepeatAnalysers.keySet().stream().sorted().collect(Collectors.toList());

        BamSlicer bamSlicer = new BamSlicer(config.MinMappingQuality, false, false, false);
        CompletableFuture<Void> bamSliceTasks = bamSlicer.queryAsync(new File(config.BamPath), readerFactory, partitions,
                false, executorService, this::processRead, this::regionComplete);
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

    public void processRead(SAMRecord read, ChrBaseRegion baseRegion)
    {
        if(read.getReadUnmappedFlag())
        {
            return;
        }

        Collection<RepeatAnalyser> repeatAnalysers = mRepeatAnalysers.get(baseRegion);

        Validate.isTrue(!repeatAnalysers.isEmpty());

        int readAlignmentStart = read.getAlignmentStart();
        int readAlignmentEnd = read.getAlignmentEnd();

        for(RepeatAnalyser analyser : repeatAnalysers)
        {
            if(BaseRegion.positionsWithin(analyser.refGenomeMicrosatellite.referenceStart(), analyser.refGenomeMicrosatellite.referenceEnd(),
                    readAlignmentStart, readAlignmentEnd))
            {
                analyser.addReadToStats(read);
            }
        }

        //GenomeRegionReadQualAnalyser analyser = getOrCreateGenomRegionData(baseRegion);
        //analyser.addReadToStats(read);


    }

    public void regionComplete(ChrBaseRegion baseRegion)
    {
        /*GenomeRegionReadQualAnalyser regionData = genomeRegionDataMap.get(baseRegion);
        if(regionData != null)
        {
            regionData.completeRegion(readTagGenerator, this::processReadProfile, this::processReadBaseSupport);

            // remove it to free up memory
            genomeRegionDataMap.remove(baseRegion);
        }*/
    }
}
