package com.hartwig.hmftools.errorprofile.microsatellite;

import static java.time.format.DateTimeFormatter.ISO_ZONED_DATE_TIME;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;

import java.io.File;
import java.io.FileNotFoundException;
import java.time.Duration;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadFactory;
import java.util.stream.Collectors;

import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.errorprofile.ErrorProfileConfig;

import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;

public class MicrosatelliteAnalyserApp
{
    public static final Logger sLogger = LogManager.getLogger(MicrosatelliteAnalyserApp.class);

    private final MicrosatelliteAnalyserConfig mConfig;

    public MicrosatelliteAnalyserApp(final ConfigBuilder configBuilder) throws ParseException
    {
        mConfig = new MicrosatelliteAnalyserConfig(configBuilder);
    }

    public int run(final String... args) throws InterruptedException, FileNotFoundException
    {
        Instant start = Instant.now();

        VersionInfo versionInfo = new VersionInfo("errorprofile.version");

        sLogger.info("ErrorProfile version: {}", versionInfo.version());

        sLogger.debug("build timestamp: {}, run args: {}",
                versionInfo.buildTime().format(ISO_ZONED_DATE_TIME), String.join(" ", args));

        if(!mConfig.isValid())
        {
            sLogger.error(" invalid config, exiting");
            return 1;
        }

        List<RefGenomeMicrosatellite> refGenomeMicrosatellites = RefGenomeMicrosatelliteFile.read(mConfig.RefGenomeMicrosatelliteFile);

        sLogger.info("loaded {} microsatellites regions", refGenomeMicrosatellites.size());

        filterSpecificRegions(refGenomeMicrosatellites);

        SampleBamProcessor sampleBamProcessor = new SampleBamProcessor(refGenomeMicrosatellites, mConfig.SamplingFraction);

        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("worker-%d").build();
        ExecutorService executorService = Executors.newFixedThreadPool(mConfig.Threads, namedThreadFactory);

        sampleBamProcessor.queryBam(mConfig, executorService);

        // now write out all the repeat stats
        MicrosatelliteSiteFile.write(MicrosatelliteSiteFile.generateFilename(mConfig.OutputDir, mConfig.SampleId),
                sampleBamProcessor.getRepeatAnalysers());

        writeMicrosatelliteStatsTable(sampleBamProcessor.getRepeatAnalysers());

        // draw a chart of the 9 ms profiles


        Instant finish = Instant.now();
        long seconds = Duration.between(start, finish).getSeconds();
        sLogger.info("run complete. Time taken: {}m {}s", seconds / 60, seconds % 60);

        return 0;
    }

    private void writeMicrosatelliteStatsTable(@NotNull final Collection<MicrosatelliteSiteAnalyser> microsatelliteSiteAnalysers)
    {
        // write two tables, one with real variant filter, one without

        List<MicrosatelliteStatsTable> msStatsTables = new ArrayList<>();

        for(MicrosatelliteSelector s : createMicrosatelliteSelectors())
        {
            MicrosatelliteStatsTable msStatsTable = new MicrosatelliteStatsTable(s.unitName());
            msStatsTable.summarise(microsatelliteSiteAnalysers.stream().filter(s::select).collect(Collectors.toList()));
            msStatsTables.add(msStatsTable);
        }

        MicrosatelliteStatsTableFile.write(MicrosatelliteStatsTableFile.generateFilename(mConfig.OutputDir, mConfig.SampleId),
                msStatsTables);
    }

    private static List<MicrosatelliteSelector> createMicrosatelliteSelectors()
    {
        // create nine summary / pivot table
        // {A/T, C/G, AT/TA, AG/GA/CT/TC, AC/CA/GT/TG, CG/GC, any 3 base, any 4 base, any 5 base}
        List<MicrosatelliteSelector> selectors = new ArrayList<>();
        selectors.add(new MicrosatelliteSelector(List.of("A", "T"), null));
        selectors.add(new MicrosatelliteSelector(List.of("C", "G"), null));
        selectors.add(new MicrosatelliteSelector(List.of("AT", "TA"), null));
        selectors.add(new MicrosatelliteSelector(List.of("AG", "GA", "CT", "TC"), null));
        selectors.add(new MicrosatelliteSelector(List.of("AC", "CA", "GT", "TG"), null));
        selectors.add(new MicrosatelliteSelector(List.of("CG", "GC"), null));
        selectors.add(new MicrosatelliteSelector(null, 3));
        selectors.add(new MicrosatelliteSelector(null, 4));
        selectors.add(new MicrosatelliteSelector(null, 5));

        return selectors;
    }

    private void filterSpecificRegions(List<RefGenomeMicrosatellite> refGenomeMicrosatellites)
    {
        /*
        if(!mConfig.SpecificRegions.isEmpty())
        {
            sLogger.info(mConfig.SpecificRegions);

            refGenomeMicrosatellites.removeIf(refGenomeHomopolymer -> mConfig.SpecificRegions.stream()
                    .noneMatch(o -> o.overlaps(refGenomeHomopolymer.genomeRegion)));
        }
         */

        // we want to get around 1000 sites for each repeat context
    }

    private static SamReaderFactory readerFactory(final ErrorProfileConfig config)
    {
        final SamReaderFactory readerFactory = SamReaderFactory.make().validationStringency(config.BamStringency);
        if(config.RefGenomeFile != null)
        {
            return readerFactory.referenceSource(new ReferenceSource(new File(config.RefGenomeFile)));
        }
        return readerFactory;
    }

    public static void main(final String... args) throws InterruptedException, FileNotFoundException, ParseException
    {
        ConfigBuilder configBuilder = new ConfigBuilder("ErrorProfile");
        MicrosatelliteAnalyserConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        setLogLevel(configBuilder);

        // set all thread exception handler
        // if we do not do this, exception thrown in other threads will not be handled and results
        // in the program hanging
        Thread.setDefaultUncaughtExceptionHandler((Thread t, Throwable e) ->
        {
            sLogger.error("[{}]: uncaught exception: {}", t, e);
            e.printStackTrace(System.err);
            System.exit(1);
        });

        System.exit(new MicrosatelliteAnalyserApp(configBuilder).run(args));
    }
}
