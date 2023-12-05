package com.hartwig.hmftools.errorprofile.repeat;

import static java.time.format.DateTimeFormatter.ISO_ZONED_DATE_TIME;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.errorprofile.repeat.RepeatProfileConstant.MIN_REPEAT_LENGTH;

import java.io.File;
import java.io.FileNotFoundException;
import java.time.Duration;
import java.time.Instant;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadFactory;

import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.errorprofile.ErrorProfileConfig;

import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class RepeatProfileApplication
{
    public static final Logger sLogger = LogManager.getLogger(RepeatProfileApplication.class);

    private final RepeatProfileConfig mConfig;

    public RepeatProfileApplication(final ConfigBuilder configBuilder) throws ParseException
    {
        mConfig = new RepeatProfileConfig(configBuilder);
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

        // find all the polymers
        IndexedFastaSequenceFile refGenome = new IndexedFastaSequenceFile(new File(mConfig.RefGenomeFile));
        List<RefGenomeHomopolymer> refGenomeHomopolymers = RefGenomeHomopolymerFinder.findHomopolymer(refGenome, MIN_REPEAT_LENGTH);

        filterSpecificRegions(refGenomeHomopolymers);

        SampleRepeatAnalyser sampleRepeatAnalyser = new SampleRepeatAnalyser(refGenomeHomopolymers, mConfig.SamplingFraction);

        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("worker-%d").build();
        ExecutorService executorService = Executors.newFixedThreadPool(mConfig.Threads, namedThreadFactory);

        sampleRepeatAnalyser.queryBam(mConfig, executorService);

        // now write out all the repeat stats
        RepeatProfileFile.write(RepeatProfileFile.generateFilename(mConfig.OutputDir, mConfig.SampleId),
                sampleRepeatAnalyser.getRepeatAnalysers());

        Instant finish = Instant.now();
        long seconds = Duration.between(start, finish).getSeconds();
        sLogger.info("run complete. Time taken: {}m {}s", seconds / 60, seconds % 60);

        return 0;
    }

    private void filterSpecificRegions(List<RefGenomeHomopolymer> refGenomeHomopolymers)
    {
        if(!mConfig.SpecificRegions.isEmpty())
        {
            sLogger.info(mConfig.SpecificRegions);

            refGenomeHomopolymers.removeIf(refGenomeHomopolymer -> mConfig.SpecificRegions.stream()
                    .noneMatch(o -> o.overlaps(refGenomeHomopolymer.genomeRegion)));
        }
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
        RepeatProfileConfig.registerConfig(configBuilder);

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

        System.exit(new RepeatProfileApplication(configBuilder).run(args));
    }
}
