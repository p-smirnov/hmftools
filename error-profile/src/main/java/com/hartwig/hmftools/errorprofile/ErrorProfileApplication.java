package com.hartwig.hmftools.errorprofile;

import static java.time.format.DateTimeFormatter.ISO_ZONED_DATE_TIME;

import java.io.File;
import java.io.FileNotFoundException;
import java.time.Duration;
import java.time.Instant;
import java.util.concurrent.atomic.AtomicLong;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.version.VersionInfo;

import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;

public class ErrorProfileApplication
{
    public static final Logger sLogger = LogManager.getLogger(ErrorProfileApplication.class);

    private final ErrorProfileConfig mConfig;

    ReadProfiler mReadProfiler = new ReadProfiler();

    ReadQualAnalyser mReadQualAnalyser;

    ReadBaseSupportFileWriter mReadBaseSupportFileWriter;

    ReadProfileFileWriter mReadProfileFileWriter;

    AtomicLong mReadsProcessed = new AtomicLong(0);

    public ErrorProfileApplication(final ConfigBuilder configBuilder) throws ParseException
    {
        mConfig = new ErrorProfileConfig(configBuilder);
        mReadBaseSupportFileWriter = new ReadBaseSupportFileWriter(
                ReadBaseSupportFileWriter.generateFilename(mConfig.OutputDir, mConfig.SampleId));
        mReadProfileFileWriter = new ReadProfileFileWriter(
                ReadProfileFileWriter.generateFilename(mConfig.OutputDir, mConfig.SampleId));
    }

    public int run(final String... args) throws InterruptedException
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

        mReadQualAnalyser = new ReadQualAnalyser();

        processBam();

        writeBaseQualityBinFile(mReadQualAnalyser.getBaseQualityBinCounter());

        // write error stats
        ErrorProfileFile.write(ErrorProfileFile.generateFilename(mConfig.OutputDir, mConfig.SampleId),
                mReadProfiler.mStats);

        for (ReadProfile errorProfile : mReadProfiler.readProfiles)
        {
            if (errorProfile == null)
            {
                sLogger.info("error profile is null");
            }
        }

        mReadBaseSupportFileWriter.close();
        mReadProfileFileWriter.close();

        Instant finish = Instant.now();
        long seconds = Duration.between(start, finish).getSeconds();
        sLogger.info("run complete. Time taken: {}m {}s", seconds / 60, seconds % 60);

        return 0;
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

    private void processBam() throws InterruptedException
    {
        mReadQualAnalyser.processBam(mConfig);
    }

    private void writeBaseQualityBinFile(BaseQualityBinCounter baseQualityBinCounter)
    {
        sLogger.info("writing base quality bin counts to output");

        // write the base qual bin counts
        BaseQualityBinCountsFile.write(BaseQualityBinCountsFile.generateFilename(mConfig.OutputDir, mConfig.SampleId),
                baseQualityBinCounter.getBaseQualityCountMap());

        TileBaseQualityBinCountsFile.write(TileBaseQualityBinCountsFile.generateFilename(mConfig.OutputDir, mConfig.SampleId),
                baseQualityBinCounter.getTileBaseQualityCountMap());
    }

    public static void main(final String... args) throws InterruptedException, FileNotFoundException, ParseException
    {
        ConfigBuilder configBuilder = new ConfigBuilder("ErrorProfile");
        ErrorProfileConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        // set all thread exception handler
        // if we do not do this, exception thrown in other threads will not be handled and results
        // in the program hanging
        Thread.setDefaultUncaughtExceptionHandler((Thread t, Throwable e) ->
        {
            sLogger.error("[{}]: uncaught exception: {}", t, e);
            e.printStackTrace(System.err);
            System.exit(1);
        });

        System.exit(new ErrorProfileApplication(configBuilder).run(args));
    }
}
