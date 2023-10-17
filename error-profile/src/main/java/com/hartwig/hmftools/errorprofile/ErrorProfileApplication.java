package com.hartwig.hmftools.errorprofile;

import java.io.File;
import java.io.FileNotFoundException;
import java.time.Duration;
import java.time.Instant;
import java.util.concurrent.atomic.AtomicLong;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.version.VersionInfo;

import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class ErrorProfileApplication
{
    public static final Logger sLogger = LogManager.getLogger(ErrorProfileApplication.class);

    private final ErrorProfileConfig mConfig;

    private VersionInfo mVersionInfo;

    ReadProfiler mReadProfiler = new ReadProfiler();

    ReadQualAnalyser mReadQualAnalyser;

    ReadBaseSupportFileWriter mReadBaseSupportFileWriter;

    ReadProfileFileWriter mReadProfileFileWriter;

    final BaseQualityBinCounter mBaseQualityBinCounter;

    AtomicLong mReadsProcessed = new AtomicLong(0);

    public ErrorProfileApplication(final ConfigBuilder configBuilder) throws ParseException
    {
        mConfig = new ErrorProfileConfig(configBuilder);
        mReadBaseSupportFileWriter = new ReadBaseSupportFileWriter(
                ReadBaseSupportFileWriter.generateFilename(mConfig.OutputDir, mConfig.SampleId));
        mReadProfileFileWriter = new ReadProfileFileWriter(
                ReadProfileFileWriter.generateFilename(mConfig.OutputDir, mConfig.SampleId));

        mBaseQualityBinCounter = new BaseQualityBinCounter();
    }

    public int run() throws InterruptedException, FileNotFoundException
    {
        Instant start = Instant.now();

        mVersionInfo = new VersionInfo("errorprofile.version");

        if(!mConfig.isValid())
        {
            sLogger.error(" invalid config, exiting");
            return 1;
        }

        mReadQualAnalyser = new ReadQualAnalyser(new IndexedFastaSequenceFile(new File(mConfig.RefGenomeFile)),
                this::processReadProfile, this::processReadBaseSupport);

        processBam();

        writeBaseQualityBinFile();

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
        long readsProcessed = mReadsProcessed.incrementAndGet();
        mReadQualAnalyser.processBam(mConfig);

        if (readsProcessed % 10_000_000 == 0)
        {
            // write the stats every 100m reads
            writeBaseQualityBinFile();
        }
    }

    private void processRead(SAMRecord read, ChrBaseRegion baseRegion)
    {
        // mReadClassifier.classifyRead(read);
        // mReadBaseSupportAnalyser.processRead(read, baseRegion);
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

    private void regionComplete(ChrBaseRegion baseRegion)
    {
        // mReadBaseSupportAnalyser.regionComplete(baseRegion);
        // mReadClassifier.classifyRead(read);
    }

    private void writeBaseQualityBinFile()
    {
        synchronized(mBaseQualityBinCounter)
        {
            sLogger.info("writing base quality bin counts to output");

            // write the base qual bin counts
            BaseQualityBinCountsFileWriter.write(BaseQualityBinCountsFileWriter.generateFilename(mConfig.OutputDir, mConfig.SampleId),
                    mBaseQualityBinCounter.getBaseQualityCountMap());

            TileBaseQualityBinCountsFileWriter.write(TileBaseQualityBinCountsFileWriter.generateFilename(mConfig.OutputDir, mConfig.SampleId),
                    mBaseQualityBinCounter.getTileBaseQualityCountMap());
        }
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

        System.exit(new ErrorProfileApplication(configBuilder).run());
    }
}
