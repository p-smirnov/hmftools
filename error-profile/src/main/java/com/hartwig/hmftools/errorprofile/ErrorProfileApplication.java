package com.hartwig.hmftools.errorprofile;

import static com.hartwig.hmftools.errorprofile.ErrorProfileUtils.createPartitions;

import java.io.File;
import java.io.FileNotFoundException;
import java.time.Duration;
import java.time.Instant;
import java.util.List;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
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

    GenomeRegionStatsCollector mGenomeRegionStatsCollector;

    public ErrorProfileApplication(final ConfigBuilder configBuilder) throws ParseException
    {
        mConfig = new ErrorProfileConfig(configBuilder);
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

        try(ReadBaseSupportFileWriter readBaseSupportFileWriter = new ReadBaseSupportFileWriter(
                ReadBaseSupportFileWriter.generateFilename(mConfig.OutputDir, mConfig.SampleId)))
        {
            try(ReadProfileFileWriter readProfileFileWriter = new ReadProfileFileWriter(
                    ReadProfileFileWriter.generateFilename(mConfig.OutputDir, mConfig.SampleId)))
            {
                mGenomeRegionStatsCollector = new GenomeRegionStatsCollector(new IndexedFastaSequenceFile(new File(mConfig.RefGenomeFile)),
                            readBaseSupportFileWriter, readProfileFileWriter);

                processBam();
            }
        }

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
        int numBamReaders = Math.max(mConfig.Threads - 1, 1);
        List<ChrBaseRegion> partitions = createPartitions(mConfig);

        AsyncBamReader.processBam(mConfig.BamPath, readerFactory(mConfig), partitions, this::processRead, this::regionComplete,
                numBamReaders, mConfig.MinMappingQuality);
    }

    private void processRead(SAMRecord read, ChrBaseRegion baseRegion)
    {
        // mReadClassifier.classifyRead(read);
        mGenomeRegionStatsCollector.processRead(read, baseRegion);
    }

    private void regionComplete(ChrBaseRegion baseRegion)
    {
        mGenomeRegionStatsCollector.regionComplete(baseRegion);
        // mReadClassifier.classifyRead(read);
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
