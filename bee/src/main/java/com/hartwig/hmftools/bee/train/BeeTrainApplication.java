package com.hartwig.hmftools.bee.train;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkCreateOutputDir;

import java.io.File;
import java.io.IOException;
import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParametersDelegate;
import com.beust.jcommander.UnixStyleUsageFormatter;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.bee.BeeUtils;
import com.hartwig.hmftools.bee.RefGenomeCompare;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.genome.gc.GCProfileFactory;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.utils.config.DeclaredOrderParameterComparator;
import com.hartwig.hmftools.common.utils.config.LoggingOptions;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.version.VersionInfo;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class BeeTrainApplication
{
    private static final Logger sLogger = LogManager.getLogger(BeeTrainApplication.class);

    // add the AmberConfig options
    @ParametersDelegate
    private final BeeTrainParams mParams = new BeeTrainParams();

    // add to the logging options
    @ParametersDelegate
    private final LoggingOptions mLoggingOptions = new LoggingOptions();

    public int run() throws IOException, InterruptedException, ExecutionException
    {
        mLoggingOptions.setLogLevel();

        VersionInfo versionInfo = new VersionInfo("bee.version");
        sLogger.info("BEE version: {}, build timestamp: {}", versionInfo.version(), versionInfo.buildTime());

        checkCreateOutputDir(mParams.outputDir);

        // most pairs will be in the same chromosome, we should search it
        // need to go by mapping start position
        processBam();

        return 0;
    }

    public void processBam() throws IOException, InterruptedException, ExecutionException
    {
        List<ChrBaseRegion> partitions;

        SamReaderFactory readerFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
        if (mParams.refGenomePath != null)
            readerFactory = readerFactory.referenceSource(new ReferenceSource(new File(mParams.refGenomePath)));

        try (SamReader samReader = readerFactory.open(new File(mParams.bamPath)))
        {
            partitions = BeeUtils.createPartitions(samReader, mParams.specificChromosomes);
        }

        // we create one less for the bam writer
        int numBamReaders = Math.max(mParams.threadCount - 1, 1);
        sLogger.info("processing bam: {} with {} bam reader threads", mParams.bamPath, numBamReaders);

        ExecutorService executorService = Executors.newFixedThreadPool(numBamReaders,
                new ThreadFactoryBuilder().setNameFormat("worker-%d").build());

        Collection<GCProfile> gcProfileList = GCProfileFactory.loadGCContent(mParams.GcProfilePath).values();

        var beeInputWriter = new BeeInputWriter(mParams.outputDir, mParams.sampleId);

        var refGenomeIndexedFasta = new IndexedFastaSequenceFile(new File(mParams.refGenomePath));
        var refGenomeCompare = new RefGenomeCompare(refGenomeIndexedFasta);

        var bqrLogic = new BqrLogic(mParams.sageBqrPath, refGenomeIndexedFasta);

        try (var readHandler = new TruthSetReadHandler(
                beeInputWriter, gcProfileList, refGenomeCompare, bqrLogic,
                mParams.readLength, mParams.minUmiGroupSize))
        {
            List<Future<?>> futures = new ArrayList<>();

            for (ChrBaseRegion chrBaseRegion : partitions)
            {
                final SamReaderFactory finalReaderFactory = readerFactory;
                futures.add(executorService.submit(() ->
                {
                    try (SamReader samReader = finalReaderFactory.open(new File(mParams.bamPath)))
                    {
                        var bamSlicer = new BamSlicer(0, true, true, true);
                        //bamSlicer.setKeepUnmapped();
                        //bamSlicer.setKeepHardClippedSecondaries();
                        bamSlicer.slice(samReader, List.of(chrBaseRegion), readHandler::handleRead);
                        readHandler.onRegionComplete(chrBaseRegion);
                    } catch (IOException e)
                    {
                        throw new RuntimeException(e);
                    }
                }));
            }

            // make sure all the tasks are complete, and also let exceptions through
            for (var future : futures)
            {
                future.get();
            }
        }
        executorService.shutdown();
    }

    public static void main(final String... args) throws IOException, InterruptedException, ExecutionException
    {
        sLogger.info("{}", LocalDateTime.now());
        sLogger.info("args: {}", java.lang.String.join(" ", args));
        BeeTrainApplication beeApp = new BeeTrainApplication();
        JCommander commander = JCommander.newBuilder()
                .addObject(beeApp)
                .build();

        // use unix style formatter
        commander.setUsageFormatter(new UnixStyleUsageFormatter(commander));
        // help message show in order parameters are declared
        commander.setParameterDescriptionComparator(new DeclaredOrderParameterComparator(BeeTrainApplication.class));

        try
        {
            commander.parse(args);
        }
        catch (com.beust.jcommander.ParameterException e)
        {
            System.out.println("Unable to parse args: " + e.getMessage());
            commander.usage();
            System.exit(1);
        }

        // set all thread exception handler
        Thread.setDefaultUncaughtExceptionHandler((Thread t, Throwable e) ->
        {
            sLogger.error("[{}]: uncaught exception: {}", t, e);
            e.printStackTrace(System.err);
            System.exit(1);
        });

        System.exit(beeApp.run());
    }

}
