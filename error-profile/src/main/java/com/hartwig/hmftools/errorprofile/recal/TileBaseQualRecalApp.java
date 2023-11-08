package com.hartwig.hmftools.errorprofile.recal;

import static java.time.format.DateTimeFormatter.ISO_ZONED_DATE_TIME;

import java.io.FileNotFoundException;
import java.time.Duration;
import java.time.Instant;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.errorprofile.BaseQualityBinCounter;
import com.hartwig.hmftools.errorprofile.TileAdjustmentKey;
import com.hartwig.hmftools.errorprofile.TileBaseQualityBin;
import com.hartwig.hmftools.errorprofile.TileBaseQualityBinCountsFile;

import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class TileBaseQualRecalApp
{
    public static final Logger sLogger = LogManager.getLogger(TileBaseQualRecalApp.class);
    private final TileBaseQualRecalConfig mConfig;

    public TileBaseQualRecalApp(final ConfigBuilder configBuilder) throws ParseException
    {
        mConfig = new TileBaseQualRecalConfig(configBuilder);
    }

    public int run(final String... args)
    {
        Instant start = Instant.now();

        VersionInfo versionInfo = new VersionInfo("errorprofile.version");

        sLogger.info("ErrorProfile version: {}", versionInfo.version());

        sLogger.debug("build timestamp: {}, run args: {}",
                versionInfo.buildTime().format(ISO_ZONED_DATE_TIME), String.join(" ", args));

        /*if(!mConfig.isValid())
        {
            sLogger.error(" invalid config, exiting");
            return 1;
        }*/

        // load up the counts
        Map<TileBaseQualityBin, BaseQualityBinCounter.Count> tileBaseQualityBinCountMap = TileBaseQualityBinCountsFile.read(
                TileBaseQualityBinCountsFile.generateFilename(mConfig.InputDir, mConfig.SampleId));

        sLogger.info("finished reading tile base qual counts");

        TileBaseQualRecal tileBaseQualRecal = new TileBaseQualRecal();
        Map<TileAdjustmentKey, TileBaseQualStats> tileStats = tileBaseQualRecal.computeTileAdjustmentStats(tileBaseQualityBinCountMap);
        Map<TileAdjustmentKey, TileBaseQualAdjustment> tileAdjustments = Maps.transformValues(tileStats, v -> v.adjustmentFunction);

        // write recal values out to file
        TileBaseQualRecalFile.write(TileBaseQualRecalFile.generateFilename(mConfig.OutputDir, mConfig.SampleId), tileAdjustments);

        List<TileBaseQualOutlier> outliers = tileStats.values().stream().flatMap(o -> o.outliers.stream())
                .collect(Collectors.toList());

        // also write the outliers
        TileBaseQualOutlierFile.write(TileBaseQualOutlierFile.generateFilename(mConfig.OutputDir, mConfig.SampleId), outliers);

        Instant finish = Instant.now();
        long seconds = Duration.between(start, finish).getSeconds();
        sLogger.info("run complete. Time taken: {}m {}s", seconds / 60, seconds % 60);

        return 0;
    }

    public static void main(final String... args) throws InterruptedException, FileNotFoundException, ParseException
    {
        ConfigBuilder configBuilder = new ConfigBuilder("ErrorProfile");
        TileBaseQualRecalConfig.registerConfig(configBuilder);

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

        System.exit(new TileBaseQualRecalApp(configBuilder).run(args));
    }
}
