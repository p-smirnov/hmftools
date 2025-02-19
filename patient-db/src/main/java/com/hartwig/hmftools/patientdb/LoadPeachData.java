package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.patientdb.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.patientdb.CommonUtils.LOGGER;
import static com.hartwig.hmftools.patientdb.CommonUtils.logVersion;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.common.peach.PeachCalls;
import com.hartwig.hmftools.common.peach.PeachCallsFile;
import com.hartwig.hmftools.common.peach.PeachGenotype;
import com.hartwig.hmftools.common.peach.PeachGenotypeFile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class LoadPeachData
{
    private static final String PEACH_GENOTYPE_TXT = "peach_genotype_txt";
    private static final String PEACH_CALLS_TXT = "peach_calls_txt";

    public static void main(@NotNull String[] args) throws ParseException, SQLException, IOException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addConfigItem(SAMPLE, SAMPLE_DESC);
        addDatabaseCmdLineArgs(configBuilder, true);
        configBuilder.addPath(PEACH_GENOTYPE_TXT, true, "Path towards the PEACH genotype txt file");
        configBuilder.addPath(PEACH_CALLS_TXT, true, "Path towards the PEACH calls txt file");

        configBuilder.checkAndParseCommandLine(args);

        String sample = configBuilder.getValue(SAMPLE);

        String peachGenotypeTxt = configBuilder.getValue(PEACH_GENOTYPE_TXT);
        String peachCallsTxt = configBuilder.getValue(PEACH_CALLS_TXT);

        try(DatabaseAccess dbWriter = databaseAccess(configBuilder))
        {
            LOGGER.info("Reading PEACH genotypes from {}", peachGenotypeTxt);
            List<PeachGenotype> peachGenotype = PeachGenotypeFile.read(peachGenotypeTxt);
            LOGGER.info(" Read {} PEACH genotypes", peachGenotype.size());

            LOGGER.info("Reading PEACH calls from {}", peachCallsTxt);
            List<PeachCalls> peachCalls = PeachCallsFile.read(peachCallsTxt);
            LOGGER.info(" Read {} PEACH calls", peachCalls.size());

            LOGGER.info("Writing PEACH into database for {}", sample);
            dbWriter.writePeach(sample, peachGenotype, peachCalls);

            LOGGER.info("Complete");
        }
        catch (Exception e)
        {
            LOGGER.error("Failed to load PEACH data", e);
            System.exit(1);
        }
    }
}
