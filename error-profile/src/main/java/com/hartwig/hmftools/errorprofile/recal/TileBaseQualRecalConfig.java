package com.hartwig.hmftools.errorprofile.recal;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeVersion;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.samtools.BamUtils.addValidationStringencyOption;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.commons.cli.ParseException;

public class TileBaseQualRecalConfig
{
    public static final String INPUT_DIR = "input_dir";

    public final String SampleId;
    public final String InputDir;
    public final String OutputDir;

    public TileBaseQualRecalConfig(final ConfigBuilder configBuilder) throws ParseException
    {
        SampleId = configBuilder.getValue(SAMPLE);
        InputDir = checkAddDirSeparator(configBuilder.getValue(INPUT_DIR));
        OutputDir = parseOutputDir(configBuilder);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, true, "sample id");
        configBuilder.addConfigItem("input_dir", true, "Input directory");
        addOutputDir(configBuilder);
        addLoggingOptions(configBuilder);
    }
}
