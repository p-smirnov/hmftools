package com.hartwig.hmftools.errorprofile;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeVersion;
import static com.hartwig.hmftools.common.samtools.BamUtils.addValidationStringencyOption;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.SPECIFIC_REGIONS;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.loadSpecificRegions;

import java.util.List;

import com.hartwig.hmftools.common.samtools.BamUtils;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.apache.commons.cli.ParseException;

import htsjdk.samtools.ValidationStringency;

public class ErrorProfileConfig
{
    public final String SampleId;
    public final String BamPath;
    public final String RefGenomeFile;
    public final String OutputDir;

    public final int MinMappingQuality;
    public final ValidationStringency BamStringency;
    public final int Threads;

    public final List<ChrBaseRegion> SpecificRegions;

    private static final String MIN_MAP_QUALITY = "min_map_quality";
    public static final int DEFAULT_MIN_MAPPING_QUALITY = 50;

    public ErrorProfileConfig(final ConfigBuilder configBuilder) throws ParseException
    {
        SampleId = configBuilder.getValue(SAMPLE);
        BamPath = configBuilder.getValue("bam");
        RefGenomeFile = configBuilder.getValue(REF_GENOME);
        OutputDir = parseOutputDir(configBuilder);
        Threads = parseThreads(configBuilder);
        MinMappingQuality = configBuilder.getInteger(MIN_MAP_QUALITY);
        BamStringency = BamUtils.validationStringency(configBuilder);
        SpecificRegions = loadSpecificRegions(configBuilder.getValue(SPECIFIC_REGIONS));
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, true, "sample id");
        configBuilder.addPath("bam", true, "path to bam file");

        addRefGenomeVersion(configBuilder);
        configBuilder.addPath(REF_GENOME, true, REF_GENOME_CFG_DESC + ", required when using CRAM files");

        addOutputDir(configBuilder);

        configBuilder.addInteger(
                MIN_MAP_QUALITY, "Minimum mapping quality for an alignment to be used", DEFAULT_MIN_MAPPING_QUALITY);

        addThreadOptions(configBuilder);
        addValidationStringencyOption(configBuilder);
        addLoggingOptions(configBuilder);

        addSpecificChromosomesRegionsConfig(configBuilder);
    }

    public boolean isValid()
    {
        checkCreateOutputDir(OutputDir);
        return true;
    }
}
