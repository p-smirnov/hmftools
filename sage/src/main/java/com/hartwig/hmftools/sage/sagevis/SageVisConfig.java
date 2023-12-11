// TODO: Not in final
//package com.hartwig.hmftools.sage.sagevis;
//
//import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
//import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
//import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
//
//import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
//import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
//
//public class SageVisConfig
//{
//    public final String InputFile;
//    public final String OutputFile;
//    public final String RefGenomeFile;
//    public final RefGenomeVersion RefGenVersion;
//
//    private static final String INPUT_FILE = "input_file";
//    private static final String OUTPUT_FILE = "output_file";
//
//    public SageVisConfig(final ConfigBuilder configBuilder)
//    {
//        InputFile = configBuilder.getValue(INPUT_FILE);
//        OutputFile = configBuilder.getValue(OUTPUT_FILE);
//
//        RefGenVersion = RefGenomeVersion.from(configBuilder);
//        RefGenomeFile = configBuilder.getValue(REF_GENOME);
//    }
//
//    public static void addConfig(final ConfigBuilder configBuilder)
//    {
//        configBuilder.addPath(INPUT_FILE, true, "Input read evidence file");
//        configBuilder.addConfigItem(OUTPUT_FILE, true, "Output HTML file");
//
//        addRefGenomeConfig(configBuilder, true);
//        addLoggingOptions(configBuilder);
//    }
//}
