package com.hartwig.hmftools.bee.train;

import static java.util.Collections.emptyList;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;

import java.util.List;

import com.beust.jcommander.Parameter;

public class BeeTrainParams
{
    public static final int DEFAULT_READ_LENGTH = 151;

    @Parameter(names = "-sample", description = "Name of sample")
    public String sampleId;

    @Parameter(names = "-bam", description = "Path to indexed bam/cram file")
    public String bamPath;

    @Parameter(names = "-" + REF_GENOME,
               required = true,
               description = "Path to the reference genome fasta file.")
    public String refGenomePath;

    @Parameter(names = "-gc_profile",
               required = true,
               description = "Location of GC Profile")
    public String GcProfilePath;

    @Parameter(names = "-output_dir",
               required = true,
               description = "Path to the output directory. "
                       + "This directory will be created if it does not already exist.")
    public String outputDir;

    @Parameter(names = "-threads",
               description = "Number of threads")
    public int threadCount = 1;

    @Parameter(names = "-read_length",
               description = "Read Length")
    public int readLength = DEFAULT_READ_LENGTH;

    @Parameter(names = "-specific_chr",
               description = "Optional: list of chromosomes separated by ,")
    public List<String> specificChromosomes = emptyList();
}
