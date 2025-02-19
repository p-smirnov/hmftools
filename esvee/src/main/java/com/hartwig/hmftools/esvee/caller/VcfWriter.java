package com.hartwig.hmftools.esvee.caller;

import static com.hartwig.hmftools.common.sv.SvVcfTags.HOTSPOT;
import static com.hartwig.hmftools.common.sv.SvVcfTags.HOTSPOT_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.PON_COUNT;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REPEAT_MASK_COVERAGE_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REPEAT_MASK_REPEAT_CLASS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REPEAT_MASK_COVERAGE;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REPEAT_MASK_REPEAT_CLASS_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REPEAT_MASK_REPEAT_TYPE;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REPEAT_MASK_REPEAT_TYPE_DESC;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.VCF_ZIP_EXTENSION;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.esvee.common.FilterType.PON;
import static com.hartwig.hmftools.esvee.common.SvConstants.ESVEE_FILE_ID;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.variant.GenotypeIds;
import com.hartwig.hmftools.esvee.caller.annotation.RepeatMaskAnnotation;
import com.hartwig.hmftools.esvee.common.FilterType;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class VcfWriter
{
    private final CallerConfig mConfig;
    private final GenotypeIds mGenotypeIds;
    private final SvDataCache mDataCache;

    private final VariantContextWriter mUnfilteredWriter;
    private final VariantContextWriter mSomaticWriter;
    private final VariantContextWriter mGermlineWriter;

    public VcfWriter(
            final CallerConfig config, final VCFHeader vcfHeader, final String gripssVersion, final GenotypeIds genotypeIds,
            final SvDataCache dataCache)
    {
        mConfig = config;
        mGenotypeIds = genotypeIds;
        mDataCache = dataCache;

        String fileSampleId = config.SampleId != null ? config.SampleId : config.ReferenceId;

        String unfilteredVcf = formVcfFilename(fileSampleId, "unfiltered");
        mUnfilteredWriter = initialiseWriter(vcfHeader, gripssVersion, unfilteredVcf);

        if(mConfig.hasTumor())
        {
            String somaticVcf = formVcfFilename(fileSampleId, "somatic");
            mSomaticWriter = initialiseWriter(vcfHeader, gripssVersion, somaticVcf);
        }
        else
        {
            mSomaticWriter = null;
        }

        if(mConfig.hasReference())
        {
            String germlineVcf = formVcfFilename(fileSampleId, "germline");
            mGermlineWriter = initialiseWriter(vcfHeader, gripssVersion, germlineVcf);
        }
        else
        {
            mGermlineWriter = null;
        }
    }

    private String formVcfFilename(final String sampleId, final String fileId)
    {
        String suffix = mConfig.OutputId != null ? "." + mConfig.OutputId + VCF_ZIP_EXTENSION : VCF_ZIP_EXTENSION;
        return mConfig.OutputDir + sampleId + "." + ESVEE_FILE_ID + "." + fileId + suffix;
    }

    private VariantContextWriter initialiseWriter(final VCFHeader vcfHeader, final String gripssVersion, final String vcfFilename)
    {
        VariantContextWriter writer = new VariantContextWriterBuilder()
                .setReferenceDictionary(vcfHeader.getSequenceDictionary())
                .setOutputFile(vcfFilename)
                .setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF)
                .build();

        // ensure genotype sample IDs match the config - also done per variant
        List<String> genotypeSampleNames = Lists.newArrayList();

        if(!mConfig.ReferenceId.isEmpty())
            genotypeSampleNames.add(mConfig.ReferenceId);

        if(!mConfig.SampleId.isEmpty())
            genotypeSampleNames.add(mConfig.SampleId);

        VCFHeader newHeader = new VCFHeader(vcfHeader.getMetaDataInInputOrder(), genotypeSampleNames);

        newHeader.addMetaDataLine(new VCFHeaderLine("gripssVersion", gripssVersion));

        for(FilterType filter : FilterType.values())
        {
            newHeader.addMetaDataLine(new VCFFilterHeaderLine(filter.vcfTag(), filter.vcfDesc()));
        }

        newHeader.addMetaDataLine(new VCFInfoHeaderLine(HOTSPOT, 1, VCFHeaderLineType.Flag, HOTSPOT_DESC));

        newHeader.addMetaDataLine(new VCFInfoHeaderLine(PON_COUNT, 1, VCFHeaderLineType.Integer, "PON count if in PON"));

        newHeader.addMetaDataLine(new VCFInfoHeaderLine(
                REPEAT_MASK_REPEAT_CLASS, 1, VCFHeaderLineType.String, REPEAT_MASK_REPEAT_CLASS_DESC));

        newHeader.addMetaDataLine(new VCFInfoHeaderLine(
                REPEAT_MASK_REPEAT_TYPE, 1, VCFHeaderLineType.String, REPEAT_MASK_REPEAT_TYPE_DESC));

        newHeader.addMetaDataLine(new VCFInfoHeaderLine(
                REPEAT_MASK_COVERAGE, 1, VCFHeaderLineType.Float, REPEAT_MASK_COVERAGE_DESC));

        writer.writeHeader(newHeader);

        return writer;
    }

    public void writeBreakends()
    {
        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chrStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

            List<Breakend> breakendList = mDataCache.getBreakendMap().get(chrStr);

            if(breakendList == null)
                continue;

            for(Breakend breakend : breakendList)
            {
                writeBreakend(breakend);
            }
        }
    }

    private void writeBreakend(final Breakend breakend)
    {
        final SvData sv = breakend.sv();

        List<Genotype> genotypes = Lists.newArrayList();

        if(mGenotypeIds.hasTumor())
            genotypes.add(new GenotypeBuilder().copy(breakend.Context.getGenotype(mGenotypeIds.TumorOrdinal)).name(mConfig.SampleId).make());

        if(mGenotypeIds.hasReference())
            genotypes.add(new GenotypeBuilder().copy(breakend.Context.getGenotype(mGenotypeIds.ReferenceOrdinal)).name(mConfig.ReferenceId).make());

        VariantContextBuilder builder = new VariantContextBuilder(breakend.Context).genotypes(genotypes).filters();

        builder.log10PError(breakend.Qual / -10.0);

        if(sv.isHotspot())
            builder.attribute(HOTSPOT, true);

        if(sv.ponCount() > 0)
            builder.attribute(PON_COUNT, breakend.sv().ponCount());

        if(sv.getRmAnnotation() != null)
        {
            final RepeatMaskAnnotation rmAnnotation = sv.getRmAnnotation();

            // remove any previously set
            builder.attribute(REPEAT_MASK_REPEAT_CLASS, rmAnnotation.RmData.ClassType);
            builder.attribute(REPEAT_MASK_REPEAT_TYPE, rmAnnotation.RmData.Repeat);
            builder.attribute(REPEAT_MASK_COVERAGE, rmAnnotation.Coverage);
        }

        // first write the unfiltered VCF with all breakends
        Set<FilterType> allFilters = breakend.sv().filters();

        Set<FilterType> germlineFilters = allFilters.stream().filter(x -> !x.germlineOnly()).collect(Collectors.toSet());

        Set<FilterType> somaticFilters = Sets.newHashSet(allFilters);

        writeBreakend(mUnfilteredWriter, builder, allFilters);

        if(mSomaticWriter != null)
        {
            if(somaticFilters.isEmpty() || (somaticFilters.size() == 1 && somaticFilters.contains(PON)))
                writeBreakend(mUnfilteredWriter, builder, somaticFilters);
        }

        if(mGermlineWriter != null && germlineFilters.isEmpty())
            writeBreakend(mGermlineWriter, builder, germlineFilters);
    }

    private static void writeBreakend(final VariantContextWriter writer, final VariantContextBuilder builder, final Set<FilterType> filters)
    {
        builder.getFilters().clear();

        if(filters.isEmpty())
            builder.filter(PASS);
        else
            filters.forEach(x -> builder.filter(x.vcfTag()));

        VariantContext variantContext = builder.make(true);

        writer.add(variantContext);
    }

    public void close()
    {
        mUnfilteredWriter.close();

        if(mSomaticWriter != null)
            mSomaticWriter.close();

        if(mGermlineWriter != null)
            mGermlineWriter.close();
    }
}
