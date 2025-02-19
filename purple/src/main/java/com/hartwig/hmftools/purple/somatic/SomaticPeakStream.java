package com.hartwig.hmftools.purple.somatic;

import static com.google.common.collect.Lists.newArrayList;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.config.PurpleConstants.CLONALITY_BIN_WIDTH;
import static com.hartwig.hmftools.purple.config.PurpleConstants.CLONALITY_MAX_PLOIDY;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.purple.fitting.ModifiableWeightedPloidy;
import com.hartwig.hmftools.purple.fitting.PeakModelData;
import com.hartwig.hmftools.purple.fitting.PeakModelFactory;

public class SomaticPeakStream
{
    public SomaticPeakStream() {}

    public List<PeakModelData> somaticPeakModel(final SomaticVariantCache somaticVariants)
    {
        if(!somaticVariants.hasData())
            return Lists.newArrayList();

        final List<ModifiableWeightedPloidy> weightedPloidies = newArrayList();

        for(SomaticVariant variant : somaticVariants.variants())
        {
            if(!variant.isPass())
                continue;

            if(variant.copyNumber() >= CLONALITY_MAX_PLOIDY)
                continue;

            if(!HumanChromosome.contains(variant.chromosome()) || !HumanChromosome.fromString(variant.chromosome()).isAutosome())
                continue;

            AllelicDepth depth = variant.tumorAlleleDepth();

            if(depth != null)
            {
                weightedPloidies.add(ModifiableWeightedPloidy.create()
                        .from(depth)
                        .setPloidy(variant.copyNumber())
                        .setWeight(1));
            }
        }

        PPL_LOGGER.debug("somatic peak uses {} variants", weightedPloidies.size());

        PeakModelFactory modelFactory = new PeakModelFactory(CLONALITY_MAX_PLOIDY, CLONALITY_BIN_WIDTH);

        List<PeakModelData> peakModels = modelFactory.model(weightedPloidies);

        /*
        for(PeakModelData peakModel : peakModels)
        {
            PPL_LOGGER.trace(format("somatic peak(%.4f wgt=%.4f) bucket(%.4f wgt=%.4f) valid(%s) subclonal(%s)",
                    peakModel.Peak, peakModel.PeakAvgWeight, peakModel.Bucket, peakModel.BucketWeight,
                    peakModel.IsValid, peakModel.IsSubclonal));
        }
        */

        return peakModels;
    }
}
