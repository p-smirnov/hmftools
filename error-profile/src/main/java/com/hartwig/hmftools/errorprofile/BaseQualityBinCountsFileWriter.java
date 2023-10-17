package com.hartwig.hmftools.errorprofile;

import java.io.File;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

import org.jetbrains.annotations.NotNull;

public class BaseQualityBinCountsFileWriter
{
    enum Column
    {
        firstOfPair,
        readPosition,
        ref,
        alt,
        trinucleotideContext,
        rawBaseQuality,
        nonVariant,
        realVariant
    }

    private static String FILE_EXTENSION = ".errorprofile.bq_bins.tsv.gz";

    public static String generateFilename(String basePath, String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    /*
            this.firstOfPair = firstOfPair;
        this.readPosition = readPosition;
        this.ref = ref;
        this.alt = alt;
        this.trinucleotideContext0 = trinucleotideContext0;
        this.trinucleotideContext1 = trinucleotideContext1;
        this.trinucleotideContext2 = trinucleotideContext2;
        this.rawBaseQuality = rawBaseQuality;

     */

    public static void write(final String filename, @NotNull final Map<BaseQualityBin, BaseQualityBinCounter.Count> baseQualityBinCountMap)
    {
        Comparator<BaseQualityBin> comparator = Comparator.comparing((BaseQualityBin o) -> !o.firstOfPair)
                .thenComparingInt((BaseQualityBin o) -> o.readPosition)
                .thenComparing((BaseQualityBin o) -> o.ref)
                .thenComparing((BaseQualityBin o) -> o.alt)
                .thenComparing((BaseQualityBin o) -> o.trinucleotideContext0)
                .thenComparing((BaseQualityBin o) -> o.trinucleotideContext1)
                .thenComparing((BaseQualityBin o) -> o.trinucleotideContext2)
                .thenComparingInt((BaseQualityBin o) -> o.rawBaseQuality);

        // sort the bins
        List<BaseQualityBin> baseQualBins = baseQualityBinCountMap.keySet().stream().sorted(comparator).collect(Collectors.toList());

        List<String> columns = Arrays.stream(Column.values()).map(Enum::name).collect(Collectors.toList());

        new DelimFileWriter().write(filename, columns, baseQualBins, (baseQualityBin, row) ->
        {
            BaseQualityBinCounter.Count count = baseQualityBinCountMap.get(baseQualityBin);
            row.set(Column.firstOfPair, baseQualityBin.firstOfPair);
            row.set(Column.readPosition, baseQualityBin.readPosition);
            row.set(Column.ref, (char)baseQualityBin.ref);
            row.set(Column.alt, (char)baseQualityBin.alt);
            row.set(Column.trinucleotideContext, baseQualityBin.getTrinucleotideContextString());
            row.set(Column.rawBaseQuality, baseQualityBin.rawBaseQuality);
            row.set(Column.nonVariant, count.nonVariantCount.get());
            row.set(Column.realVariant, count.realVariantCount.get());
        });
    }
}
