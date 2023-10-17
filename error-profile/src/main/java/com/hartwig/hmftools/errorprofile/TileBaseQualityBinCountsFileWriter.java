package com.hartwig.hmftools.errorprofile;

import java.io.File;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

import org.jetbrains.annotations.NotNull;

public class TileBaseQualityBinCountsFileWriter
{
    enum Column
    {
        flowcell,
        lane,
        tile,
        firstOfPair,
        readPosition,
        ref,
        alt,
        rawBaseQuality,
        nonVariant,
        realVariant
    }

    private static String FILE_EXTENSION = ".errorprofile.tile_bq_bins.tsv.gz";

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

    public static void write(final String filename, @NotNull final Map<TileBaseQualityBin, BaseQualityBinCounter.Count> baseQualityBinCountMap)
    {
        Comparator<TileBaseQualityBin> comparator =
                Comparator.comparing((TileBaseQualityBin o) -> o.flowcell)
                        .thenComparingInt(o -> o.lane)
                        .thenComparing(o -> o.tile)
                        .thenComparing(o -> !o.firstOfPair)
                        .thenComparingInt(o -> o.readPosition)
                        .thenComparing(o -> o.ref)
                        .thenComparing(o -> o.alt)
                        .thenComparingInt(o -> o.rawBaseQuality);

        // sort the bins
        List<TileBaseQualityBin> baseQualBins = baseQualityBinCountMap.keySet().stream().sorted(comparator).collect(Collectors.toList());

        List<String> columns = Arrays.stream(Column.values()).map(Enum::name).collect(Collectors.toList());

        new DelimFileWriter().write(filename, columns, baseQualBins, (baseQualityBin, row) ->
        {
            BaseQualityBinCounter.Count count = baseQualityBinCountMap.get(baseQualityBin);
            row.set(Column.flowcell, baseQualityBin.flowcell);
            row.set(Column.lane, baseQualityBin.lane);
            row.set(Column.tile, baseQualityBin.tile);
            row.set(Column.firstOfPair, baseQualityBin.firstOfPair);
            row.set(Column.readPosition, baseQualityBin.readPosition);
            row.set(Column.ref, (char)baseQualityBin.ref);
            row.set(Column.alt, (char)baseQualityBin.alt);
            row.set(Column.rawBaseQuality, baseQualityBin.rawBaseQuality);
            row.set(Column.nonVariant, count.nonVariantCount.get());
            row.set(Column.realVariant, count.realVariantCount.get());
        });
    }
}
