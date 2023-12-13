package com.hartwig.hmftools.errorprofile.recal;

import java.io.File;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

import org.jetbrains.annotations.NotNull;

public class TileBaseQualOutlierFile
{
    enum Column
    {
        flowcell,
        lane,
        tile,
        firstOfPair,
        rawBaseQuality,
        readPosition,
        errorCount,
        totalCount,
        recalBaseQuality
    }

    private static String FILE_EXTENSION = ".errorprofile.tile_bq_outlier.tsv.gz";

    public static String generateFilename(String basePath, String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    public static void write(final String filename, @NotNull final List<TileBaseQualOutlier> tileBaseQualOutliers)
    {
        Comparator<TileBaseQualOutlier> comparator =
                Comparator.comparing((TileBaseQualOutlier o) -> o.tileAdjustmentKey.flowcell)
                        .thenComparingInt(o -> o.tileAdjustmentKey.lane)
                        .thenComparing(o -> o.tileAdjustmentKey.tile)
                        .thenComparing(o -> !o.tileAdjustmentKey.firstOfPair)
                        .thenComparingInt(o -> o.tileAdjustmentKey.rawBaseQuality)
                        .thenComparingInt(o -> o.readPosition);

        // sort the bins
        List<TileBaseQualOutlier> sortedOutliers = tileBaseQualOutliers.stream().sorted(comparator).collect(Collectors.toList());

        List<String> columns = Arrays.stream(Column.values()).map(Enum::name).collect(Collectors.toList());

        DelimFileWriter.write(filename, columns, sortedOutliers, (outlier, row) ->
        {
            row.set(Column.flowcell, outlier.tileAdjustmentKey.flowcell);
            row.set(Column.lane, outlier.tileAdjustmentKey.lane);
            row.set(Column.tile, outlier.tileAdjustmentKey.tile);
            row.set(Column.firstOfPair, outlier.tileAdjustmentKey.firstOfPair);
            row.set(Column.rawBaseQuality, outlier.tileAdjustmentKey.rawBaseQuality);
            row.set(Column.readPosition, outlier.readPosition);
            row.set(Column.errorCount, outlier.errorCount);
            row.set(Column.totalCount, outlier.totalCount);
            row.set(Column.recalBaseQuality, outlier.empiricalBaseQuality());
        });
    }

}
