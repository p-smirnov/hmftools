package com.hartwig.hmftools.errorprofile.recal;

import java.io.File;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.file.DelimFileWriter;
import com.hartwig.hmftools.errorprofile.TileAdjustmentKey;

import org.jetbrains.annotations.NotNull;

// for each base qual recalibration, record the regression values
public class TileBaseQualRecalFile
{
    enum Column
    {
        flowcell,
        lane,
        tile,
        firstOfPair,
        rawBaseQuality,
        leftSlope,
        leftIntersect,
        midStart,
        midSlope,
        midIntersect,
        rightStart,
        rightSlope,
        rightIntersect,
        empiricalBaseQuality,
    }

    private static String FILE_EXTENSION = ".errorprofile.tile_bq_recal.tsv.gz";

    public static String generateFilename(String basePath, String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    public static void write(final String filename, @NotNull final Map<TileAdjustmentKey, TileBaseQualAdjustment> tileBaseQualAdjustmentMap)
    {
        Comparator<TileAdjustmentKey> comparator =
                Comparator.comparing((TileAdjustmentKey o) -> o.flowcell)
                        .thenComparingInt(o -> o.lane)
                        .thenComparing(o -> o.tile)
                        .thenComparing(o -> !o.firstOfPair)
                        .thenComparingInt(o -> o.rawBaseQuality);

        // sort the bins
        List<TileAdjustmentKey> baseQualBins = tileBaseQualAdjustmentMap.keySet().stream().sorted(comparator).collect(Collectors.toList());

        List<String> columns = Arrays.stream(Column.values()).map(Enum::name).collect(Collectors.toList());

        DelimFileWriter.write(filename, columns, baseQualBins, (key, row) ->
        {
            TileBaseQualAdjustment tileBaseQualAdjustment = tileBaseQualAdjustmentMap.get(key);
            row.set(Column.flowcell, key.flowcell);
            row.set(Column.lane, key.lane);
            row.set(Column.tile, key.tile);
            row.set(Column.firstOfPair, key.firstOfPair);
            row.set(Column.rawBaseQuality, key.rawBaseQuality);
            row.set(Column.leftSlope, tileBaseQualAdjustment.leftSlope);
            row.set(Column.leftIntersect, tileBaseQualAdjustment.leftIntersect);
            row.set(Column.midStart, tileBaseQualAdjustment.midStart);
            row.set(Column.midSlope, tileBaseQualAdjustment.midSlope);
            row.set(Column.midIntersect, tileBaseQualAdjustment.midIntersect);
            row.set(Column.rightStart, tileBaseQualAdjustment.rightStart);
            row.set(Column.rightSlope, tileBaseQualAdjustment.rightSlope);
            row.set(Column.rightIntersect, tileBaseQualAdjustment.rightIntersect);
            row.set(Column.empiricalBaseQuality, tileBaseQualAdjustment.empiricalBaseQuality);
        });
    }

    /*
    public static void read(final String filename, @NotNull final Map<TileBaseQualityBin, BaseQualityBinCounter.Count> baseQualityBinCountMap)
    {
        try(DelimFileReader reader = new DelimFileReader(filename))
        {
            for(DelimFileReader.Row row : reader)
            {
                String flowcell = Objects.requireNonNull(row.get(TileBaseQualityBinCountsFile.Column.flowcell)).intern();
                short lane = (short)Objects.requireNonNull(row.getInt(TileBaseQualityBinCountsFile.Column.lane)).intValue();
                short tile = (short)Objects.requireNonNull(row.getInt(TileBaseQualityBinCountsFile.Column.tile)).intValue();
                boolean firstOfPair = Objects.requireNonNull(row.getBoolean(TileBaseQualityBinCountsFile.Column.firstOfPair));
                short readPosition = (short)Objects.requireNonNull(row.getInt(TileBaseQualityBinCountsFile.Column.readPosition)).intValue();
                byte ref = (byte)Objects.requireNonNull(row.getChar(TileBaseQualityBinCountsFile.Column.ref)).charValue();
                byte alt = (byte)Objects.requireNonNull(row.getChar(TileBaseQualityBinCountsFile.Column.alt)).charValue();
                byte rawBaseQuality = (byte)Objects.requireNonNull(row.getInt(TileBaseQualityBinCountsFile.Column.rawBaseQuality)).intValue();

                TileBaseQualityBin tileBaseQualityBin = new TileBaseQualityBin(flowcell, lane, tile, firstOfPair, readPosition, ref, alt, rawBaseQuality);
                BaseQualityBinCounter.Count count = new BaseQualityBinCounter.Count();
                count.setNonVariantCount(Objects.requireNonNull(row.getInt(TileBaseQualityBinCountsFile.Column.nonVariant)));
                count.setRealVariantCount(Objects.requireNonNull(row.getInt(TileBaseQualityBinCountsFile.Column.realVariant)));

                Validate.isTrue(baseQualityBinCountMap.put(tileBaseQualityBin, count) == null, "duplicate tile bin");
            }
        }
    }
     */
}
