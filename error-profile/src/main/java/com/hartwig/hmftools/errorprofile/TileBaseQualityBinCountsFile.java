package com.hartwig.hmftools.errorprofile;

import java.io.File;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

import org.apache.commons.lang3.Validate;
import org.jetbrains.annotations.NotNull;

public class TileBaseQualityBinCountsFile
{
    enum Column
    {
        flowcell,
        lane,
        tile,
        firstOfPair,
        readPosition,
        rawBaseQuality,
        errorCount,
        totalCount
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
                        //.thenComparing(o -> o.ref)
                        //.thenComparing(o -> o.alt)
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
            //row.set(Column.ref, (char)baseQualityBin.ref);
            //row.set(Column.alt, (char)baseQualityBin.alt);
            row.set(Column.rawBaseQuality, baseQualityBin.rawBaseQuality);
            row.set(Column.errorCount, count.getErrorCount());
            row.set(Column.totalCount, count.getTotalCount());
        });
    }

    public static Map<TileBaseQualityBin, BaseQualityBinCounter.Count> read(final String filename)
    {
        try(DelimFileReader reader = new DelimFileReader(filename))
        {
            final Map<TileBaseQualityBin, BaseQualityBinCounter.Count> baseQualityBinCountMap = new HashMap<>();

            for(DelimFileReader.Row row : reader)
            {
                String flowcell = Objects.requireNonNull(row.get(Column.flowcell)).intern();
                short lane = (short)Objects.requireNonNull(row.getInt(Column.lane)).intValue();
                short tile = (short)Objects.requireNonNull(row.getInt(Column.tile)).intValue();
                boolean firstOfPair = Objects.requireNonNull(row.getBoolean(Column.firstOfPair));
                short readPosition = (short)Objects.requireNonNull(row.getInt(Column.readPosition)).intValue();
                //byte ref = (byte)Objects.requireNonNull(row.getChar(Column.ref)).charValue();
                //byte alt = (byte)Objects.requireNonNull(row.getChar(Column.alt)).charValue();
                byte rawBaseQuality = (byte)Objects.requireNonNull(row.getInt(Column.rawBaseQuality)).intValue();

                TileBaseQualityBin tileBaseQualityBin = new TileBaseQualityBin(flowcell, lane, tile, firstOfPair, readPosition,
                        // ref, alt,
                        rawBaseQuality);
                BaseQualityBinCounter.Count count = new BaseQualityBinCounter.Count();
                count.setErrorCount(Objects.requireNonNull(row.getInt(Column.errorCount)));
                count.setTotalCount(Objects.requireNonNull(row.getInt(Column.totalCount)));

                Validate.isTrue(baseQualityBinCountMap.put(tileBaseQualityBin, count) == null, "duplicate tile bin");
            }

            return baseQualityBinCountMap;
        }
    }
}
