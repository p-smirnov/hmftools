package com.hartwig.hmftools.errorprofile.repeat;

import java.io.File;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

import org.apache.commons.lang3.StringUtils;
import org.jetbrains.annotations.NotNull;

public class RepeatProfileFile
{
    enum Column
    {
        chromosome,
        start,
        end,
        base,
        numReads,
        numReadsRejected,
        readRepeatLengths,
    }

    private static String FILE_EXTENSION = ".repeat.tsv.gz";

    public static String generateFilename(String basePath, String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    public static void write(final String filename, @NotNull final Collection<RepeatAnalyser> repeatAnalysers)
    {
        Comparator<RepeatAnalyser> comparator =
                Comparator.comparing((RepeatAnalyser o) -> o.refGenomeHomopolymer.chromosome())
                        .thenComparingInt(o -> o.refGenomeHomopolymer.referenceStart())
                        .thenComparing(o -> o.refGenomeHomopolymer.referenceEnd());

        // sort the bins
        List<RepeatAnalyser> sortedAnalysers = repeatAnalysers.stream().sorted(comparator).collect(Collectors.toList());

        List<String> columns = Arrays.stream(Column.values()).map(Enum::name).collect(Collectors.toList());

        new DelimFileWriter().write(filename, columns, sortedAnalysers, (repeatAnalyser, row) ->
        {
            row.set(Column.chromosome, repeatAnalyser.refGenomeHomopolymer.chromosome());
            row.set(Column.start, repeatAnalyser.refGenomeHomopolymer.referenceStart());
            row.set(Column.end, repeatAnalyser.refGenomeHomopolymer.referenceEnd());
            row.set(Column.base, (char)repeatAnalyser.refGenomeHomopolymer.base);
            row.set(Column.numReads, repeatAnalyser.getReadRepeatMatches().size());
            row.set(Column.numReadsRejected, repeatAnalyser.getReadRepeatMatches().stream().filter(o -> o.shouldDropRead).count());
            row.set(Column.readRepeatLengths, getRepeatString(repeatAnalyser.getReadRepeatMatches()));
        });
    }

    public static String getRepeatString(final List<ReadRepeatMatch> readRepeatMatches)
    {
        List<Integer> repeatLengths = readRepeatMatches.stream().filter(o -> !o.shouldDropRead).map(ReadRepeatMatch::readRepeatLength).collect(
                Collectors.toList());
        return StringUtils.join(repeatLengths, ",");
    }
}
