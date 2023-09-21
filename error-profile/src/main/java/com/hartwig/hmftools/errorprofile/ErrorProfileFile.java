package com.hartwig.hmftools.errorprofile;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

public class ErrorProfileFile
{
    enum Column
    {
        numReads, numSplitReads, numShortFragmentReads, num5PrimeSoftClippedReads, num3PrimeSoftClippedReads
    }

    private static String FILE_EXTENSION = ".errorprofile.stats.tsv";

    public static String generateFilename(String basePath, String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    public static void write(final String filename, ErrorProfileStats errorProfileStats)
    {
        List<String> columns = Arrays.stream(Column.values()).map(Enum::name).collect(Collectors.toList());
        new DelimFileWriter().write(filename, columns, List.of(errorProfileStats), (stats, row) ->
        {
            row.set(Column.numReads, stats.numReads);
            row.set(Column.numSplitReads, stats.numSplitReads);
            row.set(Column.numShortFragmentReads, stats.numShortFragmentReads);
            row.set(Column.num5PrimeSoftClippedReads, stats.num5PrimeSoftClippedReads);
            row.set(Column.num3PrimeSoftClippedReads, stats.num3PrimeSoftClippedReads);
        });
    }
}
