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
    private static String CHROMOSOME = "chromosome";
    private static String START = "start";
    private static String END = "end";
    private static String UNIT = "unit";
    private static String NUM_READS = "numReads";
    private static String NUM_READS_REJECTED = "numReadsRejected";
    private static String COUNT_p5 = "count+5";
    private static String COUNT_p4 = "count+4";
    private static String COUNT_p3 = "count+3";
    private static String COUNT_p2 = "count+2";
    private static String COUNT_p1 = "count+1";
    private static String COUNT_p0 = "count+0";

    private static String COUNT_m5 = "count-5";
    private static String COUNT_m4 = "count-4";
    private static String COUNT_m3 = "count-3";
    private static String COUNT_m2 = "count-2";
    private static String COUNT_m1 = "count-1";

    private static String READ_REPEAT_LENGTHS = "readRepeatLengths";

    private static String FILE_EXTENSION = ".repeat.tsv.gz";

    public static String generateFilename(String basePath, String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    public static void write(final String filename, @NotNull final Collection<RepeatAnalyser> repeatAnalysers)
    {
        Comparator<RepeatAnalyser> comparator =
                Comparator.comparing((RepeatAnalyser o) -> o.refGenomeMicrosatellite.chromosome())
                        .thenComparingInt(o -> o.refGenomeMicrosatellite.referenceStart())
                        .thenComparing(o -> o.refGenomeMicrosatellite.referenceEnd());

        // sort the bins
        List<RepeatAnalyser> sortedAnalysers = repeatAnalysers.stream().sorted(comparator).collect(Collectors.toList());

        List<String> columns = List.of(CHROMOSOME, START, END, UNIT, NUM_READS, NUM_READS_REJECTED,
                COUNT_m5, COUNT_m4, COUNT_m3, COUNT_m2, COUNT_m1, COUNT_p0, COUNT_p1, COUNT_p2, COUNT_p3, COUNT_p4, COUNT_p5,
                READ_REPEAT_LENGTHS);

        // add the count columns

        DelimFileWriter.write(filename, columns, sortedAnalysers, (repeatAnalyser, row) ->
        {
            row.set(CHROMOSOME, repeatAnalyser.refGenomeMicrosatellite.chromosome());
            row.set(START, repeatAnalyser.refGenomeMicrosatellite.referenceStart());
            row.set(END, repeatAnalyser.refGenomeMicrosatellite.referenceEnd());
            row.set(UNIT,  repeatAnalyser.refGenomeMicrosatellite.unitString());
            row.set(NUM_READS, repeatAnalyser.getReadRepeatMatches().size());
            row.set(NUM_READS_REJECTED, repeatAnalyser.getReadRepeatMatches().stream().filter(o -> o.shouldDropRead).count());
            int refNumRepeat = repeatAnalyser.refGenomeMicrosatellite.numRepeat;
            row.set(COUNT_p0, getCountWithRepeatUnits(repeatAnalyser, refNumRepeat));
            row.set(COUNT_p5, getCountWithRepeatUnits(repeatAnalyser, refNumRepeat + 5));
            row.set(COUNT_p4, getCountWithRepeatUnits(repeatAnalyser, refNumRepeat + 4));
            row.set(COUNT_p3, getCountWithRepeatUnits(repeatAnalyser, refNumRepeat + 3));
            row.set(COUNT_p2, getCountWithRepeatUnits(repeatAnalyser, refNumRepeat + 2));
            row.set(COUNT_p1, getCountWithRepeatUnits(repeatAnalyser, refNumRepeat + 1));

            row.set(COUNT_m5, getCountWithRepeatUnits(repeatAnalyser, refNumRepeat - 5));
            row.set(COUNT_m4, getCountWithRepeatUnits(repeatAnalyser, refNumRepeat - 4));
            row.set(COUNT_m3, getCountWithRepeatUnits(repeatAnalyser, refNumRepeat - 3));
            row.set(COUNT_m2, getCountWithRepeatUnits(repeatAnalyser, refNumRepeat - 2));
            row.set(COUNT_m1, getCountWithRepeatUnits(repeatAnalyser, refNumRepeat - 1));

            row.set(READ_REPEAT_LENGTHS, getRepeatString(repeatAnalyser.getPassingReadRepeatMatches()));
        });
    }

    public static String getRepeatString(final List<ReadRepeatMatch> readRepeatMatches)
    {
        List<Integer> repeatLengths = readRepeatMatches.stream().filter(o -> !o.shouldDropRead).map(ReadRepeatMatch::readRepeatLength).collect(
                Collectors.toList());
        return StringUtils.join(repeatLengths, ",");
    }

    public static int getCountWithRepeatUnits(RepeatAnalyser repeatAnalyser, int numRepeatUnits)
    {
        return (int)repeatAnalyser.getPassingReadRepeatMatches().stream().filter(o -> o.numRepeatUnits() == numRepeatUnits).count();
    }
}
