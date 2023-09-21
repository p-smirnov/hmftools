package com.hartwig.hmftools.errorprofile;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.List;

import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;
import org.jetbrains.annotations.NotNull;

public class ReadProfileFileWriter implements AutoCloseable
{
    enum Column
    {
        readTag,
        readId,
        firstOfPair,
        chromosome,
        position,
        strand,
        readLength,
        cigar,
        mapQ,
        gcContent,
        averageBaseQual,
        lowQualBases,
        insertSize,
        softClipped5P,
        softClipped3P,
        polyGLength,
        homopolymerBase,
        homopolymerStart,
        homopolymerEnd,
        lowQualBeforeHomopolymer,
        lowQualAfterHomopolymer,
        tandemRepeatPattern,
        tandemRepeatStart,
        tandemRepeatEnd,
        lowQualBeforeTandemRepeat,
        lowQualAfterTandemRepeat
    }

    private static String FILE_EXTENSION = ".errorprofile.reads.tsv.gz";

    public static String generateFilename(String basePath, String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    private final CSVPrinter csvPrinter;
    public ReadProfileFileWriter(@NotNull final String filename)
    {
        CSVFormat csvFormat = CSVFormat.Builder.create()
                .setDelimiter('\t').setRecordSeparator('\n')
                .setHeader(Column.class)
                .build();

        try
        {
            csvPrinter = csvFormat.print(createBufferedWriter(filename));
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }
    }

    @Override
    public void close()
    {
        try
        {
            csvPrinter.close();
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }
    }

    synchronized public void write(List<ReadProfile> readProfileList)
    {
        try
        {
            for(ReadProfile readProfile : readProfileList)
            {
                for (Column c : Column.values())
                {
                    switch(c)
                    {
                        case readTag: csvPrinter.print(readProfile.readTag); break;
                        case readId: csvPrinter.print(readProfile.readId); break;
                        case firstOfPair: csvPrinter.print(readProfile.firstOfPair ? 1 : 2); break;
                        case chromosome: csvPrinter.print(readProfile.chromosome); break;
                        case position: csvPrinter.print(readProfile.position); break;
                        case strand: csvPrinter.print(readProfile.strand.asChar()); break;
                        case readLength: csvPrinter.print(readProfile.readLength); break;
                        case cigar: csvPrinter.print(readProfile.cigar); break;
                        case mapQ: csvPrinter.print(readProfile.mapQ); break;
                        case gcContent: csvPrinter.print(Doubles.round(readProfile.gcContent, 3)); break;
                        case averageBaseQual: csvPrinter.print(Math.round(readProfile.averageBaseQual)); break;
                        case lowQualBases: csvPrinter.print(readProfile.numLowQualBases); break;
                        case insertSize: csvPrinter.print(readProfile.insertSize); break;
                        case softClipped5P: csvPrinter.print(readProfile.num5PrimeSoftClipped); break;
                        case softClipped3P: csvPrinter.print(readProfile.num3PrimeSoftClipped); break;
                        case polyGLength: csvPrinter.print(readProfile.polyGLength); break;

                        case homopolymerBase: csvPrinter.print(readProfile.hasHomopolymer() ? readProfile.homopolymerBase : ""); break;
                        case homopolymerStart: csvPrinter.print(readProfile.hasHomopolymer() ? readProfile.homopolymerStart : ""); break;
                        case homopolymerEnd: csvPrinter.print(readProfile.hasHomopolymer() ? readProfile.homopolymerEnd : ""); break;
                        case lowQualBeforeHomopolymer: csvPrinter.print(readProfile.hasHomopolymer() ? readProfile.numLowQualBeforeHomopolymer : ""); break;
                        case lowQualAfterHomopolymer: csvPrinter.print(readProfile.hasHomopolymer() ? readProfile.numLowQualAfterHomopolymer : ""); break;

                        case tandemRepeatPattern: csvPrinter.print(readProfile.tandemRepeat != null ? readProfile.tandemRepeat.pattern : ""); break;
                        case tandemRepeatStart: csvPrinter.print(readProfile.tandemRepeat != null ? readProfile.tandemRepeat.startIndex : ""); break;
                        case tandemRepeatEnd: csvPrinter.print(readProfile.tandemRepeat != null ? readProfile.tandemRepeat.endIndex : ""); break;
                        case lowQualBeforeTandemRepeat: csvPrinter.print(readProfile.tandemRepeat != null ? readProfile.numLowQualBeforeTandemRepeat : ""); break;
                        case lowQualAfterTandemRepeat: csvPrinter.print(readProfile.tandemRepeat != null ? readProfile.numLowQualAfterTandemRepeat : ""); break;
                    }
                }
                csvPrinter.println();
            }
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }
    }
}
