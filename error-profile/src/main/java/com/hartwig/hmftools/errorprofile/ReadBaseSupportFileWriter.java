package com.hartwig.hmftools.errorprofile;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.List;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;
import org.jetbrains.annotations.NotNull;

// this file is potentially very large with millions of rows
// keep writing to disk to avoid blowing out memory
public class ReadBaseSupportFileWriter implements AutoCloseable
{
    enum Column
    {
        readTag, refPos, cigarOp, readPos53, baseQ, ref, alt, posStrandDepth,
        altPosStrandSupport, negStrandDepth, altNegStrandSupport
    }

    private static String FILE_EXTENSION = ".errorprofile.readpos.tsv.gz";

    public static String generateFilename(String basePath, String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    private final CSVPrinter csvPrinter;
    public ReadBaseSupportFileWriter(@NotNull final String filename)
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

    synchronized public void write(@NotNull final List<ReadBaseSupport> readBaseSupports)
    {
        try
        {
            for(ReadBaseSupport readBaseSupport : readBaseSupports)
            {
                for(ReadBaseSupport.PositionSupport posSupport : readBaseSupport.positionSupports)
                {
                    for(Column c : Column.values())
                    {
                        switch(c)
                        {
                            case readTag: csvPrinter.print(readBaseSupport.readTag); break;
                            //case readName: csvPrinter.print(readBaseSupport.read.getReadName()); break;
                            //case firstOfPair: csvPrinter.print(readBaseSupport.read.getFirstOfPairFlag() ? '1' : '2'); break;
                            //case strand: csvPrinter.print(readBaseSupport.read.getReadNegativeStrandFlag() ? '-' : '+'); break;
                            case refPos: csvPrinter.print(posSupport.refPosition); break;
                            //case mapQ: csvPrinter.print(readBaseSupport.read.getMappingQuality()); break;
                            //case gc: csvPrinter.print(Doubles.round(readBaseSupport.gc, 3)); break;
                            //case cigar: csvPrinter.print(readBaseSupport.read.getCigarString()); break;
                            case cigarOp: csvPrinter.print(posSupport.cigarOperator); break;
                            case readPos53: csvPrinter.print(posSupport.readPosition5To3); break;
                            case baseQ: csvPrinter.print(posSupport.baseQuality); break;
                            case ref: csvPrinter.print(posSupport.ref); break;
                            case alt: csvPrinter.print(posSupport.alt); break;
                            case posStrandDepth: csvPrinter.print(posSupport.posStrandDepth); break;
                            case altPosStrandSupport: csvPrinter.print(posSupport.posStrandSupport); break;
                            case negStrandDepth: csvPrinter.print(posSupport.negStrandDepth); break;
                            case altNegStrandSupport: csvPrinter.print(posSupport.negStrandSupport); break;
                        }
                    }
                    csvPrinter.println();
                }
            }
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }
    }
}
