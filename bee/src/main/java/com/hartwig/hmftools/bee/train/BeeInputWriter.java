package com.hartwig.hmftools.bee.train;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Arrays;

import com.hartwig.hmftools.bee.BeeInput;
import com.hartwig.hmftools.common.utils.FileWriterUtils;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;

public final class BeeInputWriter implements Closeable
{
    private enum Column
    {
        readString,
        baseQualities,
        bqrQualities,
        refGenomeBases,
        isRead1,
        isMapped,
        isMateMapped,
        gcContent,
        mappability,
        mapQuality,
        cigar,
        baseError
    }

    private static final String FILE_EXTENSION = ".bee.tsv.gz";
    private static final DecimalFormat sDecimalFormat = new DecimalFormat();

    private final CSVPrinter mCsvPrinter;

    public BeeInputWriter(String basePath, String sample) throws IOException
    {
        String filePath = basePath + File.separator + sample + FILE_EXTENSION;

        var csvFormat = CSVFormat.Builder.create()
                .setDelimiter('\t').setRecordSeparator('\n')
                .setHeader(Column.class)
                .build();

        mCsvPrinter = csvFormat.print(FileWriterUtils.createGzipBufferedWriter(filePath));
    }

    @Override
    public void close() throws IOException
    {
        mCsvPrinter.flush();
        mCsvPrinter.close();
    }

    synchronized public void write(BeeInput beeInput) throws IOException
    {
        StringBuilder stringBuilder = new StringBuilder(beeInput.baseError.length);
        for (boolean b : beeInput.baseError)
        {
            stringBuilder.append(b ? '1' : '0');
        }
        String baseErrorString = stringBuilder.toString();

        stringBuilder = new StringBuilder();

        for (int i = 0; i < beeInput.bqrBaseQualities.length; ++i)
        {
            if (i > 0)
                stringBuilder.append(',');
            stringBuilder.append(sDecimalFormat.format(beeInput.bqrBaseQualities[i]));
        }
        String bqrBaseQualityString = stringBuilder.toString();

        for (Column c : Column.values())
        {
            switch (c)
            {
                case readString:
                    mCsvPrinter.print(beeInput.readString);
                    break;
                case baseQualities:
                    mCsvPrinter.print(beeInput.baseQualityString);
                    break;
                case bqrQualities:
                    mCsvPrinter.print(bqrBaseQualityString);
                    break;
                case refGenomeBases:
                    mCsvPrinter.print(beeInput.referenceGenomeBases);
                    break;
                case isRead1:
                    mCsvPrinter.print(beeInput.isRead1 ? 1 : 0);
                    break;
                case isMapped:
                    mCsvPrinter.print(beeInput.isMapped ? 1 : 0);
                    break;
                case isMateMapped:
                    mCsvPrinter.print(beeInput.isMateMapped ? 1 : 0);
                    break;
                case gcContent:
                    mCsvPrinter.print(beeInput.gcContent);
                    break;
                case mappability:
                    mCsvPrinter.print(beeInput.mappability);
                    break;
                case mapQuality:
                    mCsvPrinter.print(beeInput.mapQuality);
                    break;
                case cigar:
                    mCsvPrinter.print(beeInput.cigar);
                    break;
                case baseError:
                    mCsvPrinter.print(baseErrorString);
                    break;
            }
        }

        mCsvPrinter.println();
    }

    private static String boolArrayToString(boolean[] boolArray)
    {
        StringBuilder stringBuilder = new StringBuilder(boolArray.length);
        for (boolean b : boolArray)
        {
            stringBuilder.append(b ? 1 : 0);
        }
        return stringBuilder.toString();
    }
}
