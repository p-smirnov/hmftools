package com.hartwig.hmftools.bamtools.compare;

import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.samtools.SupplementaryReadData;

import htsjdk.samtools.SAMRecord;

public class ReadWriter implements IDiffConsumer
{
    private final BufferedWriter mWriter;

    public ReadWriter(final CompareConfig config)
    {
        mWriter = initialiseWriter(config.OutputFile);
    }

    public boolean initialised() { return mWriter != null; }

    private BufferedWriter initialiseWriter(final String filename)
    {
        try
        {
            // write summary metrics
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("ReadId\tChromosome\tPosStart\tMismatchType\tDiff\tMateChr\tMatePos");
            writer.write("\tCigar\tFlags\tMapQual\tPaired\tIsFirst\tNegStrand\tDuplicate\tIsSupp\tSuppData");
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            BT_LOGGER.error("failed to initialise BAM comparison file: {}", e.toString());
            return null;
        }
    }

    @Override
    public synchronized void addDiff(final DiffRecord diff)
    {
        try
        {
            String diffDetails = diff.DiffList != null ? diff.DiffList.stream().collect(Collectors.joining(";")) : "";
            mWriter.write(format("%s\t%s\t%d\t%s\t%s\t%s\t%d",
                    diff.Read.getReadName(), diff.Read.getReferenceName(), diff.Read.getAlignmentStart(), diff.Type, diffDetails,
                    diff.Read.getMateReferenceName(), diff.Read.getMateAlignmentStart()));

            mWriter.write(format("\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s",
                    diff.Read.getCigarString(), diff.Read.getFlags(), diff.Read.getMappingQuality(), diff.Read.getReadPairedFlag(), diff.Read.getFirstOfPairFlag(),
                    diff.Read.getReadNegativeStrandFlag(), diff.Read.getDuplicateReadFlag(), diff.Read.getSupplementaryAlignmentFlag(),
                    diff.Read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE) ? SupplementaryReadData.extractAlignment(diff.Read).asCsv() : "N/A"));

            mWriter.newLine();
        }
        catch(IOException e)
        {
            BT_LOGGER.error("failed to write to BAM comparison file: {}", e.toString());
            System.exit(1);
        }
    }

    public void close() { closeBufferedWriter(mWriter); }

    public static String readDetails(final SAMRecord read)
    {
        boolean unmapped = read.getReadUnmappedFlag();

        return format("%s_%s_%d_%s_%s_%s",
                read.getReadName(),
                unmapped ? "unmapped" : read.getReferenceName(),
                unmapped ? 0 : read.getAlignmentStart(),
                read.getReadNegativeStrandFlag() ? "fwd" : "rev",
                read.getFirstOfPairFlag() ? "R1" : "R2",
                read.getSupplementaryAlignmentFlag() ? "supp" : "prim");
    }
}
