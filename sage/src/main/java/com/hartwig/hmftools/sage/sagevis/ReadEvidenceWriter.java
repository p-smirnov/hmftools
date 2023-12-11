// TODO: Remove in final
//package com.hartwig.hmftools.sage.sagevis;
//
//import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
//import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
//import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
//
//import java.io.BufferedWriter;
//import java.io.IOException;
//
//public class ReadEvidenceWriter
//{
//    private final String mFilename;
//    private BufferedWriter mWriter;
//
//    public ReadEvidenceWriter(final String filename)
//    {
//        mFilename = filename;
//        try
//        {
//            SG_LOGGER.info("Writing read evidence output to {}", filename);
//            mWriter = createBufferedWriter(filename, false);
//
//            // TODO: convert to tsv?
//            mWriter.write(ReadEvidenceRecord.CSV_HEADER);
//            mWriter.newLine();
//        }
//        catch(IOException e)
//        {
//            SG_LOGGER.error("Failed to initialise read evidence output file {}: {}", filename, e.toString());
//            System.exit(1);
//        }
//    }
//
//    public synchronized void write(ReadEvidenceRecord readEvidence)
//    {
//        try
//        {
//            // TODO: convert to tsv?
//            mWriter.write(readEvidence.asCsv());
//            mWriter.newLine();
//        }
//        catch(IOException e)
//        {
//            SG_LOGGER.error("Failed to write read evidence to {}: {}", mFilename, e.toString());
//            System.exit(1);
//        }
//    }
//
//    public void close()
//    {
//        closeBufferedWriter(mWriter);
//    }
//}
