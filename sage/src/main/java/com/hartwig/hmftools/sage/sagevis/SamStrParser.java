// TODO: Not in final
//package com.hartwig.hmftools.sage.sagevis;
//
//import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
//
//import java.io.File;
//
//import htsjdk.samtools.DefaultSAMRecordFactory;
//import htsjdk.samtools.SAMFileHeader;
//import htsjdk.samtools.SAMLineParser;
//import htsjdk.samtools.SAMRecord;
//import htsjdk.samtools.SamReader;
//import htsjdk.samtools.ValidationStringency;
//
//public class SamStrParser
//{
//    private static SAMLineParser mParser = null;
//
//    public static SAMRecord parse(final String samStr)
//    {
//        if (mParser == null)
//        {
//            // samFileHeader only used for setting ref index and mate ref index
//            mParser = new SAMLineParser(new DefaultSAMRecordFactory(), ValidationStringency.SILENT, new SAMFileHeader(), (SamReader)null, (File)null);
//        }
//
//        return mParser.parseLine(samStr);
//    }
//}
