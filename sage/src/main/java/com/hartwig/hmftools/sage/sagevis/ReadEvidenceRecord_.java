// TODO: Remove in final
//package com.hartwig.hmftools.sage.sagevis;
//
//import static java.lang.String.format;
//
//import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
//
//import java.util.List;
//import java.util.Map;
//import java.util.StringJoiner;
//import java.util.stream.Collectors;
//
//import com.hartwig.hmftools.sage.common.IndexedBases;
//import com.hartwig.hmftools.sage.evidence.ReadContextCounter;
//
//import org.apache.commons.compress.utils.Lists;
//
//import htsjdk.samtools.Cigar;
//import htsjdk.samtools.CigarElement;
//import htsjdk.samtools.CigarOperator;
//import htsjdk.samtools.SAMRecord;
//
//public class ReadEvidenceRecord
//{
//    public final String Sample;
//    public final String Chromosome;
//    public final String Ref;
//    public final String Alt;
//    public final ReadContextCounter.MatchType MatchType;
//    public final SAMRecord Read;
//    public final IndexedBases ReadContextBases;
//    public final int ReadIndex;
//    public final double Quality;
//
//    // TODO: Better name for position?
//    // TODO: readIndex?
//
//    public ReadEvidenceRecord(final String sample, final String chromosome, final String ref, final String alt,
//            final ReadContextCounter.MatchType matchType, final SAMRecord record, final IndexedBases readContextBases, final int readIndex,
//            final double quality)
//    {
//        Sample = sample;
//        Chromosome = chromosome;
//        Ref = ref;
//        Alt = alt;
//        MatchType = matchType;
//        Read = record;
//        ReadContextBases = readContextBases;
//        ReadIndex = readIndex;
//        Quality = quality;
//    }
//
//    // TODO: as tsv?
//    public static final List<String> HEADERS = getHeaders();
//    public static final String CSV_HEADER = HEADERS.stream().collect(Collectors.joining(","));
//
//    private static List<String> getHeaders()
//    {
//        List<String> headers = Lists.newArrayList();
//        headers.add("sample");
//        headers.add("chromosome");
//        headers.addAll(IndexedBases.HEADERS);
//        headers.add("ref");
//        headers.add("alt");
//        headers.add("matchType");
//        headers.add("readIndex");
//        headers.add("quality");
//        headers.add("samStr");
//        return headers;
//    }
//
//    public String asCsv()
//    {
//        StringJoiner csvJoiner = new StringJoiner(",");
//        csvJoiner.add(Sample);
//        csvJoiner.add(Chromosome);
//        csvJoiner.add(ReadContextBases.asCsv());
//        csvJoiner.add(Ref);
//        csvJoiner.add(Alt);
//        csvJoiner.add(MatchType.name());
//        csvJoiner.add(String.valueOf(ReadIndex));
//        csvJoiner.add(format("%.1f", Quality));
//        csvJoiner.add(Read.getSAMString().trim());
//        return csvJoiner.toString();
//    }
//
//    public static ReadEvidenceRecord fromMap(final Map<String, String> recordMap)
//    {
//        String sample = recordMap.get("sample");
//        String chromosome = recordMap.get("chromosome");
//        String ref = recordMap.get("ref");
//        String alt = recordMap.get("alt");
//        ReadContextCounter.MatchType matchType = ReadContextCounter.MatchType.valueOf(recordMap.get("matchType"));
//        IndexedBases readContextBases = IndexedBases.fromMap(recordMap);
//        int readIndex = Integer.parseInt(recordMap.get("readIndex"));
//        double quality = Double.parseDouble(recordMap.get("quality"));
//        SAMRecord read = SamStrParser.parse(recordMap.get("samStr"));
//
//        return new ReadEvidenceRecord(sample, chromosome, ref, alt, matchType, read, readContextBases, readIndex, quality);
//    }
//
//    public IndexedReadBases getIndexedReadString()
//    {
//        int unclippedStart = Read.getUnclippedStart();
//        String readString = Read.getReadString();
//        byte[] baseQuals = Read.getBaseQualities();
//
//        List<IndexedBase> indexedBases = Lists.newArrayList();
//        int baseIdx = 0;
//        for (CigarElement cigarElem : Read.getCigar().getCigarElements())
//        {
//            CigarOperator cigarOp = cigarElem.getOperator();
//            int elemLen = cigarElem.getLength();
//            for (int i = 0; i < elemLen; i++)
//            {
//                switch (cigarOp)
//                {
//                    case M:
//                    case EQ:
//                    case X:
//                    case S:
//                        indexedBases.add(new IndexedBase(readString.charAt(baseIdx), baseQuals[baseIdx]));
//                        baseIdx++;
//                        break;
//                    case I:
//                        if (!indexedBases.isEmpty())
//                            indexedBases.get(indexedBases.size() - 1).incRightInsertCount();
//
//                        baseIdx++;
//                        break;
//                    case D:
//                        indexedBases.add(IndexedBase.createDelBase());
//                        break;
//                    case H:
//                    case N:
//                        indexedBases.add(IndexedBase.createMissingBase());
//                        break;
//                }
//            }
//        }
//
//        // TODO: remove at end
//        if (baseIdx != readString.length())
//        {
//            SG_LOGGER.error("getIndexedReadString failed");
//            throw new RuntimeException("getIndexedReadString failed");
//        }
//
//        return new IndexedReadBases(indexedBases, unclippedStart);
//    }
//}
