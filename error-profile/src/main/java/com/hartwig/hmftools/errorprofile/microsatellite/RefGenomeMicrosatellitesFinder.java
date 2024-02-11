package com.hartwig.hmftools.errorprofile.microsatellite;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeVersion;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.errorprofile.microsatellite.MicrosatelliteAnalyserConstants.MAX_MICROSAT_UNIT_LENGTH;
import static com.hartwig.hmftools.errorprofile.microsatellite.MicrosatelliteAnalyserConstants.MIN_ADJACENT_MICROSAT_DISTANCE;
import static com.hartwig.hmftools.errorprofile.microsatellite.MicrosatelliteAnalyserConstants.MIN_MICROSAT_UNIT_COUNT;

import static htsjdk.samtools.util.SequenceUtil.N;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.ListIterator;
import java.util.Random;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.google.common.collect.ArrayListMultimap;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.commons.cli.ParseException;
import org.apache.commons.lang3.mutable.MutableInt;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.StringUtil;

// simple class to find microsatellites regions in reference genome
public class RefGenomeMicrosatellitesFinder
{
    public static final Logger sLogger = LogManager.getLogger(RefGenomeMicrosatellitesFinder.class);

    static class Candidate
    {
        int startIndex = 0;

        // where we are up to, it is also same as end
        int currentEndIndex = 0;
        byte[] pattern;

        String patternString;

        boolean complete = false;

        Candidate(byte[] pattern, int startIndex, int endIndex)
        {
            this.pattern = pattern;
            this.startIndex = startIndex;
            this.currentEndIndex = endIndex;
            patternString = StringUtil.bytesToString(pattern);
        }

        byte nextBase()
        {
            return pattern[(currentEndIndex - startIndex) % pattern.length];
        }

        int length()
        {
            return currentEndIndex - startIndex;
        }

        int numFullUnits()
        {
            return length() / pattern.length;
        }
    }

    // check if this is a valid repeat unit, i.e.
    // ATAT is not because it is actually 2 x AT
    static boolean isValidUnit(byte[] unit)
    {
        // check each subunit length to see if it is the same one
        for(int subunitLength = 1; subunitLength <= unit.length / 2; ++subunitLength)
        {
            boolean allMatch = true;

            if((unit.length % subunitLength) == 0)
            {
                for(int i = 0; i < subunitLength && allMatch; ++i)
                {
                    byte base = unit[i];
                    for(int j = subunitLength + i; j < unit.length; j += subunitLength)
                    {
                        if(unit[j] != base)
                        {
                            allMatch = false;
                            break;
                        }
                    }
                }
            }
            else
            {
                allMatch = false;
            }

            if (allMatch)
            {
                // found a subunit that matches the whole unit
                return false;
            }
        }
        // did not find a subunit
        return true;
    }

    public static void findMicrosatellites(ReferenceSequenceFile referenceSequenceFile, int minRepeatLength,
            Consumer<RefGenomeMicrosatellite> refGenomeMsConsumer)
    {
        int chunkSize = 100_000;
        findMicrosatellites(referenceSequenceFile, minRepeatLength, refGenomeMsConsumer, chunkSize);
    }

    // algorithm to find short tandem repeats
    // at each base, we try find longest candidate starting from this base.
    // if a microsatellite is found, we start again from the base after.
    //
    static void findMicrosatellites(ReferenceSequenceFile referenceSequenceFile, int minNumRepeats,
            Consumer<RefGenomeMicrosatellite> refGenomeMsConsumer, int chunkSize)
    {
        MutableInt microsatelliteCounter = new MutableInt(0);

        List<SAMSequenceRecord> seqRecords = referenceSequenceFile.getSequenceDictionary().getSequences()
                .stream()
                .filter(o -> HumanChromosome.contains(o.getContig()))
                .collect(Collectors.toList());

        for(SAMSequenceRecord sequenceRecord : seqRecords)
        {
            int length = sequenceRecord.getSequenceLength();

            // pending candidate, we do not accept a candidate when it is completed. We want to avoid
            // candidates that are too close to each other. They are dropped if too close.
            List<RefGenomeMicrosatellite> pendingMicrosatellies = new ArrayList<>();

            // current best candidate
            Candidate bestCandidate = null;

            // candidates
            List<Candidate> currentCandidates = new ArrayList<>();

            for(int start = 1; start <= length; start = start + chunkSize)
            {
                int endInclusive = Math.min(start + chunkSize, length) - 1;
                byte[] seq = referenceSequenceFile.getSubsequenceAt(sequenceRecord.getContig(), start, endInclusive).getBases();

                for(int i = 0; i < seq.length; ++i)
                {
                    byte base = seq[i];

                    if(base == N)
                    {
                        // if we hit an N, we clear everything
                        currentCandidates.clear();
                        bestCandidate = null;
                        continue;
                    }

                    // check all current candidates
                    ListIterator<Candidate> itr = currentCandidates.listIterator();
                    while(itr.hasNext())
                    {
                        Candidate c = itr.next();

                        // check if pattern is still valid
                        if(c.nextBase() == base)
                        {
                            // pattern continues
                            c.currentEndIndex++;
                            if(bestCandidate == null || bestCandidate.numFullUnits() < c.numFullUnits())
                            {
                                bestCandidate = c;
                            }
                        }
                        else
                        {
                            c.complete = true;
                            itr.remove();
                        }
                    }

                    // if the best candidate is completed, we want to check if it is a candidate
                    // we want to include
                    if(bestCandidate != null && bestCandidate.complete)
                    {
                        // check if the best candidate that is completed is a good one
                        int unitRepeatCount = bestCandidate.numFullUnits();

                        if(unitRepeatCount >= minNumRepeats)
                        {
                            // NOTE: for microsatellites with a pattern longer than one, we only accept full repeats
                            // i.e. ATGATGATGATGAT contains 4 full ATG plus AT at the end, even though AT is partial unit,
                            // they are excluded.

                            int baseLength = unitRepeatCount * bestCandidate.pattern.length;

                            // this is a microsatellite
                            RefGenomeMicrosatellite refGenomeMicrosatellite = new RefGenomeMicrosatellite(sequenceRecord.getContig(),
                                    bestCandidate.startIndex,
                                    bestCandidate.startIndex + baseLength - 1, // change to inclusive
                                    bestCandidate.pattern);
                            pendingMicrosatellies.add(refGenomeMicrosatellite);

                            // check the panding microsatellites, see if any can be accepted
                            checkPendingMicrosatellites(pendingMicrosatellies, refGenomeMsConsumer, microsatelliteCounter);
                        }

                        bestCandidate = null;
                    }

                    // also start new Candidates at this location
                    // all these candidates have the current base as the last base of the pattern
                    for(int j = Math.max(i - MAX_MICROSAT_UNIT_LENGTH + 1, 0); j <= i; ++j)
                    {
                        byte[] repeatUnit = Arrays.copyOfRange(seq, j, i + 1);

                        // check that this is a valid repeat unit, i.e. it is not multiple smaller unit
                        if(isValidUnit(repeatUnit))
                        {
                            // also check against existing patterns
                            if(currentCandidates.stream().noneMatch(o -> Arrays.equals(o.pattern, repeatUnit)))
                            {
                                Candidate c = new Candidate(repeatUnit, start + j, start + i + 1);
                                currentCandidates.add(c);
                            }
                        }
                    }
                }
            }

            if(pendingMicrosatellies.size() == 1)
            {
                microsatelliteCounter.increment();
                refGenomeMsConsumer.accept(pendingMicrosatellies.get(0));
            }

            sLogger.info("finished chromosome {}", sequenceRecord.getSequenceName());
        }

        sLogger.info("found {} microsatellite regions in ref genome", microsatelliteCounter);
    }

    // the aim of this code is to remove microsatellites that are too close to each other
    // We do not try to merge them for now, even though tools such as MsDetector would.
    // i.e. AAAAATAAAAAAA
    static void checkPendingMicrosatellites(List<RefGenomeMicrosatellite> pendingMicrosatellies, Consumer<RefGenomeMicrosatellite> refGenomeMsConsumer,
            MutableInt microsatelliteCounter)
    {
        int groupStart = 0;

        // check the panding microsatellites, see if any can be accepted
        for(int i = 0; i < pendingMicrosatellies.size() - 1; ++i)
        {
            RefGenomeMicrosatellite ms1 = pendingMicrosatellies.get(i);
            RefGenomeMicrosatellite ms2 = pendingMicrosatellies.get(i + 1);

            if((ms2.referenceStart() - ms1.referenceEnd()) > MIN_ADJACENT_MICROSAT_DISTANCE)
            {
                // previous group finished, if previous group only has 1 ms, we accept it, otherwise
                // remove them all from the pending list
                if(i == groupStart)
                {
                    // only 1 item, accept this
                    refGenomeMsConsumer.accept(ms1);
                    microsatelliteCounter.increment();
                    sLogger.trace("microsatellite: {}", ms1);
                }
                else
                {
                    // ms are too close to each other, remove them
                    for(int j = groupStart; j <= i; ++j)
                    {
                        sLogger.trace("reject microsatellite as too close to neighbour: {}", pendingMicrosatellies.get(j));
                    }
                }

                // update group start
                groupStart = i + 1;
            }
        }

        // we can remove anything before groupStart
        if(groupStart > 0)
        {
            pendingMicrosatellies.subList(0, groupStart).clear();
        }
    }

    // overload for easy testing
    static List<RefGenomeMicrosatellite> findMicrosatellites(ReferenceSequenceFile referenceSequenceFile, int minNumRepeats, int chunkSize)
    {
        List<RefGenomeMicrosatellite> refGenomeMicrosatellites = new ArrayList<>();
        findMicrosatellites(referenceSequenceFile, minNumRepeats, refGenomeMicrosatellites::add, chunkSize);
        return refGenomeMicrosatellites;
    }

    // filter the microsatellites such that each type of (unit, length) is approximately the target count
    static List<RefGenomeMicrosatellite> filterMicrosatellites(List<RefGenomeMicrosatellite> inputList, int targetCountPerType)
    {
        sLogger.info("filtering {} microsatellite sites, target count per type = {}", inputList.size(), targetCountPerType);

        // first put them all into a multimap
        ArrayListMultimap<Pair<String, Integer>, RefGenomeMicrosatellite> unitLengthMicrosatelliteMap = ArrayListMultimap.create();

        // should be able to use a groupby method in guava
        for(RefGenomeMicrosatellite microsatellite : inputList)
        {
            Pair<String, Integer> k = Pair.of(microsatellite.unitString(), microsatellite.numRepeat);
            unitLengthMicrosatelliteMap.put(k, microsatellite);
        }

        // use same seed for now
        Random random = new Random(0);
        List<RefGenomeMicrosatellite> filteredList = new ArrayList<>();

        for(Pair<String, Integer> msType : unitLengthMicrosatelliteMap.keySet())
        {
            List<RefGenomeMicrosatellite> l = unitLengthMicrosatelliteMap.get(msType);
            double frac = ((double)targetCountPerType) / l.size();
            if(frac < 1.0)
            {
                l.stream().filter(o -> random.nextDouble() <= frac).forEach(filteredList::add);
            }
            else
            {
                filteredList.addAll(l);
            }
        }

        sLogger.info("filtered {} microsatellite sites down to {}", inputList.size(), filteredList.size());

        return filteredList;
    }

    private static class Config
    {
        public final String refGenomeFile;
        public final String outputDir;

        public final RefGenomeVersion refGenomeVersion;

        public Config(final ConfigBuilder configBuilder) throws ParseException
        {
            refGenomeFile = configBuilder.getValue(REF_GENOME);
            outputDir = parseOutputDir(configBuilder);
            refGenomeVersion = RefGenomeVersion.from(configBuilder);
        }

        public static void registerConfig(final ConfigBuilder configBuilder)
        {
            addRefGenomeVersion(configBuilder);
            configBuilder.addPath(REF_GENOME, true, REF_GENOME_CFG_DESC + ", required when using CRAM files");

            addOutputDir(configBuilder);

            addLoggingOptions(configBuilder);

            addSpecificChromosomesRegionsConfig(configBuilder);
        }

        public boolean isValid()
        {
            checkCreateOutputDir(outputDir);
            return true;
        }
    }

    private static void findAndWriteRefGenomeMicrosatellites(Config config) throws Exception
    {
        // find all the polymers
        IndexedFastaSequenceFile refGenome = new IndexedFastaSequenceFile(new File(config.refGenomeFile));
        try (RefGenomeMicrosatelliteFile refGenomeMicrosatelliteFile = new RefGenomeMicrosatelliteFile(
                RefGenomeMicrosatelliteFile.generateFilename(config.outputDir, config.refGenomeVersion)))
        {
            RefGenomeMicrosatellitesFinder.findMicrosatellites(refGenome, MIN_MICROSAT_UNIT_COUNT, refGenomeMicrosatelliteFile::writeRow);

            //filterSpecificRegions(refGenomeMicrosatellites);
        }
    }

    public static void main(final String... args) throws Exception
    {
        ConfigBuilder configBuilder = new ConfigBuilder("ErrorProfile");
        Config.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        setLogLevel(configBuilder);

        Config config = new Config(configBuilder);

        findAndWriteRefGenomeMicrosatellites(config);

        System.exit(0);
    }


        /*
    // algorithm to find short tandem repeats
    // at each base, we try find longest candidate starting from this base.
    // if a microsatellite is found, we start again from the base after.
    //
    static void findHomopolymer(ReferenceSequenceFile referenceSequenceFile, int minNumRepeats,
            Consumer<RefGenomeMicrosatellite> refGenomeMsConsumer, int chunkSize)
    {
        int numMicrosatellites = 0;

        for(SAMSequenceRecord sequenceRecord : referenceSequenceFile.getSequenceDictionary().getSequences())
        {
            if(!HumanChromosome.contains(sequenceRecord.getContig()))
            {
                continue;
            }

            int length = sequenceRecord.getSequenceLength();

            for(int start = 1; start <= length; start = start + chunkSize)
            {
                int endInclusive = Math.min(start + chunkSize, length) - 1;
                byte[] seq = referenceSequenceFile.getSubsequenceAt(sequenceRecord.getContig(), start, endInclusive).getBases();

                for(int i = 0; i < seq.length; )
                {
                    if(seq[i] == N)
                    {
                        ++i;
                        continue;
                    }

                    // we test all candidates from this location, and choose the longest one at the end

                    // candidates
                    List<Candidate> candidates = new ArrayList<>();

                    // first step find all the candidates
                    // j is the base after pattern end
                    for(int j = i + 1; j < Math.min(i + MAX_MICROSAT_UNIT_LENGTH, seq.length); ++j)
                    {
                        int patternLength = j - i;

                        // test that this pattern does not contain N
                        if(seq[j] == N)
                        {
                            // we cannot have a pattern with N, so we can stop here
                            break;
                        }

                        // we want to weed out any cases where it is just a repeat of previous patterns,
                        // i.e. ATAT is just the same as AT
                        boolean isDuplicate = false;
                        for(Candidate c : candidates)
                        {
                            int x = patternLength / c.pattern.length;

                            if(x * c.pattern.length != patternLength)
                            {
                                // first the length must be divisible
                                continue;
                            }

                            // now test that the pattern is the same one
                            boolean isSame = true;
                            for(int k = 0; k < patternLength; ++k)
                            {
                                if(c.pattern[k % c.pattern.length] != seq[i + k])
                                {
                                    isSame = false;
                                    break;
                                }
                            }

                            if(isSame)
                            {
                                isDuplicate = true;
                                break;
                            }
                        }

                        if(!isDuplicate)
                        {
                            byte[] pattern = Arrays.copyOfRange(seq, i, j);
                            Candidate c = new Candidate(pattern);
                            candidates.add(c);
                        }
                    }

                    Candidate candidate = null;
                    int candidateEnd = -1;

                    // we got all the candidates, now try to extend them as far as we can
                    for(int j = i; j < seq.length; ++j)
                    {
                        byte base = seq[j];

                        // check all current candidates
                        ListIterator<Candidate> itr = candidates.listIterator();
                        while(itr.hasNext())
                        {
                            Candidate c = itr.next();

                            // check if pattern is still valid
                            if(c.nextBase() != base)
                            {
                                if(candidates.size() == 1)
                                {
                                    // this is last remaining one, i.e. this is the pattern
                                    candidate = c;
                                    candidateEnd = j;
                                }
                                itr.remove();
                            }
                        }
                    }

                    if(candidate != null)
                    {
                        int unitRepeatCount = (candidateEnd - i) / candidate.pattern.length;

                        if(unitRepeatCount >= minNumRepeats)
                        {
                            // this is a microsatellite
                            RefGenomeMicrosatellite refGenomeMicrosatellite = new RefGenomeMicrosatellite(sequenceRecord.getContig(),
                                    start + i,
                                    start + candidateEnd - 1,
                                    candidate.pattern,
                                    unitRepeatCount);
                            ++numMicrosatellites;
                            refGenomeMsConsumer.accept(refGenomeMicrosatellite);
                            i = candidateEnd;

                            sLogger.trace("microsatellite: {}", refGenomeMicrosatellite);

                            continue;
                        }
                    }
                    // just move to next base
                    ++i;
                }
            }

            sLogger.info("finished chromosome {}", sequenceRecord.getSequenceName());
        }

        sLogger.info("found {} microsatellite regions in ref genome", numMicrosatellites);
    }
     */
}
