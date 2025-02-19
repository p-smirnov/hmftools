package com.hartwig.hmftools.wisp.purity.cn;

import static java.lang.Math.max;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.wisp.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.wisp.purity.FileType.AMBER_LOH;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.AMBER_LOH_CN_THRESHOLD;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.AMBER_LOH_MINOR_ALLELE_THRESHOLD;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.AMBER_LOH_MIN_AF;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.AMBER_LOH_MIN_TUMOR_BAF;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.addCommonFields;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.addCommonHeaderFields;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.amber.AmberBAFFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleCopyNumberFile;
import com.hartwig.hmftools.wisp.purity.PurityConfig;
import com.hartwig.hmftools.wisp.purity.ResultsWriter;
import com.hartwig.hmftools.wisp.purity.SampleData;

import org.apache.commons.math3.distribution.PoissonDistribution;

public class AmberLohCalcs
{
    private final PurityConfig mConfig;

    private final SampleData mSample;
    private final List<PurpleCopyNumber> mCopyNumbers;

    private final BufferedWriter mWriter;

    private final Multimap<Chromosome,AmberBAF> mTumorChromosomeBafs;

    public AmberLohCalcs(final PurityConfig config, final ResultsWriter resultsWriter, final SampleData sample)
    {
        mConfig = config;
        mWriter = resultsWriter.getAmberLohWriter();
        mSample = sample;

        mCopyNumbers = Lists.newArrayList();
        mTumorChromosomeBafs = ArrayListMultimap.create();

        try
        {
            String cnFile = PurpleCopyNumberFile.generateFilenameForReading(mConfig.getPurpleDir(sample.TumorId), mSample.TumorId);
            mCopyNumbers.addAll(PurpleCopyNumberFile.read(cnFile));

            String tumorAmberFile = AmberBAFFile.generateAmberFilenameForReading(mConfig.getAmberDir(sample.TumorId), mSample.TumorId);
            mTumorChromosomeBafs.putAll(AmberBAFFile.read(tumorAmberFile, true));
        }
        catch(Exception e)
        {
            CT_LOGGER.error("sample({}) failed to load Purple and Amber data: {}", mSample.TumorId, e.toString());
        }
    }

    public AmberLohResult processSample(final String sampleId, final PurityContext purityContext)
    {
        try
        {
            String sampleAmberFile = AmberBAFFile.generateAmberFilenameForReading(mConfig.getAmberDir(sampleId), sampleId);
            Multimap<Chromosome,AmberBAF> sampleChromosomeBafs = AmberBAFFile.read(sampleAmberFile, true);

            // select only regions with an LOH
            Map<String,List<PurpleCopyNumber>> chrCopyNumbers = Maps.newHashMap();

            for(PurpleCopyNumber copyNumber : mCopyNumbers)
            {
                // purpleLOH = purpleCN %>% filter(minorAlleleCopyNumber<0.2&copyNumber>0.8)
                if(copyNumber.minorAlleleCopyNumber() > AMBER_LOH_MINOR_ALLELE_THRESHOLD
                || copyNumber.averageTumorCopyNumber() < AMBER_LOH_CN_THRESHOLD)
                {
                    continue;
                }

                List<PurpleCopyNumber> copyNumbers = chrCopyNumbers.get(copyNumber.chromosome());

                if(copyNumbers == null)
                {
                    copyNumbers = Lists.newArrayList();
                    chrCopyNumbers.put(copyNumber.chromosome(), copyNumbers);
                }

                copyNumbers.add(copyNumber);
            }

            List<RegionData> regionDataList = Lists.newArrayList();

            int totalAmberSites = 0;

            for(Map.Entry<String,List<PurpleCopyNumber>> entry : chrCopyNumbers.entrySet())
            {
                String chrStr = entry.getKey();
                Chromosome chromosome = HumanChromosome.fromString(chrStr);

                Collection<AmberBAF> tumorChrSites = mTumorChromosomeBafs.get(chromosome);
                Collection<AmberBAF> sampleChrites = sampleChromosomeBafs.get(chromosome);

                totalAmberSites += tumorChrSites.size();
                totalAmberSites += sampleChrites.size();

                for(PurpleCopyNumber copyNumber : entry.getValue())
                {
                    List<AmberBAF> tumorLohSites = tumorChrSites.stream()
                            .filter(x -> positionWithin(x.position(), copyNumber.start(), copyNumber.end())).collect(Collectors.toList());

                    List<AmberBAF> sampleLohSites = sampleChrites.stream()
                            .filter(x -> positionWithin(x.position(), copyNumber.start(), copyNumber.end())).collect(Collectors.toList());

                    RegionData regionData = buildCopyNumberRegionData(copyNumber, tumorLohSites, sampleLohSites);

                    if(regionData == null)
                        continue;

                    regionDataList.add(regionData);

                    writeLohData(mWriter, mConfig, mSample, sampleId, regionData);
                }
            }

            List<Double> lohSiteAFs = Lists.newArrayList();
            List<Double> lohSiteImpliedPurities = Lists.newArrayList();
            int totalLohSupportCount = 0;
            int totalLohSiteCount = 0;
            int totalLohRegionCount = 0;
            double totalLohCopyNumber = 0;
            double totalRegionAf = 0;

            for(RegionData regionData : regionDataList)
            {
                double regionAverageAf = regionData.averageAF();

                if(regionAverageAf >= AMBER_LOH_MIN_AF)
                    continue;

                ++totalLohRegionCount;
                double regionImpliedPurity = regionData.impliedPurity();
                double regionCopyNumber = regionData.CopyNumber.averageTumorCopyNumber();

                totalLohSiteCount += regionData.Sites.size();

                for(SiteData siteData : regionData.Sites)
                {
                    totalLohSupportCount += siteData.Support;
                    totalLohCopyNumber += regionCopyNumber;
                    totalRegionAf += regionAverageAf;

                    lohSiteAFs.add(regionAverageAf);
                    lohSiteImpliedPurities.add(regionImpliedPurity);
                }
            }

            if(lohSiteAFs.isEmpty() || totalLohSupportCount == 0)
                return AmberLohResult.INVALID_RESULT;

            double lohPercent = lohSiteAFs.size() / (double)totalAmberSites;

            Collections.sort(lohSiteAFs);
            Collections.sort(lohSiteImpliedPurities);

            int medianIndex = lohSiteAFs.size() / 2;
            double lohEstimatedPurity = lohSiteImpliedPurities.get(medianIndex);

            double lohMedianAf = lohSiteAFs.get(medianIndex);

            double lohMeanAf = totalRegionAf / (double)lohSiteAFs.size();

            double lohMeanCN = totalLohCopyNumber / totalLohSiteCount;

            PoissonDistribution poissonDistribution = new PoissonDistribution(0.5 * totalLohSupportCount);
            int observed = (int)round(lohMedianAf * totalLohSupportCount);
            double lohProbability = poissonDistribution.cumulativeProbability(observed);

            return new AmberLohResult(
                    totalLohRegionCount, totalLohSiteCount, lohEstimatedPurity, lohPercent, lohMeanCN, lohMedianAf, lohMeanAf,
                    lohProbability, totalLohSupportCount);
        }
        catch(Exception e)
        {
            CT_LOGGER.error("sample({}) failed to load Amber data: {}", sampleId, e.toString());
            e.printStackTrace();
            return null;
        }
    }

    private class SiteData
    {
        public final int Support;
        public final double TumorBaf;
        public final double SampleBaf;
        public final int TumorDepth;
        public final int SampleDepth;
        public final double AF;

        public SiteData(final double tumorBaf, final double sampleBaf, final int tumorDepth, final int sampleDepth)
        {
            TumorBaf = tumorBaf;
            SampleBaf = sampleBaf;
            TumorDepth = tumorDepth;
            SampleDepth = sampleDepth;

            double adjustedSampleBaf = TumorBaf < 0.5 ? SampleBaf : 1 - SampleBaf;
            double support = sampleDepth * adjustedSampleBaf;
            Support = (int)round(support);
            AF = SampleDepth > 0 ? support / SampleDepth : 0;
        }
    }

    private class RegionData
    {
        public final PurpleCopyNumber CopyNumber;
        public final List<SiteData> Sites;

        public RegionData(final PurpleCopyNumber copyNumber)
        {
            CopyNumber = copyNumber;
            Sites = Lists.newArrayList();
        }

        public double averageAF()
        {
            double totalAF = Sites.stream().mapToDouble(x -> x.AF).sum();
            return Sites.size() > 0 ? totalAF / Sites.size() : 0;
        }

        public double impliedPurity()
        {
            // impliedPurity=max(0,(1-2*avgAF) / (avgAF*(copyNumber-2) + 1)
            double avgAf = averageAF();
            double denom = avgAf * (CopyNumber.averageTumorCopyNumber() - 2) + 1;
            return denom > 0 ? max((1 - 2 * avgAf) / denom, 0) : 0;
        }
    }

    private RegionData buildCopyNumberRegionData(
            final PurpleCopyNumber copyNumber, final List<AmberBAF> tumorLohSites, final List<AmberBAF> sampleLohSites)
    {
        // find matching sites and accumulate VAF counts
        if(tumorLohSites.isEmpty() || sampleLohSites.isEmpty())
            return null;

        int sampleSiteIndex = 0;
        AmberBAF sampleSite = sampleLohSites.get(sampleSiteIndex);

        RegionData regionData = new RegionData(copyNumber);

        for(AmberBAF tumorSite : tumorLohSites)
        {
            if(tumorSite.tumorModifiedBAF() <= AMBER_LOH_MIN_TUMOR_BAF)
                continue;

            while(sampleSite.position() < tumorSite.position())
            {
                ++sampleSiteIndex;

                if(sampleSiteIndex >= sampleLohSites.size())
                {
                    sampleSite = null;
                    break;
                }

                sampleSite = sampleLohSites.get(sampleSiteIndex);
            }

            if(sampleSite == null)
                break;

            if(sampleSite.position() > tumorSite.position())
                continue;

            if(sampleSite.position() == tumorSite.position())
            {
                regionData.Sites.add(
                        new SiteData(tumorSite.tumorBAF(), sampleSite.tumorBAF(), tumorSite.tumorDepth(), sampleSite.tumorDepth()));
            }
        }

        return regionData;
    }

    public static BufferedWriter initialiseAmberLohWriter(final PurityConfig config)
    {
        try
        {
            String fileName = config.formFilename(AMBER_LOH);

            BufferedWriter writer = createBufferedWriter(fileName, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            addCommonHeaderFields(sj, config);

            sj.add("Chromosome").add("CnSegmentStart").add("CnSegmentEnd").add("CopyNumber");
            sj.add("SiteCount").add("SupportCount").add("AvgAF").add("ImpliedPurity");

            writer.write(sj.toString());
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            CT_LOGGER.error("failed to initialise Amber LOH file: {}", e.toString());
            return null;
        }
    }

    private static synchronized void writeLohData(
            final BufferedWriter writer, final PurityConfig config,
            final SampleData sampleData, final String sampleId, final RegionData regionData)
    {
        if(writer == null)
            return;

        try
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);
            addCommonFields(sj, config, sampleData, sampleId);

            sj.add(regionData.CopyNumber.chromosome());
            sj.add(String.valueOf(regionData.CopyNumber.start()));
            sj.add(String.valueOf(regionData.CopyNumber.end()));
            sj.add(format("%.2f", regionData.CopyNumber.averageTumorCopyNumber()));
            sj.add(String.valueOf(regionData.Sites.size()));

            int totalSupport = regionData.Sites.stream().mapToInt(x -> x.Support).sum();
            sj.add(String.valueOf(totalSupport));
            sj.add(format("%.2f", regionData.averageAF()));
            sj.add(format("%.2f", regionData.impliedPurity()));

            writer.write(sj.toString());
            writer.newLine();
        }
        catch(IOException e)
        {
            CT_LOGGER.error("failed to write Amber LOH file: {}", e.toString());
        }
    }
}