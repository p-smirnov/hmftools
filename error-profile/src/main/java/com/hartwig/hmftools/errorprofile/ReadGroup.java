package com.hartwig.hmftools.errorprofile;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

public class ReadGroup
{
    private static final org.apache.logging.log4j.Logger logger = LogManager.getLogger(ReadGroup.class);
    private static final String SUPP_ALIGNMENT_DELIM = ",";

    public static List<SupplementaryAlignment> suppAlignmentPositions(boolean firstOfPair, String suppAlignment)
    {
        if(suppAlignment == null)
        {
            return null;
        }

        List<SupplementaryAlignment> suppAignPos = new ArrayList<>();
        String[] supplementaryItems = suppAlignment.split(";");
        for(String si : supplementaryItems)
        {
            String[] items = si.split(SUPP_ALIGNMENT_DELIM);
            if(items.length < 5)
            {
                continue;
            }

            // supplementary(SA) string attribute looks like
            // 4,191039958,+,68S33M,0,0
            // the first word is the chromosome, the second is the alignment start
            SupplementaryAlignment sa = new SupplementaryAlignment(
                    firstOfPair,
                    items[0],
                    Integer.parseInt(items[1]),
                    items[2].equals("-"),
                    items[3]);
            suppAignPos.add(sa);
        }
        return suppAignPos;
    }

    public static class SupplementaryAlignment
    {
        public final boolean firstOfPair;
        public final String chromosome;
        public final int position;
        public final boolean negativeStrand;

        public final String cigar;

        public SupplementaryAlignment(boolean firstOfPair, String chromosome, int position, boolean negativeStrand, String cigar)
        {
            this.firstOfPair = firstOfPair;
            this.chromosome = chromosome;
            this.negativeStrand = negativeStrand;
            this.position = position;
            this.cigar = cigar;
        }

        public boolean isMatch(SAMRecord read)
        {
            return read.getFirstOfPairFlag() == firstOfPair && read.getAlignmentStart() == position && read.getReadNegativeStrandFlag() == negativeStrand && read.getReferenceName().equals(chromosome)
                    && cigarEqual(read.getCigarString(), cigar);
        }

        private static boolean cigarEqual(String cigar1, String cigar2)
        {
            return cigar1.replace('H', 'S').equals(cigar2.replace('H', 'S'));
        }
    }

    private final String name;
    private final List<SAMRecord> mReads = new ArrayList<>();
    private final List<SAMRecord> mSupplementaryReads = new ArrayList<>();

    public ReadGroup(String name)
    {
        this.name = name;
    }

    public List<SAMRecord> getReads()
    {
        return mReads;
    }

    public List<SAMRecord> getSupplementaryReads()
    {
        return mSupplementaryReads;
    }

    public boolean isComplete(Level logLevel)
    {
        if(getReads().isEmpty())
        {
            logger.log(logLevel, "Read is empty");
            return false;
        }

        // we check several things
        if(getReads().get(0).getReadPairedFlag() && getReads().size() != 2)
        {
            // we haven't got all the reads yet
            logger.log(logLevel, "{} missing mate pair", getReads().get(0));
            return false;
        }

        for(SAMRecord read : getReads())
        {
            String saAttribute = read.getStringAttribute(SAMTag.SA.name());
            if(saAttribute != null)
            {
                for(SupplementaryAlignment sa : suppAlignmentPositions(read.getFirstOfPairFlag(), saAttribute))
                {
                    // check if this supplementary read exists
                    if(getSupplementaryReads().stream().noneMatch(sa::isMatch))
                    {
                        logger.log(logLevel, "{} Missing supplementary read: aligned to {}:{}", getReads().get(0),
                                sa.chromosome, sa.position);
                        return false;
                    }
                }
            }
        }
        return true;
    }

    public boolean contains(SAMRecord read)
    {
        if(!read.getReadName().equals(name))
        {
            return false;
        }

        List<SAMRecord> listToLook;
        if(read.isSecondaryOrSupplementary())
        {
            listToLook = getSupplementaryReads();
        }
        else
        {
            listToLook = getReads();
        }

        // we only check chromosome and alignment start
        return listToLook.stream()
                .anyMatch(x -> x.getFirstOfPairFlag() == read.getFirstOfPairFlag() && x.getAlignmentStart() == read.getAlignmentStart()
                        && x.getReadNegativeStrandFlag() == read.getReadNegativeStrandFlag() && x.getReferenceName().equals(read.getReferenceName()) && x.getCigarString().equals(read.getCigarString()));
    }

    public List<SAMRecord> allReads()
    {
        // Implement the allReads method
        ArrayList<SAMRecord> allReads = new ArrayList<>(mReads);
        allReads.addAll(mSupplementaryReads);
        return allReads;
    }

    @Override
    public String toString()
    {
        return String.format("%s reads(%d) complete(%s)", name, mReads.size(), isComplete(null));
    }

    public boolean invariant()
    {
        // check to make sure several things:
        // 1. same record cannot appear more than once
        // 2. Reads cannot contain supplementary
        // 3. SupplementaryReads must only contain supplementary

        if(getReads().size() > 2)
        {
            return false;
        }

        for(SAMRecord read : getReads())
        {
            if(!read.getReadName().equals(name))
            {
                return false;
            }
        }

        for(SAMRecord suppRead : getSupplementaryReads())
        {
            if(!suppRead.getReadName().equals(name))
            {
                return false;
            }
        }

        boolean hasSupplementaryAlignmentFlag = getReads().stream().anyMatch(SAMRecord::getSupplementaryAlignmentFlag);
        if(hasSupplementaryAlignmentFlag)
        {
            return false;
        }
        else
        {
            return getSupplementaryReads().stream().allMatch(SAMRecord::getSupplementaryAlignmentFlag);
        }
    }

    public List<ChrBaseRegion> findMissingReadBaseRegions()
    {
        List<ChrBaseRegion> baseRegions = new ArrayList<>();
        assert invariant();

        if(getReads().size() == 1)
        {
            SAMRecord read = getReads().get(0);
            if(read.getReadPairedFlag() && read.getMateReferenceIndex() != SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX)
            {
                // note that even if the mate unmapped flag is set, we still might not be there
                // unmapped read of a mapped read pair would show up with the mate position
                if(read.getMateReferenceName().isEmpty() || read.getMateAlignmentStart() <= 0)
                {
                    // this shouldn't happen
                    logger.warn("read({}) invalid mate reference, mate ref index({})", read, read.getMateReferenceIndex());
                }
                else
                {
                    baseRegions.add(new ChrBaseRegion(read.getMateReferenceName(), read.getMateAlignmentStart(), read.getMateAlignmentStart()));
                    logger.trace("{} missing read mate: aligned to {}:{}", getReads().get(0), read.getMateReferenceName(), read.getMateAlignmentStart());
                }
            }
        }

        for(SAMRecord read : getReads())
        {
            String saAttribute = read.getStringAttribute(SAMTag.SA.name());
            if(saAttribute != null)
            {
                for(SupplementaryAlignment sa : suppAlignmentPositions(read.getFirstOfPairFlag(), saAttribute))
                {
                    // check if this supplementary read exists
                    if(getSupplementaryReads().stream().noneMatch(r -> sa.isMatch(r)))
                    {
                        baseRegions.add(new ChrBaseRegion(sa.chromosome, sa.position, sa.position));
                        logger.trace("{} Missing supplementary read: aligned to {}:{}", read, sa.chromosome, sa.position);
                    }
                }
            }
        }
        return baseRegions;
    }

    public boolean acceptRead(SAMRecord record)
    {
        if(contains(record))
        {
            return false;
        }

        if(record.isSecondaryOrSupplementary())
        {
            // bam files generated by some GATK IndelRealigner
            // uses secondary instead of supplementary alignment
            if(record.isSecondaryAlignment())
            {
                record.setSupplementaryAlignmentFlag(true);
                record.setSecondaryAlignment(false);
            }
            getSupplementaryReads().add(record);
        }
        else
        {
            getReads().add(record);
        }
        return true;
    }
}
