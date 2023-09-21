package com.hartwig.hmftools.errorprofile;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.region.Strand;

import it.unimi.dsi.fastutil.objects.Object2IntArrayMap;

// Stats for each genome position. We keep strand stats separate
public class GenomePositionStats
{
    GenomePosition mGenomePosition;

    char refBase;
    Counts posStrandCounts = new Counts();
    Counts negStrandCounts = new Counts();

    static class Counts
    {
        int totalCount = 0;
        int aCount = 0;
        int cCount = 0;
        int gCount = 0;
        int tCount = 0;
        int deleteCount = 0;

        int softClippedCount = 0;

        // we want to use normal java convention, where insert at a certain position
        // means it is inserted BEFORE
        // now the inserts of bases and counts
        private Object2IntArrayMap<String> insertCounts;

        Object2IntArrayMap<String> getInsertCounts()
        {
            if (insertCounts == null) insertCounts = new Object2IntArrayMap<>();
            return insertCounts;
        }

        int getInsertCount(String seq)
        {
            if (insertCounts == null) return 0;
            return insertCounts.getOrDefault(seq, 0);
        }

        void addToBaseCount(char base)
        {
            switch(base)
            {
                case 'A': aCount++; break;
                case 'C': cCount++; break;
                case 'G': gCount++; break;
                case 'T': tCount++; break;
                default:
                    throw new IllegalArgumentException("invalid base: " + base);
            }
        }

        int getBaseCount(char base)
        {
            switch(base)
            {
                case 'A': return aCount;
                case 'C': return cCount;
                case 'G': return gCount;
                case 'T': return tCount;
                default:
                    throw new IllegalArgumentException("invalid base: " + base);
            }
        }
    }

    Counts getStrandCounts(Strand strand)
    {
        return strand == Strand.FORWARD ? posStrandCounts : negStrandCounts;
    }

    void addInsert(Strand strand, String insertSeq)
    {
        getStrandCounts(strand).totalCount++;
        getStrandCounts(strand).getInsertCounts().compute(insertSeq, (k, count) -> (count == null) ? 1 : count + 1);
    }

    void addDelete(Strand strand)
    {
        getStrandCounts(strand).totalCount++;
        getStrandCounts(strand).deleteCount++;
    }

    void addAlignedBase(Strand strand, char base)
    {
        getStrandCounts(strand).totalCount++;
        getStrandCounts(strand).addToBaseCount(base);
    }

    void addSoftClippedBase(Strand strand)
    {
        getStrandCounts(strand).totalCount++;
        getStrandCounts(strand).softClippedCount++;
    }

    BaseSupport getAlignedBaseSupport(char base)
    {
        // we have an aligned base, lets see what it is
        // we need to check both strands
        return new BaseSupport(posStrandCounts.totalCount, negStrandCounts.totalCount,
                posStrandCounts.getBaseCount(base), negStrandCounts.getBaseCount(base));
    }

    BaseSupport getInsertSupport(String seq)
    {
        // for an inserted sequence, we find the support from the map
        // NOTE: this might not be quite correct. As a partial sequence insert could also be
        // counted together. For now we only take full sequence match
        return new BaseSupport(posStrandCounts.totalCount, negStrandCounts.totalCount,
                posStrandCounts.getInsertCount(seq), negStrandCounts.getInsertCount(seq));
    }

    BaseSupport getDeleteBaseSupport()
    {
        return new BaseSupport(posStrandCounts.totalCount, negStrandCounts.totalCount,
                posStrandCounts.deleteCount, negStrandCounts.deleteCount);
    }

    BaseSupport getSoftClippedBaseSupport(char altBase)
    {
        return new BaseSupport(posStrandCounts.totalCount, negStrandCounts.totalCount,
                posStrandCounts.softClippedCount, negStrandCounts.softClippedCount);
    }
}
