package com.hartwig.hmftools.errorprofile.repeat;

import static htsjdk.samtools.util.SequenceUtil.A;

import java.util.Arrays;

import org.junit.Assert;
import org.junit.Test;

import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;

public class RepeatAnalyserTest
{
    @Test
    public void testHandleCigarMatch()
    {
        final String chromosome = "chr1";

        // 10 repeat of As
        RefGenomeMicrosatellite refGenomeMicrosatellite = new RefGenomeMicrosatellite(chromosome, 101, 110, A);

        // test a record that matches the whole repeat
        SAMRecord record = createSamRecord(chromosome, 91, "151M", true, false);

        ReadRepeatMatch cigarHandler = ReadRepeatMatch.from(refGenomeMicrosatellite, record);

        Assert.assertEquals(10, cigarHandler.numAligned);
        Assert.assertEquals(0, cigarHandler.numDeleted);
        Assert.assertEquals(0, cigarHandler.numInserted);

        // test such that it barely matches the range
        record = createSamRecord(chromosome, 101, "10M", true, false);
        cigarHandler = ReadRepeatMatch.from(refGenomeMicrosatellite, record);

        Assert.assertEquals(10, cigarHandler.numAligned);
        Assert.assertEquals(0, cigarHandler.numDeleted);
        Assert.assertEquals(0, cigarHandler.numInserted);

        /*
        // test that it matches only 9 out of 10 bases
        record = createSamRecord(chromosome, 100, "1M10M", true, false);
        cigarHandler = ReadRepeatMatch.from(refGenomeHomopolymer, record);

        Assert.assertEquals(9, cigarHandler.numAligned);
        Assert.assertEquals(0, cigarHandler.numDeleted);
        Assert.assertEquals(0, cigarHandler.numInserted);

        // test that it matches only 9 out of 10 bases
        record = createSamRecord(chromosome, 100, "10M", true, false);
        cigarHandler = ReadRepeatMatch.from(refGenomeHomopolymer, record);

        Assert.assertEquals(9, cigarHandler.numAligned);
        Assert.assertEquals(0, cigarHandler.numDeleted);
        Assert.assertEquals(0, cigarHandler.numInserted);

         */
    }

    @Test
    public void testHandleCigarInsert()
    {
        final String chromosome = "chr1";

        // 10 repeat of As
        RefGenomeMicrosatellite refGenomeMicrosatellite = new RefGenomeMicrosatellite(chromosome, 101, 110, A);

        // test 2 base insert right before the repeat bases
        SAMRecord record = createSamRecord(chromosome, 91, "10M2I20M", true, false);

        ReadRepeatMatch cigarHandler = ReadRepeatMatch.from(refGenomeMicrosatellite, record);

        Assert.assertEquals(10, cigarHandler.numAligned);
        Assert.assertEquals(0, cigarHandler.numDeleted);
        Assert.assertEquals(2, cigarHandler.numInserted);

        // test 2 base insert right after the repeat bases
        record = createSamRecord(chromosome, 101, "10M2I20M", true, false);
        cigarHandler = ReadRepeatMatch.from(refGenomeMicrosatellite, record);

        Assert.assertEquals(10, cigarHandler.numAligned);
        Assert.assertEquals(0, cigarHandler.numDeleted);
        Assert.assertEquals(2, cigarHandler.numInserted);

        // test 2 base insert in the middle
        record = createSamRecord(chromosome, 101, "5M2I20M", true, false);
        cigarHandler = ReadRepeatMatch.from(refGenomeMicrosatellite, record);

        Assert.assertEquals(10, cigarHandler.numAligned);
        Assert.assertEquals(0, cigarHandler.numDeleted);
        Assert.assertEquals(2, cigarHandler.numInserted);

        record = createSamRecord(chromosome, 97, "5M2I20M", true, false);
        cigarHandler = ReadRepeatMatch.from(refGenomeMicrosatellite, record);

        Assert.assertEquals(10, cigarHandler.numAligned);
        Assert.assertEquals(0, cigarHandler.numDeleted);
        Assert.assertEquals(2, cigarHandler.numInserted);
    }

    @Test
    public void testHandleCigarDelete()
    {
        final String chromosome = "chr1";

        // 10 repeat of As
        RefGenomeMicrosatellite refGenomeMicrosatellite = new RefGenomeMicrosatellite(chromosome, 101, 110, A);

        // test 2 base delete right before the repeat bases, this should not count as delete of the polymer, but should disqualify the read
        SAMRecord record = createSamRecord(chromosome, 91, "8M2D20M", true, false);

        ReadRepeatMatch cigarHandler = ReadRepeatMatch.from(refGenomeMicrosatellite, record);

        Assert.assertTrue(cigarHandler.shouldDropRead);

        // test 2 base delete right after the repeat bases, this should also disqualify the read
        record = createSamRecord(chromosome, 101, "10M2D20M", true, false);
        cigarHandler = ReadRepeatMatch.from(refGenomeMicrosatellite, record);

        Assert.assertTrue(cigarHandler.shouldDropRead);

        // test 2 base delete in the middle
        record = createSamRecord(chromosome, 101, "5M2D20M", true, false);
        cigarHandler = ReadRepeatMatch.from(refGenomeMicrosatellite, record);

        Assert.assertFalse(cigarHandler.shouldDropRead);
        Assert.assertEquals(8, cigarHandler.numAligned);
        Assert.assertEquals(2, cigarHandler.numDeleted);
        Assert.assertEquals(0, cigarHandler.numInserted);

        // test 2 base delete at the start of the polymer
        record = createSamRecord(chromosome, 91, "10M2D20M", true, false);
        cigarHandler = ReadRepeatMatch.from(refGenomeMicrosatellite, record);

        Assert.assertFalse(cigarHandler.shouldDropRead);
        Assert.assertEquals(8, cigarHandler.numAligned);
        Assert.assertEquals(2, cigarHandler.numDeleted);
        Assert.assertEquals(0, cigarHandler.numInserted);

        // test 2 base delete at the end of the polymer
        record = createSamRecord(chromosome, 91, "18M2D20M", true, false);
        cigarHandler = ReadRepeatMatch.from(refGenomeMicrosatellite, record);

        Assert.assertFalse(cigarHandler.shouldDropRead);
        Assert.assertEquals(8, cigarHandler.numAligned);
        Assert.assertEquals(2, cigarHandler.numDeleted);
        Assert.assertEquals(0, cigarHandler.numInserted);

        // test 2 base delete straddling the boundary, should be dropped
        record = createSamRecord(chromosome, 91, "9M2D20M", true, false);
        cigarHandler = ReadRepeatMatch.from(refGenomeMicrosatellite, record);

        Assert.assertTrue(cigarHandler.shouldDropRead);

        // test 2 base delete straddling the boundary, should be dropped
        record = createSamRecord(chromosome, 91, "19M2D20M", true, false);
        cigarHandler = ReadRepeatMatch.from(refGenomeMicrosatellite, record);

        Assert.assertTrue(cigarHandler.shouldDropRead);

        // test 13 base delete straddling the boundary, should be dropped
        record = createSamRecord(chromosome, 91, "9M13D20M", true, false);
        cigarHandler = ReadRepeatMatch.from(refGenomeMicrosatellite, record);

        Assert.assertTrue(cigarHandler.shouldDropRead);
    }

    private static SAMRecord createSamRecord(final String chromosome, int readStart, final String cigar, boolean firstInPair, boolean isReversed)
    {
        SAMRecordSetBuilder recordBuilder = new SAMRecordSetBuilder();
        recordBuilder.setUnmappedHasBasesAndQualities(false);

        SAMRecord record = recordBuilder.addFrag(
                "read1", 1, readStart, isReversed, false, cigar, "", 10, false);

        /*
        record.setReadBases(readBases.getBytes());

        final byte[] qualities = new byte[readBases.length()];

        for(int i = 0; i < readBases.length(); ++i)
            qualities[i] = (byte)baseQual;

        record.setBaseQualities(qualities);
        */

        record.setReferenceName(chromosome);
        record.setReferenceIndex(1);

        int flags = 0;

        flags |= SAMFlag.READ_PAIRED.intValue();
        flags |= SAMFlag.PROPER_PAIR.intValue();

        if(isReversed)
            flags |= SAMFlag.READ_REVERSE_STRAND.intValue();

        if(firstInPair)
            flags |= SAMFlag.FIRST_OF_PAIR.intValue();
        else
            flags |= SAMFlag.SECOND_OF_PAIR.intValue();

        record.setFlags(flags);
        return record;
    }
}
