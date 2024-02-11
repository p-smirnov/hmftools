package com.hartwig.hmftools.errorprofile.microsatellite;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.StringUtil;

public class RefGenomeMicrosatellitesFinderTest
{
    static class TestReferenceSequenceFile implements ReferenceSequenceFile
    {
        final byte[] sequence;

        public TestReferenceSequenceFile(final String sequence)
        {
            this.sequence = StringUtil.stringToBytes(sequence);
        }

        @Override
        public SAMSequenceDictionary getSequenceDictionary()
        {
            return new SAMSequenceDictionary(List.of(new SAMSequenceRecord("chr1", sequence.length)));
        }

        @Override
        public ReferenceSequence nextSequence()
        {
            return null;
        }

        @Override
        public void reset()
        {

        }

        @Override
        public boolean isIndexed()
        {
            return false;
        }

        @Override
        public ReferenceSequence getSequence(final String s)
        {
            return null;
        }

        @Override
        public ReferenceSequence getSubsequenceAt(final String s, final long start, final long end)
        {
            return new ReferenceSequence(s, 4, Arrays.copyOfRange(sequence, (int) (start - 1), (int) end));
        }

        @Override
        public void close() throws IOException
        {

        }
    }

    @Test
    public void test()
    {
        String sequence = "ATCGATCGATCGGGGGGATCGATCAATTTTTTTCG";
        TestReferenceSequenceFile refSeqFile = new TestReferenceSequenceFile(sequence);

        // make sure the test reference sequence file is correct
        Assert.assertEquals("ATCGA", refSeqFile.getSubsequenceAt("chr1", 1, 5).getBaseString());

        Assert.assertEquals("GGGGGGA", refSeqFile.getSubsequenceAt("chr1", 12, 18).getBaseString());

        List<RefGenomeMicrosatellite> polymers = RefGenomeMicrosatellitesFinder.findMicrosatellites(refSeqFile, 4, 20);

        Assert.assertEquals(2, polymers.size());
        RefGenomeMicrosatellite polymer = polymers.get(0);
        Assert.assertEquals("chr1", polymer.genomeRegion.Chromosome);
        Assert.assertEquals(12, polymer.genomeRegion.start());
        Assert.assertEquals(17, polymer.genomeRegion.end());
        Assert.assertEquals("G", polymer.unitString());
        Assert.assertEquals(6, polymer.numRepeat);

        polymer = polymers.get(1);
        Assert.assertEquals("chr1", polymer.genomeRegion.Chromosome);
        Assert.assertEquals(27, polymer.genomeRegion.start());
        Assert.assertEquals(33, polymer.genomeRegion.end());
        Assert.assertEquals("T", polymer.unitString());
        Assert.assertEquals(7, polymer.numRepeat);

        // try different chunk size
        polymers = RefGenomeMicrosatellitesFinder.findMicrosatellites(refSeqFile, 4, 5);

        Assert.assertEquals(2, polymers.size());
        polymer = polymers.get(0);
        Assert.assertEquals("chr1", polymer.genomeRegion.Chromosome);
        Assert.assertEquals(12, polymer.genomeRegion.start());
        Assert.assertEquals(17, polymer.genomeRegion.end());
        Assert.assertEquals("G", polymer.unitString());
        Assert.assertEquals(6, polymer.numRepeat);

        polymer = polymers.get(1);
        Assert.assertEquals("chr1", polymer.genomeRegion.Chromosome);
        Assert.assertEquals(27, polymer.genomeRegion.start());
        Assert.assertEquals(33, polymer.genomeRegion.end());
        Assert.assertEquals("T", polymer.unitString());
        Assert.assertEquals(7, polymer.numRepeat);
    }

    @Test
    public void testMicrosatellites()
    {
        String sequence = "TGATCATCATCATCATCATGATCGGGGGGATCGATCAATATATATTCG";
        TestReferenceSequenceFile refSeqFile = new TestReferenceSequenceFile(sequence);

        List<RefGenomeMicrosatellite> microsatellites = RefGenomeMicrosatellitesFinder.findMicrosatellites(refSeqFile, 4, 20);

        Assert.assertEquals(3, microsatellites.size());
        RefGenomeMicrosatellite microsatellite = microsatellites.get(0);
        Assert.assertEquals("chr1", microsatellite.genomeRegion.Chromosome);
        Assert.assertEquals(3, microsatellite.genomeRegion.start());
        Assert.assertEquals(17, microsatellite.genomeRegion.end());
        Assert.assertEquals("ATC", microsatellite.unitString());
        Assert.assertEquals(5, microsatellite.numRepeat);

        microsatellite = microsatellites.get(1);
        Assert.assertEquals("chr1", microsatellite.genomeRegion.Chromosome);
        Assert.assertEquals(24, microsatellite.genomeRegion.start());
        Assert.assertEquals(29, microsatellite.genomeRegion.end());
        Assert.assertEquals("G", microsatellite.unitString());
        Assert.assertEquals(6, microsatellite.numRepeat);

        microsatellite = microsatellites.get(2);
        Assert.assertEquals("chr1", microsatellite.genomeRegion.Chromosome);
        Assert.assertEquals(38, microsatellite.genomeRegion.start());
        Assert.assertEquals(45, microsatellite.genomeRegion.end());
        Assert.assertEquals("AT", microsatellite.unitString());
        Assert.assertEquals(4, microsatellite.numRepeat);

        // try different chunk size
        microsatellites = RefGenomeMicrosatellitesFinder.findMicrosatellites(refSeqFile, 4, 5);

        Assert.assertEquals(3, microsatellites.size());
        microsatellite = microsatellites.get(0);
        Assert.assertEquals("chr1", microsatellite.genomeRegion.Chromosome);
        Assert.assertEquals(3, microsatellite.genomeRegion.start());
        Assert.assertEquals(17, microsatellite.genomeRegion.end());
        Assert.assertEquals("ATC", microsatellite.unitString());
        Assert.assertEquals(5, microsatellite.numRepeat);

        microsatellite = microsatellites.get(1);
        Assert.assertEquals("chr1", microsatellite.genomeRegion.Chromosome);
        Assert.assertEquals(24, microsatellite.genomeRegion.start());
        Assert.assertEquals(29, microsatellite.genomeRegion.end());
        Assert.assertEquals("G", microsatellite.unitString());
        Assert.assertEquals(6, microsatellite.numRepeat);

        microsatellite = microsatellites.get(2);
        Assert.assertEquals("chr1", microsatellite.genomeRegion.Chromosome);
        Assert.assertEquals(38, microsatellite.genomeRegion.start());
        Assert.assertEquals(45, microsatellite.genomeRegion.end());
        Assert.assertEquals("AT", microsatellite.unitString());
        Assert.assertEquals(4, microsatellite.numRepeat);
    }

    @Test
    public void testMicrosatellitesMinAdjacentDist()
    {
        // the first 3 microsatellites are too close to each other, they are rejected
        String sequence = "TGATCATCATCATCATCCGGGGGGAAAAAATCGATCAATATATATTCG";
        TestReferenceSequenceFile refSeqFile = new TestReferenceSequenceFile(sequence);

        List<RefGenomeMicrosatellite> microsatellites = RefGenomeMicrosatellitesFinder.findMicrosatellites(refSeqFile, 4, 20);
        Assert.assertEquals(1, microsatellites.size());
        RefGenomeMicrosatellite microsatellite = microsatellites.get(0);
        Assert.assertEquals("chr1", microsatellite.genomeRegion.Chromosome);
        Assert.assertEquals(38, microsatellite.genomeRegion.start());
        Assert.assertEquals(45, microsatellite.genomeRegion.end());
        Assert.assertEquals("AT", microsatellite.unitString());
        Assert.assertEquals(4, microsatellite.numRepeat);
    }

    @Test
    public void testIsValidUnit()
    {
        // test that we are able to weed out repeat units that are actually just multiple of smaller units
        Assert.assertTrue(RefGenomeMicrosatellitesFinder.isValidUnit(StringUtil.stringToBytes("ATA")));
        Assert.assertFalse(RefGenomeMicrosatellitesFinder.isValidUnit(StringUtil.stringToBytes("ATAT")));
        Assert.assertFalse(RefGenomeMicrosatellitesFinder.isValidUnit(StringUtil.stringToBytes("ATATAT")));
        Assert.assertFalse(RefGenomeMicrosatellitesFinder.isValidUnit(StringUtil.stringToBytes("AA")));
        Assert.assertFalse(RefGenomeMicrosatellitesFinder.isValidUnit(StringUtil.stringToBytes("TTT")));
        Assert.assertFalse(RefGenomeMicrosatellitesFinder.isValidUnit(StringUtil.stringToBytes("GGGG")));
        Assert.assertFalse(RefGenomeMicrosatellitesFinder.isValidUnit(StringUtil.stringToBytes("GGGGG")));
        Assert.assertTrue(RefGenomeMicrosatellitesFinder.isValidUnit(StringUtil.stringToBytes("ATT")));
        Assert.assertTrue(RefGenomeMicrosatellitesFinder.isValidUnit(StringUtil.stringToBytes("ATTAT")));
        Assert.assertTrue(RefGenomeMicrosatellitesFinder.isValidUnit(StringUtil.stringToBytes("ATGAT")));
        Assert.assertFalse(RefGenomeMicrosatellitesFinder.isValidUnit(StringUtil.stringToBytes("ATTATT")));
        Assert.assertFalse(RefGenomeMicrosatellitesFinder.isValidUnit(StringUtil.stringToBytes("TTATTA")));
        Assert.assertFalse(RefGenomeMicrosatellitesFinder.isValidUnit(StringUtil.stringToBytes("ATGATG")));
    }

    //@Test
    public void testReal() throws FileNotFoundException
    {
        IndexedFastaSequenceFile refGenome = new IndexedFastaSequenceFile(new File("/data/resources/bucket/reference_genome/38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"));

        String seq = refGenome.getSubsequenceAt("chr1", 5931400, 5931420).getBaseString();

        //List<RefGenomeMicrosatellite> refGenomeMicrosatellites = RefGenomeMicrosatelliesFinder.findHomopolymer(refGenome, 10);
    }
}