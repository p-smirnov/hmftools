package com.hartwig.hmftools.errorprofile.repeat;

import static htsjdk.samtools.util.SequenceUtil.G;
import static htsjdk.samtools.util.SequenceUtil.T;

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

public class RefGenomeHomopolymerFinderTest
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

        List<RefGenomeHomopolymer> polymers = RefGenomeHomopolymerFinder.findHomopolymer(refSeqFile, 4, 30);

        Assert.assertEquals(2, polymers.size());
        RefGenomeHomopolymer polymer = polymers.get(0);
        Assert.assertEquals("chr1", polymer.genomeRegion.Chromosome);
        Assert.assertEquals(12, polymer.genomeRegion.start());
        Assert.assertEquals(17, polymer.genomeRegion.end());
        Assert.assertEquals(G, polymer.base);
        Assert.assertEquals(6, polymer.numRepeat);

        polymer = polymers.get(1);
        Assert.assertEquals("chr1", polymer.genomeRegion.Chromosome);
        Assert.assertEquals(27, polymer.genomeRegion.start());
        Assert.assertEquals(33, polymer.genomeRegion.end());
        Assert.assertEquals(T, polymer.base);
        Assert.assertEquals(7, polymer.numRepeat);

        // try different chunk size
        polymers = RefGenomeHomopolymerFinder.findHomopolymer(refSeqFile, 4, 5);

        Assert.assertEquals(2, polymers.size());
        polymer = polymers.get(0);
        Assert.assertEquals("chr1", polymer.genomeRegion.Chromosome);
        Assert.assertEquals(12, polymer.genomeRegion.start());
        Assert.assertEquals(17, polymer.genomeRegion.end());
        Assert.assertEquals(G, polymer.base);
        Assert.assertEquals(6, polymer.numRepeat);

        polymer = polymers.get(1);
        Assert.assertEquals("chr1", polymer.genomeRegion.Chromosome);
        Assert.assertEquals(27, polymer.genomeRegion.start());
        Assert.assertEquals(33, polymer.genomeRegion.end());
        Assert.assertEquals(T, polymer.base);
        Assert.assertEquals(7, polymer.numRepeat);
    }

    //@Test
    public void testReal() throws FileNotFoundException
    {
        IndexedFastaSequenceFile refGenome = new IndexedFastaSequenceFile(new File("/data/resources/bucket/reference_genome/38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"));

        String seq = refGenome.getSubsequenceAt("chr1", 5931400, 5931420).getBaseString();

        List<RefGenomeHomopolymer> refGenomeHomopolymers = RefGenomeHomopolymerFinder.findHomopolymer(refGenome, 10);
    }
}