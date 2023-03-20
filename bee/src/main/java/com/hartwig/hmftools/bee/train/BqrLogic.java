package com.hartwig.hmftools.bee.train;

import java.util.List;

import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationFile;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationRecord;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

// reimplementation of the sage BQR logic
// a little simplified for this purpose
public class BqrLogic
{
    private static final Logger sLogger = LogManager.getLogger(BqrLogic.class);

    private final IndexedFastaSequenceFile mRefGenome;
    private final QualityRecalibrationMap mQualityRecalibrationMap;

    public BqrLogic(String qualityRecalibrationFilePath,
            IndexedFastaSequenceFile refGenome)
    {
        List<QualityRecalibrationRecord> bqrRecords = QualityRecalibrationFile.read(qualityRecalibrationFilePath);
        assert bqrRecords != null;
        mQualityRecalibrationMap = new QualityRecalibrationMap(bqrRecords);
        mRefGenome = refGenome;
    }

    public double[] calcBqrQualities(SAMRecord record)
    {
        sLogger.trace("read({}) cigar({})", record, record.getCigarString());

        double[] bqrQualList = new double[record.getReadLength()];

        for (int i = 0; i < record.getReadLength(); ++i)
        {
            bqrQualList[i] = record.getBaseQualities()[i];
        }

        List<AlignmentBlock> alignmentBlocks = record.getAlignmentBlocks();

        for (AlignmentBlock currentBlock : alignmentBlocks)
        {
            // 1 based to 0 based
            int readIndex = currentBlock.getReadStart() - 1;

            //
            byte[] refGenomeBases = mRefGenome.getSubsequenceAt(
                    record.getReferenceName(),
                    currentBlock.getReferenceStart() - 1, // get one more base for trinucleotide context
                    currentBlock.getReferenceStart() + currentBlock.getLength() + 1)
                    .getBases();

            byte[] trinucleotideContext = new byte[3];

            assert(currentBlock.getLength() + 2 <= refGenomeBases.length);

            // set the block
            for (int i = 0; i < currentBlock.getLength(); ++i)
            {
                byte rawQuality = record.getBaseQualities()[readIndex + i];
                byte ref = refGenomeBases[i + 1];
                byte alt = record.getReadBases()[readIndex + i];

                System.arraycopy(refGenomeBases, i, trinucleotideContext, 0, 3);
                double q = mQualityRecalibrationMap.quality(ref, alt, trinucleotideContext, rawQuality);
                bqrQualList[i] = q;

                if (!Doubles.equal(rawQuality, q))
                {
                    sLogger.trace("alt({}), ref({}), trinucleotideContext({}), raw qual({}), bqr qual({})",
                            (char) alt, (char) ref, new String(trinucleotideContext), rawQuality, q);
                }
            }
        }

        return bqrQualList;
    }
}
