package com.hartwig.hmftools.errorprofile;

public class ErrorProfileConstants
{
    public static int MIN_POLY_G_TAIL_COUNT = 4;
    public static int MIN_HOMOPOLYMER_LENGTH = 4;

    public static int MIN_TANDEM_REPEAT_COUNT = 3;

    public static int BASE_QUAL_CUTOFF = 30;

    public static double STRAND_BIAS_CUTOFF = 0.1;

    public static double READ_POS_SMOOTH_SIGMA = 4;

    // 4 points in phred space
    public static double BASEQ_OUTLIER_CUTOFF = 4;

    public static int BASEQ_OUTLIER_MEDIAN_COUNT = 5;
}
