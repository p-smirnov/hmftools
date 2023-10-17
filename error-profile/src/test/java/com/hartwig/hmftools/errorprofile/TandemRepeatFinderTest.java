package com.hartwig.hmftools.errorprofile;

import org.junit.Test;

import junit.framework.TestCase;

public class TandemRepeatFinderTest
{
    @Test
    public void testFindTandemRepeat()
    {
        TandemRepeat tandemRepeat = TandemRepeatFinder.findTandemRepeat("ATATATATCATAT");

        TestCase.assertEquals("AT", tandemRepeat.pattern);
        TestCase.assertEquals(8, tandemRepeat.length());

        tandemRepeat = TandemRepeatFinder.findTandemRepeat("ATGATGATGATCATATAT");

        TestCase.assertEquals("ATG", tandemRepeat.pattern);
        TestCase.assertEquals(11, tandemRepeat.length());

        tandemRepeat = TandemRepeatFinder.findTandemRepeat("TTTTTATGATGATGATCAGAGAGAGAGAGT");

        TestCase.assertEquals("AG", tandemRepeat.pattern);
        TestCase.assertEquals(12, tandemRepeat.length());
    }
}
