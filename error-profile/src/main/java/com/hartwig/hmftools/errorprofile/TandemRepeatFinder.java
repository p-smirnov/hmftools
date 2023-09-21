package com.hartwig.hmftools.errorprofile;

import java.util.ArrayList;
import java.util.List;
import java.util.ListIterator;

import org.jetbrains.annotations.Nullable;

public class TandemRepeatFinder
{
    static class Candidate
    {
        int startIndex = 0;
        int endIndex = 0;
        String pattern;

        Candidate(String pattern, int startIndex, int endIndex)
        {
            this.pattern = pattern;
            this.startIndex = startIndex;
            this.endIndex = endIndex;
        }

        char nextBase()
        {
            return pattern.charAt((endIndex - startIndex) % pattern.length());
        }

        int length()
        {
            return endIndex - startIndex;
        }

        int numRepeats()
        {
            return length() / pattern.length();
        }
    }

    public static @Nullable TandemRepeat findTandemRepeat(String seq)
    {
        @Nullable
        Candidate bestCandidate = null;

        // candidates
        List<Candidate> currentCandidates = new ArrayList<>();

        for(int i = 0; i < seq.length(); ++i)
        {
            char base = seq.charAt(i);

            // check all current candidates
            ListIterator<Candidate> itr = currentCandidates.listIterator();
            while(itr.hasNext())
            {
                Candidate c = itr.next();

                // check if pattern is still valid
                if(c.nextBase() == base)
                {
                    // pattern continues
                    c.endIndex++;
                    if(bestCandidate == null || bestCandidate.length() < c.length())
                    {
                        bestCandidate = c;
                    }
                }
                else
                {
                    itr.remove();
                }
            }

            // also start new Candidates at this location
            for(int j = Math.max(i - 2, 0); j < i; ++j)
            {
                char b = seq.charAt(j);

                // check that this is a potential tandem repeat, rather than homopolymer
                for(int k = j + 1; k < i + 1; ++k)
                {
                    if(seq.charAt(k) != b)
                    {
                        // this is ok as a pattern
                        String pattern = seq.substring(j, i + 1);
                        Candidate c = new Candidate(pattern, j, i + 1);
                        currentCandidates.add(c);
                    }
                }
            }
        }

        if(bestCandidate != null && bestCandidate.numRepeats() >= ErrorProfileConstants.MIN_TANDEM_REPEAT_COUNT)
        {
            TandemRepeat tandemRepeat = new TandemRepeat();
            tandemRepeat.pattern = bestCandidate.pattern;
            tandemRepeat.startIndex = bestCandidate.startIndex;
            tandemRepeat.endIndex = bestCandidate.endIndex;
            tandemRepeat.numRepeats = bestCandidate.length() / bestCandidate.pattern.length();
            return tandemRepeat;
        }

        return null;
    }
}
