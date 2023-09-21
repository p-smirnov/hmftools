package com.hartwig.hmftools.errorprofile;

import java.util.ArrayList;
import java.util.List;

import org.jetbrains.annotations.Nullable;

// a class to store information regarding tandem repeat
// in a read. Tandem repeat has to be 2 to 4 bases repeated
// more than twice
public class TandemRepeat
{
    int startIndex = -1;
    int endIndex = -1;
    String pattern = "";
    int numRepeats = 0;

    int length() { return endIndex - startIndex; }
}
