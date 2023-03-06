package com.hartwig.hmftools.bee.train;

import java.util.HashSet;
import java.util.Objects;
import java.util.Set;

// a simple class to help us avoid processing same read twice
public class ReadKeyTracker
{
    private static class ReadKey
    {
        public final String readName;
        public final boolean firstOfPair;

        public ReadKey(final String readName, boolean firstOfPair)
        {
            this.readName = readName;
            this.firstOfPair = firstOfPair;
        }

        @Override
        public boolean equals(final Object o)
        {
            if (this == o)
            {
                return true;
            }
            if (!(o instanceof ReadKey))
            {
                return false;
            }
            final ReadKey readKey = (ReadKey) o;
            return firstOfPair == readKey.firstOfPair && readName.equals(readKey.readName);
        }

        @Override
        public int hashCode()
        {
            return Objects.hash(readName, firstOfPair);
        }
    }

    private final Set<ReadKey> mReadKeys = new HashSet<>();

    public boolean tryAdd(String readName, boolean firstOfPair)
    {
        return mReadKeys.add(new ReadKey(readName, firstOfPair));
    }

    public boolean contains(ReadKey key)
    {
        return mReadKeys.contains(key);
    }
}
