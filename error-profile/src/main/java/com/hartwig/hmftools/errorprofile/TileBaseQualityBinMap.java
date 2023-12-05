package com.hartwig.hmftools.errorprofile;

import static htsjdk.samtools.util.SequenceUtil.A;
import static htsjdk.samtools.util.SequenceUtil.C;
import static htsjdk.samtools.util.SequenceUtil.G;
import static htsjdk.samtools.util.SequenceUtil.T;

import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicIntegerArray;
import java.util.concurrent.atomic.AtomicReferenceFieldUpdater;

// a map that maps the TileBaseQualityBin to values
public class TileBaseQualityBinMap
{
    static class NonDenseKey
    {
        public final String flowcell;
        public final short lane;
        public final short tile;
        public final byte rawBaseQuality;
        //public final boolean firstOfPair;

        public NonDenseKey(final String flowcell, final short lane, final short tile, final byte rawBaseQuality)//, final boolean firstOfPair)
        {
            this.flowcell = flowcell;
            this.lane = lane;
            this.tile = tile;
            this.rawBaseQuality = rawBaseQuality;
            //this.firstOfPair = firstOfPair;
        }

        @Override
        public boolean equals(final Object o)
        {
            if(this == o)
            {
                return true;
            }
            if(!(o instanceof NonDenseKey))
            {
                return false;
            }

            final NonDenseKey that = (NonDenseKey) o;

            if(lane != that.lane)
            {
                return false;
            }
            if(tile != that.tile)
            {
                return false;
            }
            if(rawBaseQuality != that.rawBaseQuality)
            {
                return false;
            }
            /*if(firstOfPair != that.firstOfPair)
            {
                return false;
            }*/
            return flowcell.equals(that.flowcell);
        }

        @Override
        public int hashCode()
        {
            // from my experiments, using 127 seems marginally better than using 33 or 37
            int result = flowcell.hashCode();
            result = 127 * result + (int) lane;
            result = 127 * result + (int) tile;
            result = 127 * result + (int) rawBaseQuality;
            //result = 127 * result + (firstOfPair ? 1 : 0);
            return result;
        }
    }

    static class DenseBins
    {
        volatile AtomicIntegerArray counts;

        private static final AtomicReferenceFieldUpdater<TileBaseQualityBinMap.DenseBins, AtomicIntegerArray> FIELD_UPDATER =
                AtomicReferenceFieldUpdater.newUpdater(TileBaseQualityBinMap.DenseBins.class, AtomicIntegerArray.class, "counts");

        public DenseBins(final int maxReadLength)
        {
            // 2 for R1 R2 and 2 for error and total count
            counts = new AtomicIntegerArray(maxReadLength * 4);
        }

        public int getTotalCount(boolean firstOfPair, int readPosition)
        {
            return counts.get(calcBinIndex(counts, firstOfPair, readPosition, false));
        }

        public int getErrorCount(boolean firstOfPair, int readPosition)
        {
            return counts.get(calcBinIndex(counts, firstOfPair, readPosition, true));
        }

        // this function must try to be thread safe while allowing read position resizing
        public void incrementCount(boolean firstOfPair, int readPosition, boolean errorCount)
        {
            addToCount(firstOfPair, readPosition, errorCount, 1);
        }

        public void addToCount(boolean firstOfPair, int readPosition, boolean errorCount, int count)
        {
            int countToAdd = count;
            counts.getAndAdd(calcBinIndex(counts, firstOfPair, readPosition, errorCount), countToAdd);
        }

        public void addToCountOld(boolean firstOfPair, int readPosition, boolean errorCount, int count)
        {
            int countToAdd = count;

            while(true)
            {
                AtomicIntegerArray countsCopy = FIELD_UPDATER.get(this);

                if(countsCopy == null)
                {
                    // another thread is resizing the array
                    continue;
                }

                int binIndex = calcBinIndex(countsCopy, firstOfPair, readPosition, errorCount);

                if(binIndex >= countsCopy.length())
                {
                    // resize array, first set the field to null
                    if(FIELD_UPDATER.compareAndSet(this, countsCopy, null))
                    {
                        // move the old counts to the new counts
                        int oldMaxReadLength = countsCopy.length() / 4;
                        int newMaxReadLength = oldMaxReadLength * 2;

                        // we managed to set it to null, we are in control
                        AtomicIntegerArray newCounts = new AtomicIntegerArray(newMaxReadLength * 4);

                        for(boolean r1 = false;; r1 = true)
                        {
                            for(boolean e = false;; e = true)
                            {
                                for(int readPos = 0; readPos < oldMaxReadLength; ++readPos)
                                {
                                    // we atomically set the count in the old array to 0, to make sure any race condition
                                    // where the count is changed after this operation can be accounted for
                                    int c = countsCopy.getAndSet(calcBinIndex(countsCopy, r1, readPos, e), 0);
                                    if(c != 0)
                                    {
                                        newCounts.getAndAdd(calcBinIndex(newCounts, r1, readPos, e), c);
                                    }
                                }
                                if(e)
                                {
                                    break;
                                }
                            }
                            if (r1)
                            {
                                break;
                            }
                        }

                        // completed the copy, make the array available for all threads
                        FIELD_UPDATER.set(this, newCounts);
                    }

                    // now should be resized, try the add again
                    continue;
                }

                countsCopy.getAndAdd(binIndex, countToAdd);

                if(countsCopy != counts)
                {
                    // the counts array has been changed, this means a resizing has occurred. We need to try to move the count to the new one
                    // if it hasn't been yet
                    countToAdd = countsCopy.getAndSet(binIndex, 0);

                    if(countToAdd != 0)
                    {
                        // the resizing thread would have used atomic compare and swap to set the moved counts to 0
                        // if there is any remaining count here, it means the count was added after the resizing thread
                        // move this value across. We have to therefore try to add this count again to the new array
                        continue;
                    }
                }
                break;
            }
        }

        public int maxReadLength()
        {
            return counts.length() / 4;
        }

        // error count is the second half, each array can store 4
        private static int calcBinIndex(AtomicIntegerArray counts, boolean firstOfPair, int readPosition, boolean errorCount)
        {
            int maxReadLength = counts.length() / 4;
            return (errorCount ? maxReadLength * 2 : 0) + (firstOfPair ? 0 : maxReadLength) + readPosition;
        }
    }

    int mMaxReadPosition = 151;
    Map<NonDenseKey, DenseBins> mTileBaseQualityCountMap = new ConcurrentHashMap<>(2000, 0.75f);

    public Map<TileBaseQualityBin, BaseQualityBinCounter.Count> toTileBaseQualityMap()
    {
        // convert it to a normal map
        Map<TileBaseQualityBin, BaseQualityBinCounter.Count> resultMap = new HashMap<>();

        for(Map.Entry<NonDenseKey, DenseBins> entry : mTileBaseQualityCountMap.entrySet())
        {
            NonDenseKey key = entry.getKey();
            DenseBins denseBins = entry.getValue();

            // get all the dense bin keys
            for(int readPosition = 0; readPosition < denseBins.maxReadLength(); ++readPosition)
            {
                for(boolean firstOfPair : new boolean[]{ true, false })
                {
                    // final String flowcell, final int lane, final int tile, final boolean firstOfPair, final int readPosition,
                    //            // final byte ref, final byte alt,
                    //            final byte rawBaseQuality
                    TileBaseQualityBin tileBaseQualityBin = new TileBaseQualityBin(key.flowcell, key.lane, key.tile, firstOfPair, readPosition, key.rawBaseQuality);

                    BaseQualityBinCounter.Count c = new BaseQualityBinCounter.Count();
                    c.setTotalCount(denseBins.getTotalCount(firstOfPair, readPosition));
                    c.setErrorCount(denseBins.getErrorCount(firstOfPair, readPosition));

                    if (c.getTotalCount() > 0)
                    {
                        // omit 0 counts
                        resultMap.put(tileBaseQualityBin, c);
                    }
                }
            }
        }

        return resultMap;
    }

    public void incrementCount(TileBaseQualityBin bin, boolean errorCount)
    {
        NonDenseKey key = toNonDenseKey(bin);

        if(bin.readPosition > mMaxReadPosition)
        {
            // we need to reorganise all the submaps
            // increaseMaxReadPosition(bin.readPosition);
        }

        DenseBins denseBins = mTileBaseQualityCountMap.computeIfAbsent(key, (k) -> new DenseBins(mMaxReadPosition));
        denseBins.incrementCount(bin.firstOfPair, bin.readPosition, errorCount);
    }

    /*
    public BaseQualityBinCounter.Count get(TileBaseQualityBin bin)
    {
        NonDenseKey key = toNonDenseKey(bin);

        if(bin.readPosition > mMaxReadPosition)
        {
            // we need to reorganise all the submaps
            // increaseMaxReadPosition(bin.readPosition);
        }

        // now get the T
        DenseBin[] denseBins = mTileBaseQualityCountMap.get(key);

        if(denseBins == null)
        {
            return null;
        }

        return denseBins[(computeIndex(bin.firstOfPair, bin.readPosition, bin.ref, bin.alt, mMaxReadPosition))].count;
    }*/

    /*
    public BaseQualityBinCounter.Count getOrCreate(TileBaseQualityBin bin)
    {
        NonDenseKey key = toNonDenseKey(bin);

        if(bin.readPosition > mMaxReadPosition)
        {
            // we need to reorganise all the submaps
            // increaseMaxReadPosition(bin.readPosition);
        }

        List<DenseBin> denseBins = mTileBaseQualityCountMap.computeIfAbsent(key, (k) -> new ArrayList<>(computeNumDenseBins(mMaxReadPosition)));

        int index = bin.readPosition;

        // TODO: fix the race condition here
        while(denseBins.size() <= index)
        {
            denseBins.add(new DenseBin((short)denseBins.size()));
        }
        return denseBins.get(index);
    }
     */

    private static NonDenseKey toNonDenseKey(TileBaseQualityBin tileBaseQualityBin)
    {
        return new NonDenseKey(tileBaseQualityBin.flowcell, tileBaseQualityBin.lane, tileBaseQualityBin.tile,
                tileBaseQualityBin.rawBaseQuality);
    }

    /*
        size_t computeIndex(alwaysT_t<Dims, std::size_t>... indexes_args) const
    {
        constexpr std::size_t dimensions[] = {Dims...};
        std::size_t indexes[] = {indexes_args...};

        size_t index = 0;
        size_t mul = 1;

        for (size_t i = 0; i != sizeof...(Dims); ++i) {
            assert(indexes[i] < dimensions[i]);
            index += indexes[i] * mul;
            mul *= dimensions[i];
        }
        assert(index < (Dims * ...));
        return index;
    }
     */

    private static int computeNumDenseBins(int maxReadPosition)
    {
        return maxReadPosition;
    }

    /*
    synchronized private void increaseMaxReadPosition(int newMaxReadPosition)
    {
        for(Map.Entry<NonDenseKey, DenseBin[]> entry : mTileBaseQualityCountMap.entrySet())
        {
            DenseBin[] oldDenseBins = entry.getValue();
            DenseBin[] newDenseBins = new DenseBin[computeNumDenseBins(newMaxReadPosition)];

            // copy old to new
            for(int firstOfPair = 0; firstOfPair <= 1; ++firstOfPair)
            {
                for(int readPosition = 0; readPosition <= mMaxReadPosition; ++readPosition)
                {
                    for(int ref = 0; ref < 4; ++ref)
                    {
                        for(int alt = 0; alt < 4; ++alt)
                        {
                            newDenseBins[computeIndex()] = oldDenseBins[computeIndex()];
                        }
                    }
                }
            }
        }
    }*/

    private static int baseToIndex(byte base)
    {
        switch(base)
        {
            case A: return 0;
            case T: return 1;
            case C: return 2;
            case G: return 3;
            default: throw new IllegalArgumentException();
        }
    }
}
