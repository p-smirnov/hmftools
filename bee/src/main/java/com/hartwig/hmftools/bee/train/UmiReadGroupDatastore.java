package com.hartwig.hmftools.bee.train;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class UmiReadGroupDatastore
{
    private static final Logger sLogger = LogManager.getLogger(UmiReadGroupDatastore.class);

    // in order to allow groups that are completed to be processed immediately, we keep them separate by
    // chromosome
    // chromosome -> position -> duplicate read groups
    Map<String, TreeMap<Integer, Collection<UmiReadGroup>>> mDupReadsMap = new HashMap<>();

    synchronized public UmiReadGroup getOrCreate(String umi, String chromosome, int startPosition, boolean firstOfPair)
    {
        TreeMap<Integer, Collection<UmiReadGroup>> positionMap = mDupReadsMap.computeIfAbsent(chromosome, k -> new TreeMap<>());

        Collection<UmiReadGroup> umiReadGroups = positionMap.computeIfAbsent(startPosition, k -> new ArrayList<>());

        // find that umi
        for (UmiReadGroup group : umiReadGroups)
        {
            if (group.firstOfPair == firstOfPair && group.umi.equals(umi))
            {
                return group;
            }
        }

        // cannot find it, make group
        UmiReadGroup group = new UmiReadGroup(umi, chromosome, startPosition, firstOfPair);
        umiReadGroups.add(group);
        return group;
    }

    synchronized public Collection<UmiReadGroup> popCompleted(String chromosome, int uptoPosition)
    {
        // get all the groups that are completed
        TreeMap<Integer, Collection<UmiReadGroup>> positionMap = mDupReadsMap.get(chromosome);

        if (positionMap == null)
        {
            return Collections.emptyList();
        }

        int previousSize = positionMap.size();

        List<UmiReadGroup> groups = new ArrayList<>();

        for (var entry : positionMap.entrySet())
        {
            if (entry.getKey() >= uptoPosition)
            {
                break;
            }

            groups.addAll(entry.getValue());
        }

        // remove all positions that are strictly less than uptoPosition
        positionMap.headMap(uptoPosition).clear();

        sLogger.debug("pop completed, previous size: {}, now size: {}, completed groups: {}", previousSize, positionMap.size(), groups.size());

        return groups;
    }
}
