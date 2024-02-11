package com.hartwig.hmftools.errorprofile.microsatellite;

import java.util.List;

public class MicrosatelliteSelector
{
    private final List<String> units;
    private final Integer unitLength;

    public MicrosatelliteSelector(List<String> units, Integer unitLength)
    {
        this.units = units;
        this.unitLength = unitLength;
    }

    public boolean select(MicrosatelliteSiteAnalyser microsatelliteSiteAnalyser)
    {
        if(units != null)
        {
            return units.contains(microsatelliteSiteAnalyser.refGenomeMicrosatellite.unitString());
        }
        else if(unitLength != null)
        {
            return microsatelliteSiteAnalyser.refGenomeMicrosatellite.unit.length == unitLength;
        }
        return false;
    }

    public String unitName()
    {
        if(units != null)
        {
            return String.join("/", units);
        }
        else if(unitLength != null)
        {
            return unitLength + "bp repeat";
        }

        return "";
    }
}
