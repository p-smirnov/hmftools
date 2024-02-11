package com.hartwig.hmftools.errorprofile.repeat;

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

    public boolean select(RepeatAnalyser repeatAnalyser)
    {
        if(units != null)
        {
            return units.contains(repeatAnalyser.refGenomeMicrosatellite.unitString());
        }
        else if(unitLength != null)
        {
            return repeatAnalyser.refGenomeMicrosatellite.unit.length == unitLength;
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
