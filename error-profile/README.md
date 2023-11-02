For homopolymer insert / delete, do we need to make sure it is added to the
start or end of the read? Right now it is done by genome position, straight from cigar,
which might or might not be the best idea.

Perhaps we need to "rejig" cigar to remove the strandedness.
Another way is to say this is a homopolymer ref length, alt length.

In a single point mutation, ref and alt would be different one base. In a delete, ref would be a sequence, alt would be nothing.
Or we could tag it on to the previous base. That however would probably make it harder for machine learning to deal with

Impplementation details

# Tile Position Adjustment

For each {tile,R1R2,BaseQual} we add an additional adjustment determined by fitting 6 parameters empirically in phred qual space: 

| Parameter         | Description                                                | Calculation                                   |     |
|-------------------|------------------------------------------------------------|-----------------------------------------------|-----|
| tileAdj           | Relative adjustment for tile as a whole                    | avgQual(tile)- avgQual(sample                 |     |
| minPos            | Base at which the error rate is assumed to be at a minimum | Base(1-2nd last base) with ~lowest error rate |     |
| minPosAdj♯♯       | Relative adjustment at last base                           | medianQual(minPos+/-2) - avgQual(tile)        |     |
| StartAdjPerBase♯♯ | Increase in adjustment per base from minPosition to base0  | Slope (1->minPos)♯                            |     |
| EndAdjPerBase♯♯   | Increase in adjustment per base from minPosition to        | Slope (minPos-> 2nd lastBase)♯                |     |
♯ Note that outliers are defined as bases with empirical base qual more than 4 points below the median of the 5 positions either side and is
ALT specific. The 1st and last base are always treated as outliers and get a base and ALT specific adjustment

♯♯ Set to 0 if any positions are missing observations.

Per tile and base qual fitted model parameters are output to a tsv file. A specific adjustment for identified outliers is stored in a
separate file keyed by {tile,R1R2,position,baseQual} and is calculated as avgQual(pos) - avgQual(tile). In SAGE the positional adjustment is
applied based on the tile id and a lookup of the parameters in addition to the standard BQR context based adjustment. For outliers, the tile
model is overridden with the specific outlier adjustment.

Perhaps better way is to use linear regression to find the regression line

See ![Input TILE_ERROR_RATE](doc/tile_empirical_bq.svg)
