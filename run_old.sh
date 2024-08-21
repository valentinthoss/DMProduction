#!/bin/bash
i=$1
m=$2
omegancdm=$3
omegacdm=$4
T=$5
if [[ $# -eq 10 ]]
then
    a=$6
    b=$7
    c=$8
    d=$9
    e=${10}
    cat PBH.ini | sed "s|MNCDM|$m|g;s|OMEGANCDM|$omegancdm|g;s|OMEGACDM|$omegacdm|g;s|TNCDM|$T|g;s|OUTPUTNAME|$i|g;s|PAR1|$a|g;s|PAR2|$b|g;s|PAR3|$c|g;s|PAR4|$d|g;s|PAR5|$e|g;s|USEFILE|0|g" > "Input/PBH_"$i".ini"

else
    psd="../Input/psd_"$i".dat"
    cat PBH.ini | sed "s|PSDFILE|$psd|g;s|MNCDM|$m|g;s|OMEGANCDM|$omegancdm|g;s|OMEGACDM|$omegacdm|g;s|TNCDM|$T|g;s|OUTPUTNAME|$i|g;s|PAR1|0|g;s|PAR2|0|g;s|PAR3|0|g;s|PAR4|0|g;s|PAR5|0|g;s|USEFILE|1|g" > "Input/PBH_"$i".ini"
fi
cd Class
./class "../Input/PBH_"$i".ini"
