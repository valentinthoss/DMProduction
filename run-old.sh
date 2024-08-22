#!/bin/bash
file=$1
i=$2
cp PBH.ini "Input/PBH_"$i".ini"
while read -r key
do
    read -r value
    sed -i "" "s|$key|$value|g;" "Input/PBH_"$i".ini"
done < $file
sed -i "" "s|OUTPUTNAME|$i|g;" "Input/PBH_"$i".ini"
cd Class
./class "../Input/PBH_"$i".ini"
