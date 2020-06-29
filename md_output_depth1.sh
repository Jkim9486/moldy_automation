#!/bin/sh
echo Run the moldy rdf,xyz script.
echo Made by Jungmin Kim, 20200620

pwd_int=$PWD

for index in $(find $pwd_int -mindepth 1 -maxdepth 1 -type d)
do
    cd $index
    pwd
    md_output_onebyall.py > /dev/null &
    echo moldy in operation ...
    echo --------------------------------------------
done
