#!/bin/sh
echo Run the moldy autorun script.
echo Made by Jungmin Kim, 20200620 
echo --------------------------------------------
pwd_int=$PWD
for index in $(find $pwd_int -mindepth 1 -maxdepth 1 -type d)
do
    cd $index
    pwd
    rm -f MDBCK*
    rm -f ah_debug
    rm -f *dump*
    g++ *.cpp -o create_file_step.out
    ./create_file_step.out
    chmod 755 execution_script
    ./execution_script &
    echo --------------------------------------------
done
sleep 1
jobs
wait
echo --------------------------------------------
echo The moldy autorun script done.
echo --------------------------------------------
echo Run the moldy rdf,xyz script.
echo Made by Jungmin Kim, 20200620
for index in $(find $pwd_int -mindepth 1 -maxdepth 1 -type d)
do
    cd $index
    pwd
    md_output_onebyall.py > /dev/null &
    echo moldy in operation ...
    echo --------------------------------------------
done
