#!/bin/sh
echo Run the moldy autorun script.
echo Made by Jungmin Kim, 20200620 
echo --------------------------------------------
for index in $(find $PWD -mindepth 1 -maxdepth 1 -type d)
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
sleep 5
jobs
echo --------------------------------------------
echo The moldy autorun script done.

