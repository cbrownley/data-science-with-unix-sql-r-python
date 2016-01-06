#!/bin/bash

START_TIME=$(date +%s)

TODAY=$(date +'%F')

echo "observation,sepal_length,sepal_width,petal_length,petal_width,species" > output_files/iris_from_sqlite_"$TODAY".csv

# create and query a sqlite database table and save the output in a file
sqlite3 input_files/DataScience.db >> output_files/iris_from_sqlite_"$TODAY".csv <<EndOfScript
	SELECT * FROM iris;
EndOfScript

perl -i -pe 's/\|/,/g' output_files/iris_from_sqlite_"$TODAY".csv

END_TIME=$(date +%s)

# $? refers to the exit status of the previous command
# 0 is success; failure statuses range from 1 to 255
# therefore the if statement says, if error code is not equal to success
# everything between <<MESSAGE and MESSAGE is a Here Document
# $(()) syntax carries out the calculation and supplies the result of the calculation
if [ $? -ne 0 ]
then
`mail -s "Data Pull Failed" clinton@clinton-mba.local << MESSAGE
The data pull failed on $TODAY.
MESSAGE`
else
`mail -s "Data Pull Success" clinton@clinton-mba.local << MESSAGE
The data pull completed successfully on $TODAY.
The process took $(($(($END_TIME - $START_TIME))/60)) minutes and $(($(($END_TIME - $START_TIME))%60)) seconds.
MESSAGE`
fi
