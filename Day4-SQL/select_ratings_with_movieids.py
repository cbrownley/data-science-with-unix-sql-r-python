#!/usr/bin/env python
# from Day3-SQL folder, run
# ./select_ratings_with_movieids.py ../Day2-Unix/output_files/movies.csv ../Day2-Unix/output_files/ratings.csv ../Day2-Unix/output_files/ratings_with_movieids.csv

import csv
import sys

input_file1 = sys.argv[1]
input_file2 = sys.argv[2]
output_file = sys.argv[3]

movie_ids = []
with open(input_file1, 'r') as in_file1:
	filereader = csv.reader(in_file1)
	for row in filereader:
		movie_ids.append(row[0])

with open(input_file2, 'r') as in_file2:
	with open(output_file, 'w') as out_file:
		filereader = csv.reader(in_file2)
		filewriter = csv.writer(out_file)
		for row in filereader:
			if row[1] in movie_ids:
				filewriter.writerow(row)
			else:
				print(row)