#!/bin/bash

destination="./"  #The user has to type their desired destination to store the raw data and processing files otherwise the script will run in the current destination

cp -R /localdisk/data/BPSM/ICA1 $destination # The first step is going to be copy the raw data files to users local machine
 
mkdir $destination/processing #Creating a directory to store the fastqc analysis files
fastqc --extract -t 10 $destination/ICA1/fastq/* #This command will run FASTQC analysis on all the files in 'fastq' directory and create 2 files per each file (html. and .zip)

find . -type f -iname "*.html" -delete # We do not need the html files, so we can delete them straight away
find . -type f -iname "*.zip" -delete
mv $destination/ICA1/fastq/*fastqc $destination/processing # Move the  folders into a separate directory


#This block will create a text file which will containt all the PASS/FAIL results from fastqc analysis and delete the source directory
cd $destination/processing
touch fastqc_summaries.txt
for dir in */;
do
	cat ./$dir/summary.txt >> ./fastqc_summaries.txt
	rm -rf ./$dir
done

#NEXT STEP, print some sort of analysis count TRUE/FAILS etc.
#There are 10lines for each sequence, we are mostly interested in "Per Base Sequence Quality" and "Overrepresented sequences"
while true; do
    read -p "Do you wish to check the 'per base sequence quality? " yn
    case $yn in
        [Yy]* ) 
		#Code to count pass/fail/warn for "per base seq quality'
			count_p=0
			count_f=0
			count_w=0
			i=0
			line_interest=2
			while read line
			do
				[ -z "$line" ] && continue
				i=$(( $i + 1 ))
				if test $i == $line_interest 
				then
					if test ${line:0:1} == "P"
					then
					count_p=$(( $count_p+1 ))
					elif test ${line:0:1} == "F"
					then
					count_f=$(( $count_f+1 ))
					elif test ${line:0:1} == "W"
					then
					count_w=$(( $count_w+1 ))
					fi
					line_interest=$(( $line_interest + 10 ))	
				fi
			done < fastqc_summaries.txt
			echo "For Per base sequence quality there are:"
			echo "$count_p Passes"
			echo "$count_f Fails"
			echo "$count_w Warnings"
		break;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac
done

while true; do
    read -p "Do you wish to check for the 'Overrepresented sequences? " yn
    case $yn in
        [Yy]* ) 
		#Code to count pass/fail/warn for "Overrepresented sequences'
			count_p=0
			count_f=0
			count_w=0
			i=0
			line_interest=9
			while read line
			do
				[ -z "$line" ] && continue
				i=$(( $i + 1 ))
				if test $i == $line_interest
				then
					if test ${line:0:1} == "P"
					then
					count_p=$(( $count_p+1 ))
					elif test ${line:0:1} == "F"
					then
					count_f=$(( $count_f+1 ))
					elif test ${line:0:1} == "W"
					then
					count_w=$(( $count_w+1 ))
					fi
					line_interest=$(( $line_interest + 10 ))	
				fi
			done < fastqc_summaries.txt
			echo "For the overrepresented sequences there are:"
			echo "$count_p Passes"
			echo "$count_f Fails"
			echo "$count_w Warnings"
		break;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac
done
