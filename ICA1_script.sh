#!/bin/bash
echo "################ Welcome #####################"
echo "Copying the raw data files to your current directory"
destination=$(pwd)  #The user has to type their desired destination to store the raw data and processing files otherwise the script will run in the current destination





# 1. 
cp -R /localdisk/data/BPSM/ICA1 $destination # The first step is going to be copy the raw data files to users local machine
 
mkdir $destination/ICA1/processing #Creating a directory to store the fastqc analysis files
fastqc --extract -t 10 $destination/ICA1/fastq/* #This command will run FASTQC analysis on all the files in 'fastq' directory and create 2 files per each file (html. and .zip)

find . -type f -iname "*.html" -delete # We do not need the html files, so we can delete them straight away
find . -type f -iname "*.zip" -delete
mv $destination/ICA1/fastq/*fastqc $destination/ICA1/processing # Move the  folders into a separate directory


#This block will create a text file which will containt all the PASS/FAIL results from fastqc analysis and delete the source directory
cd $destination/ICA1/processing
touch fastqc_summaries.txt
for dir in */;
do
	cat ./$dir/summary.txt >> ./fastqc_summaries.txt
	rm -rf ./$dir
done








# 2.
#NEXT STEP, print some sort of analysis count TRUE/FAILS etc.
#There are 10lines for each sequence, we are mostly interested in "Per Base Sequence Quality" and "Overrepresented sequences"
# I wanted to make it interactive, to the user can choose for themselves what analysis they want to run. The second line is arguably the most important as it shows the confidence in each base, but the overrepresented sequences are a good way to check for any the contaminations.
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
					line_interest=$(( $line_interest + 10 ))	# 10 lines per sequence, therefore needs to jump +10 lines to keep looking at the per base quality stuff
				fi
			done < fastqc_summaries.txt
      echo "########################################"
			echo "For Per base sequence quality there are:"
			echo "$count_p Passes"
			echo "$count_f Fails"
			echo "$count_w Warnings"
      echo "########################################"
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
      echo "########################################"
			echo "For the overrepresented sequences there are:"
			echo "$count_p Passes"
			echo "$count_f Fails"
			echo "$count_w Warnings"
      echo "########################################"
		break;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac
done
cd ..
pwd
read -p "Press enter to start the alignment process" # I need to give time for the reader to analyse the results, so the script does not continue running straight away, but because I want this just in case if the user clicks on Yes for the Overrepresented sequences, I have to include it in the case do








# 3.
# Next step is to index the Tcongo reference files so the aligning software can read it
echo "Indexing the reference file"
# MAKE SURE STARTING DIRECTORY IS CORRECT!!!!!!
bowtie2-build ./Tcongo_genome/* ./fastq/bowtie

# Now let's align the files
cd ./fastq
ct=0
no_files=$(ls -l *.fq.gz | wc -l) #Because it is a fairly long process, I want to show the user the progress and how many more files are left
for file in *gz #For every file with the ending .gz I will run this loop
  do
    if test ${file:9:10} == "1.fq.gz" # I only need to run it once per file name (Tco-XXXX), not for both file 1 and 2. Therefore I need to skip the loop for every other file ending with .fq.gz
    then
      ct=$(( $ct+1 ))
      echo "Aligning the "${file:0:8}" files ("${ct}"/"${no_files}")"  
      bowtie2 --very-fast-local -x ./bowtie -1 ${file:0:9}'1.fq.gz' -2 ${file:0:9}'2.fq.gz' | samtools view -Sb -o ${file:0:9}.bam #Pipeline to avoid storing the .sam files on local disk, instead let's keep it in RAM to save time   
    fi 
done

mv *.bam $destination/ICA1/processing #Move all the output files to a separate directory to keep it tidy and so the user can always find all the output files in one place




############ DELETE ############ 



