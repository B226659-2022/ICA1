#!/bin/bash
destination=$(pwd)  #The user can type in their desired destination to store the raw data and processing files, otherwise the script do that in the current directory

# 0. We will start by copying the raw data files to users local machine
echo "Copying the raw data files to your current directory..."
cp -R /localdisk/data/BPSM/ICA1 $destination 

echo "!Download complete!"
echo "Running FastQC Analysis..."


# 1. Run the fastqc analysis
mkdir $destination/ICA1/processing #Creating a directory to store the fastqc analysis files

no_files=$(ls -l $destination/ICA1/fastq/*.gz | wc -l) #We need to get the number of files we are working with, to know how many fastqc processes to run in parallel
find $destination/ICA1/fastq/*.gz | parallel -j $no_files "fastqc --extract -t 20" #This command will run FASTQC analysis on all the files in 'fastq' directory and create 2 files per each file (html. and .zip)

# !!!Only use the next line if you don't have Parallel/GNU installed. If that applies to you, then replace the line before with this one:
# fastqc --extract -t 2 $destination/ICA1/fastq/*gz 

#find $destination/ICA1/fastq/*.html -delete # We do not need the html files, so we can delete them straight away
#find $destination/ICA1/fastq/*.zip -delete

mv $destination/ICA1/fastq/*fastqc $destination/ICA1/processing # Move the folders into a separate directory


#This block will create a text file which will containt all the PASS/FAIL results from fastqc analysis and delete the source directory
cd $destination/ICA1/processing
touch fastqc_summaries.txt
for dir in */;
do
	cat ./$dir/summary.txt >> ./fastqc_summaries.txt
	rm -rf ./$dir &
done

echo "########################################################"
echo "!FastQC analysis complete!"








# 2. Print the analysis (PASS/FAIL/WARN) so the user can check if that satisfies him enough to continue

#There are 10lines for each sequence, we are mostly interested in "Per Base Sequence Quality" and "Overrepresented sequences". The first will show confidence in the bases retrieved and the latter will let us know if there is any possible contamination.

while true; do
    read -p "Do you wish to check the 'per base sequence quality? (y/n): " yn 
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
					line_interest=$(( $line_interest + 10 ))	# 10 lines per sequence, therefore needs to jump +10 lines to keep looking at the per base quality results.
				fi
			done < fastqc_summaries.txt
      echo ""
			echo "!For Per base sequence quality there are:"
			echo "$count_p Passes"
			echo "$count_f Fails"
			echo "$count_w Warnings"
      echo ""
		break;;
        [Nn]* ) break;;
        * ) echo "Please answer yes or no.";;
    esac
done

while true; do
    read -p "Do you wish to check for the 'Overrepresented sequences? (y/n): " yn
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
      echo ""
			echo "!For the overrepresented sequences there are:!"
			echo "$count_p Passes"
			echo "$count_f Fails"
			echo "$count_w Warnings"
      echo ""
		break;;
        [Nn]* ) break;;
        * ) echo "Please answer yes or no.";;
    esac
done
cd ..
read -p "Press enter to start the alignment process" # I need to give time for the reader to analyse the results, so the script does not continue running straight away.










# 3. Next step is to index the Tcongo reference files so the aligning software can read it
echo "Indexing the reference file..."
bowtie2-build ./Tcongo_genome/* ./fastq/bowtie
echo "!Indexing of the bowtie2 reference file complete!"

# Now let's align the files with bowtie2
cd ./fastq

echo "Aligning the sequences..."
for file in *gz #For every file with the ending .gz I will run this loop
  do
    if test ${file:9:10} == "1.fq.gz" # I only need to run it once per file name (Tco-XXXX), not for both file 1 and 2. Therefore I need to skip the loop for every other file ending with .fq.gz
    then
      bowtie2 --very-fast-local -p32 -x ./bowtie -1 ${file:0:9}'1.fq.gz' -2 ${file:0:9}'2.fq.gz' | samtools view -Sb -o ${file:0:8}.bam & #Pipeline to avoid storing the .sam files on local disk, instead let's keep it in RAM to save time 
    fi
done
wait 

no_files=$(ls -l *1.fq.gz | wc -l)
# I decided to run the indexing on a different loop, so I can run bowtie2 in parallel, which will speed up the process significantly. (couldn't use a single loop as indexing would start on a file that does not exist yet)
for file in *gz #For every file with the ending .gz I will run this loop
  do
    if test ${file:9:10} == "1.fq.gz" # I only need to run it once per file name (Tco-XXXX), not for both file 1 and 2. Therefore I need to skip the loop for every other file ending with .fq.gz
    then
      echo "Sorting and indexing the .bam files..."
      samtools sort -@$no_files ${file:0:8}.bam > ${file:0:8}_sorted.bam 
      samtools index ${file:0:8}_sorted.bam
      
      mv *sorted.bam* ${destination}/ICA1 #Move all the output files to the same directory where the .bed file is. Preparing for the step 4.
         
    fi 
done











# 4. Generating counts data
cd $destination/ICA1
echo ""
echo "Generating counts data..."
find *bam | parallel 'bedtools intersect -a *.bed -b {} -c >> {}"_bed_res.txt"' #Running bedtools and -c option adds a tab separated column with the number of reads that align to the regions of the genome



# !!! Only use the next 'for' loop if you don't have Parallel/GNU installed. If that applies to you, then replace the line before with this loop:
#for f in *.bam
#do
#  bedtools intersect -a *.bed -b $f -c >> ${f:0:8}"_bed_res.txt" & 
#done














# 5.
#Prepare a tab-delimited text file (using awk) of their mean averages of expression i.e counts generated in step4. (the last column)

#For each group:
# Clone1 induced at a particular time (3x time periods)
# Clone2 induced at a particular time (3x time periods)
# WT at a particular time             (3x time periods)

#points to do:
# get names from description file of members of a group
# look for their files in the text file generated by bedtools
# get last column of each of those lines
# get the mean of the no of lines



# for loop for each line in Tco
# grep to find a text in a file
# save names in a file [group_name.txt]


# This next part is all about dividing the file containing information about the sequences into different groups. I've done it in a more complicated way, so it can still work if new clones are added in the future
while read line
do 
  awk -v clone=2 'BEGIN{FS="\t";} 
  {
     if($4 == "0")
     {print $0 > $clone"_0.txt";}
     else if($4 == "24")
      {print $0 > $clone"_24.txt";} 
     else if($4 == "48")
      {print $0 > $clone"_48.txt";}  
  }' 
    
done<$destination/ICA1/fastq/Tco.fqfiles

mv *.txt $destination/ICA1/processing
wait #This is to let processes from step 4 to finish, we will need them completed for the next step

# for each name (line in group_name.txt) grep a last column of the line with that name from (bedfile results)
# add those numbers up divide by the number of lines in the first loop

























