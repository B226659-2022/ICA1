
destination=~/Assesments #The user has to type their desired destination to store the raw data and processing files

#cp -R /localdisk/data/BPSM/ICA1 $destination # The first step is going to be copy the raw data files to users local machine
 

mkdir $destination/processing #Creating a directory to store the fastqc files
fastqc -t 10 $destination/ICA1/fastq/* #This command will run FASTQC analysis on all the files in 'fastq' directory and create 2 files per each file (html. and .zip)

find . -type f -iname "*.html" -delete # We do not need the html files, so we can delete them straight away
mv $destination/ICA1/fastq/*fastqc.zip $destination/processing # Move the .zip files of interest into a separate directory


#The following lines will unzip the folders and delete the original .zip files
cd $destination/processing
for zip in *.zip
do
	unzip $zip
done
find . -type f -iname "*.zip" -delete 



touch fastqc_summaries.txt
for folder in */
do
cd $folder
cat summary.txt > ../fastqc_summaries.txt
cd ..
done





