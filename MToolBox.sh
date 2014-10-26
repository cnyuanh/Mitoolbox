#!/usr/bin/env bash

check_exit_status()
{
rc=$?
if [[ $rc != 0 ]]
then
	echo ""
	echo "The last process reported an error. Exit."
	exit $rc
else
	echo "Success."
	echo ""
fi
}

usage()
{
	USAGE="""
	MToolBox: a tool for heteroplasmy annotation and accurate functional analysis of mitochondrial variants from high throughput sequencing data.
	Written by Domenico Simone and Claudia Calabrese 2013-2014.
	https://sourceforge.net/projects/mtoolbox/

	You must run the MToolBox command within the folder with input files (either bam or sam or fastq).

	Input & workflow execution options (must include -i):

		-i	input file format. Mandatory argument [bam|sam|fastq|fasta]
		-m	options for mapExome script [see mapExome.py -h for details]
		-M	remove duplicate reads with PicardTools MarkDuplicates after mapExome [default: no]
		-I	perform local realignment of reads on known indels with GATK IndelRealigner [default: no]
		-a	options for assembleMTgenome script [see assembleMTgenome.py -h for details]
		-c	options for mt-classifier script [see mt-classifier.py -h for details]
		-r	reference sequence to use for read mapping (VCF output will use the same reference) [RSRS|rCRS; default: RSRS]

	Help options:

		-h	show this help message
		-v	show version

	"""
	echo "$USAGE"
}

version()
{
	VERSION=$(echo "MToolBox v0.2")
	echo $VERSION
}

# Default command lines and behaviours for scripts and programs used in the workflow
#assembleMTgenome_OPTS=""
#mt_classifier_OPTS=""
#mapExome_OPTS=""
UseMarkDuplicates=false
UseIndelRealigner=false
# export folder where MToolBox.sh is placed, it is the same folder of PicardTools and GATK jars
me=`basename $0`
export mtoolbox_folder=$(which $me | sed "s/$me//g")
export externaltoolsfolder=${mtoolbox_folder}ext_tools/

# Environment variables for executables and files required by MToolBox
export ref=RSRS
export fasta_path=/usr/local/share/genomes/
export mtdb_fasta=chrRSRS.fa
export hg19_fasta=hg19RSRS.fa
export gsnapexe=/usr/local/bin/gsnap
export gsnapdb=/usr/local/share/gmapdb/
export mtdb=chrRSRS
export humandb=hg19RSRS
export samtoolsexe=/usr/local/bin/samtools
export muscleexe=/usr/local/bin/muscle


while getopts ":hva:c:f:i:m:r:MI" opt; do
	case $opt in
		h)
			usage
			exit 1
			;;
		v)
			version
			exit 1
			;;
		a)
			assembleMTgenome_OPTS=$OPTARG
			;;
		c)
			mt_classifier_OPTS=$OPTARG
			;;
		f)
			variants_functional_annotation_OPTS=$OPTARG
			;;
		i)
			input_type=$OPTARG
			;;
		m)
			mapExome_OPTS=$OPTARG
			;;
		r)
			ref=$(echo $OPTARG | tr '[:lower:]' '[:upper:]')
			;;
		M)
			UseMarkDuplicates=true
			;;
		I)
			UseIndelRealigner=true
			;;
		\?)
			echo "Invalid option: -$OPTARG" >&2
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires an argument." >&2
			exit 1
			;;
	esac
done

# define reference
if [[ $ref == 'RCRS' ]]
then 
	export mtdb_fasta=chrRCRS.fa
	export hg19_fasta=hg19RCRS.fa
	export mtdb=chrRCRS
	export humandb=hg19RCRS
elif [[ $ref != 'RSRS' ]]
then
	echo "Reference name not valid. Abort."
	exit 1
fi


# The following lines are commented since the involved parameters
# are specified elsewhere
#
# Set thresholds for hf and tail
#export hfthreshold=0.8
#export taillength=7
#

#if (( $taillength < 5 ))
#then
#	echo "Minumum tail length required >= 5. Tail length will be set to 5."
#	export taillength=5
#i

# Check python version (2.7 required)
echo ""
echo "Check python version... (2.7 required)"
min=$(python -c "import sys; print (sys.version_info[:])[1]")
maj=$(python -c "import sys; print (sys.version_info[:])[0]")
if [[ $maj != 2 ]] || [[ $min != 7 ]]
then
echo "You need Python2.7 in order to run MToolBox. Abort."
exit 1
else
echo "OK."
echo ""
fi


# Check existence of files to be used in the pipeline; if any of them does not exist, the pipeline will be aborted.
echo "Checking files to be used in MToolBox execution..."

#-t ${mt_classifier_OPTS} \

check_files.py \
--assembleMTgenome_OPTS="${assembleMTgenome_OPTS}" \
--mapExome_OPTS="${mapExome_OPTS}" \
--mt_classifier_OPTS="${mt_classifier_OPTS}"
# Check exit status of check_files.py
rc=$?
if [[ $rc != 0 ]] ; then
	exit $rc
fi

# Function definition
fastq_input()
{ # run mapExome directly.
	# get unique list of sample IDs
	# sampleIDs=$(ls *fastq | awk 'BEGIN{FS="."}{count[$1]++}END{for (j in count) print j}')

	# map against mt genome and human genome
	# for i in $sampleIDs; do datasets=$(echo $i.*fastq); mapExome_RSRS_SamHeader.py -g ${gsnapexe} -D ${gsnapdb} -M ${mtdb} -H ${humandb} -a "${datasets}" -o OUT_${i}; done &> log_mapexome.txt
	echo ""
	echo "##### EXECUTING READ MAPPING WITH MAPEXOME..."
	echo ""
	sampleIDs=$(ls *fastq | awk 'BEGIN{FS="."}{count[$1]++}END{for (j in count) print j}')
	for i in $sampleIDs; do
		echo "mapExome for sample" ${i}", files found:" $(ls $i.*fastq)
		if (( $(ls $i.*fastq | wc -l) == 1 ))
		then
			#echo $i is 1
			mapExome.py -g ${gsnapexe} -D ${gsnapdb} -M ${mtdb} -H ${humandb} -a $i.fastq -o OUT_${i} ${mapExome_OPTS}
		elif (( $(ls $i.*fastq | wc -l) == 2 ))
		then
			#echo $i is 2
			mapExome.py -g ${gsnapexe} -D ${gsnapdb} -M ${mtdb} -H ${humandb} -a $i.R1.fastq -b $i.R2.fastq -o OUT_${i} ${mapExome_OPTS}
		elif (( $(ls $i.*fastq | wc -l) == 3 ))
		then
			if [ -s $i.fastq ]
			then 
				mapExome.py -g ${gsnapexe} -D ${gsnapdb} -M ${mtdb} -H ${humandb} -a $i.R1.fastq -b $i.R2.fastq -c $i.fastq -o OUT_${i} ${mapExome_OPTS}
			else
				rm $i.fastq
				echo "$i.fastq is an empty unpaired fastq. File has been removed."
				mapExome.py -g ${gsnapexe} -D ${gsnapdb} -M ${mtdb} -H ${humandb} -a $i.R1.fastq -b $i.R2.fastq -o OUT_${i} ${mapExome_OPTS}
			fi	
		#then
			#echo $i is 3
			#mapExome.py -g ${gsnapexe} -D ${gsnapdb} -M ${mtdb} -H ${humandb} -a $i.R1.fastq -b $i.R2.fastq -c $i.fastq -o OUT_${i} ${mapExome_OPTS}
		else (( $(ls $i.*fastq | wc -l) > 3 ))
			echo "$i not processed. Too many files."
			:
		fi
	done

	echo ""
	echo "SAM files post-processing..."
	echo ""
	# SORT SAM WITH PICARD TOOLS
	echo "##### SORTING OUT.sam FILES WITH PICARDTOOLS..."
	echo ""
	for i in $(ls -d OUT_*); do cd ${i}; java -Xmx4g \
	-Djava.io.tmpdir=`pwd`/tmp \
	-jar ${externaltoolsfolder}SortSam.jar \
	SORT_ORDER=coordinate \
	INPUT=OUT.sam \
	OUTPUT=OUT.sam.bam \
	TMP_DIR=`pwd`/tmp; cd ..; done
	check_exit_status
	# INDEXING BAM FILES WITH SAMTOOLS
	for i in $(ls -d OUT_*); do cd ${i}; ${samtoolsexe} index OUT.sam.bam; cd ..; done
	
	# REALIGN KNOWN INDELS WITH GATK
	if $UseIndelRealigner
	then
		echo ""
		echo "##### REALIGNING KNOWN INDELS WITH GATK INDELREALIGNER..."
		echo ""
		for i in $(ls -d OUT_*); do cd ${i}; \
		echo "Realigning known indels for file" ${i}"/OUT.sam.bam using" ${mtoolbox_folder}"data/MITOMAP_HMTDB_known_indels.vcf as reference..."
		java -Xmx4g \
		-Djava.io.tmpdir=`pwd`/tmp \
		-jar ${externaltoolsfolder}GenomeAnalysisTK.jar \
		-T IndelRealigner \
		-R ${mtoolbox_folder}/data/chr${ref}.fa \
		-I OUT.sam.bam \
		-o OUT.realigned.bam \
		-targetIntervals ${mtoolbox_folder}/data/intervals_file_${ref}.list  \
		-known ${mtoolbox_folder}/data/MITOMAP_HMTDB_known_indels_${ref}.vcf \
		-compress 0;
		check_exit_status; cd ..; done
	else
		for i in $(ls -d OUT_*); do cd ${i}; cat OUT.sam.bam > OUT.realigned.bam; cd ..; done
	fi
	# MARK DUPLICATES WITH PICARD TOOLS (DEFAULT: YES)
	if $UseMarkDuplicates
	then
		echo ""
		echo "##### ELIMINATING PCR DUPLICATES WITH PICARDTOOLS MARKDUPLICATES..."
		echo ""
		for i in $(ls -d OUT_*); do cd ${i}; java -Xmx4g \
		-Djava.io.tmpdir=`pwd`/tmp \
		-jar ${externaltoolsfolder}MarkDuplicates.jar \
		INPUT=OUT.realigned.bam \
		OUTPUT=OUT.sam.bam.marked.bam \
		METRICS_FILE=OUT.sam.bam.metrics.txt \
		ASSUME_SORTED=true \
		REMOVE_DUPLICATES=true \
		TMP_DIR=`pwd`/tmp; cd ..; done
	else
		for i in $(ls -d OUT_*); do cd ${i}; cat OUT.realigned.bam > OUT.sam.bam.marked.bam; cd ..; done
	fi
	# RE-CONVERT BAM OUTPUT FROM MARKDUPLICATES IN SAM.
	for i in $(ls -d OUT_*); do cd ${i}; java -Xmx4g -Djava.io.tmpdir=`pwd`/tmp -jar ${externaltoolsfolder}SamFormatConverter.jar INPUT=OUT.sam.bam.marked.bam OUTPUT=OUT.sam.bam.marked.bam.marked.sam TMP_DIR=`pwd`/tmp; cd ..; done

	for i in $(ls -d OUT_*); do cd ${i}; grep -v "^@" *marked.sam > OUT2.sam; mkdir MarkTmp; mv OUT.sam.bam MarkTmp; mv OUT.sam.bam.marked.bam MarkTmp; mv OUT.sam.bam.marked.bam.marked.sam MarkTmp; tar -czf MarkTmp.tar.gz MarkTmp; rm -R MarkTmp/; cd ..; done

	# ASSEMBLE CONTIGS, GET MT-TABLES...
	echo ""
	echo "##### ASSEMBLING MT GENOMES WITH ASSEMBLEMTGENOME..."
	echo ""
	echo "WARNING: values of tail < 5 are deprecated and will be replaced with 5"
	echo ""
	for i in $(ls -d OUT_*); do outhandle=$(echo ${i} | sed 's/OUT_//g'); cd ${i}; assembleMTgenome.py -i OUT2.sam -o ${outhandle} -r ${fasta_path} -f ${mtdb_fasta} -a ${hg19_fasta} -s ${samtoolsexe} -FCP ${assembleMTgenome_OPTS}; cd ..; done > logassemble.txt
	echo ""
	echo "##### GENERATING VCF OUTPUT..."
	# ... AND VCF OUTPUT
	VCFoutput.py -r ${ref}
}

fasta_input()
{
	if [[ $input_type = 'fasta' ]]
	then
		echo ""
		echo "##### PRE-PROCESSING OF FASTA INPUT FILES..."
		echo ""
		echo "Files to be analyzed:"
		for i in $(test_fasta.py); do bname=$(echo ${i} | awk 'BEGIN{FS="."}{print $1}'); bname_dir=OUT_${bname}; mkdir ${bname_dir}; cp ${i} ${bname_dir}/${bname}-contigs.fasta; echo ${bname}; done
		#for i in $(ls); do bname=$(echo ${i} | awk 'BEGIN{FS="."}{print $1}'); bname_dir=OUT_${bname}; mkdir ${bname_dir}; cp ${i} ${bname_dir}/${bname}-contigs.fasta; done
		check_exit_status
	fi
	echo ""
	echo "##### PREDICTING HAPLOGROUPS AND ANNOTATING/PRIORITIZING VARIANTS..."
	echo ""
	#### Haplogroup prediction and functional annotation
	# Brand new haplogroup prediction best file
	hpbest="mt_classification_best_results.csv" # change just this name for changing filename with most reliable haplogroup predictions
	echo "Best haplogroup predictions will be written in mt_classification_best_results.csv"
	echo "SampleID,Best predicted haplogroup(s)" > ${hpbest}
	for i in $(ls -d OUT_*); do inhandle=$(echo ${i} | sed 's/OUT_//g'); cd ${i}; mt-classifier.py -i ${inhandle}-contigs.fasta -s ${hpbest} -b ${inhandle} -m ${muscleexe} ${mt_classifier_OPTS}; cd ..; done

	# Functional annotation of variants
	#for i in $(ls -d OUT_*); do cd $i; variants_functional_annotation.py $hpbest ; cd ..; done
	variants_functional_annotation.py #${hpbest}
}

sam_input()
{ # convert sam to fastq and run mapExome.
	for i in $(ls *.sam); do echo "Converting sam to fastq..." ${i}; n=$(echo $i | awk 'BEGIN{FS="."}{print $1}'); java -Xmx4g \
		-Djava.io.tmpdir=`pwd`/tmp \
		-jar ${externaltoolsfolder}SamToFastq.jar \
		INPUT=$n.sam \                                              
		FASTQ=$n.R1.fastq \
		SECOND_END_FASTQ=$n.R2.fastq \
		UNPAIRED_FASTQ=$n.fastq \
		VALIDATION_STRINGENCY=SILENT \
		TMP_DIR=`pwd`/tmp; echo "Done."; done
#echo "Converting sam input(s) to fastq"
#echo "Done."
}

bam_input()
{ # convert bam to fastq and run mapExome.
	for i in $(ls *.bam); do echo "Converting bam to fastq..." ${i}; n=$(echo $i | awk 'BEGIN{FS="."}{print $1}'); java -Xmx4g \
		-Djava.io.tmpdir=`pwd`/tmp \
		-jar ${externaltoolsfolder}SamToFastq.jar \
		INPUT=$n.bam \
		FASTQ=$n.R1.fastq \
		SECOND_END_FASTQ=$n.R2.fastq \
		UNPAIRED_FASTQ=$n.fastq \
		VALIDATION_STRINGENCY=SILENT \
		TMP_DIR=`pwd`/tmp; echo "Done."; done
#echo "Converting bam input(s) to fastq"
#echo "Done."
}


if (( $# >= 1 ))
then
	if [[ $input_type = 'fasta' ]]
	then
		echo "Input type is fasta."
		fasta_input
	elif [[ $input_type = 'fastq' ]]
	then
		echo "Input type is fastq."
		fastq_input
		fasta_input
	elif [[ $input_type = 'sam' ]]
	then
		echo "Input type is sam."
		sam_input
		fastq_input
		fasta_input
	elif [[ $input_type = 'bam' ]]
	then
		echo "Input type is bam."
		bam_input
		fastq_input
		fasta_input
	else
		echo "Input format not recognized."
		exit 1
	fi
else
	echo "Input type not specified."
	exit 1
fi
