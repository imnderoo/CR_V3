#!/usr/bin/python

import os, sys, subprocess
import shutil
import argparse #AW: Package for parsing arguments
import datetime, time
import re, glob

def main():
	#AW: This block of code is easy to reuse for any script that take in arguments
	#----
	parser = argparse.ArgumentParser(description = 'Script to read VCF into DB for ALBI project')
	parser.add_argument('in_folder', help = 'Full path to folder containing input FASTQs [or input file]')
	parser.add_argument('roi_bed_path', help = 'Full path to ROI BED - Should be in resources')
	parser.add_argument('out_folder_path', help = 'Full path to output folder')
	parser.add_argument('--prefix', help = 'Prefix to add to front of sample name. Ex: "RUNID_Sample.*"') 
	parser.add_argument('--start', help = 'Stage to start in [fastq|bam] DEFAULT: fastq') # AW: Optional argument: --start option
	parser.add_argument('--end', help = 'Stage to start in [bam|vcf] DEFAULT: vcf') 
	parser.add_argument('--intermediates', help = 'Keep Intermediates [Y|N] DEFAULT: N') 

	# TO-DO: Proper handling for starting in fastq/bam or ending in bam/vcf

	#parser.add_argument('--in_file1', help = 'Name of input file 1') # AW: Each line defines the input parameter and the name to use


	#Ex: GATK_pipeline.py /path/read1.fastq.gz /path/read2/fastq.gz /path/panel.bed /path/output_folder 
	#Ex2: GATK_pipeline.py /path/read1.fastq.gz /path/read2/fastq.gz /path/panel.bed /path/output_folder --end bam
	#Ex3: GATK_pipeline.py /path/read1.fastq.gz /path/read2/fastq.gz /path/panel.bed /path/output_folder --start bam
	#Ex4: GATK_pipeline.py /path/read1.fastq.gz /path/read2/fastq.gz /path/panel.bed /path/output_folder --start bam --prefix 160428_M01234_FLOWCE

	args = parser.parse_args() # AW: This line is needed for arguments to be recognized
	#----

	if args.start:
		print("Pipeline starting at " + args.start)
	else:
		print("Pipeline starting at fastq")

	if args.end:
		print("Pipeline ending at " + args.end)
	else:	
		print("Pipeline ending at vcf")
		
	if args.prefix:
		print("Pipeline prefix: " + args.prefix)

	if not os.path.exists(args.out_folder_path):
		os.makedirs(args.out_folder_path)

	#AbsPath_to_BedFile='/media/sf_Linux/Exome/resources/BEDfiles/breast_panel.bed'
	AbsPath_to_BedFile = args.roi_bed_path # AW: Note the usage of args.parameter_namE

	# Look in the input folder and extract R1 FASTQ as set of input file.
	input_file_set = []

	for fastq in os.listdir(args.in_folder):
		#print fastq
		if "_R1_" in fastq and not "Undetermined" in fastq:
			input_file_set.append(fastq)
			
	#print input_file_set

	# Overwite default if input start is set to bam
	if args.start == "bam":
		input_file_set = []
		for bam in os.listdir(args.in_folder):
			if ".bam" in bam:
				input_file_set.append(bam)

	# SET GLOBAL VARIABLES (ex. PATHS TO RESOURCES)
	AbsPathToResources=os.path.abspath('/media/sf_resources')
	GenomeFile = "ucsc.hg19.fa" # AW genome

	AbsPathToAnnovDBFld=AbsPathToResources + '/annovar_DBs/'
	AbsPathToPubDBs=AbsPathToResources+'/gatk_DBs/'
	AbsPathToGenome=AbsPathToResources + '/genome/' + GenomeFile

	#GET ABSOLUTE PATHS OF INPUT 
	args.out_folder_path = os.path.abspath(args.out_folder_path)
	args.in_folder = os.path.abspath(args.in_folder)
	args.roi_bed_path = os.path.abspath(args.roi_bed_path)
	OutputFolder = args.out_folder_path + "/"

	#DOCKER VOLUME MOUNTS

	if not os.path.isdir(args.out_folder_path):
		os.makedirs(args.out_folder_path)
		
	dockerVolume='-v /media/sf_resources:/media/sf_resources'
	dockerVolume=dockerVolume+' '+'-v '+args.out_folder_path+':'+args.out_folder_path
	dockerVolume=dockerVolume+' '+'-v '+args.in_folder+':'+args.in_folder
	#dockerVolume=dockerVolume+' '+'-v '+os.path.dirname(args.roi_bed_path)+':'+os.path.dirname(args.roi_bed_path)

	#DOCKER COMMANDS
	bwaCMD='docker run --rm '+dockerVolume+' nderoo324/bwa'
	samtoolsCMD='docker run --rm '+dockerVolume+' nderoo324/samtools'
	picardCMD='docker run --rm '+dockerVolume+' nderoo324/picard'
	annovarCMD='docker run --rm '+dockerVolume+' nderoo324/av'
	gatkCMD='docker run --rm '+dockerVolume+' nderoo324/gatk GenomeAnalysisTK'
	exomeDepthCMD = 'docker run --rm '+dockerVolume+' nderoo324/exomedepth'
	mantaCMD = 'docker run -id --name temp_manta '+dockerVolume+' nderoo324/manta'

	# START OF VCF CALLING SCRIPT
	start_time = time.asctime( time.localtime(time.time()) )
	"""

	sample_dict = {}
	elapsed_time_dict = {}

	for inFile in input_file_set:
		######## The following quotation marks should be filled with the FastQ names
		#FastQ_R1_FullPath = args.fastq_R1_path
		#FastQ_R2_FullPath = args.fastq_R2_path
		FastQ_R1_FullPath = args.in_folder + "/" + inFile
		FastQ_R2_FullPath = args.in_folder + "/" + inFile.replace("R1", "R2")
		
		FastQ_R1 = inFile #AW: os.path.basename leaves us with just the file name
		FastQ_R2 = inFile.replace("R1", "R2")

		####### SAMPLE specific info ############
		#Sample_Number=file1[0:8]
		#Sample_Number = "_".join(FastQ_R1.split("_")[:2]) #If we want the S# straign from illumina MiSeq. It splits by _, takes first 2, joins by _
		#Sample_Number=FastQ_R1.split("_")[0] #If we want the S# straign from illumina MiSeq
		#print ("\nSample Number: ",Sample_Number)

		### Parse FASTQ file using REGEX ###

		print inFile

		matches = re.match( r'(.*?)_L(.*?)_R.*', inFile)
		Sample_Number = matches.group(1)
		lane_id = "L" + matches.group(2)
		
		Platform="Illumina"
		Library="TruSightOne"
		PU_name="1"
		
		# Add prefix if defined
		if args.prefix:
			SampleLaneIDsegment = args.prefix + "_" + Sample_Number + "-" + lane_id
			SampleIDsegment = args.prefix + "_" + Sample_Number
		else:
			SampleLaneIDsegment = Sample_Number + "-" + lane_id
			SampleIDsegment = Sample_Number
			
		sample_start_time = datetime.datetime.now()
			
		print "Output: " + OutputFolder

		OutputBamSrtFile = OutputFolder + SampleLaneIDsegment + "-aln-pe-sorted.bam"

		#Define and vopen log file for rw
		log_date = datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d')
		OutputLog = OutputFolder+SampleIDsegment+"-pipeline-"+log_date+".log"

		if args.start:
			log_write(OutputLog, "Pipeline starting at " + args.start + "\n")
		else:
			log_write(OutputLog, "Pipeline starting at fastq\n") 
		if args.end:
			log_write(OutputLog, "Pipeline ending at " + args.end + "\n")
		else:	
			log_write(OutputLog, "Pipeline ending at vcf\n")	
		if args.prefix:
			log_write(OutputLog, "Pipeline prefix: " + args.prefix + "\n")

		#	print matches.group()
		#	print (sample_id, lane_id)
			
		## BWA ############

		if args.start != 'bam':
			print ("********************* RUNNING BWA on:",FastQ_R1," and ",FastQ_R2," *******************")
			
			log_write (OutputLog, '********************* RUNNING BWA on:",FastQ_R1," and ",FastQ_R2," *******************\n')
			#print (datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S'))
			#print (AbsPathToGenome)

			# For newest samtools - pipe is tough on docker / space
			
			tmpSAM = OutputFolder + SampleLaneIDsegment + ".tmp.sam"
			
			subprocess.call([bwaCMD+" mem -t 6 \
			-R '@RG\\tID:"+SampleLaneIDsegment+"\\tSM:"+Sample_Number+"\\tPL:"+Platform+"\\tLB:"+Library+"\\tPU:"+PU_name+"' \
			"+AbsPathToGenome+" "+FastQ_R1_FullPath+" "+FastQ_R2_FullPath+" > "+tmpSAM], shell=True)
			
			tmpBAM = OutputFolder + SampleLaneIDsegment + ".tmp.bam"
			subprocess.call([samtoolsCMD+" view -Shu -@ 5 "+tmpSAM+" > "+tmpBAM], shell=True) 
			
			subprocess.call([samtoolsCMD+" sort -@ 5 -o "+OutputBamSrtFile+" "+tmpBAM], shell=True)		
			
			os.remove(tmpSAM)
			os.remove(tmpBAM)
			
			print ("********************* BWA FINISHED. OUTPUT:",OutputBamSrtFile," *******************")
			
		if SampleIDsegment not in sample_dict:
			sample_dict[SampleIDsegment] = []
			elapsed_time_dict[SampleIDsegment] = datetime.datetime.now() - sample_start_time
		else:
			elapsed_time_dict[SampleIDsegment] = elapsed_time_dict[SampleIDsegment] + (datetime.datetime.now() - sample_start_time)		
			
		sample_dict[SampleIDsegment].append(SampleLaneIDsegment)
		
		print sample_dict
		
	## Per Sample Basis ##
	for SampleIDsegment in sample_dict:
		
		sample_start_time = datetime.datetime.now()
		
		OutputBamMrgFile = OutputFolder + SampleIDsegment + "-aln-pe-sorted-merged.bam"
		OutputBamSrtDeDupFile = OutputFolder+SampleIDsegment+"-aln-pe-sorted-merged-dedup.bam"
		OutputBamSrtDeDupRecalFile = OutputFolder+SampleIDsegment+"-aln-pe-sorted-merged-dedup-recal.bam"
		RawVariantsVCF = OutputFolder+SampleIDsegment+"_raw_variants.vcf"

		#Define and vopen log file for rw
		log_date = datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d')
		OutputLog = OutputFolder+SampleIDsegment+"-pipeline-"+log_date+".log"
		
		#print sample_dict[sample]
		# Merge the multiple lane files (if they exist)
		
		print SampleIDsegment
		
		if len(sample_dict[SampleIDsegment]) == 1:
			SampleLaneIDsegment = sample_dict[SampleIDsegment][0]
			PartBamSrtFile = OutputFolder + SampleLaneIDsegment + "-aln-pe-sorted.bam"
			
			print ("os.rename(" + PartBamSrtFile + ", " + OutputBamMrgFile + ")")
			
			os.rename(PartBamSrtFile, OutputBamMrgFile)
		else:
			partBAM = []
			
			for part in sample_dict[SampleIDsegment]:
				SampleLaneIDsegment = part
				PartBamSrtFile = OutputFolder + SampleLaneIDsegment + "-aln-pe-sorted.bam"
				partBAM.append(PartBamSrtFile)
				
			subprocess.call([samtoolsCMD+" merge "+OutputBamMrgFile+" "+" ".join(partBAM)],shell=True)
		
		## Picard MarkDuplicates
		print ("******************* RUNNING Picard MarkDuplicate on "+SampleIDsegment+" *****************")
		log_write (OutputLog, "******************* RUNNING Picard MarkDuplicate on "+SampleIDsegment+" *****************\n")

		subprocess.call(\
		[picardCMD+" MarkDuplicates \
		INPUT="+OutputBamMrgFile+"  \
		OUTPUT="+OutputBamSrtDeDupFile+"  \
		METRICS_FILE="+OutputFolder+"/"+SampleIDsegment+"_PCR_duplicates  \
		CREATE_INDEX=true \
		REMOVE_DUPLICATES=true 2>&1 | tee -a "+OutputLog]
		,shell=True)

		if "bam" not in str(args.end):

			## IndelRealigner - No longer necessary with Haplotype 

			############### Base Recalibration ###########################
			print ("***************** RUNNING GATK BaseRecalibrator on "+SampleIDsegment+" *********************")
			log_write (OutputLog, "***************** RUNNING GATK BaseRecalibrator on "+SampleIDsegment+" *********************\n")
			subprocess.call(
			gatkCMD+" -T BaseRecalibrator \
			-R "+AbsPathToGenome+" \
			-I "+OutputBamSrtDeDupFile+"  \
			-L "+args.roi_bed_path+" \
			-knownSites "+AbsPathToPubDBs+"1000Genomes/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
			-knownSites "+AbsPathToPubDBs+"1000Genomes/1000G_phase1.indels.hg19.sites.vcf \
			-knownSites "+AbsPathToPubDBs+"dbSNP/dbsnp_138.hg19.vcf \
			-o  "+OutputFolder+"/"+SampleIDsegment+"_recal_data.grp"+" 2>&1 | tee -a "+OutputLog\
			,shell=True) 

			# PrintReads -Step 2 of Base Recalibration
			print ("**************** RUNNING GATK PrintReads on "+SampleIDsegment+" *****************")
			log_write (OutputLog, "**************** RUNNING GATK PrintReads on "+SampleIDsegment+" *****************\n")
			subprocess.call( \
			gatkCMD+" -T PrintReads \
			-R "+AbsPathToGenome+" \
			-I "+OutputBamSrtDeDupFile+"  \
			-L "+args.roi_bed_path+" \
			-BQSR "+OutputFolder+"/"+SampleIDsegment+"_recal_data.grp  \
			-o  "+OutputBamSrtDeDupRecalFile+" 2>&1 | tee -a "+OutputLog \
			,shell=True)

			#################  HaplotypeCaller ###########################
			print ("**************** RUNNING GATK HaplotypeCaller on "+SampleIDsegment+" ********************")
			log_write (OutputLog, "**************** RUNNING GATK HaplotypeCaller on "+SampleIDsegment+" ********************\n")
			subprocess.call( \
			gatkCMD+ " -T HaplotypeCaller \
			-R "+AbsPathToGenome+" \
			-I "+OutputBamSrtDeDupRecalFile+" \
			-L "+args.roi_bed_path+" \
			--genotyping_mode DISCOVERY \
			-stand_emit_conf 10 \
			-stand_call_conf 30 \
			-mbq 20 \
			-o "+RawVariantsVCF+" 2>&1 | tee -a "+OutputLog \
			,shell=True)

			################ ANNOVAR #################################

			print ("**************** RUNNING ANNOVAR on "+SampleIDsegment+" ********************\n")
			log_write (OutputLog, "**************** RUNNING ANNOVAR on "+SampleIDsegment+" ********************\n")
			#VCF file as input to ANNOVAR table_annovar.pl
			# -protocol refGene,avsnp144,exac03nontcga,1000g2015aug_all,esp6500siv2_all,clinvar_20160302,nci60,cosmic70,ljb26_all \
			
			subprocess.call(\
			annovarCMD+" table_annovar.pl "+RawVariantsVCF+" "+AbsPathToAnnovDBFld+" -buildver hg19 \
			-protocol refGene,avsnp144,exac03nontcga,1000g2015aug_all,esp6500siv2_all,clinvar_20160302,ljb26_all \
			-operation g,f,f,f,f,f,f \
			-argument -hgvs,,,,,, \
			-out "+OutputFolder+"/"+SampleIDsegment+"_Annov_out  \
			-remove \
			-otherinfo \
			-nastring . \
			-vcfinput "+" 2>&1 | tee -a "+OutputLog \
			,shell=True)
			
			if not args.intermediates:
				for f in os.listdir (OutputFolder):
					if  re.search("-aln-pe-sorted.ba", f) or re.search("-aln-pe-sorted-merged.ba", f) or re.search("-aln-pe-sorted-merged-dedup.ba", f) or re.search("-aln-pe-sorted-merged-dedup-realigned.ba", f):
						os.remove(os.path.join(OutputFolder, f))

		
		elapsed_time_dict[SampleIDsegment] = elapsed_time_dict[SampleIDsegment] + (datetime.datetime.now() - sample_start_time)
		
		log_write (OutputLog, "Time Elapsed: " + str(elapsed_time_dict[SampleIDsegment]) + "\n")

	"""
	# START OF CNV CALLING SCRIPT
	bam_list = glob.glob(OutputFolder+"/*-aln-pe-sorted-merged-dedup-recal.bam")
	cnv_panel_csv = os.path.join(AbsPathToResources, "bed_files", "cnv_panel_exons.csv")
	call_cnv(bam_list, OutputFolder, cnv_panel_csv, AbsPathToGenome, exomeDepthCMD, mantaCMD, dockerVolume)
		
	end_time = time.asctime( time.localtime(time.time()) )
	print ("\nStart time: " + start_time + "\n")
	print ("End time: " + end_time + "\n")

def call_cnv(bam_list, out_dir, panel_csv, genome, exomeDepthCMD, mantaCMD, dockerVolume):

	"""
	# ExomeDepth - This will run all bams at once.
	exomedepth_out = os.path.join(out_dir,"cnv","exomedepth")
	
	if not os.path.isdir(exomedepth_out):
		os.makedirs(exomedepth_out)
	
	bam_list_file = os.path.join(exomedepth_out,"bamlist.txt")
	with (open(bam_list_file, 'w')) as f:
		for in_bam in bam_list:
			f.write(in_bam + "\n")
	
	subprocess.call(" ".join([exomeDepthCMD, \
		"-b", bam_list_file, \
		"-o", out_dir, \
		"-p", panel_csv, \
		"-s", "0.2"]), shell=True)
	
	"""
	
	# Manta - Call CNV for each bam in the bam list.
	
	for in_bam in bam_list:
		in_bam = bam_list[0]
		sample_id = in_bam.split("-aln-")[0]

		manta_out = os.path.join(out_dir,"cnv","manta",sample_id)
		
		if not os.path.isdir(manta_out):
			os.makedirs(manta_out)
		
		# Manta - Builds a python workflow.py
		subprocess.call(" ".join([mantaCMD, \
			"--bam="+in_bam, \
			"--referenceFasta="+genome, \
			"--runDir="+manta_out]), shell=True)
		
		# Manta - Workflow.py performs the actual analysis
			#mantaPythonCMD = 'docker exec temp_manta python'
			#subprocess.call('docker start temp_manta', shell=True)
			#subprocess.call(" ".join([mantaPythonCMD+manta_out+"/runWorkflow.py", "-m", "local"]), shell=True)
			#subprocess.call('docker rm temp_manta', shell=True)
			

			# Gunzip and copy results to main folder
			manta_results_dir = os.path.join(manta_out,"results","variants")
			subprocess.call(" ".join(["gunzip",manta_results_dir+"/*vcf.gz"]), shell=True)	
			shutil.copy(os.path.join(manta_results_dir,"diploidSV.vcf"), os.path.join(out_dir,sample_id+"manta.cnv.vcf"))
		
		# HACK - PUT THE STUPID MANTA PYHTON LIBRARIES INTO CREATE_REPORT
		# THEN, ADD A LINE HERE THAT ADDS PYTHON PATH TO IT (RE: THE FIRST FEW LINES OF CONFIGMANTA.PY
		# THEN, WE CAN USE OUR OWN LOCAL PYTHON. SO. STUPID.
		
		###### IGNORE THE RANT. BEST SOLUTION - CREATE WRAPPER SCRIPT FOR MANTA IN THE MANTA CONTAINER SO THAT THE WHOLE THING PROCS FROM ONE RUN (CONFIG, THEN WORKFLOW FROM WITHIN)
		###### LIKE EXOMEDEPTH ####
		
	# Still need to filter resulting CNV files by the sample test manifest in the CreateReport V3 script.

#Defining internal functions:
def log_write(log_file, msg):
	log_handle=open(log_file, 'a')
	log_handle.write (msg)
	log_handle.close()

main()
	
	
