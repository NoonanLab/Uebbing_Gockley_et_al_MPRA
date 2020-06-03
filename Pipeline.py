#!/usr/bin/env python

#Written by Jake Gockley and Severin Uebbing, Yale University, 2017-2019
#Python script to run MPRA pipeline analysis on Ruddle HPC

from __future__ import division
from itertools import izip, izip_longest
import re
import sys
import os
import glob
import pysam
import argparse
import subprocess32 as subprocess
import warnings

#Read in command line arguments
def get_args():
	'''This function will parse and return command-line arguments'''
	# Assign description to the help doc
	parser = argparse.ArgumentParser(
		description='This Script Processes MPRA Read data')
	#Specify arguments
	##Run mode - NECESSARY to specify
	parser.add_argument(
		'-r', '--run', type=str, help="Specify Run Mode. Options: COMP, INERT, EXP, TAGComparison, or ANOTE", required=True)
	##Read groups Necessary in COMP, INERT, EXP, TAGComparison, and ANOTE modes
	parser.add_argument(
		'-R1', '--Read1', type=str, help='Read 1 file or comma delimited list of Read 1 files. Use a comma at the end of a single file path with * autocomplete flag', required=False)
	parser.add_argument(
		'-R2', '--Read2', type=str, help='Read 2 file or comma delimited list of Read 2 files. Use a comma at the end of a single file path with * autocomplete flag', required=False)
	##Sample name Necessary in EXP mode
	parser.add_argument(
		'-s', '--sample', type=str, help="Sample Type. Required in EXP mode ex. Rep2_1_pDNA", required=False)
	##Percent Identity Necessary in INERT Mode
	parser.add_argument(
		'-PI', '--PercentID', type=float, help="Percent ID Threshold. Must an decimal >= 0 and <= 1", required=False)
	##Anote Mode ***No two file names can be the same!!!
	parser.add_argument(
		'-I', '--INERTtag', type=str, help="Inert tag counts, required in ANOTE mode", required=False)
	parser.add_argument(
		'-C', '--COMPtag', type=str, help="Competent Library Seq Tag Counts Optional to specify in ANOTE mode", required=False)
	parser.add_argument(
		'-EP', '--PLAStag', type=str, help="Experimental Plasmid Tag Counts: can be a single file or comma delimited list of Read 1 files. Use a comma at the end of a single file path with * autocomplete flag", required=False)
	parser.add_argument(
		'-ET', '--TRANtag', type=str, help="Experimental Transcript Tag Counts: can be a single file or comma delimited list of Read 1 files. Use a comma at the end of a single file path with * autocomplete flag", required=False)
	parser.add_argument(
		'-R', '--Resource', type=str, help="File Path to Annotation Resource Files (Please End Path Directory with a Backslash!)", required=False)
	#Pipeline output file necessary in all modes			
	parser.add_argument(
		'-o', '--out', type=str, help='Specify the out file to store Pipeline Log.', required=True)
	# Barcode tag length
	parser.add_argument(
		'-t', '--taglength', type=int, help='Specify the length of the used barcode tags, required in INERT, INERT-HIQ, COMP, EXP modes.', required=False)
	#Array for all arguments passed to script
	args = parser.parse_args()
	#Assign args to variables
	Mode = args.run    	
	if args.Read1 is None:
		ROne = args.Read1
	else:
		ROne = args.Read1.split(",")
	if args.Read2 is None:
		RTwo = args.Read2
	else:
		RTwo = args.Read2.split(",")
	Samp = args.sample
	ID = args.PercentID
	out = args.out
	TagLength = args.taglength
	inertTAG = args.INERTtag
	compTAG = args.COMPtag
	if args.PLAStag is None:
		plasTAG = args.PLAStag
	else:
		plasTAG = args.PLAStag.split(",")
	if args.TRANtag is None:
		tranTAG = args.TRANtag
	else:
		tranTAG = args.TRANtag.split(",")
	resource = args.Resource
	# Return all variable values
	return Mode, ROne, RTwo, Samp, ID, out, inertTAG, compTAG, plasTAG, tranTAG, resource, TagLength

#Match return values from get_arguments() and assign to their respective variables
Mode, ROne, RTwo, Samp, ID, out, inertTAG, compTAG, plasTAG, tranTAG, resource, TagLength = get_args()

#Open the pipeline call output and print the call mode into it
PipeOUT = open(out, "w")
print >>PipeOUT, "Runing Mode: %s " % Mode

'''Define all Python functions'''
#INERT library-seq read strander function
##Original script call: python ~/Scripts/Strander.py HC_RQ5138_R1.fastq HC_RQ5138_R2.fastq
def Strander( Read1, Read2 ):
	##Counter variable to track conversation stats:
	AsIs = 0
	Reverse = 0
	Complement = 0
	RevComp = 0
	##Temp line variables
	R1_Name = ""
	R1_Seq = ""
	R1_Strand = ""
	R1_Qual = ""
	R2_Name = ""
	R2_Seq = ""
	R2_Strand = ""
	R2_Qual = ""
	#Open output files
	FivePrime = open("Inert_ReadOne.fastq", "w")
	ThreePrime = open("Inert_ReadTwo.fastq", "w")
	Garbage = open("Inert_Garbage_IDs.txt", "w")
	OUT = open("Inert_StranderOutput.txt", "w")
	#Complement dictionary
	Trans = {'A': 'T', 'T': 'A', 'C': 'G', 'G' : 'C', 'N' : 'N'}
	Assign = {1: "R1_Name = LINE1; R2_Name = LINE2", 2: "R1_Seq = LINE1; R2_Seq = LINE2", 3: "R1_Strand = LINE1; R2_Strand = LINE2", 4: "R1_Qual = LINE1; R2_Qual = LINE2"}
	FivePrimeMatch = 0
	ThreePrimeMatchR2 = 0
	ThreePrimeMatchR1 = 0
	FivePrimeMatchR2 = 0
	DoubleSeedMuts = 0
	LotsMuts = 0
	LotsMutsInverted = 0 
	ShortSeed = 0
	ShortSeedInv = 0
	TotCount = 0
	LineCounter=1
	#Loop through files simultaneously here
	#Don't need to izip longest because mate pair files are same length
	with open(Read1) as file1, open(Read2) as file2:
		for line1, line2 in izip(file1, file2):  
			line1 = line1.rstrip('\r\n')
			line2 = line2.rstrip('\r\n')
			Line1 = line1.rstrip()
			Line2 = line2.rstrip()
			#Split each line to determine value
			LINE1 = list(Line1)
			LINE2 = list(Line2)
			#Assign the temp FastQ values based on the commands stored in the array
			#This allows you to skip using
			exec(Assign[LineCounter])
			if(LineCounter < 4): #Name @
				LineCounter += 1
			else:
				#print ''.join(R1_Seq[0:35])
				LineCounter = 1
				TotCount += 1
				#Match 5' seed to read one: confirms correctly oriented fragments
				if (bool(re.search( 'CAGGTGCCAGAACATTTCTCT', ''.join(R1_Seq[0:35]))) is True):
					AsIs += 1
					FivePrimeMatch += 1     		
					print >>FivePrime, ''.join(R1_Name)+"\n",''.join(R1_Seq)+"\n",''.join(R1_Strand)+"\n",''.join(R1_Qual)
					print >>ThreePrime, ''.join(R2_Name)+"\n",''.join(R2_Seq)+"\n",''.join(R2_Strand)+"\n",''.join(R2_Qual)
				##Match 3' seed to read one: inverted oriented fragments
				elif (bool(re.search( 'CTGCTCGAAGCGGCCGGCC', ''.join(R1_Seq[0:35]))) is True):
					Reverse += 1
					ThreePrimeMatchR1 += 1
					print >>ThreePrime, ''.join(R2_Name)+"\n",''.join(R1_Seq)+"\n",''.join(R1_Strand)+"\n",''.join(R1_Qual)
					print >>FivePrime, ''.join(R1_Name)+"\n",''.join(R2_Seq)+"\n",''.join(R2_Strand)+"\n",''.join(R2_Qual)
				#Match 3' seed to read two: confirms correctly oriented fragments with mutation in 5' seed	
				elif (bool(re.search( 'CTGCTCGAAGCGGCCGGCC', ''.join(R2_Seq[0:35]))) is True):
					AsIs += 1
					ThreePrimeMatchR2 += 1
					print >>FivePrime, ''.join(R1_Name)+"\n",''.join(R1_Seq)+"\n",''.join(R1_Strand)+"\n",''.join(R1_Qual)
					print >>ThreePrime, ''.join(R2_Name)+"\n",''.join(R2_Seq)+"\n",''.join(R2_Strand)+"\n",''.join(R2_Qual)
				##Match 5' seed to read two: inverted oriented fragments with mutation in 3' seed	
				elif (bool(re.search( 'CAGGTGCCAGAACATTTCTCT', ''.join(R2_Seq[0:35]))) is True):
					Reverse += 1
					FivePrimeMatchR2 += 1
					print >>ThreePrime, ''.join(R2_Name)+"\n",''.join(R1_Seq)+"\n",''.join(R1_Strand)+"\n",''.join(R1_Qual)
					print >>FivePrime, ''.join(R1_Name)+"\n",''.join(R2_Seq)+"\n",''.join(R2_Strand)+"\n",''.join(R2_Qual)
				##Match 5' Sfi site: inverted fragments with lots of mutations 
				elif (bool(re.search( 'GGCCTAACTGGCC', ''.join(R2_Seq[20:45]))) is True):
					Reverse += 1
					LotsMutsInverted += 1
					print >>ThreePrime, ''.join(R2_Name)+"\n",''.join(R1_Seq)+"\n",''.join(R1_Strand)+"\n",''.join(R1_Qual)
					print >>FivePrime, ''.join(R1_Name)+"\n",''.join(R2_Seq)+"\n",''.join(R2_Strand)+"\n",''.join(R2_Qual)
				#Match backbone plus KpnI and XbaI seed to read one: confirms correctly oriented fragments with mutations in 5' AND 3' seeds
				elif (bool(re.search( 'CACTGCGGCTCCTGCGGTACCTCT', ''.join(R1_Seq[120:220]))) is True):
					AsIs += 1
					DoubleSeedMuts += 1
					print >>FivePrime, ''.join(R1_Name)+"\n",''.join(R1_Seq)+"\n",''.join(R1_Strand)+"\n",''.join(R1_Qual)
					print >>ThreePrime, ''.join(R2_Name)+"\n",''.join(R2_Seq)+"\n",''.join(R2_Strand)+"\n",''.join(R2_Qual)
				#Match 5' Sfi site: confirms correctly oriented fragments with lots of mutations
				elif (bool(re.search( 'GGCCTAACTGGCC', ''.join(R1_Seq[20:45]))) is True):
					AsIs += 1
					LotsMuts += 1
					print >>FivePrime,''.join(R1_Name)+"\n",''.join(R1_Seq)+"\n",''.join(R1_Strand)+"\n" ,''.join(R1_Qual)
					print >>ThreePrime,''.join(R2_Name)+"\n",''.join(R2_Seq)+"\n" ,''.join(R2_Strand)+"\n" ,''.join(R2_Qual)
				#Match shorter 5'seed: confirms correctly oriented fragments with lots of mutations
				elif (bool(re.search( 'CAGGTGCCAGAACA', ''.join(R1_Seq[3:23]))) is True):
					AsIs += 1
					ShortSeed += 1
					print >>FivePrime, ''.join(R1_Name)+"\n",''.join(R1_Seq)+"\n",''.join(R1_Strand)+"\n",''.join(R1_Qual)
					print >>ThreePrime, ''.join(R2_Name)+"\n",''.join(R2_Seq)+"\n",''.join(R2_Strand)+"\n",''.join(R2_Qual)
				##Match shorter 5'seed: confirms inverted fragments with lots of mutations
				elif (bool(re.search( 'CAGGTGCCAGAACA', ''.join(R2_Seq[3:23]))) is True):
					Reverse += 1
					ShortSeedInv += 1
					print >>ThreePrime, ''.join(R2_Name)+"\n",''.join(R1_Seq)+"\n",''.join(R1_Strand)+"\n",''.join(R1_Qual)
					print >>FivePrime, ''.join(R1_Name)+"\n",''.join(R2_Seq)+"\n",''.join(R2_Strand)+"\n",''.join(R2_Qual)
				else:
					print >>Garbage, ''.join(R1_Name)+"\t"+''.join(R1_Qual)+"\t"+''.join(R2_Qual)
	TotalPerRead = 	100*((int(AsIs)+int(Reverse))/int(TotCount))			
	print >>OUT, "Total Correct (5' ends Read 1 and 3' ends Read2): %s" % AsIs
	print >>OUT, "Total Swapped (5' ends Read 2 and 3' ends Read1): %s" % Reverse	
	print >>OUT, "Percent of Reads Retained: %s" % TotalPerRead	
	print >>OUT, "Total Correct 5' Seed Matches to R1 (Correct): %s" % FivePrimeMatch
	print >>OUT, "Total Correct 3' Seed Matches to R1 (Inverted): %s" % ThreePrimeMatchR1
	print >>OUT, "Total Correct 3' Seed Matches to R2 (Correct): %s" % ThreePrimeMatchR2 
	print >>OUT, "Total Correct 5' Seed Matches to R2 (Inverted): %s" % FivePrimeMatchR2
	print >>OUT, "Total R2 Align Upstream SfiI (Inverted): %s" % LotsMutsInverted
	print >>OUT, "Total Align KpnI, XbaI, and backbone (Correct): %s" % DoubleSeedMuts
	print >>OUT, "Total R1 Align 5' SfiI (Correct): %s" % LotsMuts
	print >>OUT, "Total Short Seed Correct Orientation: %s" % ShortSeed
	print >>OUT, "Total Short Seed Inverted Orientation: %s" % ShortSeedInv

# Performs a couple of quality checks on tags and their sequences
# Originally written in perl as tag-check.pl and later translated to python by SU
def tag_check(Tags, CREs, TagLength):
	tag = Tags.split('.')
	tag[-2] = tag[-2] + '_qual'
	Tagsnew = '.'.join(tag)
	cre = CREs.split('.')
	cre[-2] = cre[-2] + '_qual'
	CREsnew = '.'.join(cre)
	good = {}
	tagcnt = 1
	name = ""
	seq = ""
	m = re.compile('[B-I]{%d}' % TagLength)
	with open(Tags) as ti:
		to = open(Tagsnew,"w")
		for line in ti:
			LINE = line.rstrip('\r\n')
			if tagcnt == 1 and not LINE.startswith("@"):
				raise Exception('Error: Row count error.')
			elif tagcnt == 1:
				name = LINE
			elif tagcnt == 2:
				seq = LINE
			elif tagcnt == 3 and LINE is not '+':
				raise Exception('Error: Row count error.')
			elif tagcnt == 4:
				if m.match(LINE):
					print >>to, name + "\n" + seq + "\n" + "+\n" + LINE
					good[name] = 1
				tagcnt = 0
			tagcnt = tagcnt +1
		to.close()
	tagcnt = 1
	with open(CREs) as ci:
		goodcre = 0
		co = open(CREsnew,"w")
		for line in ci:
			LINE = line.rstrip('\r\n')
			if tagcnt == 1 and not LINE.startswith("@"):
				raise Exception('Error: Row count error.')
			elif tagcnt == 1 and LINE in good:
				print >>co, LINE
				goodcre = 1
			elif tagcnt == 2 and goodcre == 1:
				print >>co, LINE + '\n+'
			elif tagcnt == 4:
				if goodcre == 1:
					print >>co, LINE
				tagcnt = 0
				goodcre = 0
			tagcnt = tagcnt +1
		co.close()

#Define the function that matches fragments to tags
##Original Matcher call: python ~/Scripts/Matcher.py Cres.fastq Tags.fastq
def Matcher(Read1, Read2):
	##Temp line variables
	R1_Name = ""
	R1_Seq = ""
	R1_Strand = ""
	R1_Qual = ""
	R2_Name = ""
	R2_Seq = ""
	R2_Strand = ""
	R2_Qual = ""
	#Open output files
	FivePrime = open("Inert_Compiled_Cre_tags.txt", "w")
	Cres = open("Inert_Final_Cres.fastq", "w")
	Tags = open("Inert_Final_Tags.fastq", "w")
	OUT = open("Inert_MiscCut_Sequences.txt", "w")
	Assign = {1: "R1_Name = LINE1; R2_Name = LINE2", 2: "R1_Seq = LINE1; R2_Seq = LINE2", 3: "R1_Strand = LINE1; R2_Strand = LINE2", 4: "R1_Qual = LINE1; R2_Qual = LINE2"}
	LineCounter=1
	#Loop through files simultaneously here
	#Don't need to izip longest because mate pair files are same length
	UnCut=0
	with open(Read1) as file1, open(Read2) as file2:
		for line1, line2 in izip(file1, file2):  
			LINE1 = line1.rstrip('\r\n')
			LINE2 = line2.rstrip('\r\n')
			exec(Assign[LineCounter])
			if(LineCounter < 4): 
				LineCounter += 1
			else:
				LineCounter = 1
				#Split names to exclude read number
				L1 = re.split('\s+', R1_Name)
				L2 = re.split('\s+', R2_Name)
				#Match read ID's
				if (L1[0] == L2[0]):				
					#If no cuts were able to be made then register in Uncut variable
					if (R1_Seq == R2_Seq):
						UnCut += 1
					#Else print out cres (fastq) and tags (fastq) and cre-tag (txt) pairing files
					else:
						print >>FivePrime, L1[0]+"\t"+R1_Seq+"\t"+R2_Seq
						print >>Cres, ''.join(R1_Name)+"\n",''.join(R1_Seq)+"\n",''.join(R1_Strand)+"\n",''.join(R1_Qual)
						print >>Tags, ''.join(R2_Name)+"\n",''.join(R2_Seq)+"\n",''.join(R2_Strand)+"\n",''.join(R2_Qual)
				else:
					print "ERROR"
	print >>OUT, "Uncut sequences: %s" % UnCut

# Pads FastQ file with adapter sequences to ensure correct mapping of similar fragments
# Originally written in perl as pad-fastq.pl and then translated to python by SU
def pad_fastq(infile, outfile):
	with open(infile) as f:
		padcnt = 1
		o = open(outfile,"w")
		for line in f:
			LINE = line.rstrip('\r\n')
			if padcnt == 2:
				LINE = 'ACTGGCCGCTTGACG' + LINE + 'CACTGCGGCTCCTGC'
			elif padcnt == 4:
				LINE = 'IIIIIIIIIIIIIII' + LINE + 'IIIIIIIIIIIIIII'
				padcnt = 0
			padcnt = padcnt +1
			print >>o, LINE
		o.close()

#Replaces Contig based alignment coordinates with MPRA fragment names
#Specify the genomic parsing function to score the alignments
##OG Call: python FINAL_Genomic_Parser.py <INUPUT.bam> ~/Scripts/LibTranslation.txt 
####Changed: ~/Scripts/LibTranslation.txt to ~/Genomes/MPRA/Lib_Align_Files/Condensed_HC_Lib/Masked_LibTranslation.txt
####Changed: print >>Parsed, name,"\t",UniqID -> print >>Parsed, name,"\t",KEY[UniqID]
#### Skips importing genomic align MPRA translator script
####Changed: Cigar_Tag_Anote.py rolled into this def too. Will print out tags in final file
####Compiled_Cre_tags = Inert_Compiled_Cre_tags.txt
def GenomicParse(bamFile, Translator, Compiled_Cre_tags, TagLength):
	KEY = {}
	F1= open(Translator)
	#Load in alignment coordinate to CRE translation file
	for line in F1:
		LINE = line.rstrip('\r\n')
		lst = LINE.split('\t')
		KEY[lst[1]] = lst[0]
	F1.close()
	#Load in CRE - tag translation file
	F2= open(Compiled_Cre_tags)
	ID_Tag = {}
	for line in F2:
		LINE = line.rstrip('\r\n')
		Lst = LINE.split('\t')
		Id = Lst[0].strip(' ')
		ID = Lst[0].strip('@')
		Tag = Lst[2].strip(' ')
		if len(Tag) == TagLength:
			ID_Tag[ID] = Tag
		else:
			pass
	F2.close()
	#Open output files
	Parsed = open("CigarParsedMatched.txt","w")
	print >>Parsed,"Read_Name","\t","Tag","\t","Reference","\t","Cigar","\t","NM","\t","Match+Mismatch","\t","Matches","\t","MisMatches","\t","Insertion","\t","Deletion","\t","Skip","\t","SoftClip","\t","HardClip","\t","Padding"
	NoTAG = open("CigarParsedNoTag.txt","w")
	print >>NoTAG,"Read_Name","\t","Reference","\t","Cigar","\t","NM","\t","Match+Mismatch","\t","Matches","\t","MisMatches","\t","Insertion","\t","Deletion","\t","Skip","\t","SoftClip","\t","HardClip","\t","Padding"
	Failed = open("CigarParsedFailed.txt","w")
	print >>Failed,"Read_Name","\t","Reference","\t","Cigar","\t","NM","\t","Match+Mismatch","\t","Matches","\t","MisMatches","\t","Insertion","\t","Deletion","\t","Skip","\t","SoftClip","\t","HardClip","\t","Padding"
	#Load in aligned reads
	bamFP = pysam.Samfile(bamFile,"rb");
	for read in bamFP:
		if( not( read.is_unmapped ) ): #if it's mapped
			cigarLine=read.cigar;
			Tags = read.tags
			cigarTwo=read.cigarstring;
			name=read.qname
			FAName=bamFP.getrname(read.tid)
			Legth_on_Ref = read.alen
			Start = read.reference_start
			End = read.reference_end
			Flag = read.flag
			UniqID = '_'.join([FAName,str(Start),str(End),str(Legth_on_Ref)])
			Seq=read.query
			MM=0
			Insertion = 0
			Deletion = 0
			Skip = 0 
			SoftClip = 0
			HardClip = 0
			Padding = 0
			NM = 0
			#Get edit distance tag
			for (TagType,TagLength) in Tags:
				if(  TagType == "NM"): #NM number
					NM = TagLength     
				else:
					pass
			#Pull alignment info
			for (cigarType,cigarLength) in cigarLine:
				if(  cigarType == 0): #match/mismatch
					MM = MM+cigarLength               
				elif(cigarType == 1): #insertions
					Insertion = Insertion+cigarLength
				elif(cigarType == 2): #deletion
					Deletion = Deletion+cigarLength
				elif(cigarType == 3): #skip
					Skip = Skip+cigarLength
				elif(cigarType == 4): #soft clipping
					SoftClip = SoftClip+cigarLength
				elif(cigarType == 5): #hard clipping
					HardClip = HardClip+cigarLength
				elif(cigarType == 6): #padding
					Padding = Padding+cigarLength
				else:
					print "Wrong CIGAR number";
					sys.exit(1);
			#Calculate the amount of mismatches and matches
			Mismatches = NM-(Insertion+Deletion)
			Matches = MM-Mismatches
			#Print out the properly and improperly aligned reads to output
			if UniqID in KEY:
				if name in ID_Tag:
					print >>Parsed, name,"\t",ID_Tag[name],"\t",KEY[UniqID],"\t",cigarTwo,"\t",NM,"\t",MM,"\t",Matches,"\t",Mismatches,"\t",Insertion,"\t",Deletion,"\t",Skip,"\t",SoftClip,"\t",HardClip,"\t",Padding
				else:
					print >>NoTAG, name,"\t",KEY[UniqID],"\t",cigarTwo,"\t",NM,"\t",MM,"\t",Matches,"\t",Mismatches,"\t",Insertion,"\t",Deletion,"\t",Skip,"\t",SoftClip,"\t",HardClip,"\t",Padding
			else:
				print >>Failed, name,"\t",UniqID,"\t",cigarTwo,"\t",NM,"\t",MM,"\t",Matches,"\t",Mismatches,"\t",Insertion,"\t",Deletion,"\t",Skip,"\t",SoftClip,"\t",HardClip,"\t",Padding

#Filters for unique tag-Cre matches only
##RepeatFilter.py script rolled in to print Unique and repeat tags files at once
def MultipleTagged( File ):
	##Original script call: python ~/Scripts/MultipleTagged.py Translated_Parsed_Cigar_withTAG_TrimedPI_gt88_temp.txt
	F1 = open(File)
	TotalFile = open("Inert_RepeatedTags.txt", "w")
	Tags = {}
	RepTags = {}
	for line in F1:
		LINE = line.rstrip('\r\n')
		Lst = LINE.split('\t')
		if Lst[1].strip(' ') in Tags:
			if Lst[2].strip(' ') == Tags[Lst[1].strip(' ')]:
				pass
			else:
				RepTags[Lst[1].strip(' ')] = Lst[1].strip(' ')
		else:
			Tags[Lst[1].strip(' ')] = Lst[2].strip(' ')
	F1.close()
	#Print out tags for repeated tag file
	for key, value in RepTags.iteritems():
		print >>TotalFile, key
	##Save the full file entry of repeated tags for optional analysis
	F2 = open(File)
	TotalFil = open("Inert_SamEntries_RepeatedTags.txt", "w") 
	FinalFile = open("UniqTags_Translated_Parsed_Cigar_withTAG_TrimedPI_gt88.txt", "w")
	for line in F2:    
		LINE = line.rstrip('\r\n')
		Lst = LINE.split('\t')
		if Lst[1] in RepTags:
			print >>TotalFil, LINE
		else:
			print >>FinalFile, LINE
	F2.close()

#Tabulates tag counts		
def Tabulator( File1 ):
	TagsCount = {}
	TagsLibFrag = {}
	
	F1= open(File1)
	for line in F1:    
		LINE = line.rstrip('\r\n')
		Lst = LINE.split('\t')
		if Lst[1] in TagsCount:
			TagsCount[Lst[1]] += 1
		else:
			TagsCount[Lst[1]]  = 1
			TagsLibFrag[Lst[1]] = Lst[2]
	F1.close()
	TagFile = open("Inert_Tag_Counts.txt", "w")
	LibCounts = {}
	for key, value in TagsCount.iteritems():
		Entry = '\t'.join([key,str(TagsCount[key]),TagsLibFrag[key]])
		print >>TagFile, Entry
		if TagsLibFrag[key] in LibCounts:
			LibCounts[TagsLibFrag[key]] += 1
		else:
			LibCounts[TagsLibFrag[key]] = 1
	LibFile = open("Inert_Lib_Counts.txt", "w")
	for key, value in LibCounts.iteritems():
		Entry = '\t'.join([key,str(LibCounts[key])])
		print >>LibFile, Entry

def ProcessR_sorted(infile):
	##Original script call: python ~/Scripts/ProcessR_sorted.py CigarParsedMatched_sort.txt
	tag = ''
	cre = ''
	ccount = int(0)
	o = open('R_Processed.tsv', "w")
	print >>o, "Tag\tAlignment\tmaxID\tn"
	with open(infile, "r") as f:
	        next(f)
	        for line in f:
	                line = line.rstrip()
	                T = re.split('\t| \t',line)
	                if T[1] == 'Tag':
	                        continue
	                else:
	                        if T[1] == tag and T[2] == cre:
	                                ccount +=1
	                                T[5:13] = map(int,T[5:13])
	                                ID = (T[5]-sum(T[7:13]))/T[5]
	                                if maxID < ID:
	                                        maxID = ID
	                        else:
	                                if ccount >0:
	                                        entry = [tag, cre, str(maxID), str(ccount)]
	                                        print >>o, '\t'.join(entry)
	                                tag = T[1]
	                                cre = T[2]
	                                ccount = 1
	                                T[5:13] = map(int,T[5:13])
	                                ID = (T[5]-sum(T[7:13]))/T[5]
	                                maxID = ID
	entry = [tag, cre, str(maxID), str(ccount)]
	print >>o, '\t'.join(entry)
	o.close()

def ProcessR_sorted_v2(infile, usefile, garbagefile):
	##Original script call: python ~/Scripts/ProcessR_sorted_v2.py R_Processed.tsv R_Processed_use.tsv R_Processed_trash.tsv
	tag = ''
	cre = ''
	garbage = []
	with open(infile, "r") as f:
	        for line in f:
	                line = line.rstrip()
	                T = re.split('\t| \t',line)
	                if T[0] == 'Tag':
	                        continue
	                else:
	                        U = T[1].split('_')
	                        if T[0] == tag and U[0] != cre:
	                                garbage.append(tag)
	                        elif T[0] != tag:
	                                tag = T[0]
	                                cre = U[0]
	o = open(usefile, "w")
	wo = open(garbagefile, 'w')
	with open(infile, "r") as f:
	        for line in f:
	                line = line.rstrip()
	                T = re.split('\t| \t',line)
	                if T[0] == 'Tag':
	                        print >>o, line
	                        print >>wo, line
	                else:
	                        if T[0] in garbage:
	                                print >>wo, line
	                        else:
	                                print >>o, line
	o.close()
	wo.close()

#Competent Library-seq read strander function
##Original script call: python ~/Scripts/TAG_Strander_V2016_10_31.py R1_Trimmed.fastq R2_Trimmed.fastq
def Comp_TAG_Strander( Read1, Read2 ):
	##Counter variable to track conversation stats:
	AsIs = 0
	Reverse = 0
	Complement = 0
	RevComp = 0
	##Temp line variables
	R1_Name = ""
	R1_Seq = ""
	R1_Strand = ""
	R1_Qual = ""
	R2_Name = ""
	R2_Seq = ""
	R2_Strand = ""
	R2_Qual = ""
	#Open output files
	FivePrime = open("Comp_ReadOne.fastq", "w")
	ThreePrime = open("Comp_ReadTwo.fastq", "w")
	Garbage = open("Comp_Garbage_IDs.txt", "w")
	OUT = open("Comp_TAG_StranderOutput.txt", "w")
	#Complement dictionary
	Trans = {'A': 'T', 'T': 'A', 'C': 'G', 'G' : 'C', 'N' : 'N'}
	Assign = {1: "R1_Name = LINE1; R2_Name = LINE2", 2: "R1_Seq = LINE1; R2_Seq = LINE2", 3: "R1_Strand = LINE1; R2_Strand = LINE2", 4: "R1_Qual = LINE1; R2_Qual = LINE2"}
	FivePrimeMatch = 0
	ThreePrimeMatchR2 = 0
	ThreePrimeMatchR1 = 0
	FivePrimeMatchR2 = 0
	DoubleSeedMuts = 0
	LotsMuts = 0
	LotsMutsInverted = 0 
	ShortSeed = 0
	ShortSeedInv = 0
	TotCount = 0
	LineCounter=1
	#Loop through files simultaneously here
	#Don't need to izip longest because mate pair files are same length
	with open(Read1) as file1, open(Read2) as file2:
		for line1, line2 in izip(file1, file2):  
			line1 = line1.rstrip('\r\n')
			line2 = line2.rstrip('\r\n')
			Line1 = line1.rstrip()
			Line2 = line2.rstrip()
			#Split each line to determine value
			LINE1 = list(Line1)
			LINE2 = list(Line2)
			#Assign the temp FastQ values based on the commands stored in the array
			#This allows you to skip using
			exec(Assign[LineCounter])
			if(LineCounter < 4): #Name @	
				LineCounter += 1
			else:
				#print ''.join(R1_Seq[0:35])
				LineCounter = 1
				TotCount += 1
				#Match 5' seed to read one: confirms correctly oriented fragments
				if (bool(re.search( 'CAAGAAGGGCGGCAAGAT', ''.join(R1_Seq[0:30]))) is True):
					AsIs += 1
					FivePrimeMatch += 1     		
					print >>FivePrime, ''.join(R1_Name)+"\n",''.join(R1_Seq)+"\n",''.join(R1_Strand)+"\n",''.join(R1_Qual)
					print >>ThreePrime, ''.join(R2_Name)+"\n",''.join(R2_Seq)+"\n",''.join(R2_Strand)+"\n",''.join(R2_Qual)
				##Match 3' seed to read one: inverted oriented fragments
				elif (bool(re.search( 'TGTCTGCTCGAAGCGGCC', ''.join(R1_Seq[0:30]))) is True):
					Reverse += 1
					ThreePrimeMatchR1 += 1
					print >>ThreePrime, ''.join(R2_Name)+"\n",''.join(R1_Seq)+"\n",''.join(R1_Strand)+"\n",''.join(R1_Qual)
					print >>FivePrime, ''.join(R1_Name)+"\n",''.join(R2_Seq)+"\n",''.join(R2_Strand)+"\n",''.join(R2_Qual)
				#Match 3' seed to read two: confirms correctly oriented fragments with mutation in 5' seed	
				elif (bool(re.search( 'CATGTCTGCTCGAAGCGG', ''.join(R2_Seq[0:30]))) is True):
					AsIs += 1
					ThreePrimeMatchR2 += 1
					print >>FivePrime, ''.join(R1_Name)+"\n",''.join(R1_Seq)+"\n",''.join(R1_Strand)+"\n",''.join(R1_Qual)
					print >>ThreePrime, ''.join(R2_Name)+"\n",''.join(R2_Seq)+"\n",''.join(R2_Strand)+"\n",''.join(R2_Qual)
				#Match backbone plus XbaI seed to read two: confirms correctly oriented fragments with mutations in 5' AND 3' seeds	
				elif (bool(re.search( 'TGTAATAATTCTAGA', ''.join(R2_Seq[0:50]))) is True):
					Reverse += 1
					FivePrimeMatchR2 += 1
					print >>ThreePrime, ''.join(R2_Name)+"\n",''.join(R1_Seq)+"\n",''.join(R1_Strand)+"\n",''.join(R1_Qual)
					print >>FivePrime, ''.join(R1_Name)+"\n",''.join(R2_Seq)+"\n",''.join(R2_Strand)+"\n",''.join(R2_Qual)
				elif (bool(re.search( 'TGTAATAATTCTAGA', ''.join(R1_Seq[0:50]))) is True):
					AsIs += 1
					DoubleSeedMuts += 1
					print >>FivePrime, ''.join(R1_Name)+"\n",''.join(R1_Seq)+"\n",''.join(R1_Strand)+"\n",''.join(R1_Qual)
					print >>ThreePrime, ''.join(R2_Name)+"\n",''.join(R2_Seq)+"\n",''.join(R2_Strand)+"\n",''.join(R2_Qual)
				#Match 5' Sfi site: confirms correctly oriented fragments with lots of mutations 
				elif (bool(re.search( 'AAGCTAGTCGGGGC', ''.join(R1_Seq[70:110]))) is True):
					AsIs += 1
					LotsMuts += 1
					print >>FivePrime,''.join(R1_Name)+"\n",''.join(R1_Seq)+"\n",''.join(R1_Strand)+"\n" ,''.join(R1_Qual)
					print >>ThreePrime,''.join(R2_Name)+"\n",''.join(R2_Seq)+"\n" ,''.join(R2_Strand)+"\n" ,''.join(R2_Qual)
				else:
					print >>Garbage, ''.join(R1_Name)+"\t"+''.join(R1_Qual)+"\t"+''.join(R2_Qual)
	TotalPerRead = 	100*((int(AsIs)+int(Reverse))/int(TotCount))		
	print >>OUT, "Total Correct (5' ends Read 1 and 3' ends Read2): %s" % AsIs
	print  >>OUT, "Total Swapped (5' ends Read 2 and 3' ends Read1): %s" % Reverse	
	print  >>OUT, "Percent of Reads Retained: %s" % TotalPerRead	
	print  >>OUT, "Total Correct 5' Seed Matches to R1 (Correct): %s" % FivePrimeMatch
	print  >>OUT, "Total Correct 3' Seed Matches to R1 (Inverted): %s" % ThreePrimeMatchR1	
	print  >>OUT, "Total Correct 3' Seed Matches to R2 (Correct): %s" % ThreePrimeMatchR2 
	print  >>OUT, "Total Correct 5' Seed Matches to R2 (Inverted): %s" % FivePrimeMatchR2
	print  >>OUT, "Total Align KpnI, XbaI, and backbone (Correct): %s" % DoubleSeedMuts
	print  >>OUT, "Total R1 Align 5' SfiI (Correct): %s" % LotsMuts

#Tag tabulator script: used for both competent library and experiment sequencing call on "Tags.fastq"
##Original script call: time python Tag_Tabulator.py Tags.fastq
def TagTabulator(File, Sample, TagLength):
	F1= open(File)
	name = '_'.join([Sample,"Tag_Counts.txt"])
	TotalFile = open(name, "w")
	Tags = {}
	counter = 0
	for line in F1:
			LINE = line.rstrip('\r\n')
			Lst = LINE.split('\t')
			ele = list(Lst[0])
			counter += 1
			if counter == 2:
				if len(ele) == TagLength: # Barcode length needs to be encoded variably!
					if Lst[0] in Tags:
						Tags[Lst[0]] += 1
					else:
						Tags[Lst[0]] = 1
				else:
					pass
			elif counter == 4:
				counter = 0
			else:
				pass
	for key, value in Tags.iteritems():
		Entry = '\t'.join([key,str(value)])
		print >>TotalFile, Entry

######Filter competent tags for those only in inert library
####Original script call: python ~/Scripts/Common_Tags.py Competant_Counts.txt Inert_Counts.txt
def TagFilter(Inert,Competant):
	F1 = open(Competant)
	F2 = open(Inert)
	CompTags = {}
	for line in F1:
		LINE = line.rstrip('\r\n')
		lst = LINE.split('\t')
		CompTags[lst[0]] = lst[1]
	F1.close();
	#Open Output Files
	OUT1 = open("Inert_TagCounts_In_Comp.txt", "w")
	OUT2 = open("Comp_TagCounts_In_Inert.txt", "w")
	for line in F2:
		LINE = line.rstrip('\r\n')
		lst = LINE.split('\t')
		if lst[0] in CompTags:
			print >>OUT1, LINE
			print >>OUT2, lst[0],"\t",CompTags[lst[0]]	
		else:
			pass

#Experimental library-seq read strander function
#Library strand orienting function for experimental MPRA reads (need separate function for different primers)
##Original script call: time python ~/Scripts/Tag_Strander_2016_11_21.py R1_pDNA_Trimmed.fastq R2_pDNA_Trimmed.fastq
def Exp_TAG_Strander( Read1, Read2, Sample ):
	##Counter variable to track conversation stats:
	AsIs = 0
	Reverse = 0
	Complement = 0
	RevComp = 0
	##Temp line variables
	R1_Name = ""
	R1_Seq = ""
	R1_Strand = ""
	R1_Qual = ""
	R2_Name = ""
	R2_Seq = ""
	R2_Strand = ""
	R2_Qual = ""
	#Open output files
	FivePrime = open("_".join([Sample,"_ReadOne.fastq"]), "w")
	ThreePrime = open("_".join([Sample,"_ReadTwo.fastq"]), "w")
	Garbage = open("_".join([Sample,"_Garbage_IDs.txt"]), "w")
	OUT = open("_".join([Sample,"_StranderOutput.txt"]), "w")
	#Complement dictionary
	Trans = {'A': 'T', 'T': 'A', 'C': 'G', 'G' : 'C', 'N' : 'N'}
	Assign = {1: "R1_Name = LINE1; R2_Name = LINE2", 2: "R1_Seq = LINE1; R2_Seq = LINE2", 3: "R1_Strand = LINE1; R2_Strand = LINE2", 4: "R1_Qual = LINE1; R2_Qual = LINE2"}
	FivePrimeMatch = 0
	ThreePrimeMatchR2 = 0
	ThreePrimeMatchR1 = 0
	FivePrimeMatchR2 = 0
	MidMatchCor = 0
	MidMatchRev = 0
	LastDitchRev = 0
	LastDitchCor = 0
	TopThirdMatchCor = 0
	TopThirdMatchRev = 0
	DoubleSeedMuts = 0
	LotsMuts = 0
	LotsMutsInverted = 0
	ShortSeed = 0
	ShortSeedInv = 0
	TotCount = 0
	LineCounter=1
	#Loop through files simultaneously here
	#Don't need to izip longest because mate pair files are same length
	with open(Read1) as file1, open(Read2) as file2:
		for line1, line2 in izip(file1, file2):
			line1 = line1.rstrip('\r\n')
			line2 = line2.rstrip('\r\n')
			Line1 = line1.rstrip()
			Line2 = line2.rstrip()
			#Split each line to determine value
			LINE1 = list(Line1)
			LINE2 = list(Line2)
			#Assign the temp FastQ values based on the commands stored in the array
			#This allows you to skip using
			exec(Assign[LineCounter])
			if(LineCounter < 4): #Name @
					LineCounter += 1
			else:
				#print ''.join(R1_Seq[0:35])
				LineCounter = 1
				TotCount += 1
				#Match 5' seed to read one: confirms correctly oriented fragments
				if (bool(re.search( 'CAAGAAGGGCGGCAAGAT', ''.join(R1_Seq[0:30]))) is True):
					AsIs += 1
					FivePrimeMatch += 1
					print >>FivePrime, ''.join(R1_Name)+"\n",''.join(R1_Seq)+"\n",''.join(R1_Strand)+"\n",''.join(R1_Qual)
					print >>ThreePrime, ''.join(R2_Name)+"\n",''.join(R2_Seq)+"\n",''.join(R2_Strand)+"\n",''.join(R2_Qual)
				##Match 3' seed to read one: inverted oriented fragments
				elif (bool(re.search( 'ACGCTCTTCCGATCT', ''.join(R1_Seq[0:30]))) is True):
					Reverse += 1
					ThreePrimeMatchR1 += 1
					print >>ThreePrime, ''.join(R2_Name)+"\n",''.join(R1_Seq)+"\n",''.join(R1_Strand)+"\n",''.join(R1_Qual)
					print >>FivePrime, ''.join(R1_Name)+"\n",''.join(R2_Seq)+"\n",''.join(R2_Strand)+"\n",''.join(R2_Qual)
				#Match 3' seed to read two: confirms correctly oriented fragments with mutation in 5' seed      
				elif (bool(re.search( 'ACGCTCTTCCGATCT', ''.join(R2_Seq[0:30]))) is True):
					AsIs += 1
					ThreePrimeMatchR2 += 1
					print >>FivePrime, ''.join(R1_Name)+"\n",''.join(R1_Seq)+"\n",''.join(R1_Strand)+"\n",''.join(R1_Qual)
					print >>ThreePrime, ''.join(R2_Name)+"\n",''.join(R2_Seq)+"\n",''.join(R2_Strand)+"\n",''.join(R2_Qual)
				#Match backbone plus XbaI seed to read two: confirms correctly oriented fragments with mutations in 5' AND 3' seeds 
				elif (bool(re.search( 'TCTAGAATTATTACA', ''.join(R2_Seq[25:65]))) is True):
					Reverse += 1
					FivePrimeMatchR2 += 1
					print >>ThreePrime, ''.join(R2_Name)+"\n",''.join(R1_Seq)+"\n",''.join(R1_Strand)+"\n",''.join(R1_Qual)
					print >>FivePrime, ''.join(R1_Name)+"\n",''.join(R2_Seq)+"\n",''.join(R2_Strand)+"\n",''.join(R2_Qual)
				elif (bool(re.search( 'TGTAATAATTCTAGA', ''.join(R1_Seq[90:110]))) is True):
					AsIs += 1
					DoubleSeedMuts += 1
					print >>FivePrime, ''.join(R1_Name)+"\n",''.join(R1_Seq)+"\n",''.join(R1_Strand)+"\n",''.join(R1_Qual)
					print >>ThreePrime, ''.join(R2_Name)+"\n",''.join(R2_Seq)+"\n",''.join(R2_Strand)+"\n",''.join(R2_Qual)
				#Match 5' Sfi site: confirms correctly oriented fragments with lots of mutations
				elif (bool(re.search( 'GATTCTCATTAAGGCCA', ''.join(R1_Seq[45:90]))) is True):
					AsIs += 1
					LotsMuts += 1
					print >>FivePrime,''.join(R1_Name)+"\n",''.join(R1_Seq)+"\n",''.join(R1_Strand)+"\n" ,''.join(R1_Qual)
					print >>ThreePrime,''.join(R2_Name)+"\n",''.join(R2_Seq)+"\n" ,''.join(R2_Strand)+"\n" ,''.join(R2_Qual)
				#Match Middle  to Read one: confirms correctly oriented fragments
				elif (bool(re.search( 'CCAAGAAGGGCGGCAA', ''.join(R1_Seq[65:95]))) is True):
					AsIs += 1
					MidMatchCor += 1                     
					print >>FivePrime, ''.join(R1_Name)+"\n",''.join(R1_Seq)+"\n",''.join(R1_Strand)+"\n",''.join(R1_Qual)
					print >>ThreePrime, ''.join(R2_Name)+"\n",''.join(R2_Seq)+"\n",''.join(R2_Strand)+"\n",''.join(R2_Qual)
				##Match Middle seed to Read two: inverted oriented fragments
				elif (bool(re.search( 'CCAAGAAGGGCGGCAA', ''.join(R2_Seq[65:95]))) is True):
					Reverse += 1
					MidMatchRev += 1
					print >>ThreePrime, ''.join(R2_Name)+"\n",''.join(R1_Seq)+"\n",''.join(R1_Strand)+"\n",''.join(R1_Qual)
					print >>FivePrime, ''.join(R1_Name)+"\n",''.join(R2_Seq)+"\n",''.join(R2_Strand)+"\n",''.join(R2_Qual)
				#Match Middle Inv to Read one: confirms correctly oriented fragments
				elif (bool(re.search( 'TTGCCGCCCTTCTTGG', ''.join(R2_Seq[50:80]))) is True):
					AsIs += 1
					MidMatchCor += 1                     
					print >>FivePrime, ''.join(R1_Name)+"\n",''.join(R1_Seq)+"\n",''.join(R1_Strand)+"\n",''.join(R1_Qual)
					print >>ThreePrime, ''.join(R2_Name)+"\n",''.join(R2_Seq)+"\n",''.join(R2_Strand)+"\n",''.join(R2_Qual)
				##Match 3' Inv seed to Read two: inverted oriented fragments
				elif (bool(re.search( 'TTGCCGCCCTTCTTGG', ''.join(R1_Seq[50:80]))) is True):
					Reverse += 1
					MidMatchRev += 1
					print >>ThreePrime, ''.join(R2_Name)+"\n",''.join(R1_Seq)+"\n",''.join(R1_Strand)+"\n",''.join(R1_Qual)
					print >>FivePrime, ''.join(R1_Name)+"\n",''.join(R2_Seq)+"\n",''.join(R2_Strand)+"\n",''.join(R2_Qual)
				#Match TopThird  to Read one: confirms correctly oriented fragments
				elif (bool(re.search( 'ACTGACCGGCAAGTTG', ''.join(R1_Seq[15:45]))) is True):
					AsIs += 1
					MidMatchCor += 1                     
					print >>FivePrime, ''.join(R1_Name)+"\n",''.join(R1_Seq)+"\n",''.join(R1_Strand)+"\n",''.join(R1_Qual)
					print >>ThreePrime, ''.join(R2_Name)+"\n",''.join(R2_Seq)+"\n",''.join(R2_Strand)+"\n",''.join(R2_Qual)
				##Match TopThird seed to Read two: inverted oriented fragments
				elif (bool(re.search( 'ACTGACCGGCAAGTTG', ''.join(R2_Seq[15:45]))) is True):
					Reverse += 1
					MidMatchRev += 1
					print >>ThreePrime, ''.join(R2_Name)+"\n",''.join(R1_Seq)+"\n",''.join(R1_Strand)+"\n",''.join(R1_Qual)
					print >>FivePrime, ''.join(R1_Name)+"\n",''.join(R2_Seq)+"\n",''.join(R2_Strand)+"\n",''.join(R2_Qual)
				##MatchTopThird Inv seed to Read two: inverted oriented fragments
				elif (bool(re.search( 'CAACTTGCCGGTCAGT', ''.join(R1_Seq[80:110]))) is True):
					Reverse += 1
					TopThirdMatchRev += 1
					print >>ThreePrime, ''.join(R2_Name)+"\n",''.join(R1_Seq)+"\n",''.join(R1_Strand)+"\n",''.join(R1_Qual)
					print >>FivePrime, ''.join(R1_Name)+"\n",''.join(R2_Seq)+"\n",''.join(R2_Strand)+"\n",''.join(R2_Qual)
				#MatchTopThird Inv to Read one: confirms correctly oriented fragments
				elif (bool(re.search( 'CAACTTGCCGGTCAGT', ''.join(R2_Seq[80:110]))) is True):
					AsIs += 1
					TopThirdMatchCor += 1                     
					print >>FivePrime, ''.join(R1_Name)+"\n",''.join(R1_Seq)+"\n",''.join(R1_Strand)+"\n",''.join(R1_Qual)
					print >>ThreePrime, ''.join(R2_Name)+"\n",''.join(R2_Seq)+"\n",''.join(R2_Strand)+"\n",''.join(R2_Qual)
				#Last Ditch to Read one: confirms correctly oriented fragments
				elif (bool(re.search( 'TCCGCGAGA', ''.join(R1_Seq[45:65]))) is True):
					AsIs += 1
					LastDitchCor += 1                     
					print >>FivePrime, ''.join(R1_Name)+"\n",''.join(R1_Seq)+"\n",''.join(R1_Strand)+"\n",''.join(R1_Qual)
					print >>ThreePrime, ''.join(R2_Name)+"\n",''.join(R2_Seq)+"\n",''.join(R2_Strand)+"\n",''.join(R2_Qual)
				##Last Ditch to Read two: inverted oriented fragments
				elif (bool(re.search( 'TCCGCGAGA', ''.join(R2_Seq[45:65]))) is True):
					Reverse += 1
					LastDitchRev += 1
					print >>ThreePrime, ''.join(R2_Name)+"\n",''.join(R1_Seq)+"\n",''.join(R1_Strand)+"\n",''.join(R1_Qual)
					print >>FivePrime, ''.join(R1_Name)+"\n",''.join(R2_Seq)+"\n",''.join(R2_Strand)+"\n",''.join(R2_Qual)
				##MatchTopThird Inv seed to Read two: inverted oriented fragments
				elif (bool(re.search( 'TCTCGCGGA', ''.join(R1_Seq[65:80]))) is True):
					Reverse += 1
					LastDitchRev += 1
					print >>ThreePrime, ''.join(R2_Name)+"\n",''.join(R1_Seq)+"\n",''.join(R1_Strand)+"\n",''.join(R1_Qual)
					print >>FivePrime, ''.join(R1_Name)+"\n",''.join(R2_Seq)+"\n",''.join(R2_Strand)+"\n",''.join(R2_Qual)
				#MatchTopThird Inv to Read one: confirms correctly oriented fragments
				elif (bool(re.search( 'TCTCGCGGA', ''.join(R2_Seq[65:80]))) is True):
					AsIs += 1
					LastDitchCor += 1                     
					print >>FivePrime, ''.join(R1_Name)+"\n",''.join(R1_Seq)+"\n",''.join(R1_Strand)+"\n",''.join(R1_Qual)
					print >>ThreePrime, ''.join(R2_Name)+"\n",''.join(R2_Seq)+"\n",''.join(R2_Strand)+"\n",''.join(R2_Qual)    
				else:
					print >>Garbage, ''.join(R1_Name)+"\t"+''.join(R1_Qual)+"\t"+''.join(R2_Qual)
	TotalPerRead = 	100*((int(AsIs)+int(Reverse))/int(TotCount))
	print >>OUT, "Total Correct (5' ends Read 1 and 3' ends Read2): %s" % AsIs
	print >>OUT, "Total Swapped (5' ends Read 2 and 3' ends Read1): %s" % Reverse
	print >>OUT, "Percent of Reads Retained: %s " % TotalPerRead
	print >>OUT, "Total Correct 5' Seed Matches to R1 (Correct): %s" % FivePrimeMatch
	print >>OUT, "Total Correct 3' Seed Matches to R1 (Inverted): %s" % ThreePrimeMatchR1
	print >>OUT, "Total Correct 3' Seed Matches to R2 (Correct): %s" % ThreePrimeMatchR2
	print >>OUT, "Total Correct 5' Seed Matches to R2 (Inverted): %s" % FivePrimeMatchR2
	print >>OUT, "Total Correct Middle Matches (Correct): %s" % MidMatchCor
	print >>OUT, "Total Correct Middle Matches (Inverted): %s" % MidMatchRev
	print >>OUT, "Total Correct Last Ditch 9bp Match (Correct): %s" % LastDitchCor
	print >>OUT, "Total Correct Last Ditch 9bp Match (Inverted): %s" % LastDitchRev
	print >>OUT, "Total Correct Top Third Matches (Correct): %s" % TopThirdMatchCor
	print >>OUT, "Total Correct Top Third (Inverted): %s" % TopThirdMatchRev
	print >>OUT, "Total Align KpnI, XbaI, and backbone (Correct): %s" % DoubleSeedMuts
	print >>OUT, "Total R1 Align 5' SfiI (Correct): %s" % LotsMuts

#Compile Count files
#OG Call: #python BasicMerge.py pDNA_Filtered_Rep_1_2_Tags.txt pDNA_Filtered_Rep_2_2_Tags.txt cDNA_Filtered_Rep_1_2_Tags.txt cDNA_Filtered_Rep_2_2_Tags.txt Merged_CountsForR.txt
def Compiler( Inert, Comp, pDNA, cDNA ):
	#Load Inert Counts
	Inertcnts = {}
	Align = {}
	F0 = open(Inert)
	for line in F0:
		LINE = line.rstrip('\r\n')
		lst = LINE.split('\t')
		Inertcnts[lst[0]] = lst[1]
		Align[lst[0]] = lst[2]
	F0.close();
	#Load Competent Counts
	if Comp is None:
		pass
	else:
		Compcnts = {}
		F1 = open(Comp)

		for line in F1:
			LINE = line.rstrip('\r\n')
			lst = LINE.split('\t')
			Compcnts[lst[0]] = lst[1]
		F1.close();
	#Make Header
	if Comp is None:
		header = ["Tag","Alignment","Inert_Counts"]
	else:
		header = ["Tag","Alignment","Inert_Counts", "Comp_Counts"]
	#Load cDNA File Counts
	cDNA_Dicts = [] 
	for name in cDNA:		
		if (bool(re.search( '/', name) is True)):
			DictName = (name.split('/')[-1]).split('.')[0]
			header.append(DictName)
			exec( 'DictName = {} ')			
			exec('cDNA_Dicts.append( DictName )')
		else:
			DictName = (name.split('/')[-1]).split('.')[0]
			header.append(DictName)
			exec( 'DictName = {} ')
			exec('cDNA_Dicts.append( DictName )')
		#Loop Through File and load counts hash
		FP = open(name)
		for line in FP:
			LINE = line.rstrip('\r\n')
			lst = LINE.split('\t')
			#Load To File Specific dict
			exec('DictName[lst[0]] = lst[1]')
		FP.close();
	#Load pDNA Files Counts
	pDNA_Dicts = []
	for name in pDNA:
		if (bool(re.search( '/', name) is True)):	
			DictName = (name.split('/')[-1]).split('.')[0]
			header.append(DictName)
			exec( 'DictName = {} ')
			exec('pDNA_Dicts.append( DictName )')
		else:
			DictName = (name.split('/')[-1]).split('.')[0]
			header.append(DictName)
			exec( 'DictName = {} ')
			exec('pDNA_Dicts.append( DictName )')
		#Loop Through and Load Each Counts Hash
		FP = open(name)
		for line in FP:
			LINE = line.rstrip('\r\n')
			lst = LINE.split('\t')		
			#Load File Specific 
			exec('DictName[lst[0]] = str(lst[1])')
		FP.close();
	###Print out all values by tags based on Inert Tag hash entries
	#Print File Header
	OUT = open("Merged_CountsForR.txt", "w")
	print >>OUT, '\t'.join(header)
	#Get the Inert Tag Counts Data
	for Tag in Inertcnts:
		entry =[Tag]
		entry.append(Align[Tag])
		entry.append(Inertcnts[Tag])
		#Get the Competent Tag Counts if exists
		if Comp is None:
			pass
		else:
			if Tag in Compcnts:
				entry.append(Compcnts[Tag])
			else:
				entry.append(0)
		#Print the cDNA Tag counts
		for cDNA_File in cDNA_Dicts:
			exec('entry.append(cDNA_File[Tag]) if Tag in cDNA_File else entry.append(0)')
			#exec(statement.format('if Tag in cDNA_File:, entry.append(cDNA_File[Tag]), else:, entry.append(0)'))
		#Print the pDNA Tag counts
		for pDNA_File in pDNA_Dicts:
			exec('entry.append(pDNA_File[Tag]) if Tag in pDNA_File else entry.append(0)')
			#exec(statement.format('if Tag in pDNA_File:, entry.append(pDNA_File[Tag]), else:, entry.append(0)')) 
		ENTRY = '\t'.join(map(str, entry))
		print >>OUT, ENTRY

#Annotation Sub Routines
#Anotion Subroutine A
##python AnnoterA.py 		
## Files needed CHIMP_LIB_intial_File Merged_CountsForR.txt
def AnnoterA( resourceDir ):
	#Translates Between Ortholog names
	OrthologTranslator = {}
	#Stores the coordinates of every fragment
	Coordlibrary = {}
	#refs = resource.split('/')
	#refs.append("Master_Ortholog_TranslatorFin")
	Refs = resourceDir+"Master_Ortholog_TranslatorFin.txt"
	F0 = open(Refs)
	for line in F0:
		LINE = line.rstrip('\r\n')
		lst = LINE.split('\t')
		#Load the Ortholog Translation file
		OrthologTranslator[lst[0]] = lst[2]
		OrthologTranslator[lst[2]] = lst[0]
		#Load the Coordinate Library 
		Coordlibrary[lst[0]] = lst[1]
		Coordlibrary[lst[2]] = lst[3]	
	F0.close();
	OUT = open("Annotated_with_coords.txt", "w")
	Files = "Merged_CountsForR.txt"
	F1 = open(Files)
	#Stores every represented fragment
	InSequencedLib = {}
	fileLength = None
	for line in F1:
		LINE = line.rstrip('\r\n')
		lst = LINE.split('\t')
		fileLength = len(lst)
		#Ignore Header
		#####Adjust for number of samples!!
		if lst[0] == "Tag":
			print >>OUT, '\t'.join([lst[0],lst[1], "Ortholog_Seq","Chr","Start","Stop", "Species",'\t'.join(map(str,lst[2:]))])
		else:
			#Load into dict to compare to OrthologTranslator{} later
			InSequencedLib[lst[1]] = lst[1]
			#Test Fragment Species
			if (bool(re.search( 'Chimp', lst[1])) is True):
				species = "PanTro2"
			else:
				species = "hg19"
			#Parse Coordinates in order to pint in bed format
			Coords = re.split('[- :]', Coordlibrary[lst[1]])
			chr = Coords[0]
			Start = Coords[1]
			End = Coords[2].split('[_]')[0]
			Ortholog = OrthologTranslator[lst[1]]
			Entry = '\t'.join([lst[0],lst[1],Ortholog,chr,str(Start),str(End),species,'\t'.join(map(str,lst[2:]))])
			print >>OUT, Entry
	fileLength = fileLength-2
	for entry in OrthologTranslator:
		if entry in InSequencedLib:
			pass
		else:
			#Test Fragment Species
			if (bool(re.search( 'Chimp', entry)) is True):
				species = "PanTro2"
			else:
				species = "hg19"
			#Parse Coordinates in order to pint in bed format
			Coords = re.split('[- :]', Coordlibrary[entry])
			chr = Coords[0]
			Start = Coords[1]
			End = Coords[2].split('[_]')[0]
			Ortholog = OrthologTranslator[entry]
			Zeroed_Entry = [0]
			for x in range(1, fileLength):
				Zeroed_Entry.append(0)
			Entry = '\t'.join(["NA",entry,Ortholog,chr,str(Start),str(End),species,'\t'.join(map(str,Zeroed_Entry))])	
			print >>OUT, Entry	
	F1.close();
	
#Anntotion Sub-Routine B
def AnnoterB( resouceDir ):
	ChimpAllele = {}
	HumanAllele = {}
	F0 = open(resouceDir+"Allele_Master.txt")
	for line in F0:
		LINE = line.rstrip('\r\n')
		lst = LINE.split('\t')
		if lst[0] == "PanTro2_Coord":
			pass
		else:
			#Load the allele into each species hash
			Allele = '->'.join([lst[1],lst[3]])
			ChimpAllele[lst[0]] = Allele
			HumanAllele[lst[2]] = Allele
	F0.close();
	#Annotate Human
	OUT = open("Hg_Subs_IntersectAnottate.txt", "w")
	F1 = open(resouceDir+"Human_Subs_Intersect.txt")
	for line in F1:
		LINE = line.rstrip('\r\n')
		lst = LINE.split('\t')
		if (bool(re.search( 'NoSUBsCntrl', lst[3])) is True):
			Entry = '\t'.join([lst[0],lst[1],lst[2],"NA","NA","NA","NA","NA"])
		else:
			coord = lst[4]+':'+lst[5]+'-'+lst[6]
			#Human fragments are zero based named
			Pos = (int(lst[5])-int(lst[1]))+1
			Entry = '\t'.join([lst[0],lst[1],lst[2],lst[4],lst[5],lst[6],str(Pos),HumanAllele[coord]])
		print >>OUT, Entry
	F1.close();
	OUT2 = open("PanTro2_Subs_IntersectAnottate.txt", "w")
	F2 = open(resouceDir+"PanTro2_Subs_Intersect.txt")
	for line in F2:
		LINE = line.rstrip('\r\n')
		lst = LINE.split('\t')
		if (bool(re.search( 'NoSUBsCntrl', lst[3])) is True):
			Entry = '\t'.join([lst[0],lst[1],lst[2],"NA","NA","NA","NA","NA"])
		else:
			coord = lst[4]+':'+lst[5]+'-'+lst[6]
			#Chimp fragments are one based named
			Pos = (int(lst[5])+1-int(lst[1]))+1
			Entry = '\t'.join([lst[0],lst[1],lst[2],lst[4],lst[5],lst[6],str(Pos),ChimpAllele[coord]])
		print >>OUT2, Entry
	F2.close();

#Anntotion Sub-Routine C
##This Subroutine takes The master Allele file and annotates the intersected human and chimp fragment fiels with position and allele info
#Original Script Call: python AnnoteC.py 
def AnnoterC( resouceDir ):
	HumanPos = {}
	HumanAllele = {}
	#Load in human data
	F0 = open(resouceDir+"Hg_Subs_IntersectAnottate.txt")
	for line in F0:
		LINE = line.rstrip('\r\n')
		lst = LINE.split('\t')
		coord = lst[0]+':'+lst[1]+'-'+lst[2]
		#Add Position
		if coord in HumanPos:
			HumanPos[coord].append(lst[6])
		else:
			HumanPos[coord] = [lst[6]]
		#Add Allele
		if coord in HumanAllele:
			HumanAllele[coord].append(lst[7])
		else:
			HumanAllele[coord] = [lst[7]]
	F0.close();
	ChimpPos = {}
	ChimpAllele = {}
	#Load in chimp data
	F1 = open(resouceDir+"PanTro2_Subs_IntersectAnottate.txt")
	for line in F1:
		LINE = line.rstrip('\r\n')
		lst = LINE.split('\t')
		coord = lst[0]+':'+lst[1]+'-'+lst[2]
		#Add Position
		if coord in ChimpPos:
			ChimpPos[coord].append(lst[6])
		else:
			ChimpPos[coord] = [lst[6]]
		#Add Allele
		if coord in ChimpAllele:
			ChimpAllele[coord].append(lst[7])
		else:
			ChimpAllele[coord] = [lst[7]]
	F1.close();
	#Print out human
	HUMAN = open("Human_Frag_BPandAllele_Anonted.txt", "w")
	F2 = open(resouceDir+"Hg19_Fragments_Srt_Fin.bed")
	for line in F2:
		LINE = line.rstrip('\r\n')
		lst = LINE.split('\t')
		coord = lst[0]+':'+lst[1]+'-'+lst[2]
		POS = ';'.join(HumanPos[coord])
		Allele = ';'.join(HumanAllele[coord])
		Entry = '\t'.join([lst[0],lst[1],lst[2],lst[3],Allele,POS])
		print >>HUMAN, Entry	
	F2.close();
	#Print out chimp
	CHIMP = open("Chimp_Frag_BPandAllele_Anonted.txt", "w")
	F3 = open(resouceDir+"PanTro2_Fragments_srt.bed")
	for line in F3:
		LINE = line.rstrip('\r\n')
		lst = LINE.split('\t')
		coord = lst[0]+':'+lst[1]+'-'+lst[2]
		POS = ';'.join(ChimpPos[coord])
		Allele = ';'.join(ChimpAllele[coord])
		Entry = '\t'.join([lst[0],lst[1],lst[2],lst[3],Allele,POS])
		print >>CHIMP, Entry	
	F3.close();

#Anntotion Sub-Routine D
###This SubRoutine annotates the master counts file with the allele changes and fragment positions of those changes
#Original Script Call: python AnnoterD.py 
def AnnoterD( resouceDir ):
	MasterPos = {}
	MasterAllele = {}
	#Load in human data
	F0 = open(resouceDir+"Human_Frag_BPandAllele_Anonted.txt")
	for line in F0:
		LINE = line.rstrip('\r\n')
		lst = LINE.split('\t')
		#Adjust No Subs control label to be a match to aligned file
		if (bool(re.search( 'NoSUBsCntrl', lst[3])) is True):
			AdjName = lst[3].replace("NoSUBsCntrl","NoSUBsCnTrl")
			MasterPos[AdjName] = lst[4]
			MasterAllele[AdjName] = lst[5]
		else:
			MasterPos[lst[3]] = lst[4]
			MasterAllele[lst[3]] = lst[5]
	F0.close();
	#Load in chimp data
	F1 = open(resouceDir+"Chimp_Frag_BPandAllele_Anonted.txt")
	for line in F1:
		LINE = line.rstrip('\r\n')
		lst = LINE.split('\t')
		if (bool(re.search( 'NoSUBsCntrl', lst[3])) is True):
			AdjName = lst[3].replace("NoSUBsCntrl","NoSUBsCnTrl_Chimp")
			MasterPos[AdjName] = lst[4]
			MasterAllele[AdjName] = lst[5]
		else:
			MasterPos[lst[3]] = lst[4]
			MasterAllele[lst[3]] = lst[5]
	F1.close();
	#Open out and move through main counts file to annotate
	OUT = open("COUNTS_ANNOTATED.txt", "w")
	F2 = open("Annotated_with_coords.txt")
	for line in F2:
		LINE = line.rstrip('\n\r')
		lst = LINE.split('\t')
		if lst[0] == "Tag":
			Entry = '\t'.join(["Tag","Alignment","Ortholog_Seq","Chr","Start","Stop","Species","eVar_BPChange","eVar_Pos",'\t'.join(map( str, lst[7:] ))])
			print >>OUT, Entry
		else:
			if len(lst) < 11:
				#Catches Errors
				print LINE
			else:
				if (bool(re.search( 'Cntrl', lst[1])) is True):
					frag = lst[1].replace("Cntrl","CnTrl")
				else:
					frag = lst[1]
				Entry = '\t'.join([lst[0],lst[1],lst[2],lst[3],lst[4],lst[5],lst[6],MasterPos[frag],MasterAllele[frag],'\t'.join(map( str, lst[7:] ))])
				print >>OUT, Entry

#######Run Pipeline
####Inert Library Processing pipeline
if (Mode == "INERT"):
	##Example Call:
	#python Pipeline.py -r INERT -PI .88 -R1 Inert_R1.fastq.gz -R2 Inert_R2.fastq.gz -t 17 -o Inert_Pipeline_OUT.txt
	if ID is None:
		raise Exception('No Percent ID is specified please specify a Percent ID: 1 > P_ID > 0')
	elif ROne is None:
		raise Exception('Read one file or file paths not specified, please specify ReadOne.fastq.gz with the -R1 flag')
	elif RTwo is None:
		raise Exception('Read two file or file paths not specified, please specify ReadTwo.fastq.gz with the -R2 flag')
	elif (resource is None):
			raise Exception('Must specify path to resource files include trailing back slash!')
	else:
		print >>PipeOUT, ''.join(["\nPercent ID set at: ",str(ID)])
		if (ID <= 0 or ID > 1):
			raise Exception('Percent ID is out of Specified Range Please Specify as: 1 > P_ID > 0')
		elif len(RTwo) != len(RTwo):
			warnings.warn("You Specified an Unequal Number of Read One and Read Two Files or Auto Completed File Paths!")
			print >>PipeOUT, "You Specified an Unequal Number of Read One and Read Two Files or File Paths!"			
		else:
			#Print out Called files to pipeline log file
			print >>PipeOUT, "\nRead One File(s)"
			for r in ROne:
				print >>PipeOUT, "%s" % r
		
			print >>PipeOUT, "\nRead Two File(s)"
			for r in RTwo:
				print >>PipeOUT, "%s" % r
			print >>PipeOUT, "Pipeline Command Line Calls:"
			#Pull Reads
			cmd = ''.join(["zcat ", ' '.join(ROne), " > Inert_ReadONE.fastq"])
			print >>PipeOUT, cmd
			subprocess.Popen(cmd, shell=True).wait()
			cmd = ''.join(["zcat ", ' '.join(RTwo), " > Inert_ReadTWO.fastq"])
			print >>PipeOUT, cmd
			subprocess.Popen(cmd, shell=True).wait()
			#Orient Reads so that 5 + Strand is Read One and 3' - Strand is Read Two:
			Strander( "Inert_ReadONE.fastq", "Inert_ReadTWO.fastq" )
			#Trim Low quality ends of reads 
			cmd = "time fastx_trimmer -Q 33 -f 10 -i Inert_ReadOne.fastq -o Inert_R1_Trimmed.fastq"
			print >>PipeOUT, cmd
			subprocess.Popen(cmd, shell=True).wait()
			#Align the mate pairs into a single contig
			cmd = "time pear -j 20 -q 33 -f Inert_R1_Trimmed.fastq -r Inert_ReadTwo.fastq -o Inert_AlignedReads >Inert_Pear_Log.txt"
			print >>PipeOUT, cmd
			subprocess.Popen(cmd, shell=True).wait()
			#Trim Sequences into tags and CREs
			cmd = "cutadapt -e 0.16 -g ACTGGCCGCTTGACG -o Inert_5P_Cut_Seqs.fastq Inert_AlignedReads.assembled.fastq > Inert_5PrimeTrim_Log.txt"
			print >>PipeOUT, cmd
			subprocess.Popen(cmd, shell=True).wait()
			cmd = "cutadapt -e 0.16 -a AGATCGGAAGAGCGTC -o Inert_3PCut_Seqs.fastq Inert_5P_Cut_Seqs.fastq > Inert_3PrimeTrim_Log.txt"
			print >>PipeOUT, cmd
			subprocess.Popen(cmd, shell=True).wait()
			cmd = "cutadapt -e 0.16 -g CACTGCGGCTCCTGCGGTACCTCTAGA -o Inert_Tags.fastq Inert_3PCut_Seqs.fastq > Inert_Tag_Trim_Log.txt"
			print >>PipeOUT, cmd
			subprocess.Popen(cmd, shell=True).wait()
			cmd = "cutadapt -e 0.16 -a CACTGCGGCTCCTGCGGTACCTCTAGA -o Inert_Cres.fastq Inert_3PCut_Seqs.fastq  > Inert_Cre_Trim_Log.txt"
			print >>PipeOUT, cmd
			subprocess.Popen(cmd, shell=True).wait()
			#Match Cres to Tags for Translation file
			Matcher( "Inert_Cres.fastq", "Inert_Tags.fastq" )
			#Align Cres to library:
			BTtmp = []
			for file in os.listdir(resource):
				if file.endswith('.rev.1.bt2'):
					BTtmp.append(str(file[:-10]))
			if(len(BTtmp) >1):
				raise Exception('More than one Bowtie2 index present in resource folder!')
			else:
				IndexName = ''.join([resource,BTtmp[0]])
			print IndexName
			cmd = "bowtie2 -S Inert.sam --end-to-end -x "+IndexName+" -U Inert_Final_Cres.fastq"
			print >>PipeOUT, cmd
			subprocess.Popen(cmd, shell=True).wait()
			cmd = "samtools view -Sb  Inert.sam  > Inert.bam"
			print >>PipeOUT, cmd
			subprocess.Popen(cmd, shell=True).wait()
			#Annotate each alignment with alignment deatail, Tags, and translate alignments to CREs
			## Resource file path integrate here!!
			TransName = ''.join([resource,"Masked_LibTranslation.txt"])
			GenomicParse("Inert.bam", TransName, "Inert_Compiled_Cre_tags.txt", TagLength)
			#Score Alingments with Percent Identity
			###LOOK AT R SCRIPT DIFFS
			cmd = "cp "+resource+"ProcessR_Pipeline.R ."
			print >>PipeOUT, cmd
			subprocess.Popen(cmd, shell=True).wait()
			cmd = ''.join(["R CMD BATCH \'--args ",str(ID),"\' ProcessR_Pipeline.R"])
			print >>PipeOUT, cmd
			subprocess.Popen(cmd, shell=True).wait()
			#Cleanup files 
			cmd = "sed 's/ //g' Parsed_Cigar_withTAG_TrimedPI_gt88_temp.txt > temp"
			print >>PipeOUT, cmd
			subprocess.Popen(cmd, shell=True).wait()
			cmd = 'sed \'s/"//g\' temp > temp2'
			print >>PipeOUT, cmd
			subprocess.Popen(cmd, shell=True).wait()
			cmd = "mv temp2 Parsed_Cigar_withTAG_TrimedPI_gt88_temp.txt"
			print >>PipeOUT, cmd
			subprocess.Popen(cmd, shell=True).wait()
			cmd = "rm temp"
			print >>PipeOUT, cmd
			subprocess.Popen(cmd, shell=True).wait()
			#Filter out CREs with tags that represent multiple CREs
			##LOOK FOR DIFFS HERRE TOO
			MultipleTagged( "Parsed_Cigar_withTAG_TrimedPI_gt88_temp.txt" )
			Tabulator( "UniqTags_Translated_Parsed_Cigar_withTAG_TrimedPI_gt88.txt" )
####Inert Library Processing pipeline
elif (Mode == "INERT-HIQ"):
	##Example Call:
	#python Pipeline.py -r INERT-HIQ -R1 Inert_R1.fastq.gz -R2 Inert_R2.fastq.gz -t 17 -o Inert_Pipeline_OUT.txt
	if ROne is None:
		raise Exception('Read one file or file paths not specified, please specify ReadOne.fastq.gz with the -R1 flag')
	elif RTwo is None:
		raise Exception('Read two file or file paths not specified, please specify ReadTwo.fastq.gz with the -R2 flag')
	elif (resource is None):
			raise Exception('Must specify path to resource files include trailing back slash!')
	else:
		if len(RTwo) != len(RTwo):
			warnings.warn("You Specified an Unequal Number of Read One and Read Two Files or Auto Completed File Paths!")
			print >>PipeOUT, "You Specified an Unequal Number of Read One and Read Two Files or File Paths!"			
		else:
			#Print out Called files to pipeline log file
			print >>PipeOUT, "\nRead One File(s)"
			for r in ROne:
				print >>PipeOUT, "%s" % r
			print >>PipeOUT, "\nRead Two File(s)"
			for r in RTwo:
				print >>PipeOUT, "%s" % r
			print >>PipeOUT, "Pipeline Command Line Calls:"
			#Pull Reads
			cmd = ''.join(["zcat ", ' '.join(ROne), " >Inert_ReadONE.fastq"])
			print >>PipeOUT, cmd
			subprocess.Popen(cmd, shell=True).wait()
			cmd = ''.join(["zcat ", ' '.join(RTwo), " >Inert_ReadTWO.fastq"])
			print >>PipeOUT, cmd
			subprocess.Popen(cmd, shell=True).wait()
			#Orient Reads so that 5 + Strand is Read One and 3' - Strand is Read Two:
			Strander("Inert_ReadONE.fastq", "Inert_ReadTWO.fastq")
			#Trim Low quality ends of reads 
			cmd = "time fastx_trimmer -Q 33 -f 10 -i Inert_ReadOne.fastq -o Inert_R1_Trimmed.fastq"
			print >>PipeOUT, cmd
			subprocess.Popen(cmd, shell=True).wait()
			#Align the mate pairs into a single contig
			cmd = "time pear -j 20 -q 33 -f Inert_R1_Trimmed.fastq -r Inert_ReadTwo.fastq -o Inert_AlignedReads >Inert_Pear_Log.txt"
			print >>PipeOUT, cmd
			subprocess.Popen(cmd, shell=True).wait()
			#Trim Sequences into tags and CREs
			cmd = "cutadapt -e 0.16 -g ACTGGCCGCTTGACG -o Inert_5P_Cut_Seqs.fastq Inert_AlignedReads.assembled.fastq >Inert_5PrimeTrim_Log.txt"
			print >>PipeOUT, cmd
			subprocess.Popen(cmd, shell=True).wait()
			cmd = "cutadapt -e 0.16 -a AGATCGGAAGAGCGTC -o Inert_3PCut_Seqs.fastq Inert_5P_Cut_Seqs.fastq >Inert_3PrimeTrim_Log.txt"
			print >>PipeOUT, cmd
			subprocess.Popen(cmd, shell=True).wait()
			cmd = "cutadapt -e 0.16 -g CACTGCGGCTCCTGCGGTACCTCTAGA -o Inert_Tags.fastq Inert_3PCut_Seqs.fastq >Inert_Tag_Trim_Log.txt"
			print >>PipeOUT, cmd
			subprocess.Popen(cmd, shell=True).wait()
			cmd = "cutadapt -e 0.16 -a CACTGCGGCTCCTGCGGTACCTCTAGA -o Inert_Cres.fastq Inert_3PCut_Seqs.fastq  >Inert_Cre_Trim_Log.txt"
			print >>PipeOUT, cmd
			subprocess.Popen(cmd, shell=True).wait()
			#Some tag quality checks
			tag_check("Inert_Tags.fastq", "Inert_Cres.fastq", TagLength)
			#Match Cres to Tags for Translation file
			Matcher("Inert_Cres_qual.fastq", "Inert_Tags_qual.fastq")
			#Pad FastQ file with adapter sequences to ensure correct mapping of similar frag sequences
			pad_fastq("Inert_Final_Cres.fastq", "Inert_Final_Cres_with-adapters.fastq")			
			#Align Cres to library:
			BTtmp = []
			for file in os.listdir(resource):
				if file.endswith('_with_adapters.rev.1.bt2'):
					BTtmp.append(str(file[:-10]))
			if(len(BTtmp) >1):
				raise Exception('More than one Bowtie2 index present in resource folder!')
			else:
				IndexName = ''.join([resource,BTtmp[0]])
			print IndexName
			cmd = "bowtie2 --rg-id inert -S Inert.sam -x "+IndexName+" -U Inert_Final_Cres_with-adapters.fastq"
			print >>PipeOUT, cmd
			subprocess.Popen(cmd, shell=True).wait()
			cmd = "samtools view -Sb Inert.sam >Inert.bam"
			print >>PipeOUT, cmd
			subprocess.Popen(cmd, shell=True).wait()
			#Annotate each alignment with alignment deatail, Tags, and translate alignments to CREs
			## Resource file path integrate here!!
			TransName = ''.join([resource,"Masked_LibTranslation.txt"])
			GenomicParse("Inert.bam", TransName, "Inert_Compiled_Cre_tags.txt", TagLength)
			cmd = "sort -k2,2 -k3,3 CigarParsedMatched.txt >CigarParsedMatched_sort.txt"
			print >>PipeOUT, cmd
			subprocess.Popen(cmd, shell=True).wait()
			ProcessR_sorted("CigarParsedMatched_sort.txt")
			ProcessR_sorted_v2("R_Processed.tsv","R_Processed_use.tsv","R_Processed_trash.tsv")
			#Cleanup files 
			# cmd = "sed 's/ //g' Parsed_Cigar_withTAG_TrimedPI_gt88_temp.txt > temp"
			# print >>PipeOUT, cmd
			# subprocess.Popen(cmd, shell=True).wait()
			# cmd = 'sed \'s/"//g\' temp > temp2'
			# print >>PipeOUT, cmd
			# subprocess.Popen(cmd, shell=True).wait()
			# cmd = "mv temp2 Parsed_Cigar_withTAG_TrimedPI_gt88_temp.txt"
			# print >>PipeOUT, cmd
			# subprocess.Popen(cmd, shell=True).wait()
			# cmd = "rm temp"
			# print >>PipeOUT, cmd
			# subprocess.Popen(cmd, shell=True).wait()
			# #Filter out CREs with tags that represent multiple CREs
			# ##LOOK FOR DIFFS HERRE TOO
			# MultipleTagged("Parsed_Cigar_withTAG_TrimedPI_gt88_temp.txt")
			# Tabulator("UniqTags_Translated_Parsed_Cigar_withTAG_TrimedPI_gt88.txt")
####Competent Library Processing pipeline
elif (Mode == "COMP"):
	##Example Call:
	#python Pipeline.py -r COMP -R1 Comp_R1.fastq.gz -R2 Comp_R2.fastq.gz -t 17 -o Comp_Pipeline_OUT.txt
	#test if the sample designation was set
	if ROne is None:
		raise Exception('Read one file or file paths not specified, please specify ReadOne.fastq.gz with the -R1 flag')
	elif RTwo is None:
		raise Exception('Read two file or file paths not specified, please specify ReadTwo.fastq.gz with the -R2 flag')
	elif len(RTwo) != len(RTwo):
		warnings.warn("You Specified an Unequal Number of Read One and Read Two Files or Auto Completed File Paths!")
		print >>PipeOUT, "You Specified an Unequal Number of Read One and Read Two Files or File Paths!"
	else:
		#Print out Called files to pipeline log file
		print >>PipeOUT, "\nRead One File(s)"
		for r in ROne:
			print >>PipeOUT, "%s" % r
		
		print >>PipeOUT, "\nRead Two File(s)"
		for r in RTwo:
			print >>PipeOUT, "%s" % r
		#Pull Reads
		cmd = ''.join(["zcat ", ' '.join(ROne), " > Comp_ReadONE.fastq"])
		print >>PipeOUT, cmd
		subprocess.Popen(cmd, shell=True).wait()
		cmd = ''.join(["zcat ", ' '.join(RTwo), " > Comp_ReadTWO.fastq"])
		print >>PipeOUT, cmd
		subprocess.Popen(cmd, shell=True).wait()
		#Trim poor seq quality at end of reads
		cmd = "time fastx_trimmer -Q 33 -l 131 -i Comp_ReadONE.fastq -o Comp_R1_Trimmed.fastq"
		print >>PipeOUT, cmd
		subprocess.Popen(cmd, shell=True).wait()
		cmd = "time fastx_trimmer -Q 33 -l 126 -i Comp_ReadTWO.fastq -o Comp_R2_Trimmed.fastq"
		print >>PipeOUT, cmd
		subprocess.Popen(cmd, shell=True).wait()
		#Orient Reads so that 5 + Strand is Read One and 3' - Strand is Read Two
		Comp_TAG_Strander( "Comp_R1_Trimmed.fastq", "Comp_R2_Trimmed.fastq" )
		#Align the mate pairs into a single contig
		cmd = "time pear -j 20 -q 33 -f Comp_ReadOne.fastq -r Comp_ReadTwo.fastq -o Comp_AlignedReads > Comp_Pear_Log.txt"
		print >>PipeOUT, cmd
		subprocess.Popen(cmd, shell=True).wait()
		cmd = "cutadapt -e 0.16 -g GTGTAATAATTCTAGA -o Comp_5P_Cut_Seqs.fastq Comp_AlignedReads.assembled.fastq > Comp_5PrimeTrim_Log.txt"
		print >>PipeOUT, cmd
		subprocess.Popen(cmd, shell=True).wait()
		cmd = "cutadapt -e 0.16 -a AGATCGGAAGAGCGTC -o Comp_Tags.fastq Comp_5P_Cut_Seqs.fastq > Comp_3PrimeTrim_Log.txt"
		print >>PipeOUT, cmd
		subprocess.Popen(cmd, shell=True).wait()
		#Tabulate Tag Counts - Outputs: CompLib_Tag_Counts.txt
		TagTabulator("Comp_Tags.fastq", "CompLib", TagLength)
####Experimental Read Processing Pipeline
elif (Mode == "EXP"):
	##Example call:
	#python Pipeline.py -r EXP -s pDNA_Test -R1 pDNA_R1.fastq.gz -R2 pDNA_R2.fastq.gz -t 17 -o pNDA_EXP_Pipeline_OUT.txt
	#python Pipeline.py -r EXP -s cDNA_Test -R1 cDNA_R1.fastq.gz -R2 cDNA_R2.fastq.gz -t 17 -o cNDA_EXP_Pipeline_OUT.txt
	#Test if the sample designation was set
	if Samp is None:
		raise Exception('Sample Name Not Specified please specify ex: Rep2_1_cDNA')
	elif ROne is None:
		raise Exception('Read one file or file paths not specified, please specify ReadOne.fastq.gz with the -R1 flag')
	elif RTwo is None:
		raise Exception('Read two file or file paths not specified, please specify ReadTwo.fastq.gz with the -R2 flag')
	elif len(RTwo) != len(RTwo):
		warnings.warn("You Specified an Unequal Number of Read One and Read Two Files or Auto Completed File Paths!")
		print >>PipeOUT, "You Specified an Unequal Number of Read One and Read Two Files or File Paths!"
	else:
		#Print out called files to pipeline log file
		print >>PipeOUT, "\nRead One File(s)"
		for r in ROne:
			print >>PipeOUT, "%s" % r
		print >>PipeOUT, "\nRead Two File(s)"
		for r in RTwo:
			print >>PipeOUT, "%s" % r
		print >>PipeOUT, ''.join(["\nSample Name Specified As: ",Samp])
		#Pull reads
		cmd = ''.join(["zcat ", ' '.join(ROne), " > ",Samp,"_ReadONE.fastq"])
		print >>PipeOUT, cmd
		subprocess.Popen(cmd, shell=True).wait()
		cmd = ''.join(["zcat ", ' '.join(RTwo), " > ",Samp,"_ReadTWO.fastq"])
		print >>PipeOUT, cmd
		subprocess.Popen(cmd, shell=True).wait()
		#Trim poor seq quality at end of reads
		cmd = ''.join(["time fastx_trimmer -Q 33 -l 131 -i ",Samp,"_ReadONE.fastq -o ",Samp,"_ReadONE_Trimmed.fastq"])
		print >>PipeOUT, cmd
		subprocess.Popen(cmd, shell=True).wait()
		cmd = ''.join(["time fastx_trimmer -Q 33 -l 126 -i ",Samp,"_ReadTWO.fastq -o ",Samp,"_ReadTWO_Trimmed.fastq"])
		print >>PipeOUT, cmd
		subprocess.Popen(cmd, shell=True).wait()
		#Orient reads so that 5'-strand is read one and 3'-strand is read two
		Exp_TAG_Strander( ''.join([Samp,"_ReadONE_Trimmed.fastq"]), ''.join([Samp,"_ReadTWO_Trimmed.fastq"]), Samp )
		#Align the mate pairs into a single contig
		cmd = ''.join(["time pear -j 20 -q 33 -f ",Samp,"__ReadOne.fastq"," -r ",Samp,"__ReadTwo.fastq"," -o ",Samp,"_AlignedReads > Pear_",Samp,"_Log.txt"])
		print >>PipeOUT, cmd
		subprocess.Popen(cmd, shell=True).wait()
		#Trim to get tag sequences
		cmd = ''.join(["cutadapt -e 0.16 -g GTGTAATAATTCTAGA -o ",Samp,"_5P_Cut_Seqs.fastq ",Samp,"_AlignedReads.assembled.fastq > ",Samp,"_5PrimeTrim_Log.txt"])
		print >>PipeOUT, cmd
		subprocess.Popen(cmd, shell=True).wait()
		cmd = ''.join(["cutadapt -e 0.16 -a AGATCGGAAGAGCGTC -o ",Samp,"_Tags.fastq ",Samp,"_5P_Cut_Seqs.fastq > ",Samp,"_3PrimeTrim_Log.txt"])
		print >>PipeOUT, cmd
		subprocess.Popen(cmd, shell=True).wait()
		#Tabulate tag counts
		TagTabulator(''.join([Samp,"_Tags.fastq"]), Samp, TagLength)
####Filter tags present in competent library which are present in the inert library
elif (Mode == "TAGComparison"):
	##Example call:
	#python Pipeline.py -r TAGComparison -I Inert_Tag_Counts.txt -C CompLib_Tag_Counts.txt -o ComparisonTags_Pipeline_OUT.txt
	if (inertTAG is None or compTAG is None):
		raise Exception('Must specify Inert and Competent Library Tag Count Files')
	else:
		TagFilter( inertTAG,compTAG )
####Create a master annotation file from all tag calls across replicates and competent library sequencing
elif (Mode == "ANOTE"):
	##Example call:
	##python Pipeline.py -r ANOTE -R /File/Path/to/Resource/Files/ -I Inert_Tag_Counts.txt -C CompLib_Tag_Counts.txt -ET Rep1_Trans_Tag_Counts.txt,Rep2_Trans_Tag_Counts.txt -EP Rep1_Plas_Tag_Counts.txt,Rep2_Plas_Tag_Counts.txt -o ANNOTATION_Pipeline_OUT.txt
	if inertTAG is None:
		raise Exception('Must specify Inert Library Tag Counts File')
	elif resource is None:
		raise Exception('Must Specify a Path to MPRA Anotation Resource Files')
	elif compTAG is None:
		warnings.warn("Warning: You Did Not Specify a Competent Libray Sequencing Tag Counts File!!!")
		print >>PipeOUT, "Warning: You Did Not Specify a Competent Libray Sequencing Tag Counts File!!!"
	elif len(plasTAG) != len(tranTAG):
		warnings.warn("You Specified an Unequal Number of Plasmid and Transcript Tag-Count Files!!!!")
		print >>PipeOUT, "You Specified an Unequal Number of Plasmid and Transcript Tag-Count Files!!!!"	
	else:
		Compiler( inertTAG, compTAG, plasTAG, tranTAG )
		AnnoterA( resource )
		AnnoterB( resource )
		AnnoterC( resource )
		AnnoterD( resource )
else:
	raise Exception('Run Mode Improperly Specified please specify as COMP  INERT')
#Close the pipeline call output
PipeOUT.close()
