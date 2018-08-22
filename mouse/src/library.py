#-*- coding: utf-8 -*-

# Librairie partagée (ScriptA and ScriptB)

# Import
import subprocess
import re
import math
import pandas as pd



#***************************************
#Sequence class and Functions Associated
#***************************************
class Seq(object):
	# Sequence define by its header and sequence and functions associated
	def __init__(self, seqrecord_object):
		self.header = seqrecord_object.id
		self.seq = seqrecord_object.seq
		#getAttributes()


	def getAttributes(self):
		# TODO get the attibute names directly from the fasta header. Minimise hard coding.
		attributes = self.header.split("|")
		self.geneID = attributes[0].split(":")[1]
		self.transcriptID = attributes[1].split(":")[1]
		self.geneName = attributes[2].split(":")[1]
		self.UTR5PEnd = attributes[3].split(":")[1]
		self.UTR5PStart = attributes[4].split(":")[1]
		self.UTR3PEnd = attributes[5].split(":")[1]
		self.UTR3PStart = attributes[6].split(":")[1]
		self.cDNAStart = attributes[7].split(":")[1]
		self.cDNAEnd = attributes[8].split(":")[1]


	def writeP3UTR(self):
		#Write every 3'UTR sequence in a fasta file named "p3utr.fasta".
		with open("P3UTR", "a") as p3utr_fasta:
			p3UTR = str(self.seq)[int(self.cDNAEnd):]
			if len(p3UTR) > 0:
				p3utr_fasta.write("{}\n".format(">" + str(self.transcriptID)))
				p3utr_fasta.write("{}\n".format(p3UTR))
		return(len(p3UTR))


	def writeP5UTR(self):
		#Write every 5'UTR sequence in a fasta file named "p3utr.fasta".
		with open("P5UTR", "a") as p5utr_fasta:
			p5UTR = str(self.seq)[0:int(self.cDNAStart)]
			if len(p5UTR) > 0:
				p5utr_fasta.write("{}\n".format(">" + str(self.transcriptID)))
				p5utr_fasta.write("{}\n".format(p5UTR))
		return(len(p5UTR))


	def getKozak(self, k, j):
		'''
		Input:
		k: default, k = 10: number of nucleotide before START-1 (A)
		and after START+1 (G) of the START codon (ATG)
		>> k---START---k is the Kozak sequence
		j: default, j = 20: number of nucleotide to select around
		the kozak sequence
		Output:
		kozak: which is the sequence of the context of the Kozak
		'''
		startM1 = int(self.cDNAStart) - 1  # position start-1
		startP1 = int(self.cDNAStart) + 1  # position start+1

		contextM1 = (startM1 - k) - j
		contextP1 = (startP1 + k) + j

		# Cas ou cadre du context_kozak dépasse l'indice de début de sequence
		# Est_ce Normal ou pas ??
		contextM1 = abs(contextM1)

		contextKozak = str(self.seq)[contextM1:contextP1]
		seqKozak = str(self.seq)[(startM1-k):(startP1+k)]

		return([contextKozak, seqKozak])

	def get_uORF(self):
		self.getAttributes()
		ideal_size = 300
		ORF_indices = []
		codons_stops_visite = []
		codon_start = re.compile('ATG')
		codon_stop = re.compile('T(AG|AA|GA)')

		if codon_start.search(str(self.seq)[0:int(self.cDNAStart)]) is None:
			return('not uORF detected')

		for start in codon_start.finditer(str(self.seq)[0:int(self.cDNAStart)]):
			for stop in codon_stop.finditer(str(self.seq)[start.start():]):
				if(stop.start() % 3 == 0 and stop.start() not in codons_stops_visite):
					ORF = str(self.seq)[start.start():start.start()+stop.end()]
					if(len(ORF)>ideal_size):
						#print str(self.seq)[start.start():start.start()+stop.end()]
						ORF_indices.append("{}:{};".format(start.start(), stop.end()))
						codons_stops_visite.append(stop.start())
						break
		return(["".join(ORF_indices), len(ORF_indices)])

	def get_dORF(self):
		self.getAttributes()
		ideal_size = 300
		ORF_indices = []
		codons_stops_visite = []
		codon_start = re.compile('ATG')
		codon_stop = re.compile('T(AG|AA|GA)')

		if codon_start.search(str(self.seq)[int(self.cDNAStart):]) is None:
			return('not dORF detected')

		for start in codon_start.finditer(str(self.seq)[int(self.cDNAStart):]):
			for stop in codon_stop.finditer(str(self.seq)[start.start():]):
				if(stop.start() % 3 == 0):
					if stop.start() + start.start() not in codons_stops_visite:
						ORF = str(self.seq)[int(self.cDNAStart) + start.start():start.start()+stop.end()]
						if(len(ORF) > ideal_size):
							#print "\nORF",ORF, start.start() + stop.start()
							ORF_indices.append("{}:{};".format(int(self.cDNAStart) + start.start(), start.start()+stop.end()))
							codons_stops_visite.append(start.start() + stop.start())
							break
		return(["".join(ORF_indices), len(ORF_indices)])

	def getFeatures(self):
		self.getAttributes()
		tab = pd.DataFrame()
		tab["ensembl_gene_id"] = pd.Series(self.geneID).astype(str)
		tab["ensembl_transcript_id"] = pd.Series(self.transcriptID).astype(str)
		tab["uORF"] = pd.Series(self.get_uORF()[1])
		tab["5PLen"] = pd.Series(self.writeP5UTR()).astype(float)
		tab["dORF"] = pd.Series(self.get_dORF()[1])
		tab["3PLen"] = pd.Series(self.writeP3UTR()).astype(float)
		tab["Kozak_Context"] = pd.Series(self.getKozak(10,20)[0]).astype(str)
		tab["Kozak_Sequence"] = pd.Series(self.getKozak(10,20)[1]).astype(str)
		return(tab)


#****************
#Others functions
#****************
#For ScriptA
#***********
def get_utr5MAX(cdna_feat_row):
	#print("getting UtrMax...")
	# take the cdna_feat-row with multiples utrs

	utr5s_start = cdna_feat_row["5' UTR start"].values[0].split(";")
	utr5s_end = cdna_feat_row["5' UTR end"].values[0].split(";")

	size_liste = []
	for i in range(len(utr5s_start)):
		size = int(utr5s_end[i])-int(utr5s_start[i])
		size_liste.append(size)
		indice_max = size_liste.index(max(size_liste))
		max_utr5_start = int(utr5s_start[indice_max])
		max_utr5_end = int(utr5s_end[indice_max])
	return([max_utr5_start, max_utr5_end])


def get_utr3MAX(cdna_feat_row):
	#print("getting UtrMax...")
	utr3s_start = cdna_feat_row["3' UTR start"].values[0].split(";")
	utr3s_end = cdna_feat_row["3' UTR end"].values[0].split(";")

	size_liste = []
	for i in range(len(utr3s_start)):
		size = int(utr3s_end[i])-int(utr3s_start[i])
		size_liste.append(size)
		indice_max = size_liste.index(max(size_liste))
		max_utr3_start = int(utr3s_start[indice_max])
		max_utr3_end = int(utr3s_end[indice_max])

	return([max_utr3_start, max_utr3_end])


def get_cDNAMIN(cdna_feat_row):
	#print("getting cDNA_Start_Min...")
	# take the cdna_feat-row with multiples utrs
	cDNA_start = cdna_feat_row["cDNA coding start"].values[0].split(";")
	cDNA_end = cdna_feat_row["cDNA coding end"].values[0].split(";")

	min_cDNA_start = min(map(int, cDNA_start))
	min_CDNA_end = max(map(int, cDNA_end))

	return([min_cDNA_start, min_CDNA_end])

def txt2fasta(cdna_feat_table, fastaOut):
	for i in range(cdna_feat_table.shape[0]):
		ligne = pd.DataFrame(cdna_feat_table.loc[i, :]).transpose()
		fastaOut.write(">GeneID:{}|TranscriptID:{}|GeneName:{}|5P_UTR_end:{}|5P_UTR_start:{}|3P_UTR_end:{}|3P_UTR_end:{}|cDNAstart:{}|cDNAend:{}\n".format(
		ligne["Gene stable ID"].values[0], ligne["Transcript stable ID"].values[0],
		ligne["Gene name"].values[0], ligne["5' UTR end"].values[0], ligne["5' UTR start"].values[0],
		ligne["3' UTR end"].values[0], ligne["3' UTR start"].values[0],
		ligne["cDNA coding start"].values[0], ligne["cDNA coding end"].values[0]
		))
		fastaOut.write("{}\n".format(ligne["cDNA sequences"].values[0]))
	print("txt to Fasta conversion done!")


#***********
#For ScriptB
#***********
def RNAfold_calcul(input):
	with open(input, "r") as input:
		sortie = subprocess.check_output('RNAfold --noPS --verbose --jobs', stdin=input, shell=True)
		foldrex = re.compile('(-[0-9]+\.[0-9]+|\s0\.0)')
		foldinf = re.compile('ENS[0-9A-Z]*')
		mfe = pd.Series(foldrex.findall(sortie), index = foldinf.findall(sortie))
		return(mfe)
