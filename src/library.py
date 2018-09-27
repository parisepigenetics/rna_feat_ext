#-*- coding: utf-8 -*-

# Librairie partagée (ScriptA and ScriptB)

# Import
import subprocess
import re
import pandas as pd
import numpy as np

#FIXME we need to do a bit of refactoring in names and objects!

#***************************************
# Sequence class and Functions Associated
#***************************************
class Seq(object):
    # Sequence define by its header and sequence and functions associated
    def __init__(self, seqrecord_object):
        self.header = seqrecord_object.id
        self.seq = seqrecord_object.seq
        # getAttributes()

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
        # Write every 3'UTR sequence in a fasta file named "p3utr.fasta".
        with open("P3UTR", "a") as p3utr_fasta:
            p3UTR = str(self.seq)[int(self.cDNAEnd):]
            if len(p3UTR) > 0:
                p3utr_fasta.write("{}\n".format(">" + str(self.transcriptID)))
                p3utr_fasta.write("{}\n".format(p3UTR))
        return(len(p3UTR))

    def writeP5UTR(self):
        # Write every 5'UTR sequence in a fasta file named "p3utr.fasta".
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
                    if(len(ORF) > ideal_size):
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
                        ORF = str(self.seq)[int(self.cDNAStart) +
                                            start.start():start.start()+stop.end()]
                        if(len(ORF) > ideal_size):
                            #print "\nORF",ORF, start.start() + stop.start()
                            ORF_indices.append("{}:{};".format(
                                int(self.cDNAStart) + start.start(), start.start()+stop.end()))
                            codons_stops_visite.append(start.start() + stop.start())
                            break
        return(["".join(ORF_indices), len(ORF_indices)])

    def getFeatures(self):
        self.getAttributes()
        tab = pd.DataFrame()
        tab["ensembl_gene_id"] = pd.Series(self.geneID).astype(str)
        tab["ensembl_transcript_id"] = pd.Series(self.transcriptID).astype(str)
        #tab["uORF"] = pd.Series(self.get_uORF()[1])
        tab["5PLen"] = pd.Series(self.writeP5UTR()).astype(float)
        #tab["dORF"] = pd.Series(self.get_dORF()[1])
        tab["3PLen"] = pd.Series(self.writeP3UTR()).astype(float)
        tab["Kozak_Context"] = pd.Series(self.getKozak(10, 20)[0]).astype(str)
        tab["Kozak_Sequence"] = pd.Series(self.getKozak(10, 20)[1]).astype(str)
        return(tab)


def RNAfold_calcul(inputA):
    with open(inputA, "r") as input:
        sortie = subprocess.check_output('RNAfold --noPS --verbose --jobs', stdin=input, shell=True)
        foldrex = re.compile('(-[0-9]+\.[0-9]+|\s0\.0)')
        foldinf = re.compile('ENS[0-9A-Z]*')
        mfe = pd.Series(foldrex.findall(sortie), index=foldinf.findall(sortie))
        return(mfe)


def RBPs_Motif(fimo_tab):
    dic_motif = {}
    for i in range(len(fimo_tab)):
        id_info = fimo_tab.iloc[i:i+1, 2:3].values[0][0]
        if id_info is not np.NaN:
            regex = re.compile('ENS[A-Z]*T[0-9]*')
            transcript = regex.findall(id_info)[0]
            if transcript not in dic_motif:
                dic_motif[transcript] = [fimo_tab.iloc[i:i+1, 0:1].values[0][0]]
            else:
                mtf_id = fimo_tab.iloc[i:i+1, 0:1].values[0][0]
                mtif_split = dic_motif[transcript][0].split(";")
                if mtf_id not in mtif_split:
                    dic_motif[transcript][0] = dic_motif[transcript][0] + '***' + mtf_id
    tab = pd.DataFrame(dic_motif, index=range(1)).transpose()
    tab.columns = ['motif_ID']

    return tab
