# -*- coding: utf-8 *-*
"""Python module for RNA fetures extraction from ENSEMBL derived fasta files.
"""

__version__ = "0.2a02"

import os
import re
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import GC
import tempfile
import shlex


# CLASSES Interface.
class ENSEMBLSeqs(object):
    """Class to represent RNA related sequence features from ENSEMBLself.

    Needs a Bio.Seq.Record generator object to initialise."""

    def __init__(self, bioSeqRecsGen):
        self.gen = bioSeqRecsGen
        self.bioSeqRecs = self.getbioSeqRecs()

    def getbioSeqRecs(self):
        """Expands the Bio.Seq.Rec generator to a list of Bio.Record objects.

        Also puts some SeqIO.record member variables in place, namely the id, the gene name the description and the features."""
        recs = []
        for rec in self.gen:
            # Extract the gene name.
            rec.name = rec.description.split("|")[2].split(":")[1]
            descr = rec.description.split("|")[-1]
            # Extract all the features as a dictionary (apart from the last one which was the description.).
            feat = dict(item.split(":") for item in rec.description.split("|")[1:-1])
            rec.features = feat
            #TODO for the moment all the features are stored as a dictionary and not as proper SeqFeature objects. (perhaps we can stick with that and there is no need to change it.)
            rec.description = descr
            recs.append(rec)
        return recs


class FeaturesExtract(object):
    """Claas to extract features."""

    def __init__(self, bioSeqRecs):
        """Initialise with a list of SeqIO records."""
        self.bioSeqRecs = bioSeqRecs
        # Temporary files of the UTRs
        self.tf5p = tempfile.NamedTemporaryFile(mode="a", delete=False)
        self.tf3p = tempfile.NamedTemporaryFile(mode="a", delete=False)

    def collectFeatures(self):
        """Collect the features that do not need external computations."""
        # Initialise a pabdas data frame.
        pdf = pd.DataFrame(columns = ['ensembl_gene_id', 'gene_name', 'coding_len', '5pUTR_len', '5pUTR_GC', '3pUTR_len', '3pUTR_GC', 'Kozak_Sequence', 'Kozak_Context'])
        for rec in self.bioSeqRecs:
            # Get Koxaks
            seqKozak, contKozak = self.getKozak(rec)
            # Fetch UTRs, lenghts and GCs
            utr3len, utr3gc = self.get3PUTR(rec, self.tf3p)
            utr5len, utr5gc = self.get5PUTR(rec, self.tf5p)
            # Length of coding region.
            codeLen = int(rec.features["cDNAend"]) - int(rec.features["cDNAstart"])
            # Add to pandas data frame.
            pdfe = [rec.features["GeneID"], rec.name, codeLen, utr5len, "{0:.2f}".format(utr5gc), utr3len, "{0:.2f}".format(utr3gc), seqKozak, contKozak]
            pdf.loc[rec.id] = pdfe
        self.tf5p.close()
        self.tf3p.close()
        return pdf

    def getKozak(self, rec, s = 10, c = 20):
        """Extract both Kozak sequence and context from a SeqIO record.
        (s and c define the extremeties of a Koxak sequence and are chosen by convention.)
        return a tuple of Kozak sequence and Kozak context.

        If ATG is located near the 5'UTR start, it is most likely that either seqs will not be retrieved as seld.base[-5:10] returns blank.

        In this case, we test whether the value left to ':' is negative or not.
        If it is < 0, we simply take the seq from 0 as : self.bases[0:self.cDNAStart) + 2) + s]"""
        feat = rec.features
        seq = rec.seq
        # Kozak sequence
        if int(feat["cDNAstart"]) - 1 - s < 0:
            kozakSeq = seq[0:(int(feat["cDNAstart"]) + 2) + s]
        else:
            kozakSeq = seq[(int(feat["cDNAstart"]) - 1 - s):(int(feat["cDNAstart"]) + 2)+ s]
        # Kozak context
        if int(feat["cDNAstart"]) -1 - s - c < 0:
            kozakContext = seq[0:(int(feat["cDNAstart"]) + 2)+ s + c]
        else:
            kozakContext = seq[((int(feat["cDNAstart"]) - 1 - s) - c):(int(feat["cDNAstart"]) + 2) + s + c]
        return(str(kozakSeq), str(kozakContext))

    def get3PUTR(self, rec, tf3p):
        """Fetch 3PUTR sequence, append it to fasta file and return its length and GC content."""
        utr3p = rec.seq[int(rec.features["cDNAend"]):]
        tf3p.write(">{}_3PUTR\n{}\n".format(rec.id, utr3p))
        return(len(utr3p), GC(utr3p))

    def get5PUTR(self, rec, tf5p):
        """Fetch 5PUTR sequence, append it to fasta file and return its length and GC content."""
        utr5p = rec.seq[0:int(rec.features["cDNAstart"])-1]
        tf5p.write(">{}_5PUTR\n{}\n".format(rec.id, utr5p))
        return(len(utr5p), GC(utr5p))

    def calculateFeatures(self):
        """Method to perfom all feature calculation.

        Basically invokes external software to make calculations."""
        # Calculates 5pUTR folding etc.
        fe5p = self.calculateFreeEnergy(self.tf5p.name, "5pUTR")
        fe3p = self.calculateFreeEnergy(self.tf3p.name, "3pUTR")
        #motifs = self.predictBinding()
        # Merge data frames and return.
        return pd.concat([fe5p, fe3p], axis=1, sort=False)

    def calculateFreeEnergy(self, ffile, col):
        """Method to perform the free energy calculation by RNAfold and parsing of the results.

        Needs the RNA Vienna package to be installed."""
        pdf = pd.DataFrame(columns = ['{}_MFE'.format(col), '{}_MfeBP'.format(col)])
        cmd = 'RNAfold --verbose --noPS --jobs -i {}'.format(ffile)
        proc = subprocess.check_output(shlex.split(cmd))
        proc = proc.decode() # Convert it to proper string.
        outList = proc.splitlines() # Retrieve STDOUT.
        for i in range(0, len(outList), 3):
            idt = outList[i][1:-6] # To exclude the _UTR suffix.
            seq = outList[i+1]
            fold = outList[i+2]
            mfeRE = re.compile('-[0-9]+\.[0-9]+')
            mfe = float(mfeRE.search(fold).group())
            mfeBp = mfe / float(len(seq))
            pdf.loc[idt] = [mfe, mfeBp]
        return pdf

    def predictBinding(self, ffile, motifs):
        """Method to run a FIMO search and parse and collect the results.

        Requires the intallation of the MEME suite."""
        #TODO implement the method!!!
        subprocess.call('fimo --verbosity 1 ' + motifs + ' ' + ffile.name, shell=True)
        fimo_tab = pd.read_csv("fimo_out/fimo.tsv", sep="\t")
        fimo_tab = fimo_tab.reset_index(drop = True)
        return None

    def __del__(self):
        """Cleanup the temp files"""
        os.remove(self.tf5p.name)
        os.remove(self.tf3p.name)



# FUNCTIONS
def geneIDs2Fasta(listID, out_fasta, dataset):
    print("Connection to server....")
    server = BiomartServer( "http://www.ensembl.org/biomart/")
    print("....done!")
    dt = server.datasets[dataset]
    print("Connexion to ENSEMBL dataset.")

    #TODO provide the required attributes in an external parameters file.
    listAttrib = ['ensembl_gene_id', 'ensembl_transcript_id',
                  'external_gene_name', 'transcript_start', 'transcript_end',
                  '5_utr_end', '5_utr_start', '3_utr_end', '3_utr_start',
                  'transcription_start_site', 'transcript_biotype',
                  'cdna_coding_start', 'cdna_coding_end', 'cdna', 'description']
    print("Fetching data.....")
    response = dt.search({'filters': {'ensembl_gene_id': listID}, 'attributes' : listAttrib}, header = 1)
    print("...fetching done!")

    # Convert stringIO to pandas data frame.
    print("...Convert string to PD.data.frame.")
    stream = StringIO(response.text)
    cdna = pd.read_table(stream, sep="\t")

    # remove NA and reset index
    cdna = cdna.dropna()
    cdna = cdna.reset_index()

    #TODO TO TEST all these before publishing anything!
    # select only longest 5'UTR and 3'UTR
    print("...Select the longest UTRs!")
    for i in range(cdna.shape[0]):
        ligne = pd.DataFrame(cdna.loc[i,:]).transpose()
        #For 5_UTR_start
        indices = rnaFeaturesLib.get_utr5MAX(ligne)
        if(len(cdna.iloc[i:i+1,:]["5' UTR start"].values[0].split(";")) > 1):
            cdna.iloc[i:i+1,7:8] = indices[1]
            cdna.iloc[i:i+1,8:9] = indices[0]
        if(len(cdna.iloc[i:i+1,:]["3' UTR start"].values[0].split(";")) > 1):
            cdna.iloc[i:i+1,9:10] = indices[1]
            cdna.iloc[i:i+1,10:11] = indices[0]

    # select longest cDNA
    print("...Select the longest cDNA")
    for i in range(cdna.shape[0]):
        ligne = pd.DataFrame(cdna.loc[i,:]).transpose()
        # smallest cDNA start value
        min_cdna_start = rnaFeaturesLib.get_cDNAstartMIN(ligne)
        # therefore: biggest cDNA end value
        max_cdna_end = rnaFeaturesLib.get_cDNAendMAX(ligne)
        cdna.iloc[i, cdna.columns.get_loc('cDNA coding start')] = min_cdna_start
        cdna.iloc[i, cdna.columns.get_loc('cDNA coding end')] = max_cdna_end

    # Exportat to FASTA format
    print("Return FASTA file!")
    rnaFeaturesLib.txt2fasta(cdna, out_fasta)


def get_utr5MAX(cdna_feat_row):
    #take the cdna_feat-row with multiples utrs
    #print("getting UtrMax...")
    utr5s_start = cdna_feat_row["5' UTR start"].values[0].split(";")
    utr5s_end = cdna_feat_row["5' UTR end"].values[0].split(";")
    size_liste = []
    for i in range(len(utr5s_start)):
        size = int(utr5s_end[i])-int(utr5s_start[i])
        size_liste.append(size)
    indice_max = size_liste.index(max(size_liste))
    max_utr5_start = int(utr5s_start[indice_max])
    max_utr5_end = int(utr5s_end[indice_max])
    return([max_utr5_start,max_utr5_end])


def get_utr3MAX(cdna_feat_row):
    #take the cdna_feat-row with multiples utrs
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
    return([max_utr3_start,max_utr3_end])


def get_cDNAstartMIN(cdna_feat_row):
    #take the cdna_feat-row with multiples utrs
    #print("getting cDNA_Start_Min...")
    cDNA_start = cdna_feat_row["cDNA coding start"].values[0].split(";")
    min_cDNA_start = min(map(int, cDNA_start))
    return(min_cDNA_start)


def get_cDNAendMAX(cdna_feat_row):
    #take the cdna_feat-row with multiples utrs
    #print("getting cDNA_End_Max...")
    cDNA_end = cdna_feat_row["cDNA coding end"].values[0].split(";")
    max_cDNA_end = max(map(int, cDNA_end))
    return(max_cDNA_end)


def txt2fasta(cdna_feat_table, fastaOut):
    with open(fastaOut + ".fasta", "w+") as ff:
        for i in range(cdna_feat_table.shape[0]):
            ligne = pd.DataFrame(cdna_feat_table.loc[i,:]).transpose()
            ff.write(">{} |GeneID:{}|GeneName:{}|5P_UTR_end:{}|5P_UTR_start:{}|3P_UTR_end:{}|3P_UTR_end:{}|cDNAstart:{}|cDNAend:{}|{}\n".format(ligne["Transcript stable ID"].values[0], ligne["Gene stable ID"].values[0], ligne["Gene name"].values[0], ligne["5' UTR end"].values[0], ligne["5' UTR start"].values[0], ligne["3' UTR end"].values[0], ligne["3' UTR start"].values[0], ligne["cDNA coding start"].values[0], ligne["cDNA coding end"].values[0], ligne["Gene description"].values[0]))
            ff.write("{}\n".format(ligne["cDNA sequences"].values[0]))
    print("txt to Fasta conversion done!")


def RNAfold_calcul(utr_fasta, out_mfe):
    #print("RNAfold_Calcul ...")
    with open(utr_fasta, "r") as inputfile, open(out_mfe + ".mfe", "w+") as output:
        subprocess.call("RNAfold --noPS --jobs", stdin = inputfile, stdout = output, shell = True)
    print("RNAfold_Calcul ... Done !")
    return(output)


def getFoldingEnergy(input_mfe):
    with open(str(input_mfe), "r") as rnafoldfile:
        tot = rnafoldfile.readlines()
    #TODO use the delimiter of the RNAFold output file to extract the folding energy (the split function). Remove the REGEXPs.
    foldrex = re.compile('(-[0-9]+\.[0-9]+|\s0\.0)')
    strtot = ' '.join(tot) # convertion liste en string
    mfe = foldrex.findall(strtot)
    return(float(mfe[0]))
