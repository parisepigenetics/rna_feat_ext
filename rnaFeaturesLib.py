# -*- coding: utf-8 *-*
"""Python module for RNA fetures extraction from ENSEMBL derived fasta files.
"""
__version__ = "0.3a05"

import os
import sys
import io
import re
import tempfile
import shlex
import shutil
import subprocess
from biomart import BiomartServer
import pandas as pd
from Bio.SeqUtils import GC


# CLASSES Interface.
class ENSEMBLSeqs(object):
    """Class to represent RNA related sequence features from ENSEMBLself.

    Needs a Bio.Seq.Record generator object to initialise."""

    def __init__(self, bioSeqRecsGen, exprTrans):
        self.gen = bioSeqRecsGen
        if exprTrans:
            self.exprTranscripts = self.parse_expressed_transcripts(exprTrans)
        else:
            self.exprTranscripts = None
        self.bioSeqRecs = self.get_bio_seqrec()

    def get_bio_seqrec(self):
        """Expands the Bio.Seq.Rec generator to a list of Bio.Record objects.

        Also puts some SeqIO.record member variables in place, namely the id, the gene name the description and the features."""

        def _construct_bio_seq(recs, rec):
            """Auxiiary function to construct a Bio.Seq object."""
            # Extract the gene name.
            rec.name = rec.description.split("|")[2].split(":")[1]
            descr = rec.description.split("|")[-1]
            # Extract all the features as a dictionary (apart from the last one which was the description.).
            feat = dict(item.split(":") for item in rec.description.split("|")[1:-1])
            rec.features = feat  # TODO for the moment all the features are stored as a dictionary and not as proper SeqFeature objects. (perhaps we can stick with that and there is no need to change it.)
            rec.description = descr
            recs.append(rec)

        recs = []
        for rec in self.gen:
            if self.exprTranscripts is None:
                _construct_bio_seq(recs, rec)
            if self.exprTranscripts and (rec.id in self.exprTranscripts):
                _construct_bio_seq(recs, rec)
        return recs

    def parse_expressed_transcripts(self, exprTrans):
        """Parses the given transcript expression file and returns a set of the expressed transcripts for quicker lookup."""
        trDf = pd.read_table(exprTrans, header=None)  # We do not really care about headers we only need the first column.
        return set(trDf[0])


class FeaturesExtract(object):
    """Claas to extract features."""

    def __init__(self, bioSeqRecs, utr3len, utrFiles):
        """Initialise with a list of SeqIO records."""
        self.bioSeqRecs = bioSeqRecs
        # Temporary files of the UTRs
        self.tf5p = tempfile.NamedTemporaryFile(mode="a", delete=False)
        self.tf3p = tempfile.NamedTemporaryFile(mode="a", delete=False)
        self.utr3len = utr3len
        self.utrFiles = utrFiles

    def collect_features(self):
        """Collect the features that do not need external computations."""
        # Initialise a pabdas data frame.
        pdf = pd.DataFrame(columns=['ensembl_gene_id', 'gene_name', 'coding_len', '5pUTR_len', '5pUTR_GC', '3pUTR_len', '3pUTR_GC', 'Kozak_Sequence', 'Kozak_Context'])
        for rec in self.bioSeqRecs:
            # Fetch UTRs, lenghts and GCs
            res3 = get_3utr(rec, self.tf3p, self.utr3len)
            if res3:  # Conditon for 3pUTR length.
                utr3len, utr3gc = res3
                utr5len, utr5gc = get_5utr(rec, self.tf5p)
                # Get Kozaks
                seqKozak, contKozak = get_kozak(rec)
                # Length of coding region.
                codeLen = int(rec.features["cDNA_end"]) - int(rec.features["cDNA_start"])
                # Add to pandas data frame.
                pdfe = [rec.features["GeneID"], rec.name, codeLen, utr5len, "{0:.2f}".format(utr5gc), utr3len, "{0:.2f}".format(utr3gc), seqKozak, contKozak]
                pdf.loc[rec.id] = pdfe
            else:
                continue
        self.tf5p.close()
        self.tf3p.close()
        return pdf

    def calculate_features(self):
        """Method to perfom all feature calculation.

        Basically invokes external software to make calculations."""
        # Calculates 5pUTR folding etc.
        fe5p = calculate_free_energy(self.tf5p.name, "5pUTR")
        fe3p = calculate_free_energy(self.tf3p.name, "3pUTR")
        # motifs = self.predictBinding()
        # Merge data frames and return.
        return pd.concat([fe5p, fe3p], axis=1, sort=False)

    def __del__(self):
        """Cleanup the temp files"""
        if self.utrFiles:
            shutil.move(self.tf5p.name, self.utrFiles[0])
            shutil.move(self.tf3p.name, self.utrFiles[1])
        else:
            os.remove(self.tf5p.name)
            os.remove(self.tf3p.name)


# FUNCTIONS
def get_kozak(rec, s=10, c=20):
    """Extract both Kozak sequence and context from a SeqIO record.
    (s and c define the extremeties of a Koxak sequence and are chosen by convention.)
    return a tuple of Kozak sequence and Kozak context.

    If ATG is located near the 5'UTR start, it is most likely that either seqs will not be retrieved as seld.base[-5:10] returns blank.

    In this case, we test whether the value left to ':' is negative or not.
    If it is < 0, we simply take the seq from 0 as : self.bases[0:self.cDNA_start) + 2) + s]"""
    feat = rec.features
    seq = rec.seq
    # Kozak sequence
    if int(feat["cDNA_start"]) - 1 - s < 0:
        kozakSeq = seq[0:(int(feat["cDNA_start"]) + 2) + s]
    else:
        kozakSeq = seq[(int(feat["cDNA_start"]) - 1 - s):(int(feat["cDNA_start"]) + 2) + s]
    # Kozak context
    if int(feat["cDNA_start"]) - 1 - s - c < 0:
        kozakContext = seq[0:(int(feat["cDNA_start"]) + 2) + s + c]
    else:
        kozakContext = seq[((int(feat["cDNA_start"]) - 1 - s) - c):(int(feat["cDNA_start"]) + 2) + s + c]
    return(str(kozakSeq), str(kozakContext))


def get_3utr(rec, tf3p, lim):
    """Fetch 3PUTR sequence, append it to fasta file and return its length and GC content.

    Impose limits to UTRs to upto <lim> nucleotides."""
    utr3p = rec.seq[int(rec.features["cDNA_end"]):]
    utr3plen = len(utr3p)
    if utr3plen >= lim:
        return None
    if utr3plen:
        tf3p.write(">{}_3PUTR\n{}\n".format(rec.id, utr3p))
    else:
        tf3p.write(">{}_3PUTR\n{}\n".format(rec.id, "N"))
    return(utr3plen, GC(utr3p))


def get_5utr(rec, tf5p):
    """Fetch 5PUTR sequence, append it to fasta file and return its length and GC content."""
    utr5p = rec.seq[0:int(rec.features["cDNA_start"])-1]
    if utr5p:
        tf5p.write(">{}_5PUTR\n{}\n".format(rec.id, utr5p))
    else:
        tf5p.write(">{}_5PUTR\n{}\n".format(rec.id, "N"))
    return(len(utr5p), GC(utr5p))


def calculate_free_energy(ffile, col):
    """Method to perform the free energy calculation by RNAfold and parsing of the results.

    Needs the RNA Vienna package to be installed."""
    pdf = pd.DataFrame(columns=['{}_MFE'.format(col), '{}_MfeBP'.format(col)])
    cmd = 'RNAfold --verbose --noPS --jobs -i {}'.format(ffile)
    proc = subprocess.check_output(shlex.split(cmd))
    proc = proc.decode()  # Convert it to proper string.
    outList = proc.splitlines()  # Retrieve STDOUT.
    for i in range(0, len(outList), 3):
        idt = outList[i][1:-6]  # To exclude the _UTR suffix.
        seq = outList[i+1]
        fold = outList[i+2]
        mfeRE = re.compile(r"[-+]?\d*\.\d+|\d+")
        mfe = float(mfeRE.search(fold).group())
        mfeBp = mfe / float(len(seq))
        pdf.loc[idt] = [mfe, mfeBp]
    return pdf


def predict_binding(ffile, motifs):
    """Method to run a FIMO search and parse and collect the results.

    Requires the intallation of the MEME suite."""
    # TODO implement the method!!!
    subprocess.call('fimo --verbosity 1 ' + motifs + ' ' + ffile.name, shell=True)
    fimo_tab = pd.read_csv("fimo_out/fimo.tsv", sep="\t")
    fimo_tab = fimo_tab.reset_index(drop=True)
    return None


def get_gene_ids(listID, out_fasta, dataset):
    """Function to connect to ENSEBL and retrieve data."""
    print("Connection to server.", file=sys.stderr)
    server = BiomartServer("http://www.ensembl.org/biomart/")
    dt = server.datasets[dataset]
    print("Connexion to ENSEMBL dataset...", file=sys.stderr)
    listAttrib = ['ensembl_gene_id', 'ensembl_transcript_id',
                  'external_gene_name', 'transcript_start', 'transcript_end',
                  '5_utr_end', '5_utr_start', '3_utr_end', '3_utr_start',
                  'transcription_start_site', 'transcript_biotype',
                  'cdna_coding_start', 'cdna_coding_end', 'cdna', 'description']
    print("Fetching data...", file=sys.stderr)
    # Collect data from the ENSEMBL datasets.
    dataf = pd.DataFrame()
    for chunk in chunks(listID, 100):
        res = dt.search({'filters': {'ensembl_gene_id': chunk}, 'attributes': listAttrib}, header=1)
        # Reading stream to a pandas data frame.
        dfTMP = pd.read_table(io.StringIO(res.text), sep='\t', encoding='utf-8')
        # Concatenate data frames.
        dataf = pd.concat([dataf, dfTMP], axis=0, sort=False)
        print('Fetching...', file=sys.stderr)
    print("...fetch done!", file=sys.stderr)
    # Cleanup data frame lines that do not correspond to protein coding genes.
    cdna = dataf[dataf['Transcript type'] == 'protein_coding']
    # remove NA from DF and reset index
    cdna.dropna(inplace=True)
    # FIXME TODO CAREFULL with the commenting of the following line as we need to keep that and remove the above if we do a transcription selection together with the geneID.
    #cdna = cdna.dropna(subset=['cDNA coding start', 'cDNA coding end'])
    cdna.reset_index(drop=True, inplace=True)
    # select only longest 5'UTR and 3'UTR
    # FIXME also this must be integrated with the fixme before.
    for i in range(cdna.shape[0]):
        ligne = pd.DataFrame(cdna.loc[i, :]).transpose()
        # For 5_UTR_start
        indices = get_utr5MAX(ligne)
        if len(cdna.iloc[i:i+1, :]["5' UTR start"].values[0].split(";")) > 1:
            cdna.iloc[i:i+1, 7:8] = indices[1]
            cdna.iloc[i:i+1, 8:9] = indices[0]
        if len(cdna.iloc[i:i+1, :]["3' UTR start"].values[0].split(";")) > 1:
            cdna.iloc[i:i+1, 9:10] = indices[1]
            cdna.iloc[i:i+1, 10:11] = indices[0]
    # select longest cDNA
    for index, row in cdna.iterrows():
        # If there is multiple annotation for cDNA start/end we take the longest possible.!!!
        min_cdna_start = min(map(int, row['cDNA coding start'].split(";")))
        max_cdna_end = max(map(int, row['cDNA coding end'].split(";")))
        # Put them back to the original data frame.
        cdna.at[index, 'cDNA coding start'] = min_cdna_start
        cdna.at[index, 'cDNA coding end'] = max_cdna_end
    # Exportat to FASTA format
    txt2fasta(cdna, out_fasta)


def get_utr5MAX(cdna_feat_row):
    """Select the longest 5'UTR from a cdna_feat-row with multiples utrs."""
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
    """Select the longest 3'UTR from a cdna_feat-row with multiples utrs."""
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


def get_cDNAstartMIN_OBSOLETE(cdna_feat_row):
    """Select the minimum cdna_start feature from multiples cDNA starts."""
    cDNA_start = cdna_feat_row["cDNA coding start"].values[0].split(";")
    min_cDNA_start = min(map(int, cDNA_start))
    return min_cDNA_start


def get_cDNAendMAX_OBSOLETE(cdna_feat_row):
    """Select the maximum cdna_end feature from multiples cDNA ends."""
    cDNA_end = cdna_feat_row["cDNA coding end"].values[0].split(";")
    max_cDNA_end = max(map(int, cDNA_end))
    return max_cDNA_end


def txt2fasta(cdna_feat_table, fastaOut):
    """Write the ENSEMBL features to a fasta file with the appropriate first fasta line."""
    with fastaOut as ff:
        for i, r in cdna_feat_table.iterrows():
            ff.write(">{} |GeneID:{}|GeneName:{}|5pUTR_start:{}|5pUTR_end:{}|3pUTR_start:{}|3pUTR_end:{}|cDNA_start:{}|cDNA_end:{}|{}\n".format(r["Transcript stable ID"], r["Gene stable ID"], r["Gene name"], r["5' UTR start"], r["5' UTR end"], r["3' UTR start"], r["3' UTR end"], r["cDNA coding start"], r["cDNA coding end"], r["Gene description"]))
            ff.write("{}\n".format(r["cDNA sequences"]))


def chunks(l, n):
    """Yield successive n-sized chunks from a list l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]
