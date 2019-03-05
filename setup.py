#!/usr/bin/python3

from distutils.core import setup

config = {
  'name': 'rnaFeaturesLib',
  'version': rnaFeaturesLib.__version__,
  'author': 'Costas Bouyioukos',
  'author_email': 'costas.bouyioukos@univ-paris-diderot.fr',
  'url': 'github.com/parisepigenetics/rna_feat_ext',
  'description': 'A software package for the Computer Aided Design of NEW To Nature chromosomes.',
  'long_description': open("README.md").read(),
  'download_url': 'https://github.com/parisepigenetics/rna_feat_ext.git',
  'py_modules': ['rnaFeaturesLib'],
  'scripts': ['bin/fasta2table.py', 'bin/geneIDs2fasta.py'],
  'requires': ['biomart', 'biopython', 'pandas', 'prettytable'],
  #'data_files': [('data', ['testRNAfeatExt_IDs.txt'])],
  'license': 'GPL v3.0 or later',
  'classifiers': ['Programming Language :: Python', 'Topic :: Science :: Computational Biology'],
}

setup(**config)
