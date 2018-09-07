#!/usr/bin/python
from optparse import OptionParser
#import xml.etree.ElementTree as etree  
import re
import os
from Bio import SeqIO
import sys
import util

def get_options():
    parser = OptionParser()
    parser.add_option("-f", "--fasta_file", dest="filename",
                  help="input fasta file of protein sequences", metavar="FILE")
    (options, args) = parser.parse_args()
    return options


#  MAIN
options = get_options()
WINDOW_SIZE=40
#options.filename = "uniprot-organism%3A-caenorhabditis+elegans-.xml"
print("Filename: ", options.filename)
with open(options.filename, "r") as filehandle:
    for record in SeqIO.parse(filehandle, "fasta"):
        seq =str(record.upper().seq)
        seq = util.remove_whitespace( seq )
        if util.is_fg(seq):
            seq_emph = re.sub("FG", "[***FG***]", seq)
            print(record.id, seq_emph)
            # Print disordered scores for FGs
            fi_list= util.get_fold_index(seq, w=WINDOW_SIZE)
            seq_id=1
            for aa,aa_plus_1, fi in zip(seq[:], 
                                        seq[1:], 
                                        fi_list[int(WINDOW_SIZE/2):]):
                if(aa=='F' and aa_plus_1=='G'):
                    print("FG {:4d} fold-index {:+.2f}".format(seq_id, fi))
                seq_id= seq_id+1
