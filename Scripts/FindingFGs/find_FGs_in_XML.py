#!/usr/bin/python
from optparse import OptionParser
import xml.etree.ElementTree as etree  
import re
import os
import util

def get_options():
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                  help="input xml file", metavar="FILE")
    (options, args) = parser.parse_args()
    return options


#  MAIN
options = get_options()
WINDOW_SIZE=40
#options.filename = "uniprot-organism%3A-caenorhabditis+elegans-.xml"
print("Filename: ", options.filename)
context = etree.iterparse(options.filename, events = ('start','end') )
root = None
for event, elem in context:
    # get hold of root
    if event == "start" and root is None:
        root = elem 
    # extract sequence of each entry and print fg sequences info
    elif elem.tag == "{http://uniprot.org/uniprot}entry" and event=="end":
            name = elem.findtext("{http://uniprot.org/uniprot}name")
            # find full protein name
            full_prot_name = "unknown"
            protein = elem.find("{http://uniprot.org/uniprot}protein") 
            if protein is not None: 
                for prot_elem in protein.iter():
                    if prot_elem.tag == "{http://uniprot.org/uniprot}fullName":
                        full_prot_name = prot_elem.text
            prim_gene_name = "unknown"
            gene = elem.find("{http://uniprot.org/uniprot}gene")
            if gene is not None: 
                for gene_elem in gene.iter():
                    if gene_elem.tag == "{http://uniprot.org/uniprot}name" and gene_elem.get("type") == "primary":
                        prim_gene_name = gene_elem.text
            sequence = elem.findtext("{http://uniprot.org/uniprot}sequence")
            sequence = util.remove_whitespace( sequence )
            if util.is_fg(sequence) \
             and not re.search("nup|nuclear pore complex|nucleoporin", full_prot_name, re.IGNORECASE):
                sequence_emph = re.sub("FG", "[***FG***]", sequence)
#                print(name, full_prot_name, prim_gene_name)
                print(name, full_prot_name, prim_gene_name, sequence_emph)
                # Print disordered scores for FGs
                fi_list= util.get_fold_index(sequence, w=WINDOW_SIZE)
                seq_id=1
                for aa,aa_plus_1, fi in zip(sequence[:],
                                            sequence[1:], 
                                            fi_list[int(WINDOW_SIZE/2):]):
                    if(aa=='F' and aa_plus_1=='G'):
                        print("FG {:4d} fold-index {:+.2f}".format(seq_id, fi))
                    seq_id= seq_id+1
                # Print GO annotations
                dbRefs = elem.findall("{http://uniprot.org/uniprot}dbReference")
                for dbRef in dbRefs:
                    if(dbRef.get("type") == "GO"):
                        for property in dbRef.findall("{http://uniprot.org/uniprot}property"):
                            if(property.get("type") == "term"):
                                print("GO: ", name, "[" + full_prot_name + " ; " + prim_gene_name + "]", property.get("value"))
            root.clear()
#    else:
#        print elem.tag, len(elem), event


