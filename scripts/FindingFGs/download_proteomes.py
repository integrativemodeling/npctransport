import requests
import pandas as pd
import time
from get_proteome_list import get_proteome_list

INPUT_FOLDER='Input/'

def download_one_proteome(uni_prot_id, save_in='../Input/proteomes_archae/'):
    """
    Takes in a Proteome ID (unit_prot_id), downloads the proteome in
    fasta format, saves the file in save_in folder.
    """
    # Construct an URL to access
    url_start = 'https://www.uniprot.org/uniprot/?query=proteome:'
    url_end = '&compress=no&force=true&format=fasta'
    url = url_start + str(uni_prot_id) + url_end

    # Access the URL, retrieve the fasta, which comes as an attachment,
    # then saves in the save_in folder
    fasta = requests.get(url)
    file_name =  save_in + str(uni_prot_id) + '_proteome.fasta'
    with open(file_name, 'wb') as infile:
        infile.write(fasta.content)

# When run, download multiple proteomes
if __name__ == "__main__":
    # download archae proteomes
    proteome_list_archae = get_proteome_list(INPUT_FOLDER+'archae_list.tab')
    for uni_prot_id in proteome_list_archae[:30]: #[100:150]:
        print('Downloading {} ...'.format(uni_prot_id))
        download_one_proteome(uni_prot_id, save_in=INPUT_FOLDER+'proteomes_archae/')
        time.sleep(.1)   # courtesy for the server

    # download bacteria proteomes
    proteome_list_bacteria = get_proteome_list(INPUT_FOLDER+'bacteria_list.tab')
    for uni_prot_id in proteome_list_bacteria[:30]: # [100:150]:
        print('Downloading {} ...'.format(uni_prot_id))
        download_one_proteome(uni_prot_id, save_in=INPUT_FOLDER+'proteomes_bacteria/')
        time.sleep(.1)   # courtesy for the server
