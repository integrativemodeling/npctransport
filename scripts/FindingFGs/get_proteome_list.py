import pandas as pd

def get_proteome_list(list_file):
    """
    Takes a .tab file downloaded manually from Uniprot.
    It should contain a Proteome ID as its first column.

    This function just reads the .tab file and returns the
    Proteome ID column
    """
    df = pd.read_csv(list_file, sep='\t')
    proteome_list = df['Proteome ID']
    return proteome_list


if __name__ == "__main__":
    archae_file = 'archae_list.tab'
    bacteria_file = 'bacteria_list.tab'

    proteome_list_archae = get_proteome_list(archae_file)
    proteome_list_bacteria = get_proteome_list(bacteria_file)

    print('===== Archae =====')
    print(proteome_list_archae[:10])
    print('===== Bacteria =====')
    print(proteome_list_bacteria[:10])