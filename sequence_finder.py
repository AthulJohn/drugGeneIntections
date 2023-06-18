from Bio import SeqIO
from Bio import Entrez

def get_protein_sequence(protein_symbol):

    # Set up the Entrez API
    Entrez.email = "johnpc716@gmail.com"  

    # Search for the protein symbol in UniProt using Entrez
    handle = Entrez.esearch(db="protein", term=protein_symbol)
    record = Entrez.read(handle)
    handle.close()

    # Get the protein accession number from the search results
    protein_accession = record['IdList'][0]

    # Fetch the protein record from UniProt using the protein accession number
    handle = Entrez.efetch(db="protein", id=protein_accession, rettype="fasta", retmode="text")
    protein_record = SeqIO.read(handle, "fasta")
    handle.close()

    # Extract the protein sequence from the protein record
    protein_sequence = protein_record.seq
    return protein_sequence