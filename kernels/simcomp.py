

# import pubchempy as pcp

# def simcomp_kernel(cid1, cid2, threshold=0.8):
#     sids1 = pcp.get_sids(cid1)
#     sids2 = pcp.get_sids(cid2)
#     n_common = 0
#     print(sids1,sids2)
#     for sid1 in sids1:
#         for sid2 in sids2:
#             sim = pcp.similarity(sid1, sid2)
#             print(sim)
#             if sim >= threshold:
#                 n_common += 1
#     return n_common

# print(simcomp_kernel(11342,14564))

# import pubchempy as pcp
# from rdkit import Chem
# from rdkit.Chem import AllChem
# from rdkit.Chem.Fingerprints import FingerprintMols
# from rdkit.SimDivFilters import rdSimDivPickers

# # define a list of drug names
# drug_names = ['aspirin', 'ibuprofen', 'acetaminophen']

# # get the PubChem CIDs for the drugs
# cids = [pcp.get_compounds(name, 'name')[0].cid for name in drug_names]

# # get the SMILES for the drugs
# print(cids)
# smiles = [pcp.Compound.from_cid(cid).isomeric_smiles for cid in cids]

# # # convert the SMILES to RDKit moleculesfrom rdkit import Ch

# smiles_list = ['CCO', 'CCC', 'CCCN', 'CCNC']
# mols = [Chem.MolFromSmiles(smi) for smi in smiles_list]
# fps = [AllChem.GetMorganFingerprint(mol, 2) for mol in mols]
# similarity_matrix = [[AllChem.DataStructs.TanimotoSimilarity(fp1, fp2) for fp2 in fps] for fp1 in fps]

# print(similarity_matrix)


from rdkit import Chem
import pubchempy as pcp
from rdkit.Chem import AllChem
from rdkit import DataStructs

def calculate_simcomp_similarity(drugs):
    similarity_matrix = []
    for i in range(len(drugs)):
        similarity_row = []
        for j in range(len(drugs)):
            if i == j:
                similarity_row.append(1.0)  # Similarity of a drug to itself is 1.0
            else:
                if(drugs[i]==None or drugs[j]==None):
                    similarity_row.append(0.0)
                    continue
                mol1 = Chem.MolFromSmiles(drugs[i])
                mol2 = Chem.MolFromSmiles(drugs[j])
                fp1 = AllChem.GetMorganFingerprint(mol1, 2)  # Morgan fingerprint with radius 2
                fp2 = AllChem.GetMorganFingerprint(mol2, 2)
                similarity = DataStructs.TanimotoSimilarity(fp1, fp2)  # Tanimoto coefficient
                similarity_row.append(similarity)
        similarity_matrix.append(similarity_row)
    return similarity_matrix

