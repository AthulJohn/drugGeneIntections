import numpy as np
from rdkit import Chem
import pubchempy as pcp
from rdkit.Chem import GraphDescriptors
from scipy.spatial.distance import pdist, squareform

def calculate_lambdak_similarity(drugs):
    num_drugs = len(drugs)
    similarity_matrix = np.zeros((num_drugs, num_drugs))
    
    for i in range(num_drugs):
        for j in range(i, num_drugs):
            if(drugs[i]==None or drugs[j]==None):
                similarity_matrix[i,j]=0.0
                similarity_matrix[j,i]=0.0
                continue
            mol1 = Chem.MolFromSmiles(drugs[i])
            mol2 = Chem.MolFromSmiles(drugs[j])
            
            # Calculate the graph-level descriptor (lambda)
            lambda1 = GraphDescriptors.Chi0n(mol1)
            lambda2 = GraphDescriptors.Chi0n(mol2)
            
            # Calculate the λ-Kernel similarity
            similarity = np.exp(-((lambda1 - lambda2) ** 2))
            
            similarity_matrix[i, j] = similarity
            similarity_matrix[j, i] = similarity
    
    return similarity_matrix.tolist()

