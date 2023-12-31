from ordered_set import OrderedSet
import requests
import json
import pandas as pd

import pubchempy as pcp
import math
import numpy as np
import matplotlib.pyplot as plt
from drug_class import DrugClass
from kernels.smith_waterman_score import smith_waterman_score
from kernels.spectrum_kernel import compute_spectrum_matrix
from kronrls_mkl.kronrls_mkl_fun import kronrls_mkl
from kernels.lambda_kernel import calculate_lambdak_similarity
from sequence_finder import get_protein_sequence
from kernels.simcomp import calculate_simcomp_similarity
from similar_drugs import findSimilarDrugs

TOKEN = '43303cb2339290ef71c054f1371f3e6136b86a91'


def findaminoseqs(genedatas):
    aminoseqs=[]
    notableGenes=[]
    totallength=0
    for(i,j) in genedatas.iterrows():
        totallength+=1
        if(j['score']>0.1):
            aminoseqs.append(get_protein_sequence(j['gene_symbol']))
            notableGenes.append(j)
    
    print("\nFound ",totallength,"  total genes and ",len(notableGenes)," notable genes")
    print("Notable Genes: ",end=' ')
    for gene in notableGenes:
        print(gene['gene_symbol'],"(",gene['score'],")",end=', ')
    return aminoseqs,notableGenes


diseaseCode=input("Enter Input Disease Code (Eg: C0005426):") #Eg: C0000737, C0000744, C0000768, C0004903, C0003873, Code:C0005426, C0007286
if(diseaseCode==''):
    diseaseCode='C0005426'# Default Disease: Biliary Tract Neoplasm
res=requests.get('https://www.disgenet.org/api/disease/'+diseaseCode,headers = {'Authorization': 'Bearer '+TOKEN})  
diseaseData=pd.DataFrame(res.json())
print('Input Disease: ',diseaseData['disease_name'][0],'      |       Disease Code: ',diseaseCode)
res=requests.get('https://www.disgenet.org/api/gda/disease/'+diseaseCode,headers = {'Authorization': 'Bearer '+TOKEN})  
genedatas=pd.DataFrame(res.json())
drugKernels=[]
geneKernels=[]

#finding Gene Kernels

aminoseqs,notableGenes=findaminoseqs(genedatas)
geneKernels.append(smith_waterman_score(aminoseqs))
geneKernels.append(compute_spectrum_matrix(aminoseqs,3))
geneKernels.append(compute_spectrum_matrix(aminoseqs,4))
print("\nCalculated Smith Waterman, Spectrum Kernels for Genes\n")
numGenes=len(aminoseqs)



drugsList =OrderedSet([])
interactions = []

# finding Known Interactions
ind1=0
for gene in notableGenes:
    drug_gene = requests.get('https://dgidb.org/api/v2/interactions?genes='+gene['gene_symbol'])
    drugs=drug_gene.json()['matchedTerms']
    ind2=0
    for i in drugs[0]['interactions']:
        index=drugsList.add(DrugClass(i['drugName']))
        if(index<len(interactions)):
            interactions[index][ind1]=1
        else:
            interactions.append([0]*numGenes)
            interactions[index][ind1]=1
        ind2+=1
        break
    ind1+=1
print("\nFound ",len(drugsList)," drugs from ",len(interactions)," interactions")
drugsListN = findSimilarDrugs(drugsList)
print("Found ",len(drugsListN)-len(drugsList)," similar drugs")
print("Thus total number of drugs: ",len(drugsListN))
print("Drugs: ",end=' ')
for drug in drugsListN:
    print(drug.name,end=', ')
interactions.extend([[0]*numGenes]*( len(drugsListN)-len(drugsList)))
cids = [drug.cid for drug in drugsListN]

# get the SMILES for the drugs
smiles = []
for cid in cids:
    if(cid==None):
        smiles.append(None) 
    else:
        smiles.append(pcp.Compound.from_cid(cid).isomeric_smiles)
drugKernels.append(calculate_simcomp_similarity(smiles))
drugKernels.append(calculate_lambdak_similarity(smiles))
print("\nCalculated Simcomp, Lambda Kernels for Drugs\n")

# For Testing purposes, use the below values as kernels. 
# geneKernels=[[[1.0000000000000002, 0.19518001458970666, 0.19518001458970666, 0.16903085094570333, 0.16903085094570333, 0.19518001458970666], [0.07575540190785703, 1.0, 0.08469711416439278, 0.08469711416439278, 0.07575540190785703, 0.07575540190785703], [0.10992994198586253, 0.12290541152149842, 0.9999999999999998, 0.10992994198586253, 0.10283004486475836, 0.10992994198586253], [0.08512565307587487, 0.1098967455659645, 0.0982946374365981, 1.0000000000000002, 0.665760253376428, 0.0982946374365981], [0.08660254037844388, 0.1, 0.09354143466934853, 0.6773108592072034, 1.0, 0.1], [0.10511766624552735, 0.10511766624552735, 0.10511766624552735, 0.10511766624552735, 0.10511766624552735, 1.0000000000000002]],[[1.0000000000000002, 0.19518001458970666, 0.19518001458970666, 0.16903085094570333, 0.16903085094570333, 0.19518001458970666], [0.07575540190785703, 1.0, 0.08469711416439278, 0.08469711416439278, 0.07575540190785703, 0.07575540190785703], [0.10992994198586253, 0.12290541152149842, 0.9999999999999998, 0.10992994198586253, 0.10283004486475836, 0.10992994198586253], [0.08512565307587487, 0.1098967455659645, 0.0982946374365981, 1.0000000000000002, 0.665760253376428, 0.0982946374365981], [0.08660254037844388, 0.1, 0.09354143466934853, 0.6773108592072034, 1.0, 0.1], [0.10511766624552735, 0.10511766624552735, 0.10511766624552735, 0.10511766624552735, 0.10511766624552735, 1.0000000000000002]]]
# drugKernels=[[[1.0, 0.0, 0.13636363636363635, 0.1171875, 0.8791208791208791, 0.17088607594936708, 0.11627906976744186], [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.13636363636363635, 0.0, 1.0, 0.12121212121212122, 0.1282051282051282, 0.9, 0.12030075187969924], [0.1171875, 0.0, 0.12121212121212122, 1.0, 0.11627906976744186, 0.15328467153284672, 0.9830508474576272], [0.8791208791208791, 0.0, 0.1282051282051282, 0.11627906976744186, 1.0, 0.1625, 0.11538461538461539], [0.17088607594936708, 0.0, 0.9, 0.15328467153284672, 0.1625, 1.0, 0.15217391304347827], [0.11627906976744186, 0.0, 0.12030075187969924, 0.9830508474576272, 0.11538461538461539, 0.15217391304347827, 1.0]],[[1.0, 0.0, 0.13636363636363635, 0.1171875, 0.8791208791208791, 0.17088607594936708, 0.11627906976744186], [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.13636363636363635, 0.0, 1.0, 0.12121212121212122, 0.1282051282051282, 0.9, 0.12030075187969924], [0.1171875, 0.0, 0.12121212121212122, 1.0, 0.11627906976744186, 0.15328467153284672, 0.9830508474576272], [0.8791208791208791, 0.0, 0.1282051282051282, 0.11627906976744186, 1.0, 0.1625, 0.11538461538461539], [0.17088607594936708, 0.0, 0.9, 0.15328467153284672, 0.1625, 1.0, 0.15217391304347827], [0.11627906976744186, 0.0, 0.12030075187969924, 0.9830508474576272, 0.11538461538461539, 0.15217391304347827, 1.0]]]
# interactions=[[1, 0, 0, 0, 0, 0], [0, 0, 1, 1, 0, 0], [0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]]
print("\nPredicting new interactions using KronRLS-MKL...")
A=kronrls_mkl(np.array(drugKernels),np.array(geneKernels),np.array(interactions),1)


for gene in notableGenes:
    print(gene['gene_symbol'],"(",gene['score'],")",end=' ')

for i in range(len(drugsListN)):
    print(drugsListN[i].name,": ",end=' ')
    drug_score=0
    for j in range(len(A[i])):
        print(A[i][j],end=' ')
        drug_score+=abs(A[i][j]*notableGenes[j]['score'])
    print("SCORE: ",drug_score)




