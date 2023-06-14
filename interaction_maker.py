from ordered_set import OrderedSet
import requests
import json
import pandas as pd

import pubchempy as pcp
import math
import numpy as np
import matplotlib.pyplot as plt
from drug_class import DrugClass
from kronrls_mkl.kronrls_mkl_fun import kronrls_mkl
from kernels.lambda_kernel import calculate_lambdak_similarity
from sequence_finder import get_protein_sequence
from kernels.simcomp import calculate_simcomp_similarity
from similar_drugs import findSimilarDrugs

from kernels.smith_waterman_score import smith_waterman_score
from kernels.spectrum_kernel import compute_spectrum_matrix
TOKEN = '43303cb2339290ef71c054f1371f3e6136b86a91'
# disease - gene association
# res=requests.get('https://www.disgenet.org/api/gda/source/ALL?limit=10',headers = {'Authorization': 'Bearer 43303cb2339290ef71c054f1371f3e6136b86a91'})
# #print(res.json())
# pd.set_option('display.max_columns', None)
# pands=pd.DataFrame(res.json())
# print(pands[['geneid','disease_name']])

##From here uncomment CODENEW


def findaminoseqs(genedatas):
    details = genedatas.apply(lambda x : True
            if x['score']>0.1 else False, axis = 1)
    geneLength=len(details[details == True].index)
    swkernel=[[None for i in range(geneLength)] for j in range(geneLength)]
    #Finding aminoacid sequences of each gene
    aminoseqs=[]
    notableGenes=[]
    for(i,j) in genedatas.iterrows():
        if(j['score']>0.1):
            print("Finding sequence of ",j['gene_symbol'])
            aminoseqs.append(get_protein_sequence(j['gene_symbol']))
            notableGenes.append(j)
    return aminoseqs,notableGenes


    # ind1,ind2=0,0
    # for(i,j) in genedatas.iterrows():
    #     ind2=0
    #     if(j['score']>0.1):
    #         for(k,l) in genedatas.iterrows():
    #             if(l['score']>0.1):
    #                 swkernel[ind1][ind2]=smith_waterman_score(aminoseqs[ind1],aminoseqs[ind2])
    #                 print("Calculating SWKernel for ",ind1,j['gene_symbol']," and ",ind2,l['gene_symbol']," : ",swkernel[ind1][ind2])
    #                 # print(j['geneid'],j['gene_symbol'],j['score'])
    #                 ind2+=1
    #         ind1+=1
    #     if(ind1==2):
    #         break
    return swkernel

diseaseCode=input("Enter Input Disease Code:") #Eg: C0000737, C0000744, C0000768, C0004903, C0003873, Code:C0005426, C0007286
res=requests.get('https://www.disgenet.org/api/disease/'+diseaseCode,headers = {'Authorization': 'Bearer '+TOKEN})  
diseaseData=pd.DataFrame(res.json())
print('Disease: ',diseaseData['disease_name'][0],'  |  Disease Class: ',diseaseData['disease_class'][0])
res=requests.get('https://www.disgenet.org/api/gda/disease/'+diseaseCode,headers = {'Authorization': 'Bearer '+TOKEN})  
genedatas=pd.DataFrame(res.json())


drugKernels=[]
geneKernels=[]

#finding Gene Kernels

aminoseqs,notableGenes=findaminoseqs(genedatas)
geneKernels.append(smith_waterman_score(aminoseqs))
# geneKernels.append(compute_spectrum_matrix(aminoseqs,3))
# geneKernels.append(compute_spectrum_matrix(aminoseqs,4))

numGenes=len(aminoseqs)


drugsList =OrderedSet([])
interactions = []
# finding Known Interactions


ind1=0
for gene in notableGenes:
    drug_gene = requests.get('https://dgidb.org/api/v2/interactions?genes='+gene['gene_symbol'])
    drugs=drug_gene.json()['matchedTerms']
    # print('Gene '+drugs[0]['geneName'],' affected by Drugs: ',end='')
    ind2=0
    # drugsList.append([])
    for i in drugs[0]['interactions']:
        index=drugsList.add(DrugClass(i['drugName']))
        if(index<len(interactions)):
            interactions[index][ind1]=1
        else:
            interactions.append([0]*numGenes)
            interactions[index][ind1]=1
        # print(drugsList[ind1][ind2].drug_id)
        ind2+=1
        break
    ind1+=1

print(interactions)
print(" with the drugs")
drugsListN = findSimilarDrugs(drugsList)
interactions.extend([[0]*numGenes]*( len(drugsListN)-len(drugsList)))
cids = [drug.cid for drug in drugsListN]

# get the SMILES for the drugs
print(cids)
smiles = []
for cid in cids:
    if(cid==None):
        smiles.append(None) 
    else:
        smiles.append(pcp.Compound.from_cid(cid).isomeric_smiles)
drugKernels.append(calculate_simcomp_similarity(smiles))
drugKernels.append(calculate_lambdak_similarity(smiles))
print("Unique drugs: ",len(drugsListN))


# geneKernels=[[[1.0000000000000002, 0.19518001458970666, 0.19518001458970666, 0.16903085094570333, 0.16903085094570333, 0.19518001458970666], [0.07575540190785703, 1.0, 0.08469711416439278, 0.08469711416439278, 0.07575540190785703, 0.07575540190785703], [0.10992994198586253, 0.12290541152149842, 0.9999999999999998, 0.10992994198586253, 0.10283004486475836, 0.10992994198586253], [0.08512565307587487, 0.1098967455659645, 0.0982946374365981, 1.0000000000000002, 0.665760253376428, 0.0982946374365981], [0.08660254037844388, 0.1, 0.09354143466934853, 0.6773108592072034, 1.0, 0.1], [0.10511766624552735, 0.10511766624552735, 0.10511766624552735, 0.10511766624552735, 0.10511766624552735, 1.0000000000000002]]]
# drugKernels=[[[1.0, 0.0, 0.13636363636363635, 0.1171875, 0.8791208791208791, 0.17088607594936708, 0.11627906976744186], [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.13636363636363635, 0.0, 1.0, 0.12121212121212122, 0.1282051282051282, 0.9, 0.12030075187969924], [0.1171875, 0.0, 0.12121212121212122, 1.0, 0.11627906976744186, 0.15328467153284672, 0.9830508474576272], [0.8791208791208791, 0.0, 0.1282051282051282, 0.11627906976744186, 1.0, 0.1625, 0.11538461538461539], [0.17088607594936708, 0.0, 0.9, 0.15328467153284672, 0.1625, 1.0, 0.15217391304347827], [0.11627906976744186, 0.0, 0.12030075187969924, 0.9830508474576272, 0.11538461538461539, 0.15217391304347827, 1.0]]]
# interactions=[[1, 0, 0, 0, 0, 0], [0, 0, 1, 1, 0, 0], [0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]]

A=kronrls_mkl(np.array(drugKernels),np.array(geneKernels),np.array(interactions),1)

for drug in drugsListN:
    print(drug.name)

for i in range(len(notableGenes)):
    print(notableGenes[i]['gene_symbol'],": ")
    for val in A[i]:
        print(val,end=' ')
    print("\n")





#to here






# disease => gene 
# --------------------
# a) List of Drugs
# print("Fetching genes")
# pd.set_option('display.max_columns', None)
# pands=pd.DataFrame(res.json())

# print("Genes Count: ",len(pands))
# # print(pands['disease_name'])
# # b) accessing each disease=> gene
# for i in pands['symbol']:
#     print("For ",i)
#     dis_gene=requests.get('https://www.disgenet.org/api/gda/gene/'+str(i),headers = {'Authorization': 'Bearer 43303cb2339290ef71c054f1371f3e6136b86a91'})
#     # print(dis_gene.json())
#     print('Gene '+i+' associated with Diseases: ',len(dis_gene.json()),end='')
#     for j in dis_gene.json():
#             if(j['score']>0.3):
#                 print(j['diseaseid'],j['disease_name'],end=' ')
        
#     print("\n")
#         # print('Disease '+i+' associated with Genes: ',end='')

# # # Api for fetching Drug-Gene data
# # for gen in pands['gene_symbol']:
# #     drug_gene = requests.get('https://dgidb.org/api/v2/interactions?genes='+gen)
# #     drugs=drug_gene.json()['matchedTerms']
# #     print('Gene '+drugs[0]['geneName'],' affected by Drugs: ',end='')
# #     for i in drugs[0]['interactions']:
# #         print(i['drugName'],end=' ')
# #     print("\n")

# # Api for fetching Drug-Gene data

# # drug_gene = requests.get('https://dgidb.org/api/v2/interactions?count=10&page=1')
# # drugs=drug_gene.json()['records']
# # for i in drugs:
# #     print('Gene '+i['gene_name']+' affected by Drug '+i['drug_name'])
    
# # print("\n")
# print("\n\n\n")

# # drug => Gene
# # a) fetching Drug list
# # drugss=[]
# # druglist = requests.get('https://dgidb.org/api/v2/drugs?count=10&page=1')
# # # print(druglist.json())
# # drug_list=druglist.json()['records']
# # for i in drug_list:
# #     # print('Drug '+i['name'])
# #     drugss.append(i['name'])

# # # b)fetching Drug=> gene data

# # for j in drugss:
# #     drug_gene = requests.get('https://dgidb.org/api/v2/interactions?drugs='+j)
# #     # print(drug_gene.json())
# #     drugs=drug_gene.json()['matchedTerms']
# #     print('Drug '+drugs[0]['drugName'],' affects Gene: ',end='')
# #     for i in drugs[0]['interactions']:
# #         print(i['geneName'],end=' ')
# #     print("\n")

