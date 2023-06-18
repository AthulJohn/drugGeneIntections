# from ordered_set import OrderedSet
# import requests
# import json
# import pandas as pd

# import pubchempy as pcp
# import math
# import numpy as np
# import matplotlib.pyplot as plt
# from drug_class import DrugClass
# from kronrls_mkl.kronrls_mkl_fun import kronrls_mkl
# from kernels.lambda_kernel import calculate_lambdak_similarity
# from sequence_finder import get_protein_sequence
# from kernels.simcomp import calculate_simcomp_similarity
# from similar_drugs import findSimilarDrugs

# from kernels.smith_waterman_score import smith_waterman_score
# from kernels.spectrum_kernel import compute_spectrum_matrix
# TOKEN = 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'



# def findaminoseqs(genedatas):
#     details = genedatas.apply(lambda x : True
#             if x['score']>0.1 else False, axis = 1)
#     geneLength=len(details[details == True].index)
#     swkernel=[[None for i in range(geneLength)] for j in range(geneLength)]
#     #Finding aminoacid sequences of each gene
#     aminoseqs=[]
#     notableGenes=[]
#     for(i,j) in genedatas.iterrows():
#         if(j['score']>0.1):
#             print("Finding sequence of ",j['gene_symbol'])
#             aminoseqs.append(get_protein_sequence(j['gene_symbol']))
#             notableGenes.append(j)
#     return aminoseqs,notableGenes

# diseaseCode=input("Enter Input Disease Code:") #Eg: C0000737, C0000744, C0000768, C0004903, C0003873, Code:C0005426, C0007286
# res=requests.get('https://www.disgenet.org/api/disease/'+diseaseCode,headers = {'Authorization': 'Bearer '+TOKEN})  
# diseaseData=pd.DataFrame(res.json())
# res=requests.get('https://www.disgenet.org/api/gda/disease/'+diseaseCode,headers = {'Authorization': 'Bearer '+TOKEN})  
# genedatas=pd.DataFrame(res.json())

# drugKernels=[]
# geneKernels=[]
# #finding Gene Kernels
# aminoseqs,notableGenes=findaminoseqs(genedatas)
# geneKernels.append(smith_waterman_score(aminoseqs))
# geneKernels.append(compute_spectrum_matrix(aminoseqs,3))
# geneKernels.append(compute_spectrum_matrix(aminoseqs,4))
# numGenes=len(aminoseqs)

# # finding Known Interactions
# drugsList =OrderedSet([])
# interactions = []
# ind1=0
# for gene in notableGenes:
#     drug_gene = requests.get('https://dgidb.org/api/v2/interactions?genes='+gene['gene_symbol'])
#     drugs=drug_gene.json()['matchedTerms']
#     ind2=0
#     for i in drugs[0]['interactions']:
#         index=drugsList.add(DrugClass(i['drugName']))
#         if(index<len(interactions)):
#             interactions[index][ind1]=1
#         else:
#             interactions.append([0]*numGenes)
#             interactions[index][ind1]=1
#         ind2+=1
#         break
#     ind1+=1
# drugsListN = findSimilarDrugs(drugsList)
# interactions.extend([[0]*numGenes]*( len(drugsListN)-len(drugsList)))
# cids = [drug.cid for drug in drugsListN]

# # get the SMILES for the drugs
# smiles = []
# for cid in cids:
#     if(cid==None):
#         smiles.append(None) 
#     else:
#         smiles.append(pcp.Compound.from_cid(cid).isomeric_smiles)
# drugKernels.append(calculate_simcomp_similarity(smiles))
# drugKernels.append(calculate_lambdak_similarity(smiles))
# print("Unique drugs: ",len(drugsListN))

# A=kronrls_mkl(np.array(drugKernels),np.array(geneKernels),np.array(interactions),1)
# for drug in drugsListN:
#     print(drug.name)
# for i in range(len(notableGenes)):
#     print(notableGenes[i]['gene_symbol'],": ")
#     for val in A[i]:
#         print(val,end=' ')
#     print("\n")





                CDKN2A ( 0.3 ) ELF3 ( 0.3 ) ERBB2 ( 0.3 ) IDH1 ( 0.3 ) IDH2 ( 0.3 ) PRKACB ( 0.3 ) 
RIGOSERTIB :  0.27655251124838637 -0.07003171111108354 0.08202873370926878 -0.008278800232630905 0.025103543222849304 -0.03506861215409118 SCORE:  0.14911917350349302


ERTUMAXOMAB :  0.09560856811569955 0.1949351098733744 1.0037970893873422 0.30463656149918195 0.09135416698996852 0.09698762089851129 SCORE:  0.5361957350292234


ENASIDENIB :  0.029849801128002604 -0.016627458480976125 0.04351223999860901 0.02156695361975296 0.2173398140977029 -0.0007036878882332393 SCORE:  0.09887998656398303


FASUDIL :  -0.032402763111455846 -0.037463036344157874 0.06560194003416597 -0.009915297554884335 0.002987772882353528 0.27694931210028084 SCORE:  0.12759603660818952


Rigosertib sodium :  0.17801357188310413 -0.04943838810254842 0.03766946865955807 -0.013778713457584324 -0.0007740853199459903 -0.022488543571885462 SCORE:  0.09064883129838792      


Enasidenib mesylate :  0.00021770882279466101 -0.01201274092402998 0.0032989201088918048 -0.008241066677596987 0.13082777943763524 0.008549432929748902 SCORE:  0.04894429467020927   


Fasudil hydrochloride :  -0.030242298858326275 -0.03588694962098199 0.05590565906050173 -0.010444160686190095 0.002283105799127035 0.2532840591661265 SCORE:  0.1164138699573761  