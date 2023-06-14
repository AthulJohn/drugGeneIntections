# import requests

# def find_similar_drugs(drug_name, threshold=90):
#     # Search for reference drug by name
#     url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{drug_name}/cids/JSON"
#     response = requests.get(url)
#     if response.status_code != 200:
#         return None
#     cids = response.json()["IdentifierList"]["CID"]
#     if not cids:
#         return None
    
#     # Get substructures of reference drug
#     # url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cids[0]}/substructure/similarity/cids/JSON?Threshold={threshold}"
#     url="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/c00035/similarity/JSON"
#     payload = {
#     "identity_type": "name",
#     "identity": "acetaminophen",
#     "threshold": 90,
#     "max_records": 10,
#     "list_return": "name"
#     }
#     response = requests.get(url,params=payload)
    
#     print(response.reason, response.text)
#     if response.status_code != 200:
#         return None
#     similar_cids = response.json()["IdentifierList"]["CID"]
    
#     # Get drug names of similar drugs
#     similar_drugs = []
#     for cid in similar_cids:
#         url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalName/JSON"
#         response = requests.get(url)
#         if response.status_code == 200:
#             name = response.json()["PropertyTable"]["Properties"][0]["CanonicalName"]
#             similar_drugs.append(name)
#     return similar_drugs

# print(find_similar_drugs('Acetaminophen'))

from ordered_set import OrderedSet
import pubchempy as pcp

from drug_class import DrugClass
def findSimilarDrugs(DrugList):
# Search for the compound by CID
    ordersetList=OrderedSet(DrugList)
    # ordersetList.
    print("Length of ORderSet:",len(ordersetList))
    for drug in DrugList:
        c = pcp.get_compounds(drug.name, 'name')

        # Use the CID to search for similar compounds
        if not c:
            continue
        similar_compounds = pcp.get_compounds(identifier=c[0].cid, searchtype='similarity', threshold=90, listkey_count=5)

        for comp in similar_compounds:
            if(len(comp.synonyms)>0):
                ordersetList.add(DrugClass(comp.synonyms[0]))
        print("Last LEngth/",len(ordersetList))
    print("Final length: ",len(ordersetList))
    return ordersetList
