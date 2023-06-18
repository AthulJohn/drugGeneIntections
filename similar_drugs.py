

from ordered_set import OrderedSet
import pubchempy as pcp

from drug_class import DrugClass
def findSimilarDrugs(DrugList):
# Search for the compound by CID
    ordersetList=OrderedSet(DrugList)
    for drug in DrugList:
        c = pcp.get_compounds(drug.name, 'name')

        # Use the CID to search for similar compounds
        if not c:
            continue
        similar_compounds = pcp.get_compounds(identifier=c[0].cid, searchtype='similarity', threshold=90, listkey_count=5)

        for comp in similar_compounds:
            if(len(comp.synonyms)>0):
                ordersetList.add(DrugClass(comp.synonyms[0]))
    return ordersetList
