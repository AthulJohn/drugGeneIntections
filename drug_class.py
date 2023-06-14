import re
import requests
import pubchempy as pcp


class DrugClass:
    name=""
    drug_id=None
    cid=None

    def __init__(self, name, drug_class_id=None):
        self.name = name
        self.drug_id=self.findDrugID(name)
        self.cid=self.findDrugCID(name)

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return NotImplemented
        return self.drug_id==other.drug_id
    
    def __hash__(self):
        return hash(self.drug_id)

    def findDrugID(self,name):
        drug_id=None
        druglist = requests.get('https://rest.kegg.jp/find/drug/'+name)
        content=druglist.content.decode("utf-8")
        # initializing substrings
        sub1 = "dr:"
        sub2 = "\t"
        s=str(re.escape(sub1))
        e=str(re.escape(sub2))
        res=re.findall(s+"(.*)"+e,content)
        suba1="\t"
        suba2=" ("
        s=str(re.escape(suba1))
        e=str(re.escape(suba2))
        drugnames=re.findall(s+"(.*?)"+e,content)
        try:
            return res[drugnames.index(name.capitalize())]
        except:
            return ""
    
    def findDrugCID(self,name):
        if(len(pcp.get_compounds(name, 'name'))==0):
            return None
        return pcp.get_compounds(name, 'name')[0].cid 
       
