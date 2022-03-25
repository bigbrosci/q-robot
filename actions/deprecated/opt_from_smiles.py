#!/usr/bin/env python3

from rdkit import Chem
from rdkit.Chem import AllChem

m = Chem.MolFromSmiles('Cc1ccccc1')
m2 = Chem.AddHs(m)
AllChem.EmbedMolecule(m2)
AllChem.MMFFOptimizeMolecule(m2)
#print(Chem.MolToMolBlock(m2))

cids = AllChem.EmbedMultipleConfs(m2, numConfs=10)

rmslist = []
AllChem.AlignMolConformers(m2, RMSlist=rmslist)
#for i in rmslist:
#    print(i)

res = AllChem.MMFFOptimizeMoleculeConfs(m2, maxIters=600)
print(res)
#print(Chem.MolToMolBlock(res))
