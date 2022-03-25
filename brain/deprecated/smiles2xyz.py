#!/usr/bin/env python3
import sys, os
from openbabel import openbabel
from openbabel import pybel 

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Lipinski
from rdkit.Chem import AllChem

from  xyz2mol import *
from data import *
from gpm import *
from smiles2xyz import * 

import itertools
from itertools import combinations
import numpy as np
import matplotlib.pyplot as plt
from statistics import mean
from scipy import stats


import platform, yaml, time
import pandas as pd
from pgradd.GroupAdd.Library import GroupLibrary
from pgradd.Error import GroupMissingDataError
import pgradd.ThermoChem

#lib = GroupLibrary.Load('GRWSurface2018')
lib = GroupLibrary.Load('BensonGA')

def save_xyz(coords, name):
    out_name = name + '.xyz'
    file_out = open(out_name, 'w')
    file_out.write(str(len(coords)) + '\n')
    file_out.write(name + '\n')
    file_out.writelines(coords)
    file_out.close()
    
def canonical_hack(mol):
    new_smiles = Chem.MolToSmiles(mol, isomericSmiles=False)    
    m = Chem.MolFromSmiles(new_smiles)
    new_smiles = Chem.MolToSmiles(m, isomericSmiles=False)    
    return new_smiles
def smiles_file_name(smiles):
    '''Avoid using the "(", ")", "["  and "]" for the file names'''
    file_name = smiles.replace('(', 'L').replace(')', 'R').replace('[', 'l').replace(']','r').replace('=', 'e') 
    return file_name

def get_smiles_m1(xyz_file):
    '''Method using the xyz2mol.py script '''
    charged_fragments = False
    quick = True
    atomicNumList, charge, xyz_coordinates = read_xyz_file(xyz_file)
    mol = xyz2mol(atomicNumList, charge, xyz_coordinates, charged_fragments, quick)
    formula = CalcMolFormula(mol)
    new_smiles = canonical_hack(mol)
    return new_smiles, formula

def get_smiles_m1_log(log_file):
    '''Method using the xyz2mol.py script '''
    charged_fragments = False
    quick = True
    charge = 0
    coords = file_analyzer(log_file)
    atomicNumList = []
    xyz_coordinates = []
    for coord in coords:
        infor = coord.rstrip().split()
        atom_num = dict_element_2.get(infor[0])
        atomicNumList.append(int(atom_num))
        xyz_coordinates.append([float(i) for i in infor[1:]])
    mol = xyz2mol(atomicNumList, charge, xyz_coordinates, charged_fragments, quick)
    formula = CalcMolFormula(mol)
    new_smiles = canonical_hack(mol)
    return new_smiles, formula

def get_smiles_m2(xyz_file):
    '''Method using the openbabel and RDkit '''
    mols = pybel.readfile('xyz', xyz_file)
    for mol in mols:
        formula = mol.formula
        smiles = mol.write('smi').split()[0]
        mol = Chem.MolFromSmiles(smiles)
        new_smiles = canonical_hack(mol)
        return new_smiles, formula

def get_smiles_m2_log(log_file):
    '''Method using the openbabel and RDkit '''
    mols = pybel.readfile('g09', log_file)
    for mol in mols:
        formula = mol.formula
        smiles = mol.write('smi').split()[0]
        mol = Chem.MolFromSmiles(smiles)
        new_smiles = canonical_hack(mol)
        return new_smiles, formula
    
def get_smiles_file(xyz_file):
    '''Method using the openbabel and RDkit '''
    smiles, formula  = get_smiles_m1(xyz_file)
    if '=' in smiles and ']' in smiles:
        smiles, formula = get_smiles_m2(xyz_file)
    smiles_file = smiles_file_name(smiles)
#    os.rename(xyz_file, smiles_file+'.xyz')
    #print(smiles)
    return smiles
    

def get_pdb_from_xyz(xyz_file):
    ''' Use pybel to convert the xyz file to pdb file'''
    mols = pybel.readfile('xyz', xyz_file)
    for mol in mols:
        mol.write('pdb', 'test.pdb',overwrite=True)

def get_mol_from_xyz(xyz_file):
    ''' Use pybel to convert the xyz file to mol file'''
    mols = pybel.readfile('xyz', xyz_file)
    name = xyz_file.split('.')[0] + '.mol'
    for mol in mols:
        mol.write(format='mol', filename=name, overwrite=True)
        mol.write(format='mpc', filename=name, overwrite=True)

def get_mol_from_babel(smiles):
    name = smiles_file_name(smiles) + '.mol'
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("smi", "mol")
    mol = openbabel.OBMol()
    mol = obConversion.ReadString(mol, smiles)   
    mol.AddHydrogens()
    obConversion.WriteFile(mol, name)
#    mol = readstring("smi", smiles)
#    mol.make3D(forcefield='mmff94',steps=10)
#    mol.write(format='mol', filename=name, overwrite=True)

def get_mol_from_rdkit(smiles):
    m0 = Chem.MolFromSmiles(smiles)
    m = Chem.AddHs(m0)
    AllChem.EmbedMolecule(m)
    name = smiles_file_name(smiles) + '_rd.mol'
    print(Chem.MolToMolBlock(m),file=open(name,'w+'))
    #    Chem.rdmolfiles.MolToMolFile(m, name)
    
### Part-1: Smiles and xyz conversions
def get_xyz_from_openbabel(smiles):
    '''use openbabel to convert the SMILES to xyz file '''
    mol = pybel.readstring("smi", smiles)
#    obConversion = openbabel.OBConversion()
#    obConversion.SetInAndOutFormats("smi", "mol")
#    mol = openbabel.OBMol()
#    obConversion.ReadString(mol, smiles)   
#    mol.AddHydrogens()
#    obConversion.WriteFile(mol, smiles_file_name(smiles) + '.xyz')
    name = smiles_file_name(smiles) 
    mol.make3D(forcefield='mmff94',steps=10)  #print(pybel.forcefields)
#    name = 'test_pybel' #smiles_file_name(smiles) 
    mol.write(format='xyz', filename=name+'.xyz', overwrite=True)
#    mol.write(format='mpc', filename=name+'.mop', overwrite=True)
#    mol.write(format='mol', filename=name+'.mol', overwrite=True)
    
def smiles_convert(smiles):
    '''Get thee Canonical SMILES 
    1) use openbabel to covert the initial smiles to xyz 
    2) read xyz, and convert to the canonical smiles
    '''     
    get_xyz_from_openbabel(smiles)
    file_xyz = smiles_file_name(smiles) + '.xyz'
    charged_fragments = False 
    quick = True
    atomicNumList, charge, xyz_coordinates = read_xyz_file(file_xyz)
    mol = xyz2mol(atomicNumList, charge, xyz_coordinates, charged_fragments, quick)
    new_smiles = canonical_hack(mol)
    
    return new_smiles
  
def get_coords_from_chemml(smiles):
    '''Use Chemml to convert the SMILES to xyz file, the same as using RDkit '''
    mol = Molecule(smiles, input_type='smiles')
    mol.hydrogens('add')
    mol.to_xyz(optimizer='MMFF', mmffVariant='MMFF94s', maxIters=100)
    
    xyz = mol.xyz.geometry
    symbols=mol.xyz.atomic_symbols
    N_num=len(symbols)
    
    coords = []
    for num, ele in enumerate(xyz):
        line = '%s\t%s\n' %(symbols[num][0], ' '.join([str(round(i,5)) for i in ele]))
        coords.append(line)
        
    name =  'test_chemml'#smiles_file_name(smiles)
    save_xyz(coords, name)
#    return coords

def get_coords_from_rdkit(smiles):
#    m2 = Chem.AddHs(Chem.MolFromSmiles('Oc1c(O)c(O)ccc1'))
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    energy = 100000000
    confIds = AllChem.EmbedMultipleConfs(m, 1000)
    for confId in confIds:
        m_i = Chem.Mol(m, confId)
        AllChem.EmbedMolecule(m_i, useRandomCoords=True)
        name = 'test_rdkit.mol' #smiles_file_name(smiles) + '_rd.mol'
        print(Chem.MolToMolBlock(m_i),file=open(name,'w+'))
        prop = AllChem.MMFFGetMoleculeProperties(m_i, mmffVariant="MMFF94") 
        ff =AllChem.MMFFGetMoleculeForceField(m_i, prop)
        if ff.CalcEnergy() < energy:
            energy = ff.CalcEnergy()
            for c in m_i.GetConformers():
                N_num = c.GetNumAtoms()
                xyz = c.GetPositions()
                
            symbols=[atom.GetSymbol() for atom in m_i.GetAtoms()] 
    
            coords = []       
            name = smiles_file_name(smiles)
            for num, ele in enumerate(xyz):
                line = '%s  %s\n' %(symbols[num][0], '  '.join(['%.5f' %(i) for i in ele]))
                coords.append(line)
#            save_gjf(coords, name)
            save_xyz(coords, name)
#            save_poscar()
#            return coords

def save_out(file_name, radical_list):
    f_out = open(file_name, 'w')
    f_out.write('\n'.join(radical_list))
    f_out.close()

####Part-2 Get ALL the SMILES from RING output of one reaction
def get_radicals(reaction):
    infor = reaction.rstrip().split()[0].split('>>')
    reactant, product_1, product_2 = ['','','']
    if len(infor[0]) > len(infor[1]):
        radicals = infor[0].replace('.]', 'D]')
        reactant = infor[1]
    else:    
        reactant = infor[0]
        radicals = infor[1].replace('.]', 'D]')

    products = [ i.replace('D]', '.]') for i in radicals.split('.')]
    return reactant, products

def get_smiles_from_ring(ring_out):
    ''' ring_out is the file reactions_SMILES.txt'''
    f = open(ring_out, 'r')
    lines = f.readlines()
    f.close()
    
    reaction_dict = {}
    reaction_name_list =[]
    reaction_num_list = []
    molecules = []
    radicals = []
    for num, line in enumerate(lines):
        if 'Reactions of rule' in line:
            reaction_name = line.strip().split()[-1][:-1]
            reaction_name_list.append(reaction_name)
            reaction_num_list.append(num)
    reaction_num_list.append(len(lines))        
    
    for num, reaction in enumerate(reaction_name_list):
        '''Count the slices for one type of reaction '''
        start_num = reaction_num_list[num] + 1
        end_num   = reaction_num_list[num+1]
        key = reaction
        value = lines[start_num:end_num]
        reaction_dict.update({key:value})
        
    reaction_out = open('reactions.out', 'w')    
    for key, value in reaction_dict.items():
        for reaction in value:
            reactant, products = get_radicals(reaction)
            try: 
                reactant, products = get_radicals(reaction)
                if '.' in reactant:
                    if reactant not in radicals:
                        radicals.append(reactant)
                else:
                    if reactant not in molecules:
                        molecules.append(reactant)
                for i in products:
                    if i not in radicals:
                        radicals.append(i)
                reaction_out.write('%s\t%s\t%s\n' %(key, reactant, '\t'.join(products) ))
            except:
                 print(reaction)
    reaction_out.close()        
#    
    mol_ring = list(dict.fromkeys(molecules))
    radicals_ring = list(dict.fromkeys(radicals))
    
    save_out('mol_ring.out', mol_ring)
    save_out('rad_ring.out', radicals_ring)

#### Part-3: Analyze the SMILES from RING 
def analyze_one_radical_smiles(radical):

    def try_to_convert_one_radical_smiles(radical):
        out='Yes'
        try:
            new_smiles = smiles_convert(radical)
            if '[' not in new_smiles:
                out='No'
        except ValueError:
            out='No'
        return out

    new_radical = radical 

    if '=[C.]' in radical or '[C.]=' in radical:
        new_radical = radical
    elif '[O.]' in radical:
        new_radical = radical
    else:
        if '[C.]' in radical:
            if radical[-4:] == '[C.]':
                '''  XXX-CH2*  '''
                new_radical = radical.replace('[C.]', '[CH2]')
            elif '[C.])' in radical:
                ''' XXX-(XXX-CH2*)-XXX '''
                new_radical = radical.replace('[C.]', '[CH2]')
            elif radical[:5] == '[C.](':
                if try_to_convert_one_radical_smiles(radical.replace('[C.]', '[CH]')) =='Yes':
                    ''' XXX-CH*-XXX '''
                    new_radical = radical.replace('[C.]', '[CH]')
                else:
                    ''' XXX-C*-XXX  '''
                    new_radical = radical.replace('[C.]', '[C]')
            elif '[C.](' in radical:
                '''  XXX-C*(XX)-XX '''
                new_radical = radical.replace('[C.]', '[C]')
            else:
                new_radical = radical.replace('[C.]', '[CH]')
    return new_radical

def analyze_ring_molecule_smiles(ring_molecule_out):
    file_in = ring_molecule_out
    f = open(file_in, 'r')
    lines = f.readlines()
    f.close()
    
    mol_rdkit = []
    for line in lines:
        smiles = line.rstrip() 
        mol = Chem.MolFromSmiles(smiles)
        try: 
            smiles_rdkit  = Chem.MolToSmiles(mol,isomericSmiles=False)
            mol_rdkit.append(smiles_rdkit)
        except:
            print(smiles)
    save_out('mol_final_rdkit.out', mol_rdkit)    
    
def analyze_ring_radical_smiles(ring_radical_out):
    file_in = ring_radical_out
    f = open(file_in, 'r')
    lines = f.readlines()
    f.close()

    radical_total = []
#    radical_ring  = []
#    radical_chain_t = []
#    radical_chain_d = []
#    radical_chain_m = []
#    radical_o     = []
#    radical_others = []

    for line in lines:
        radical = line.rstrip()
        try: 
            new_radical = analyze_one_radical_smiles(radical)
            final_smiles = smiles_convert(new_radical)
            radical_total.append(final_smiles)
        except:
            print('Failed to convert %s.' %(radical))
        
#        print(radical, new_radical)
#        radical_total.append(new_radical)
#        if '=[C.]' in radical or '[C.]=' in new_radical:
#            radical_ring.append(new_radical)
#        elif '[O.]' in new_radical:
#            radical_o.append(new_radical)
#        elif '[CH2]' in new_radical:
#            radical_chain_t.append(new_radical)
#        elif '[CH]' in new_radical:
#            radical_chain_d.append(new_radical)
#        elif '[C]' in new_radical:
#            radical_chain_m.append(new_radical)
#        else:
#            radical_others.append(new_radical) 

    save_out('rad_final_rdkit.out', radical_total)
#    save_out('radical_ring.out', radical_ring)
#    save_out('radical_o.out', radical_o)
#    save_out('radical_chain_t.out', radical_chain_t)
#    save_out('radical_chain_d.out', radical_chain_d)
#    save_out('radical_chain_m.out', radical_chain_m)
#    save_out('radical_others.out', radical_others)

### Hbonds
#from __future__ import print_function
#from rdkit import Chem
#from rdkit.Chem import ChemicalFeatures
#from rdkit import RDConfig
#import os
#from rdkit.Chem import Lipinski,rdMolDescriptors
#
#
#m = Chem.MolFromSmiles('Cc1ccc(OC(O)Cc2ccc(O)cc2)cc1')
#for bond in m.GetBonds():
#    print(bond.GetBondType())
##print(rdMolDescriptors.CalcNumHBA(m))
##print(rdMolDescriptors.CalcNumHBD(m))
#
##print(Lipinski.NHOHCount(m))
##print(Lipinski.NumAromaticRings(m))
##print(Lipinski.RingCount(m))
##print(Lipinski.NumHAcceptors(m))
##print(Lipinski.NumHDonors(m))
#
#
##fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
##factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
##
##feats = factory.GetFeaturesForMol(m)
##for i in feats:
##    if i.GetType() in ['SingleAtomDonor', 'SingleAtomAcceptor']:
##        print(i.GetType())
##print(i.GetType() for i in feats)



def get_patt_num(smiles, patt):
    m = Chem.MolFromSmiles(smiles)
    try:
        N_patt = len(m.GetSubstructMatches(Chem.MolFromSmarts(patt)))
    except:
        N_patt = 0
    return N_patt


def replace_patt(smiles, patt, repl):
    '''Replace one pattern with another one'''
    patt = Chem.MolFromSmarts(patt)
    repl = Chem.MolFromSmarts(repl)
    m = Chem.MolFromSmiles(smiles)
    rms = AllChem.ReplaceSubstructs(m,patt,repl,replaceAll=True, replacementConnectionPoint=0)
    new_smiles = Chem.MolToSmiles(rms[0],isomericSmiles=False, kekuleSmiles=False)
    return new_smiles


# In[ ]:


def pre_treat_smiles(smiles):
    '''Pretreat the SMILES'''
    smiles = smiles.replace('[CH2]', 'C').replace('[CH]', 'C').replace('[C]', 'C').replace('[c]', 'c')
    smiles = smiles.replace('(C)', '').replace('(CC)', '').replace('(CCC)', '').replace('(CCCC)', '')
    smiles = smiles.replace('(OCC)', '(OC)').replace('(OCC)', 'OC')    
    smiles = smiles.replace('(c1ccccc1)', '').replace('(Cc1ccccc1)', '').replace('(CCc1ccccc1)', '').replace('(c2ccccc2)', '')  
    return smiles


# In[ ]:


def update_arOH_alone(smiles):   
    '''There is one Ar-OH but no H-bond formed in the ring.'''
    pattern = '[cX3;!$(cO)][cX3;$(c[O;H1])][cX3;!$(cO)]'
    return get_patt_num(smiles, pattern)


# In[ ]:


def update_arO_alone(smiles):  
    '''Ar-O* radical alone will increase the system energy '''   
    pattern = '[cX3;!$(c[O;H1])]c([O;X1v1+0])[cX3;!$(c[O;H1])]'
    return get_patt_num(smiles, pattern)


# In[ ]:


def update_arOR_crowd(smiles): 
    ''' c(R)c(OC)c(R) (R is not H) 
    the OC in the middle will perpendicular to the Ar ring, 
    and this will increase the energy'''
    # pattern = '[cX3;!$(c[H])][cX3;$(cO[!H])][cX3;!$(c[H])]'
    pattern = '[cX3;!$([c;H1])][cX3;$(cO[!H])][cX3;!$([c;H1])]'
    return get_patt_num(smiles, pattern)


# In[ ]:


def update_arO_hbond(smiles):  
    '''Hbonds formed by Ar-O* and adjacent Ar-OH in the same ring'''
    pattern = '[cX3;$(c[O;H1])]c([O;X1v1+0])'
    return get_patt_num(smiles, pattern)


# In[ ]:


def update_ar_hbond(smiles):  
    ''' Hbonds formed by two Ar-OH groups in the same ring'''
    pattern_1 = '[cX3;$(c[O;H1])][cX3;$(c[O;H1])]'
    pattern_2 = '[cX3;$(c[O;H1])][cX3;$(cO[!H])]' # Hbonds formed by Ar-OH and Ar-OC
    list_p = [pattern_1, pattern_2]
    return sum([get_patt_num(smiles, p) for p in list_p])


# In[ ]:


def update_ar_ali_hbond(smiles):  
    '''Calculate the H bonds formed by the OH in the ring and O from Aliphatic chain'''
    pattern_1 = '[cX3;$(c[O;H1])][cX3;$(cC)]C([OH])'
    pattern_2 = '[cX3;$(c[O;H1])][cX3;$(cC)][CX4;!$(C[O;H1])]C([OH])'
    list_p = [pattern_1, pattern_2]
    return sum([get_patt_num(smiles, p) for p in list_p])

# In[ ]:


def update_ar_ar_hbond(smiles): 
    '''Calculate the H bond numbers from two OH in two rings: (O)cc-cc(O)'''
    pattern = '[cX3;$(c[O;H1])][cX3;$(cc)][cX3;$(cc)][cX3;$(c[O;H1])]'
    return get_patt_num(smiles, pattern)


# In[ ]:


def update_ali_hbond(smiles):  
    '''Calculate the H bond numbers for aliphatic part'''
    pattern_1 = 'C([OH])C([OH])' # C(O)C(O)
    pattern_2 = 'C([OH])[CX4;!$(C[O;H1])]C([OH])'  #C(O)CC(O)
    pattern_3 = 'C([OH])[CX4;!$(C[O;H1])][CX4;!$(C[O;H1])]C([OH])' #C(O)CCC(O)
    pattern_4 = 'C([OH])[CX4;!$(C[O;H1])][CX4;!$(C[O;H1])][CX4;!$(C[O;H1])]C([OH])' #C(O)CCCC(O)
    
    pattern_r1 = 'C([OH])C([O;X1v1+0])' # C(O)C([O])
    pattern_r2 = 'C([OH])[CX4;!$(C[O;H1])]C([O;X1v1+0])'  #C(O)CC([O])
    pattern_r3 = 'C([OH])[CX4;!$(C[O;H1])][CX4;!$(C[O;H1])]C([O;X1v1+0])' #C(O)CCC([O])
    pattern_r4 = 'C([OH])[CX4;!$(C[O;H1])][CX4;!$(C[O;H1])][CX4;!$(C[O;H1])]C([O;X1v1+0])' #C(O)CCCC([O])
    
    list_p = [pattern_1, pattern_2, pattern_3, pattern_4, pattern_r1, pattern_r2, pattern_r3, pattern_r4]
    return sum([get_patt_num(smiles, p) for p in list_p])


# In[ ]:


def update_ether_hbond(smiles): 
    '''Some H bond patterns in ether linkages like: Ar-O-CH2OH, CH3-O-CH2OH'''
    pattern_r1 = 'cO[CX4;$(C[O;H1])]' # Ar-O-CH(R)OH
    pattern_r2 = 'cO[CX4;!$(C[O;H1])][CX4;$(C[O;H1])]' # Ar-O-CH2CH(R)OH
    pattern_r3 = 'cO[CX4;!$(C[O;H1])][CX4;!$(CO)][CX4;$(C[O;H1])]'# Ar-O-CH2CH2-CH(R)-OH
    
    pattern_a1 = '[CX4;!$(C[O;H1])]O[CX4;$(C[O;H1])]' # C-O-CH(R)OH
    pattern_a2 = '[CX4;!$(C[O;H1])]O[CX4;!$(C[O;H1])][CX4;$(C[O;H1])]' # C-O-CH2CH(R)OH
    pattern_a3 = '[CX4;!$(C[O;H1])]O[CX4;!$(C[O;H1])][CX4;!$(C[O;H1])][CX4;$(C[O;H1])]'# C-O-CH2CH2-CH(R)-OH
    
    list_p = [pattern_r1, pattern_r2, pattern_r3, pattern_a1, pattern_a2, pattern_a3]
    return sum([get_patt_num(smiles, p) for p in list_p])


# In[ ]:


def update_none_Arrings(smiles):  
    '''Calcuclate the intra-molecule rings like: c1-C-C-C-c1, 
    this kind of rings will increase system energy'''
    mol = Chem.MolFromSmiles(smiles)
    return rdMolDescriptors.CalcNumAliphaticRings(mol)


# In[9]:


def analyze_smiles(smiles):
    
    smiles = pre_treat_smiles(smiles) 
    NH_ar = update_ar_hbond(smiles)
    NH_ar_ar = update_ar_ar_hbond(smiles)
    NH_ar_ali = update_ar_ali_hbond(smiles)
    NH_ali = update_ali_hbond(smiles)
    NH_ether = update_ether_hbond(smiles)
    NH_total = NH_ar + NH_ar_ar + NH_ar_ali + NH_ali + NH_ether 
    
    num_arO_hbond  = update_arO_hbond(smiles)  # Separate the stronger H bonds
    
    num_arO_alone = update_arO_alone(smiles)
    num_arOR_crowd = update_arOR_crowd(smiles)
    num_none_Arrings = update_none_Arrings(smiles)
 #   print(NH_total, num_arO_alone, num_arOR_crowd, num_arO_hbond, num_none_Arrings )
    return NH_total, num_arO_alone, num_arOR_crowd, num_arO_hbond, num_none_Arrings 
