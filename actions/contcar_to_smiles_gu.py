#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Qiang , Vlachos Group, University of Delaware, 2021
"""
import os,sys
from ase.io import read
import numpy as np
from rdkit import Chem
from ase import Atoms as ase_Atoms

from rdkit import Chem
from rdkit.Chem import RWMol
from collections import defaultdict

class adsorbate(object):
    """
    This is an adsorbate graph class that converts atomic coordinates to rdkit 
    molecular graph object, Mol. Use "LoadByCovalentRadius" to initialize.
    
    Class Variables
    soan: selected organic atomic number. These atoms are considered adosbates
    rcov: covalent radius. Info available in wikipedia.
    
    Class Attributes
    ASEAtoms:                   ASE Atoms object.
    RdkitMol:                   Rdkit Mol object.
    SurfaceAtomSymbols:         List of symbols of surface atoms.
    ASEAtomIndex2RdKitAtomIndex: Index mapping from ASE atoms to Rdkit Mol
    RdKitAtomIndex2ASEAtomIndex: Index mapping from Rdkit Mol to ASE Atoms.
    """
    # selected organic atomic numbers
    soan = [1,6,7,8] 
    # atomic number -> covalent radius
    rcov = {1:0.31, 6:0.76, 7:0.71, 8:0.66, 26:1.26, 27:1.21, 28:1.21, 29:1.21,\
            44:1.16, 45:1.21, 46:1.26, 47:1.46, 75:1.21, 77:1.21 ,78:1.21, 79:1.21}
    
    def __init__(self,ASEAtoms,RdkitMol,SurfaceAtomSymbols, \
                 ASEAtomIndex2RdKitAtomIndex, RdKitAtomIndex2ASEAtomIndex):
        
        assert isinstance(ASEAtoms,ase_Atoms)
        assert isinstance(RdkitMol,Chem.Mol)
        assert isinstance(ASEAtomIndex2RdKitAtomIndex,dict)
        assert isinstance(RdKitAtomIndex2ASEAtomIndex,dict)
        if isinstance(SurfaceAtomSymbols,str):
            SurfaceAtomSymbols = [SurfaceAtomSymbols]
        else:
            assert isinstance(SurfaceAtomSymbols,list)
        self.ASEAtoms = ASEAtoms
        self.RdkitMol = RdkitMol
        self.SurfaceAtomSymbols = SurfaceAtomSymbols
        self.ASEAtomIndex2RdKitAtomIndex = ASEAtomIndex2RdKitAtomIndex
        self.RdKitAtomIndex2ASEAtomIndex = RdKitAtomIndex2ASEAtomIndex
    
    @classmethod
    def LoadByCovalentRadius(cls,CoordinateFPath, SurfaceAtomSymbols, \
        rfacup = 1.35,rfacdown = 0.6, z_vector = 2):
        """ 
        This function reads file using ASE read, and construts molecular graph
        in rdkit object, Mol. See manuscript for overall algorithm.
        
        
        Input List
        CoordinateFPath:    path to ASE readable coordinate file.
        SurfaceAtomSymbols: List of atomic symbols of surface atoms.
        rfacup:             Upper percentage limit for determining connectivity.
        rfacdown:           Lower percentage limit for determining connectivity.
        z_vector:           index of cell basis vector that is orthogonal to surface.
        
        Output List
        adsorbate class
        """
        
        # initialize
        ASEAtomIndex2RdKitAtomIndex = dict()
        RdKitAtomIndex2ASEAtomIndex = dict()
        if isinstance(SurfaceAtomSymbols,str):
            SurfaceAtomSymbols = [SurfaceAtomSymbols]
        else:
            assert isinstance(SurfaceAtomSymbols,list)
        # load POSCAR
        AseAtoms = read(CoordinateFPath)
        # if none given for surface layer z coordinate, average the top layer atomic coordinate
        _, SurfaceAtomIndex = DetermineSurfaceLayerZ(AseAtoms, SurfaceAtomSymbols, ZVecIndex = z_vector)

        # (p)eriodic (b)oundary (c)ondition(s)
        PBCs = [[0,0,0]]
        if AseAtoms.pbc[0]:
            temp = np.add(PBCs,[1,0,0])
            temp = np.concatenate((temp,np.add(PBCs,[-1,0,0])))
            PBCs = np.concatenate((PBCs,temp))
        if AseAtoms.pbc[1]:
            temp = np.add(PBCs,[0,1,0])
            temp = np.concatenate((temp,np.add(PBCs,[0,-1,0])))
            PBCs = np.concatenate((PBCs,temp))
        if AseAtoms.pbc[2]:
            temp = np.add(PBCs,[0,0,1])
            temp = np.concatenate((temp,np.add(PBCs,[0,0,-1])))
            PBCs = np.concatenate((PBCs,temp))
                    
        # Get organic atoms from the DFT calculations (their index and atomic number)
        ans = AseAtoms.get_atomic_numbers() # (a)tomic (n)umber(s)
        oai = list() #organic atom index in the atoms object
        oan = list() #organic atomic number
        for i in range(0,AseAtoms.__len__()):
            if ans[i] in cls.soan:
                oai.append(i)
                oan.append(ans[i])
        
        # Determine connectivity of the organic atoms
        adj_mat = np.zeros((oai.__len__(),oai.__len__())) # adjacency matrix
        for i in range(0,oai.__len__()):
            for j in range(i+1,oai.__len__()):
                if cls._DetermineConnectivity(AseAtoms,oai[i],oai[j],PBCs,rfacup,rfacdown):
                    adj_mat[i,j] = 1
        
        # construct mol object
        RdkitMol = Chem.Mol()
        RdkitMol = Chem.RWMol(RdkitMol)
        
        #print(oan)
        ## add atom
        ### organic atoms
        for i in range(0,len(oan)):
            atom = Chem.Atom(int(oan[i]))
            atom.SetNoImplicit(True) # this allows molecule to have radical atoms
            atom.SetBoolProp('Adsorbed',False)
            RdkitMol.AddAtom(atom)
            ASEAtomIndex2RdKitAtomIndex[oai[i]] = i
            RdKitAtomIndex2ASEAtomIndex[i] = oai[i]
        ### surface atoms
        for index in SurfaceAtomIndex:
            atom = Chem.Atom(AseAtoms[index].symbol)
            atom.SetBoolProp('SurfaceAtom',True)
            atom.SetBoolProp('Occupied',False)
            i = RdkitMol.AddAtom(atom)
            ASEAtomIndex2RdKitAtomIndex[index] = i
            RdKitAtomIndex2ASEAtomIndex[i] = index
        
        ## add bond
        ### between organic atoms
        for i in range(0,oai.__len__()):
            for j in range(i+1,oai.__len__()):
                if adj_mat[i,j] == 1:
                    RdkitMol.AddBond(i,j,order=Chem.rdchem.BondType.SINGLE)
                    
        ### between surface atoms
        for i in range(0,len(SurfaceAtomIndex)):
            for j in range(i+1,len(SurfaceAtomIndex)):
                if cls._DetermineConnectivity(AseAtoms,SurfaceAtomIndex[i],SurfaceAtomIndex[j],PBCs,rfacup,rfacdown):
                    RdkitMol.AddBond(ASEAtomIndex2RdKitAtomIndex[SurfaceAtomIndex[i]],ASEAtomIndex2RdKitAtomIndex[SurfaceAtomIndex[j]],order=Chem.rdchem.BondType.ZERO)
                    
        ## assign radicals
        Chem.AssignRadicals(RdkitMol)
        
        ## set smilesSymbol
        for atom in RdkitMol.GetAtoms():
            if atom.GetSymbol() in ['C','O'] and atom.GetNumRadicalElectrons() == 0:
                atom.SetProp("smilesSymbol",'[' + atom.GetSymbol() + str(atom.GetNumRadicalElectrons())+ ']')
            elif atom.GetNumRadicalElectrons() > 0:
                atom.SetProp("smilesSymbol",atom.GetSymbol() + str(atom.GetNumRadicalElectrons()))
            
        # Find surface binding atom. This is done by finding all the radical atoms
        rai_rdkit = list() # radical atom index for rdkit mol
        rai_ase = list() # radical atom index for rdkit ase atoms object
        for atom in RdkitMol.GetAtoms():
            if atom.GetNumRadicalElectrons() > 0:
                rai_rdkit.append(atom.GetIdx())
                rai_ase.append(oai[atom.GetIdx()])
        
        # Surface connectivity
        for i in range(0,len(rai_ase)):
            for j in range(0,len(SurfaceAtomIndex)):
                if cls._DetermineConnectivity(AseAtoms,rai_ase[i],SurfaceAtomIndex[j],PBCs,rfacup,rfacdown):
                    RdkitMol.AddBond(rai_rdkit[i],ASEAtomIndex2RdKitAtomIndex[SurfaceAtomIndex[j]],order=Chem.rdchem.BondType.ZERO)
                    RdkitMol.GetAtomWithIdx(ASEAtomIndex2RdKitAtomIndex[SurfaceAtomIndex[j]]).SetBoolProp('Occupied',True)
                    RdkitMol.GetAtomWithIdx(rai_rdkit[i]).SetBoolProp('Adsorbed',True)
        
        # assign binding site.
        for i in range(0,len(rai_rdkit)):
            a = RdkitMol.GetAtomWithIdx(rai_rdkit[i])
            nsurf = 0
            for neighbor_atom in a.GetNeighbors():
                if neighbor_atom.GetSymbol() in SurfaceAtomSymbols:
                    nsurf += 1
            a.SetProp("smilesSymbol",a.GetProp("smilesSymbol") + '_' + str(nsurf) + 'fold')
            
        adsorbate = cls(AseAtoms,RdkitMol,SurfaceAtomSymbols, \
                 ASEAtomIndex2RdKitAtomIndex, RdKitAtomIndex2ASEAtomIndex)
        
        return adsorbate 

    @classmethod
    def _DetermineConnectivity(cls,AseAtoms,i,j,PBCs,rfacup,rfacdown):
        """
        Determine connectivity between atom i and j. See equation (1) in the 
        manuscript.
        
        Input List
        ASEAtoms:           ASE atoms containing adsorbate/surface system
        PBCs:               Periodic Boundary Conditions. e.g., (1,0,0) means 
                            cell repeats in first basis vector but not others.
        rfacup:             upper tolerance factor
        rfacdown:           lower tolerance factor
        
        Output List
        Bool:               True if connected, false if not.
        """
        xyz1 = AseAtoms[i].position
        # loop over periodic cells
        for PBC in PBCs:
            xyz2 = AseAtoms[j].position + np.dot(PBC,AseAtoms.cell)
            # Criteria:
            # TolFaclower * ideal_distance < distance < TolFacupper * ideal_distance 
            # ideal ideal_distance = Rcov(Atom1) + Rcov(Atom2)
            d = np.linalg.norm(xyz1-xyz2) # distance
            i_d = cls.rcov[AseAtoms[i].number] + cls.rcov[AseAtoms[j].number] # ideal distance
            if d <= i_d*rfacup and d >= i_d*rfacdown:
                return True
        return False
    
    
def DetermineSurfaceLayerZ(ASEAtoms, SurfaceAtomSymbols, ZVecIndex = 2, ztol = 0.5):
    """
    Find top layer surface atom z coordinates by averaging
    atoms within ztol (angstrom) of the top most atoms are selected for averaging
    
    Input List
    ASEAtoms:           ASE atoms containing adsorbate/surface system.
    SurfaceAtomSymbols: Symbol of surface atoms.
    ZVecIndex:          index of cell basis vector that is orthogonal to surface.
    ztol:               Atoms within ztol(angstrom) of the top most atoms are selected as 
                        surface atoms.
    Output List
    SurfaceLayerZ:      z coordinate of surface layer.
    SurfaceAtomIndex:   Index of surface atoms.
    """
    assert isinstance(ASEAtoms,ase_Atoms)
    # get highest surface atom coordinate
    zmax = 0
    zs = ASEAtoms.get_scaled_positions()[:,2]
    for i in range(0,len(ASEAtoms)):
        if ASEAtoms[i].symbol in SurfaceAtomSymbols and zmax < zs[i]:
            zmax = zs[i]
            
    # determine z coordinate. average out top layer
    ztol = ztol/np.linalg.norm(ASEAtoms.cell[2,:])
    SurfaceAtomIndex = list()
    SurfZs = list()
    for i in range(0,len(ASEAtoms)):
        if ASEAtoms[i].symbol in SurfaceAtomSymbols and zmax - ztol < zs[i]:
            SurfZs.append(zs[i])
            SurfaceAtomIndex.append(i)
    SurfaceLayerZ = np.array(SurfZs).mean()

    return SurfaceLayerZ, SurfaceAtomIndex


##################

def _GetSMILES(mol,idxlist):
    tmol = mol.__copy__() #(t)emporary
    tmol = RWMol(tmol)
    for AtomIdx in range(tmol.GetNumAtoms()-1,-1,-1):
        if AtomIdx not in idxlist:
            tmol.RemoveAtom(AtomIdx)
    return Chem.MolToSmiles(tmol)

def LumpH(molecule):
    """
    Lump hydrogen atoms as a single atom. Note that Si, Al, Mg, Na are used as 
    pseudoatoms. However, this does not affect printing SMILES, as smilesSymbol
    are appropriately set.
    """
    molecule = Chem.RWMol(molecule)
    Hidx = list()
    for i in range(0,molecule.GetNumAtoms()):
        atom = molecule.GetAtomWithIdx(i)
        if atom.GetSymbol() != 'H':
            NumH = 0
            for neighbor_atom in atom.GetNeighbors():
                if neighbor_atom.GetSymbol() == 'H':
                    NumH += 1
                    Hidx.append(neighbor_atom.GetIdx())
            if NumH == 4:
                a = Chem.Atom('Si')
                a.SetProp('smilesSymbol','H4')
                idx = molecule.AddAtom(a)
                molecule.AddBond(atom.GetIdx(),idx,Chem.rdchem.BondType.QUADRUPLE)
                molecule.GetAtomWithIdx(idx).SetNoImplicit(True)
            elif NumH == 3:
                a = Chem.Atom('Al')
                a.SetProp('smilesSymbol','H3')
                idx = molecule.AddAtom(a)
                molecule.AddBond(atom.GetIdx(),idx,Chem.rdchem.BondType.TRIPLE)
                molecule.GetAtomWithIdx(idx).SetNoImplicit(True)
            elif NumH == 2:
                a = Chem.Atom('Mg')
                a.SetProp('smilesSymbol','H2')
                idx = molecule.AddAtom(a)
                molecule.AddBond(atom.GetIdx(),idx,Chem.rdchem.BondType.DOUBLE)
                molecule.GetAtomWithIdx(idx).SetNoImplicit(True)
            elif NumH == 1:
                a = Chem.Atom('Na')
                a.SetProp('smilesSymbol','H')
                idx = molecule.AddAtom(a)
                molecule.AddBond(atom.GetIdx(),idx,Chem.rdchem.BondType.SINGLE)
                molecule.GetAtomWithIdx(idx).SetNoImplicit(True)
    Hidx.sort(reverse=True)
    for i in Hidx:
        molecule.RemoveAtom(i)
    return molecule

def GetAdsorbateGraphsOfRadius(mol,SurfaceAtomSymbols,radius):
    """
    Adsorbate Graph mining Tool. Molecular graph of the surface and adsorbate 
    are inputed. Graphs radially increases similar to Morgan Fingerprint.
    
    Input List
    mol:                RdKit Mol object of surface and adsorbate
    SurfaceAtomSymbols: Surface atom symbols
    radius:             Desired radius
    
    Output List
    SMILES_count:       dictionary of SMILES -> count
    """
    # check if mol is rdkit mol    
    # each smiles contains indexes, which are later converted to mol.
    assert isinstance(mol,Chem.rdchem.Mol)
    assert radius >= 1
    # Go through the molecule, add non-surface atoms to relevant atom list
    NSAL = list() # (N)on(S)urface (A)tom (L)ist
    SAL = list() # (S)urface (A)tom (L)ist
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in SurfaceAtomSymbols:
            NSAL.append(atom.GetIdx())
        else:
            SAL.append(atom.GetIdx())
    # Get finger prints
    smiless = list() # (F)inger(p)rint(s). It's in atom index, and will be
                    # converted to SMILES later. This is done to remove duplicates
    ## Radius = 1 case
    for idx in NSAL:
        smiless.append([idx])
    ## Radius > 1 cases
    for RA in NSAL:
        OAL = list() # (O)ld (A)tom (L)ist
        NAL = list() # (N)ew (A)tom (L)ist
        NAL.append(RA)
        for i in range(1,radius):
            ACL = list(NAL) # (A)tom to (C)heck (L)ist
            NAL = list()
            for AC in ACL:
                OAL.append(AC)
                for NeighborAtom in mol.GetAtomWithIdx(AC).GetNeighbors():
                    NeighborAtomIdx = NeighborAtom.GetIdx()
                    if NeighborAtomIdx not in ACL+OAL+NAL and \
                        NeighborAtom.GetSymbol() not in SurfaceAtomSymbols:
                        NAL.append(NeighborAtomIdx)
            smiles = sorted(OAL+NAL)
            if smiles not in smiless:
                smiless.append(smiles)
    
    # Add surface atom that is neighbor to SMILES
    for smiles in smiless:
        RSAL = list() # (R)elevant (S)urface (A)tom (L)ist
        for SA in SAL:
            for NeighborAtom in mol.GetAtomWithIdx(SA).GetNeighbors():
                NeighborAtomIdx = NeighborAtom.GetIdx()
                if NeighborAtomIdx in smiles:
                    RSAL.append(SA)
        smiles += RSAL
    # convert atom indexes to SMILES
    SMILES_count = defaultdict(int)
    for smiles in smiless:
        SMILES_count[_GetSMILES(mol,smiles)] += 1
    return SMILES_count
    
def FindAllAdsorbateGraphsOfLengthN(mol1,SurfaceAtomSymbols,NMaxEdge,NMinEdge=0,valence_based = False, debug = False):
    """
    Adsorbate Graph mining Tool. Molecular graph of the surface and adsorbate 
    are inputed. Graphs of any shape with lengths specified are obtained
    
    Input List
    mol1:               RdKit Mol object of surface and adsorbate
    SurfaceAtomSymbols: Surface atom symbols
    NMaxEdge:           Maximum number of edges. if -1, find everything
    NMinEdge:           Minimum number of edges
    valence_based:      Valence based graphs obtained
    debug:              Debugging function
    
    Output List
    SMILES_count:       dictionary of SMILES -> count
    """
    # AdsorbateSMILES of length N is found. Not radial.
    # if NMaxEdge = -1, find everything
    
    # mol2 is copied, then surface atoms are removed, and the enumeration is performed
    # check if mol is rdkit mol    
    assert isinstance(mol1,Chem.rdchem.Mol)
    mol1 = LumpH(mol1)
    print(Chem.MolToSmiles(mol1, isomericSmiles=False))
    mol2 = mol1.__copy__() # mol without surface atoms
    mol2 = Chem.RWMol(mol2)
    print(Chem.MolToSmiles(mol2, isomericSmiles=False))
    # Go through the molecule, add non-surface atoms to relevant atom list
    Mol2toMol1Map = dict()
    NSAL = list() # (N)on(S)urface (A)tom (L)ist
    SAL = list() # (S)urface (A)tom (L)ist
    mol2i = 0
    for mol1i in range (0,mol1.GetNumAtoms()):
        atom = mol1.GetAtomWithIdx(mol1i)
        if atom.GetSymbol() not in SurfaceAtomSymbols:
            NSAL.append(mol1i)
            Mol2toMol1Map[mol2i] = mol1i
            mol2i += 1
        else:
            SAL.append(mol1i)
            mol2.RemoveAtom(mol2i)
    # check N
    if NMaxEdge >= mol2.GetNumBonds():
        NMaxEdge = mol2.GetNumBonds()
    if NMinEdge >= mol2.GetNumBonds():
        NMinEdge = mol2.GetNumBonds()
    if NMaxEdge == -1:
        NMaxEdge = mol2.GetNumBonds()
    # (S)ubgraph (B)ond (I)ndex (L)ist(s)
    SBILs = list()
    if NMaxEdge > 0:
        SBILs = Chem.FindAllSubgraphsOfLengthMToN(mol2,NMinEdge,NMaxEdge)
    # valence_based
    if valence_based:
        bonds = mol1.GetBonds()
        # remove surface atom - surface atom bond
        for i in range(len(bonds)-1,-1,-1):
            if bonds[i].GetBeginAtom().GetSymbol() in SurfaceAtomSymbols and \
                bonds[i].GetEndAtom().GetSymbol() in SurfaceAtomSymbols:
                mol1.RemoveBond(bonds[i].GetBeginAtomIdx(),bonds[i].GetEndAtomIdx())
        atoms = mol1.GetAtoms()
        # Separate binding sites
        for i in range(len(atoms)-1,-1,-1):
            if atoms[i].GetSymbol() in SurfaceAtomSymbols:
                neighbors = atoms[i].GetNeighbors()
                neighboridx = list()
                for neighbor_atom in neighbors:
                    neighboridx.append(neighbor_atom.GetIdx())
                    
                if len(neighboridx) > 1: # it has connection to organic atoms
                    for j in range(1,len(neighboridx)):
                        # add valence
                        mol1.RemoveBond(atoms[i].GetIdx(),neighboridx[j])
                        new_surf_atom = Chem.Atom(SurfaceAtomSymbols)
                        new_surf_atom_idx = mol1.AddAtom(new_surf_atom)
                        mol1.AddBond(new_surf_atom_idx,neighboridx[j],Chem.rdchem.BondType.ZERO)
                    
    # debug
    if debug:
        print(Chem.MolToSmiles(mol1))
    # Get finger prints
    smiless = list() # (F)inger(p)rint(s). It's in atom index, and will be
                    # converted to SMILES later. This is done to remove duplicates
    ## length = 0 case
    if NMinEdge ==0:
        for idx in NSAL:
            smiless.append([idx])
    ## Radius > 1 cases
    # add organic atoms
    for SBIL in SBILs:
        for SBI in SBIL:
            AL = set() # (A)tom (L)ist
            for BI in SBI:
                B = mol2.GetBondWithIdx(BI)
                AL.add(Mol2toMol1Map[B.GetBeginAtomIdx()])
                AL.add(Mol2toMol1Map[B.GetEndAtomIdx()])
            smiless.append(list(AL))
    
    # Add surface atom that is neighbor to SMILES
    # add surface atoms
    for smiles in smiless:
        RSAL = list() # (R)elevant (S)urface (A)tom (L)ist
        for SA in SAL:
            for NeighborAtom in mol1.GetAtomWithIdx(SA).GetNeighbors():
                NeighborAtomIdx = NeighborAtom.GetIdx()
                if NeighborAtomIdx in smiles:
                    RSAL.append(SA)
        smiles += RSAL
    # convert atom indexes to SMILES
    SMILES_count = defaultdict(int)
    for smiles in smiless:
        SMILES_count[_GetSMILES(mol1,smiles)] += 1
    return SMILES_count



SurfaceAtomSymbols = ['Ru']
################################ USER INPUT ###################################
# Get SMILES
#fpath = './CONTCAR'
fpath = sys.argv[1]
# Convert atomic coordinates (VASP) to molecular graph.
mole = adsorbate.LoadByCovalentRadius(fpath, SurfaceAtomSymbols = SurfaceAtomSymbols)
# Mine adsorbate graphs and print
print(FindAllAdsorbateGraphsOfLengthN(mole.RdkitMol,SurfaceAtomSymbols,1))
    
