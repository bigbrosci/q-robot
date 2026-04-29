#!/usr/bin/env python3
import os
dict_element = {'1':'H','2':'He','3':'Li','4':'Be','5':'B','6':'C','7':'N','8':'O','9':'F','10':'Ne','11':'Na','12':'Mg','13':'Al','14':'Si','15':'P','16':'S','17':'Cl','18':'Ar','19':'K','20':'Ca','21':'Sc','22':'Ti','23':'V','24':'Cr','25':'Mn','26':'Fe','27':'Co','28':'Ni','29':'Cu','30':'Zn','31':'Ga','32':'Ge','33':'As','34':'Se','35':'Br','36':'Kr','37':'Rb','38':'Sr','39':'Y','40':'Zr','41':'Nb','42':'Mo','43':'Tc','44':'Ru','45':'Rh','46':'Pd','47':'Ag','48':'Cd','49':'In','50':'Sn','51':'Sb','52':'Te','53':'I','54':'Xe','55':'Cs','56':'Ba','57':'La','58':'Ce','59':'Pr','60':'Nd','61':'Pm','62':'Sm','63':'Eu','64':'Gd','65':'Tb','66':'Dy','67':'Ho','68':'Er','69':'Tm','70':'Yb','71':'Lu','72':'Hf','73':'Ta','74':'W','75':'Re','76':'Os','77':'Ir','78':'Pt','79':'Au','80':'Hg','81':'Tl','82':'Pb','83':'Bi','84':'Po','85':'At','86':'Rn','87':'Fr','88':'Ra','89':'Ac','90':'Th','91':'Pa','92':'U','93':'Np','94':'Pu','95':'Am','96':'Cm','97':'Bk','98':'Cf','99':'Es','100':'Fm','101':'Md','102':'No','103':'Lr','104':'Rf','105':'Db','106':'Sg','107':'Bh','108':'Hs','109':'Mt',}

#### Part-0 Dictionaries
atomic_mass = dict(H=1.01, He=4.00, Li=6.94, Be=9.01, B=10.81, C=12.01,
                   N=14.01, O=16.00, F=19.00, Ne=20.18, Na=22.99, Mg=24.31,
                   Al=26.98, Si=28.09, P=30.97, S=32.07, Cl=35.45, Ar=39.95,
                   K=39.10, Ca=40.08, Sc=44.96, Ti=47.87, V=50.94, Cr=52.00,
                   Mn=54.94, Fe=55.85, Co=58.93, Ni=58.69, Cu=63.55, Zn=65.39,
                   Ga=69.72, Ge=72.61, As=74.92, Se=78.96, Br=79.90, Kr=83.80,
                   Rb=85.47, Sr=87.62, Y=88.91, Zr=91.22, Nb=92.91, Mo=95.94,
                   Tc=98.00, Ru=101.07, Rh=102.91, Pd=106.42, Ag=107.87,
                   Cd=112.41, In=114.82, Sn=118.71, Sb=121.76, Te=127.60,
                   I=126.90, Xe=131.29, Cs=132.91, Ba=137.33, La=138.91,
                   Ce=140.12, Pr=140.91, Nd=144.24, Pm=145.00, Sm=150.36,
                   Eu=151.96, Gd=157.25, Tb=158.93, Dy=162.50, Ho=164.93,
                   Er=167.26, Tm=168.93, Yb=173.04, Lu=174.97, Hf=178.49,
                   Ta=180.95, W=183.84, Re=186.21, Os=190.23, Ir=192.22,
                   Pt=195.08, Au=196.97, Hg=200.59, Tl=204.38, Pb=207.2,
                   Bi=208.98, Po=209.00, At=210.00, Rn=222.00, Fr=223.00,
                   Ra=226.00, Ac=227.00, Th=232.04, Pa=231.04, U=238.03,
                   Np=237.00, Pu=244.00, Am=243.00, Cm=247.00, Bk=247.00,
                   Cf=251.00, Es=252.00, Fm=257.00, Md=258.00, No=259.00,
                   Lr=262.00, Rf=261.00, Db=262.00, Sg=266.00, Bh=264.00,
                   Hs=269.00, Mt=268.00)

### VDW Parameters fro DFT_D2 method
c6dict = {
'Au': '7.308', 'Ni': '2.6263', 'Ru': '4.1678',
'Cu': '2.740', 'O' : '0.7000', 'Si': '9.2300',
'H' : '0.140', 'C' : '1.7500', 'B' : '3.130',
'N' : '1.230', 'Br': '12.470', 'S' : '5.570',
'Pd': '5.510', 'Pt': '7.0000', 'Cl': '6.070',
'Ag': '5.481', 'Co': '2.5650', 'Rh': '4.364',
'Ir': '6.163', 'Os': '5.8780', 'P': '7.84',
'Fe': '2.600',
} 
     
r0dict = {
'Au': '1.823', 'Ni': '1.562', 'Ru': '1.639',
'Cu': '1.562', 'O' : '1.342', 'Si': '1.716',
'H': '1.001',  'C' : '1.452', 'B': '1.485',
'N': '1.397',  'Br': '1.749', 'S': '1.683',
'Pd': '1.69',  'Pt': '1.75',  'Cl': '1.639',
'Ag': '1.819', 'Co': '1.349', 'Rh': '1.677',
'Ir': '1.698', 'Os': '1.504', 'P': '1.705',
'Fe':'1.40',
} 

### Convert CH3OH to 1412 
dict_l = {'1312' : 'HOCH2', '1300' : 'CH3', '1212' : 'HOCH', '1211' : 'OCH2',
          '1200' : 'CH2', '1112' : 'HOC', '1111' : 'OCH', '1100' : 'CH', 
          '1011' : 'OC', '1000' : 'C', '0000' : '0000'}
dict_r = {'1312' : 'CH2OH', '1300' : 'CH3', '1212' : 'CHOH', '1211' : 'CH2O', 
          '1200' : 'CH2', '1112' : 'COH', '1111' : 'CHO', '1100' : 'CH', 
          '1011' : 'CO', '1000' : 'C', '0000' : '0000'}

#dict_r_r = {y:x for x,y in dict_r.items()}

## Crystal structure of elements: from https://en.wikipedia.org/wiki/Periodic_table_(crystal_structure)
bcc = ['V',  'Cr', 'Mn', 'Fe', 'Nb', 'Pb']
hcp = ['Mg', 'Sc', 'Ti', 'Co', 'Zn', 'Y', 'Zr', 'Tc', 'Ru', 'Cd', 'Hf', 'Re', 'Os']
fcc = ['Al', 'Ca', 'Ni', 'Cu', 'Rh', 'Pd', 'Ag', 'Ir', 'Pt', 'Au']

### Metal  Bulk structures from DFT_D2 
### More information: data_base/metal_bulks 
dict_metals_vdWNo           = {
'Ag':(-10.87635094,4,4.1446,4.1446),
'Co':(-14.07386001,2,2.4924,4.0221),
'Cu':(-14.9138106,4,3.6361,3.6361),
'Fe':(-16.47275657,2,2.8317,2.8317),
'Ir':(-35.40904766,4,3.8713,3.8713),
'Ni':(-21.8703144,4,3.5167,3.5167),
'Os':(-22.50851382,2,2.7528,4.3496),
'Pd':(-20.86566531,4,3.9337,3.9337),
'Pt':(-24.40310988,4,3.9656,3.9656),
'Rh':(-29.11414683,4,3.8227,3.8227),
'Ru':(-18.51185829,2,2.7142,4.2856),
}

dict_metals_vdWD2              = {
'Ag':(-11.23613381,4,4.126,4.126),
'Co':(-14.4547676,2,2.4968,4.0249),
'Cu':(-15.32521363,4,3.6202,3.6202),
'Fe':(-16.82023298,2,2.8323,2.8323),
'Ir':(-35.40904766,4,3.8713,3.8713),
'Ni':(-22.32245708,4,3.5004,3.5004),
'Os':(-22.8652351,2,2.7516,4.3534),
'Pd':(-21.37910959,4,3.9262,3.9262),
'Pt':(-24.99443446,4,3.9533,3.9533),
'Rh':(-29.5771939,4,3.818,3.818),
'Ru':(-18.7500579,2,2.7086,4.2772),
}

dict_metals_vdWD3zero           = {
'Ag':(-12.83721745,4,4.0711,4.0711),
'Co':(-14.76776158,2,2.4729,3.9814),
'Cu':(-16.95007138,4,3.5683,3.5683),
'Fe':(-17.09658076,2,2.8068,2.8068),
'Ir':(-38.24513552,4,3.8377,3.8377),
'Ni':(-23.4791717,4,3.4744,3.4744),
'Os':(-23.64822617,2,2.7331,4.3278),
'Pd':(-23.17971002,4,3.8822,3.8822),
'Pt':(-27.46886041,4,3.9171,3.9171),
'Rh':(-31.41993461,4,3.7845,3.7845),
'Ru':(-19.62094079,2,2.6903,4.2535),
} 

dict_metals_vdWD3BJ           = {
'Ag':(-13.19451205,4,4.0705,4.0705),
'Co':(-14.92272042,2,2.4681,3.9901),
'Cu':(-17.27528277,4,3.5675,3.5675),
'Fe':(-17.19762711,2,2.8079,2.8079),
'Ir':(-38.30745938,4,3.8429,3.8429),
'Ni':(-23.75578398,4,3.475,3.475),
'Os':(-23.76187666,2,2.7346,4.3281),
'Pd':(-23.50021344,4,3.8853,3.8853),
'Pt':(-27.46653366,4,3.9264,3.9264),
'Rh':(-31.69505713,4,3.7863,3.7863),
'Ru':(-19.75899154,2,2.6903,4.2547),
}

dict_metals_vdWoptB86b           = {
'Ag':(0.84421148,4,4.0901,4.0901),
'Co':(-9.30800993,2,2.4691,3.9848),
'Cu':(-3.71482413,4,3.6011,3.6011),
'Fe':(-11.88398153,2,2.807,2.807),
'Ir':(-26.41771245,4,3.8606,3.8606),
'Ni':(-11.5764003,4,3.4893,3.4893),
'Os':(-18.08068787,2,2.7463,4.3393),
'Pd':(-9.75241654,4,3.9032,3.9032),
'Pt':(-15.32228447,4,3.9468,3.9468),
'Rh':(-19.68788189,4,3.8044,3.8044),
'Ru':(-13.98875677,2,2.7016,4.2659),
}

dict_metals_vdWoptB88           = {
'Ag':(1.40470723,4,4.1283,4.1283),
'Co':(-8.77101385,2,2.4848,4.0093),
'Cu':(-2.95489179,4,3.6271,3.6271),
'Fe':(-11.29492318,2,2.8219,2.8219),
'Ir':(-24.71850215,4,3.8838,3.8838),
'Ni':(-10.64867601,4,3.5095,3.5095),
'Os':(-17.10995344,2,2.7612,4.3628),
'Pd':(-8.88228611,4,3.9337,3.9337),
'Pt':(-14.03984804,4,3.9779,3.9779),
'Rh':(-18.431108,4,3.8292,3.8292),
'Ru':(-13.22362589,2,2.7186,4.2924),
}

dict_metals_vdWoptPBE          = { 
'Ag':(2.34929913,4,4.1635,4.1635),
'Co':(-8.31527839,2,2.497,4.0299),
'Cu':(-2.11142609,4,3.6523,3.6523),
'Fe':(-10.88720327,2,2.8362,2.8362),
'Ir':(-23.65016927,4,3.8895,3.8895),
'Ni':(-9.71082594,4,3.5305,3.5305),
'Os':(-16.60709913,2,2.7643,4.3684),
'Pd':(-7.78402982,4,3.9514,3.9514),
'Pt':(-12.97109645,4,3.9896,3.9896),
'Rh':(-17.26723459,4,3.8384,3.8384),
'Ru':(-12.64376043,2,2.7247,4.303),
}

dict_metals_vdWDF           = {
'Ag':(4.0353394,4,4.2423,4.2423),
'Co':(-7.25596362,2,2.5256,4.0764),
'Cu':(-0.37256033,4,3.7055,3.7055),
'Fe':(-9.84560669,2,2.8679,2.8679),
'Ir':(-20.72902652,4,3.9242,3.9242),
'Ni':(-7.72319436,4,3.5737,3.5737),
'Os':(-15.07593176,2,2.7827,4.3982),
'Pd':(-5.56339295,4,4.0099,4.0099),
'Pt':(-10.42459996,4,4.0316,4.0316),
'Rh':(-14.70909191,4,3.8802,3.8802),
'Ru':(-11.2466035,2,2.7463,4.3373),
}

dict_metals_vdWDF2           = {
'Ag':(2.46205389,4,4.3095,4.3095),
'Co':(-7.41526434,2,2.5486,4.1115),
'Cu':(-1.50771134,4,3.7476,3.7476),
'Fe':(-9.84028843,2,2.8894,2.8894),
'Ir':(-19.86383932,4,3.9863,3.9863),
'Ni':(-8.48079358,4,3.6087,3.6087),
'Os':(-14.23269719,2,2.8207,4.4615),
'Pd':(-6.59190857,4,4.0843,4.0843),
'Pt':(-10.49169121,4,4.1092,4.1092),
'Rh':(-14.73865299,4,3.941,3.941),
'Ru':(-10.93941608,2,2.7854,4.4018),
}

dict_metals_vdWrevDF2           = { 
'Ag':(-0.27916318,4,4.0971,4.0971),
'Co':(-9.76906381,2,2.471,3.9878),
'Cu':(-4.85812653,4,3.6066,3.6066),
'Fe':(-12.3206222,2,2.8075,2.8075),
'Ir':(-27.26723233,4,3.8611,3.8611),
'Ni':(-12.61884988,4,3.4901,3.4901),
'Os':(-18.48565924,2,2.7464,4.3394),
'Pd':(-10.78990441,4,3.9073,3.9073),
'Pt':(-16.20871775,4,3.9471,3.9471),
'Rh':(-20.59614884,4,3.8063,3.8063),
'Ru':(-14.41463483,2,2.7022,4.2667),
}

dict_metals_vdWSCAN          = {
'Ag':(-124.65340237,4,4.0546,4.0546),
'Co':(-34.15820236,2,2.441038,4.00374),
'Cu':(-59.77914967,4,3.544,3.544),
'Fe':(-35.94045354,2,2.83662,2.83662),
'Ir':(-284.84888428,4,3.78252,3.78252),
'Ni':(-63.82323938,4,3.456518,3.456518),
'Os':(-143.51450486,2,2.713,4.2823),
'Pd':(-130.41067436,4,3.8763,3.8763),
'Pt':(-280.65462449,4,3.8881,3.8881),
'Rh':(-134.82854838,4,3.7675,3.7675),
'Ru':(-69.65946843,2,2.6718,4.2149),
}

dict_metals_vdWTS2            = {
'Ag':(-12.66940228,4,4.0668,4.0668),
'Co':(-15.71926358,2,2.4658,3.9443),
'Cu':(-17.34038071,4,3.5446,3.5446),
'Fe':(-18.29812446,2,2.7774,2.7774),
'Ir':(-37.63483518,4,3.8423,3.8423),
'Ni':(-25.8389368,4,3.423,3.423),
'Os':(-22.50646554,2,2.7528,4.3496),
'Pd':(-21.93008933,4,3.9098,3.9098),
'Pt':(-26.40654184,4,3.9656,3.9656),
'Rh':(-32.29610409,4,3.7663,3.7663),
'Ru':(-20.63176818,2,2.662,4.2324),
}

dict_metals_vdWTS21           = {
'Ag':(-12.63336341,4,4.0702,4.0702),
'Co':(-15.63272259,2,2.4641,3.9621),
'Cu':(-17.28859311,4,3.5471,3.5471),
'Fe':(-18.23393479,2,2.7822,2.7822),
'Ir':(-37.58869156,4,3.8427,3.8427),
'Pd':(-21.90946451,4,3.9112,3.9112),
'Pt':(-26.39090267,4,3.9335,3.9335),
'Rh':(-32.22941813,4,3.7671,3.7671),
'Ru':(-20.58541811,2,2.6631,4.2355),
}

dict_metals_vdWTS4            = {
'Ag':(-12.19633179,4,4.0967,4.0967),
'Co':(-14.78174354,2,2.4718,3.9912),
'Cu':(-16.17623115,4,3.5939,3.5939),
'Fe':(-17.19163894,2,2.806,2.806),
'Ir':(-36.86481076,4,3.8551,3.8551),
'Ni':(-23.21815361,4,3.4858,3.4858),
'Pd':(-22.0857718,4,3.912,3.912),
'Pt':(-25.88436299,4,3.944,3.944),
'Rh':(-30.53944471,4,3.8012,3.8012),
'Ru':(-19.2539375,2,2.6999,4.2669),
}

######  DFT + U parameters 
## The user need to modify it by him/herself  <<<  VERT IMPORTTANT!!! 
u_value = {
'Ag': 5.0,
'Cu': 5.0,
'Fe': 5.0,
'Ir': 5.0,
'Ni': 5.0,
'Pd': 5.0,
'Pt': 5.0,
'Rh': 5.0,
'Co': 5.0,
'Ru': 5.0,
'Os': 5.0,        
'Au': 5.0,        
'Ti': 5.1,
'Zn': 5.0,
'Sn': 5.0,
}

j_value = {
'Ag': 1.0,
'Cu': 1.0,
'Fe': 1.0,
'Ir': 1.0,
'Ni': 1.0,
'Pd': 1.0,
'Pt': 1.0,
'Rh': 1.0,
'Co': 1.0,
'Ru': 1.0,
'Os': 1.0,
'Au': 1.0,
'Ti': 1.0,        
'Zn': 1.0,        
'Sn': 1.0,        
}

#### Information about the potcar databse 
def get_potcar_data():
    home = os.path.expanduser('~')
    data_potcars_file = home + '/bin/q-robot/books/potpaw_PBE.52/data_potcars'
    file_in = open(data_potcars_file, 'r')
    data = file_in.read()
    file_in.close()
    return eval(data)    

### Group User ID in BSC
dict_id_nl_bsc = {
'Qiang_Li':'iciq72010',
}

list_id_nl_tekla = ['qli']


ja = {
'Accounts of Chemical Research':'Acc. Chem. Res.',
'ACS Applied Materials and Interfaces':'ACS Appl. Mater. Interfaces',
'ACS Applied Materials and Interfaces':'ACS Appl. Mater. Interfaces',
'Advances in Catalysis':'Adv. Catal.',
'AIP Conference Proceedings':'AIP Conf. Proc.',
'Angewandte Chemie International Edition':'Angew. Chem. Int. Ed.',
'Annual Review of Plant Biology':'Annu. Rev. Plant Biol.',
'Applied Surface Science':'Appl. Surf. Sci.',
'Catalysis Communications':'Catal. Commun.',
'Catalysis Letters':'Catal. Lett.',
'Catalysis Today':'Catal. Today',
'Chemical Communications':'Chem. Commun.',
'Chemical Engineering and Technology':'Chem. Eng. Technol.',
'Chemical Engineering Science':'Chem. Eng. Sci.',
'Chemical Physics Letters':'Chem. Phys. Lett.',
'Chemical Reviews':'Chem. Rev.',
'Chemical Science':'Chem. Sci.',
'Chemical Society Reviews':'Chem. Soc. Rev.',
'ChemPhysChem':'ChemPhysChem',
'ChemSusChem':'ChemSusChem',
'Computational Materials Science':'Comput. Mater. Sci.',
'Dalton Transactions':'Dalton Trans.',
'Energy and Environmental Science':'Energy Environ. Sci.',
'Energy and Fuels':'Energy Fuels',
'Faraday Discussions':'Faraday Discuss.',
'Green Chemistry':'Green Chem.',
'Industrial and Engineering Chemistry Research':'Ind. Eng. Chem. Res.',
'International Journal of Hydrogen Energy':'Int. J. Hydrogen Energy',
'Journal of physics. Condensed matter : an Institute of Physics journal':'J. Phys. Condens. Matter',
'Journal of Alloys and Compounds':'J. Alloys Compd.',
'Journal of Applied Polymer Science':'J. Appl. Polym. Sci.',
'Journal of Catalysis':'J. Catal.',
'Journal of Chemical Information and Modeling':'J. Chem. Inf. Model.',
'Journal of Chemical Physics':'J. Chem. Phys.',
'Journal of Chemical Theory and Computation':'J. Chem. Theory Comput.',
'Journal of Computational Chemistry':'J. Comput. Chem.',
'Journal of electroanalytical chemistry and interfacial electrochemistry':'J. Electroanal. Chem.',
'Journal of Electroanalytical Chemistry':'J. Electroanal. Chem.',
'Journal of Molecular Catalysis A: Chemical':'J. Mol. Catal. A: Chem.',
'Journal of Physical Chemistry A':'J. Phys. Chem. A',
'Journal of Physical Chemistry B':'J. Phys. Chem. B',
'Journal of Physical Chemistry C':'J. Phys. Chem. C',
'Journal of Physical Chemistry':'J. Phys. Chem.',
'The Journal of Physical Chemistry':'J. Phys. Chem.',
'Journal of Physical Chemistry. A':'J. Phys. Chem. A',
'Journal of Physical Chemistry. B':'J. Phys. Chem. B',
'Journal of Physics and Chemistry of Solids':'J. Phys. Chem. Solids',
'Journal of Power Sources':'J. Power Sources',
'Journal of the American Chemical Society':'J. Am. Chem. Soc.',
'Journal of Molecular Modeling':'J. Mol. Model.',
'Langmuir':'Langmuir',
'Nanotechnology':'Nanotechnology',
'Nature Communications':'Nat. Commun.',
'Nature Materials':'Nat. Mater.',
'Nature':'Nature',
'Nature Catalysis': 'Nat. Catal.',
'Organic and Biomolecular Chemistry':'Org. Biomol. Chem.',
'Organic and Biomolecular Chemistry':'Org. Biomol. Chem.',
'Physical Chemistry Chemical Physics':'Phys. Chem. Chem. Phys.',
'Physical Review B':'Phys. Rev. B',
'Physical Review Letters':'Phys. Rev. Lett.',
'RSC Advances':'RSC Adv.',
'Science':'Science',
'Surface Science Reports':'Surf. Sci. Rep.',
'Surface Science':'Surf. Sci.',
'Tetrahedron':'Tetrahedron',
'The Journal of Chemical Physics':'J. Chem. Phys.',
'The Journal of Physical Chemistry B':'J. Phys. Chem. B',
'The Journal of Physical Chemistry C':'J. Phys. Chem. C',
'THEOCHEM':'THEOCHEM',
'Topics in Catalysis':'Top. Catal.',
'Zeitschrift fuer Anorganische und Allgemeine Chemie':'Z. Anorg. Allg. Chem.',
'Platinum Metals Review':'Platin. Met. Rev.',
'Theochem':'Theochem',
'ACS Central Science':'ACS Cent. Sci.',
'Physical Review':'Phys. Rev.',
'Applied catalysis B Environmental':'Appl. Catal. B',
'{Applied Catalysis A: General': 'Appl. Catal. A',
'Science Advances':'Sci. Adv.',
'Sensors and Actuators B: Chemical':'Sens. Actuator B-Chem.',
'Catalysis Science and Technology':'Catal. Sci. Tech.',
'Zeitschrift für Physikalische Chemie':'Z. Phys. Chem.-Stoch. Ve.',
'Journal of Industrial and Engineering Chemistry':'J. Ind. Eng. Chem.',
'Renewable and Sustainable Energy Reviews':'Rev. Mod. Phys.',
'Journal of Physics F: Metal Physics':'J. Phys. F: Metal Phys.',
'Transactions of the Faraday Society':'Trans. Faraday Soc.',
'The Chemical Educator':'Chem. Educator',
'Biofuels, Bioproducts and Biorefining':'Biofuel. Bioprod. Bioref.',
'ACS Catalysis':'ACS Catal.',
'Proceedings of the Royal Society of London. Series A, Mathematical and physical sciences':'Proc. R. Soc. Lond. Math. Phys. Sci.',
'Industrial crops and products':'Ind. Crops Prod.',
'Catalysts':'Catalysts',
'Nature Chemistry':'Nat. Chem.',
'ChemCatChem':'ChemCatChem',
'Applied catalysis. A General':'Appl. Catal. A Gen.',
'Chemistry: A European Journal':'Chem.: Eur. J.',
} 

#f_out = open('/home/robot/bin/Q_robot/books/books_have_read/my_ja_dict.txt', 'w')
#f_out.write(str(ja))
#f_out.close()

       
BibEntries = [
'article',
'book',
'booklet',
'inbook',
'incollection',
'inproceedings',
'manual',
'mastersthesis',
'misc',
'phdthesis',
'proceedings',
'techreport',
'unpublished',
]

BibFields = [
'address',
'annote',
'author',
'booktitle',
'chapter',
'crossref',
'edition',
'editor',
'howpublished',
'institution',
'journal',
'key',
'month',
'note',
'number',
'organization',
'pages',
'publisher',
'school',
'series',
'title',
'type',
'volume',
'year',
'abstract',
'doi',
'keywords',
'ISSN',
]

def save_dict_csv(dict_in, name):
    import csv
    out_name = name + '.csv'
    w = csv.writer(open(out_name, 'w'))
    for key, val in dict.items():
        w.writerow([key,val]) 

def save_dict_json(dict_in, name):
    import json 
    out_name = name + '.json'
    json = json.dumps(dict_in)
    f = open(out_name, 'w')
    f.write(json)
    f.close()

def save_dict_txt(dict_in, name):
    out_name = name + '.txt'
    f = open(out_name, 'w')
    f.write(str(dict_in))
    f.close()


def eval_dict_txt(file_in):
    dict_txt  = eval(open(file_in).read())
    return dict_txt

def eval_dict_csv(file_in):
    import csv
    dict_csv = {}
    with open('example.csv') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
            dict_csv.update({row[0]:row[1:]})
    return dict_csv

def eval_dict_json(file_in):
    import json
    with open(file_in) as f:
      dict_json = json.load(f)
    return dict_json


