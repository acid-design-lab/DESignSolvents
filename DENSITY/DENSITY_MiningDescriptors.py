# Importing the necessary libraries
import pandas as pd
import numpy as np
import re
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Descriptors3D
from rdkit.Chem import AllChem
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem.rdchem import GetPeriodicTable
from rdkit import Chem
from rdkit.Chem import AllChem



density_df = pd.read_csv('DENSITY/density_without_rep.csv') #Reading the database without repetition

# Loading the database of individual components
ind_comp_df = pd.read_csv('DENSITY/individual_compounds_df_new_names.csv')

def get_molecular_formula(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    formula = Chem.rdMolDescriptors.CalcMolFormula(molecule)
    return formula

### Molecular weight
# Adding a line with the molecular weight of substances
f_get_MW = lambda x: Descriptors.MolWt(Chem.MolFromSmiles(x))
ind_comp_df['Molecular weight'] = ind_comp_df['CanonicalSMILES'].apply(f_get_MW)

### 3D Descriptors
#Moments of inertia
#PMI1
def f_get_PMI1(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    ps = AllChem.EmbedParameters() 
    ps.embedFragmentsSeparately = False
    ps.useRandomCoords = True
    ps.useExpTorsionAnglePrefs=True
    ps.useBasicKnowledge=True
    ps.randomSeed = 42
    flag = AllChem.EmbedMultipleConfs(mol, 10, ps)
    AllChem.MMFFOptimizeMoleculeConfs(mol)
    try:
        res = Descriptors3D.PMI1(mol) #Calculation of descriptors
    except:
        return None
    return res

#PMI2
def f_get_PMI2(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    ps = AllChem.EmbedParameters() 
    ps.embedFragmentsSeparately = False
    ps.useRandomCoords = True
    ps.useExpTorsionAnglePrefs=True
    ps.useBasicKnowledge=True
    ps.randomSeed = 42
    flag = AllChem.EmbedMultipleConfs(mol, 10, ps)
    AllChem.MMFFOptimizeMoleculeConfs(mol)
    try:
        res = Descriptors3D.PMI2(mol) #Calculation of descriptors
    except:
        return None
    return res


#PMI3
def f_get_PMI3(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    ps = AllChem.EmbedParameters() 
    ps.embedFragmentsSeparately = False
    ps.useRandomCoords = True
    ps.useExpTorsionAnglePrefs=True
    ps.useBasicKnowledge=True
    ps.randomSeed = 42
    flag = AllChem.EmbedMultipleConfs(mol, 10, ps)
    AllChem.MMFFOptimizeMoleculeConfs(mol)
    try:
        res = Descriptors3D.PMI3(mol) #Calculation of descriptors
    except:
        return None
    return res


#Applying functions and adding columns
ind_comp_df['PMI1'] = ind_comp_df['IsomericSMILES'].apply(f_get_PMI1)
ind_comp_df['PMI2'] = ind_comp_df['IsomericSMILES'].apply(f_get_PMI2)
ind_comp_df['PMI3'] = ind_comp_df['IsomericSMILES'].apply(f_get_PMI3)


### Structural fragments
# Adding structural fragments
def f_get_fragments(smile):
  mols = [Chem.MolFromSmiles(i) for i in smile] #Getting a list of molecules
  calc = MoleculeDescriptors.MolecularDescriptorCalculator(x[0] for x in Descriptors._descList if re.search(r'fr_', x[0]) is not None) #Selection of descriptors starting with fr_
  desc_names = calc.GetDescriptorNames() #Getting descriptor names
  Mol_descriptors = [] #Blank list to fill with descriptor values
  for mol in mols:
    mol = Chem.AddHs(mol) 
    descriptors = calc.CalcDescriptors(mol) #Calculation of descriptors
    Mol_descriptors.append(descriptors) #Adding descriptors to a list
  return Mol_descriptors, desc_names

Mol_descriptors, desc_names = f_get_fragments(ind_comp_df['IsomericSMILES']) #Getting the result of the function
ind_comp_df = ind_comp_df.join(pd.DataFrame(Mol_descriptors, columns = desc_names)) #Adding to the table


# Function for adding molecular descriptors
def RDKit_descriptors(smiles):
  mols = [Chem.MolFromSmiles(i) for i in smiles] #Getting a list of molecules
  calc = MoleculeDescriptors.MolecularDescriptorCalculator(x[0] for x in Descriptors._descList if x[0] in ['Topological Torsions', 'HeavyAtomCount', 'NumHAcceptors', 'NumHDonors', 'NumHeteroatoms', 'NumRotatableBonds', 'NumValenceElectrons', 'RingCount']) #Selecting descriptors from a sheet_
  desc_names = calc.GetDescriptorNames() #Getting Descriptor Names
  Mol_descriptors = [] #Blank list to fill with descriptor values
  for mol in mols:
    mol = Chem.AddHs(mol) 
    descriptors = calc.CalcDescriptors(mol) #Calculation of descriptors
    Mol_descriptors.append(descriptors) #Adding descriptors to a list
  return Mol_descriptors, desc_names

Mol_descriptors, desc_names = RDKit_descriptors(ind_comp_df['IsomericSMILES'])
ind_comp_df = ind_comp_df.join(pd.DataFrame(Mol_descriptors, columns = desc_names)) #Adding to the table


def Get_Numb_Atoms(smile):
    # We get list with elements and their number
    MolecularFormula = get_molecular_formula(smile) #Finding the molecular formula from smiles
    chars = re.findall(r'[a-zA-Z]+', MolecularFormula) #We select all letter symbols (element designations)
    nums = re.findall(r'\d+', MolecularFormula) #We select all the digits (indexes)
    new_chars = [] #To record elements
    new_nums = [] #For writing indexes
    # Since 1 are not put in molecular formulas, but they need to be taken into account, we will highlight them
    for elem in chars:
        if (sum(i.isupper() for i in elem) >= 2): # Check if there are combinations of several elements, this is indicated by capital letters
            for i in re.sub(r'([A-Z])', r' \1', elem).split(): #We divide the string by capital letters
                new_chars.append(i) #Adding new values
            for _ in range(len(re.sub(r'([A-Z])', r' \1', elem).split()) - 1):
                new_nums.append(1) #Adding the required number of 1
            try:
                new_nums.append(int(nums[chars.index(elem)])) #Adding a number from the main list
            except:
                new_nums.append(1)
        elif len(chars) != len(nums):
            new_chars.append(chars[chars.index(elem)])
            try:
                new_nums.append(int(nums[chars.index(elem)]))
            except:
                new_nums.append(1)
        else:
    # If there are no two large letters, just return the element and the index
            new_chars.append(chars[chars.index(elem)])
            new_nums.append(int(nums[chars.index(elem)]))
    
    #Let's assemble a list of elements
    Elements_List = ['Li', 'C', 'N', 'O', 'F', 'Na', 'Mg', 'Al', 'P', 'S', 'Cl', 'K', 'Ca', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Br']
    Numb_of_elem = []
    for elem in Elements_List:
        elem_ind = new_chars.index(elem) if elem in new_chars else -1
        if elem_ind != -1:
            Numb_of_elem.append(new_nums[elem_ind])
        else:
            Numb_of_elem.append(0)
    
    return Numb_of_elem


#Adding the number of elements to a table with individual components
Elements_List1 = ['Li', 'C', 'N', 'O', 'F', 'Na', 'Mg', 'Al', 'P', 'S', 'Cl', 'K', 'Ca', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Br']
ind_comp_df[Elements_List1] = None
for ind in ind_comp_df.index:
    Numb_of_elem = Get_Numb_Atoms(ind_comp_df.loc[ind, 'IsomericSMILES'])
    for elem in Elements_List1:
        ind_comp_df.loc[ind, elem] = Numb_of_elem[Elements_List1.index(elem)]
ind_comp_df


## Van der Waals radii and volumes
def Get_VdWVolume(smile):
    per_table = GetPeriodicTable() #Loading the periodic system
    mol = Chem.MolFromSmiles(smile) #We get the molecule
    RA = Descriptors.NumAromaticRings(mol) #Number of aromatic rings
    RNA = Descriptors.RingCount(mol) - Descriptors.NumAromaticRings(mol) #Number of non-aromatic rings
    mol_h = Chem.AddHs(mol) #Molecule taking into account H
    ps = AllChem.EmbedParameters() 
    ps.embedFragmentsSeparately = False
    ps.useRandomCoords = True
    flag = AllChem.EmbedMultipleConfs(mol_h, 10, ps)
    AllChem.MMFFOptimizeMolecule(mol_h)
    NB = mol_h.GetNumAtoms() - 1 + RA +  RNA #Number of bonds in a molecule
    # We get lists with elements and their number
    MolecularFormula = get_molecular_formula(smile) #Finding the molecular formula from smiles
    chars = re.findall(r'[a-zA-Z]+', MolecularFormula) #We select all letter symbols (element designations)
    nums = re.findall(r'\d+', MolecularFormula) #We select all the digits (indexes)
    new_chars = [] #To record elements
    new_nums = [] #To record indexes
    # Since 1 are not put in molecular formulas, but they need to be taken into account, we will highlight them
    for elem in chars:
        if (sum(i.isupper() for i in elem) >= 2): # Check if there are combinations of several elements, this is indicated by capital letters
            for i in re.sub(r'([A-Z])', r' \1', elem).split(): #We divide the string by capital letters
                new_chars.append(i) #Adding new values
            for _ in range(len(re.sub(r'([A-Z])', r' \1', elem).split()) - 1):
                new_nums.append(1) #Adding the required number of 1
            try:
                new_nums.append(int(nums[chars.index(elem)])) #Adding a number from the main list
            except:
                new_nums.append(1)
        elif len(chars) != len(nums):
            new_chars.append(chars[chars.index(elem)])
            try:
                new_nums.append(int(nums[chars.index(elem)]))
            except:
                new_nums.append(1)
        else:
    # If there are no two large letters, just return the element and the index
            new_chars.append(chars[chars.index(elem)])
            new_nums.append(int(nums[chars.index(elem)]))
  # Find the sum of the contributions of atoms
        sum_atom_contributions = 0
        for i in range(len(new_chars)):
            radius = per_table.GetRvdw(new_chars[i]) #V-D-W radius
            sum_atom_contributions += new_nums[i] * 4 / 3 * np.pi * radius**3
        # The basic formula for the calculation
        VdWVolume = sum_atom_contributions - 5.98*NB - 14.7*RA - 3.8*RNA
        return VdWVolume
    

ind_comp_df['VdWVolume, A^3'] = ind_comp_df['IsomericSMILES'].apply(Get_VdWVolume) #Applying a function to a dataframe


density_df = density_df.drop(['Unnamed: 0', 'Unnamed: 0.1'], axis = 1)
Elements_List1 = ['Li', 'C', 'N', 'O', 'F', 'Na', 'Mg', 'Al', 'P', 'S', 'Cl', 'K', 'Ca', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Br']
ind_comp_df = ind_comp_df.drop_duplicates(['IsomericSMILES']) #Removing duplicate rows from the database if there are any
ind_comp_df.index = ind_comp_df['IsomericSMILES'] #We use IsomericSMILES as an index
Descr_list =  ['Molecular weight',
       'NumValenceElectrons', 'HeavyAtomCount', 'NumHAcceptors', 'NumHDonors',
       'NumHeteroatoms', 'NumRotatableBonds', 'RingCount', 'VdWVolume, A^3', 'PMI1', 'PMI2', 'PMI3'] + Elements_List1 #List of descriptors to be added
for desc in Descr_list:
  f_get_desc = lambda x: ind_comp_df.loc[x][desc] if isinstance(x, str) else 0 #Function for getting descriptors from the table for individual substances
  #There is a separate column for each component of the system
  for num_comp in range(3):
    name_new_column = desc + '#' + str(num_comp + 1)
    name_old_column = 'isomer_smiles' + '#' + str(num_comp + 1)
    density_df[name_new_column] = density_df[name_old_column].apply(f_get_desc)


# Replace None in column X#3 (molar fraction) with 0 for correct calculations
density_df['X#3 (molar fraction)'] = density_df['X#3 (molar fraction)'].replace(np.nan, 0)
Desc_list_new = ['Molecular weight', 'HeavyAtomCount', 'NumHAcceptors', 'NumHDonors', 'NumHeteroatoms', 'NumRotatableBonds', 'RingCount']+ Elements_List1 #List of descriptors
#Averaging descriptors
for desc in Desc_list_new:
  density_df[desc] = density_df[desc + '#' + '1'] * density_df['X#1 (molar fraction)'] + density_df[desc + '#' + '2'] * density_df['X#2 (molar fraction)'] + density_df[desc + '#' + '3'] * density_df['X#3 (molar fraction)']

density_df


#InertialShapeFactor
def f_get_InertialShapeFactor_gen():
    #Expressions for moments of inertia
    pm1 = density_df['PMI1#1'] * density_df['X#1 (molar fraction)'] + density_df['PMI1#2'] * density_df['X#2 (molar fraction)'] + density_df['PMI1#3'] * density_df['X#3 (molar fraction)']
    pm2 = density_df['PMI2#1'] * density_df['X#1 (molar fraction)'] + density_df['PMI2#2'] * density_df['X#2 (molar fraction)'] + density_df['PMI2#3'] * density_df['X#3 (molar fraction)']
    pm3 = density_df['PMI3#1'] * density_df['X#1 (molar fraction)'] + density_df['PMI3#2'] * density_df['X#2 (molar fraction)'] + density_df['PMI3#3'] * density_df['X#3 (molar fraction)']
    #Formula for InertialShapeFactor
    return pm2 / (pm1*pm3)


#SpherocityIndex
def f_get_SpherocityIndex_gen():
    #Expressions for moments of inertia
    pm1 = density_df['PMI1#1'] * density_df['X#1 (molar fraction)'] + density_df['PMI1#2'] * density_df['X#2 (molar fraction)'] + density_df['PMI1#3'] * density_df['X#3 (molar fraction)']
    pm2 = density_df['PMI2#1'] * density_df['X#1 (molar fraction)'] + density_df['PMI2#2'] * density_df['X#2 (molar fraction)'] + density_df['PMI2#3'] * density_df['X#3 (molar fraction)']
    pm3 = density_df['PMI3#1'] * density_df['X#1 (molar fraction)'] + density_df['PMI3#2'] * density_df['X#2 (molar fraction)'] + density_df['PMI3#3'] * density_df['X#3 (molar fraction)']
    #Formula for SpherocityIndex
    return 3 * pm1 / (pm1+pm2+pm3)



density_df['InertialShapeFactor'] = f_get_InertialShapeFactor_gen()
density_df['SpherocityIndex'] = f_get_SpherocityIndex_gen()



# Adding structural descriptors to the main table
Descr_list =  [x for x in ind_comp_df.columns if re.search(r'fr_', x) is not None] #List of descriptors
for desc in Descr_list:
  f_get_desc = lambda x: ind_comp_df.loc[x][desc] if isinstance(x, str) else 0 #The function of adding a descriptor
  #Adding to the table with the design of the name
  for num_comp in range(3):
    name_new_column = desc + '#' + str(num_comp + 1)
    name_old_column = 'isomer_smiles' + '#' + str(num_comp + 1)
    density_df[name_new_column] = density_df[name_old_column].apply(f_get_desc)


# Let's calculate the average by descriptors
density_df['X#3 (molar fraction)'] = density_df['X#3 (molar fraction)'].replace(np.nan, 0)
Desc_list_new = [x for x in ind_comp_df.columns if re.search(r'fr_', x) is not None]
for desc in Desc_list_new:
  density_df[desc] = density_df[desc + '#' + '1'] * density_df['X#1 (molar fraction)'] + density_df[desc + '#' + '2'] * density_df['X#2 (molar fraction)'] + density_df[desc + '#' + '3'] * density_df['X#3 (molar fraction)']


#Function for getting a list of elements and their number in a molecule
def Get_Numb_Atoms(smile):
    # We get sheets with elements and their number
    MolecularFormula = get_molecular_formula(smile) #Finding the molecular formula from smiles
    chars = re.findall(r'[a-zA-Z]+', MolecularFormula) #We select all letter symbols (element designations)
    nums = re.findall(r'\d+', MolecularFormula) #We select all the digits (indexes)
    new_chars = [] #To record elements
    new_nums = [] #To writing indexes
    # Since 1 are not put in molecular formulas, but they need to be taken into account, we will highlight them
    for elem in chars:
        if (sum(i.isupper() for i in elem) >= 2): # Check if there are combinations of several elements, this is indicated by capital letters
            for i in re.sub(r'([A-Z])', r' \1', elem).split(): #We divide the string by capital letters
                new_chars.append(i) #Adding new values
            for _ in range(len(re.sub(r'([A-Z])', r' \1', elem).split()) - 1):
                new_nums.append(1) #Adding the required number of 1
            try:
                new_nums.append(int(nums[chars.index(elem)])) #Adding a number from the main list
            except:
                new_nums.append(1)
        elif len(chars) != len(nums):
            new_chars.append(chars[chars.index(elem)])
            try:
                new_nums.append(int(nums[chars.index(elem)]))
            except:
                new_nums.append(1)
        else:
    # If there are no two large letters, just return the element and the index
            new_chars.append(chars[chars.index(elem)])
            new_nums.append(int(nums[chars.index(elem)]))
    return new_chars, new_nums


def Get_mass_fraction_metal(smiles):
    try:
        per_table = GetPeriodicTable()
        
        Metal_list = ['Li', 'Be', 'Na', 'Mg', 'Al', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 
                        'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ha', 'Rb', 'Sr', 'Y', 'Zr', 
                        'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Cs', 
                        'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 
                        'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 
                        'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi']

        MetAt = 0
        new_chars, new_nums = Get_Numb_Atoms(smiles)
        for elem in new_chars:
            if elem in Metal_list:
                MetAt += new_nums[new_chars.index(elem)] * per_table.GetAtomicWeight(elem)
    except:
        MetAt = 0

    return MetAt


density_df['Metal_frac_gen'] = (density_df['X#1 (molar fraction)']*density_df['isomer_smiles#1'].apply(Get_mass_fraction_metal)+density_df['X#2 (molar fraction)']*density_df['isomer_smiles#2'].apply(Get_mass_fraction_metal)+density_df['X#3 (molar fraction)']*density_df['isomer_smiles#3'].apply(Get_mass_fraction_metal))/(density_df['X#1 (molar fraction)']*density_df['Molecular weight#1']+density_df['X#2 (molar fraction)']*density_df['Molecular weight#2']+density_df['X#3 (molar fraction)']*density_df['Molecular weight#3'])
    

density_df_ML = density_df[['Component#1', 'Component#2', 'Component#3', 'X#1 (molar fraction)', 'X#2 (molar fraction)', 'X#3 (molar fraction)','I', 'II', 'III', 'IV', 'V', 
                                     'Temperature, K', 'Density, g/cm^3','isomer_smiles#1','isomer_smiles#2','isomer_smiles#3','VdWVolume, A^3#1','VdWVolume, A^3#2','VdWVolume, A^3#3',
                                     'NumHeteroatoms','RingCount', 'InertialShapeFactor', 'SpherocityIndex', 'Metal_frac_gen'] + Elements_List1]

density_df_ML.to_csv('density_df_with_elem.csv')
density_df_ML