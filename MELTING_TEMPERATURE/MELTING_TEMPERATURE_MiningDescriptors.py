import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit import Chem
from rdkit.Chem import AllChem

df = pd.read_csv('Melting_temperature.csv') #Reading the database without repetition

# Loading the database of individual components
individual_compounds_df = pd.unique(df[['Smiles#1', 'Smiles#2']].values.ravel()) #List of unique substances
individual_compounds_df = pd.DataFrame(individual_compounds_df) #Transfer to dataframe
individual_compounds_df.columns = ['Compound'] # Renaming a column


def get_mols(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    ps = AllChem.EmbedParameters()
    ps.embedFragmentsSeparately = False
    ps.useRandomCoords = True
    ps.useExpTorsionAnglePrefs=True
    ps.useBasicKnowledge=True
    ps.randomSeed = 42
    flag = AllChem.EmbedMultipleConfs(mol, 10, ps)
    AllChem.MMFFOptimizeMolecule(mol,maxIters=5000)
    return mol


 
 # Function for adding molecular descriptors
def RDKit_descriptors(smiles):
  mols = [Chem.MolFromSmiles(i) for i in smiles] #Getting a list of molecules
  calc = MoleculeDescriptors.MolecularDescriptorCalculator(x[0] for x in Descriptors._descList if x[0] in ['fr_Al_COO', 'fr_Ar_COO', 'fr_Ar_N','fr_Ar_OH', 'fr_NH0', 'fr_NH1', 'fr_amide']) #Selecting descriptors from a sheet
  desc_names = calc.GetDescriptorNames() #Getting Descriptor Names
  Mol_descriptors = [] #Blank list to fill with descriptor values
  for mol in mols:
    mol = Chem.AddHs(mol) 
    descriptors = calc.CalcDescriptors(mol) #Calculation of descriptors
    Mol_descriptors.append(descriptors) #Adding descriptors to a list
  return Mol_descriptors, desc_names

Mol_descriptors, desc_names = RDKit_descriptors(individual_compounds_df['Compound'])
individual_compounds_df = individual_compounds_df.join(pd.DataFrame(Mol_descriptors, columns = desc_names)) #Adding to the table


def calc_QED(smiles):
    mols = [Chem.MolFromSmiles(i) for i in smiles] #Getting a list of molecules
    rdkit_qed = [list(Chem.QED.properties(mol)) for mol in mols]
    header = ['MW', 'ALOGP', 'HBA', 'HBD', 'PSA', 'ROTB', 'AROM', 'ALERTS']
    df = pd.DataFrame(rdkit_qed,columns=header)
    return df

df_QED = calc_QED(individual_compounds_df['Compound'])
individual_compounds_df = individual_compounds_df.join(df_QED[['MW','HBD', 'AROM', 'ALERTS']]) #Adding to the table


def get_nHM(smiles):
    atomic_list = [3,11,19,37,55,87,
                   4,12,20,38,56,88,
                   21,22,23,24,25,26,27,28,29,30,
                   39,40,41,42,43,44,45,46,47,48,
                   72,73,74,75,76,77,78,79,80,
                   13,31,49,50,81,82,83]

    n_HM = 0
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.RemoveHs(mol)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in atomic_list:
            n_HM += 1
    return n_HM


individual_compounds_df['n_HM'] = individual_compounds_df['Compound'].apply(get_nHM)

#individual_compounds_df.to_csv('ind_comp_desc.csv')
individual_compounds_df.index = individual_compounds_df['Compound']


# Adding structural descriptors to the main table
Descr_list =  ['MW', 'HBD', 'fr_Al_COO', 'fr_Ar_COO',
       'fr_Ar_N', 'fr_Ar_OH', 'fr_NH0', 'fr_NH1', 'fr_amide', 'AROM', 'ALERTS',
       'n_HM'] #List of descriptors
for desc in Descr_list:
  f_get_desc = lambda x: individual_compounds_df.loc[x][desc] if isinstance(x, str) else 0 #The function of adding a descriptor
  #Adding to the table with the design of the name
  for num_comp in range(2):
    name_new_column = desc + '#' + str(num_comp + 1)
    name_old_column = 'Smiles' + '#' + str(num_comp + 1)
    df[name_new_column] = df[name_old_column].apply(f_get_desc)


# Let's calculate the average by descriptors
Desc_list_new = ['MW', 'HBD', 'fr_Al_COO', 'fr_Ar_COO',
       'fr_Ar_N', 'fr_Ar_OH', 'fr_NH0', 'fr_NH1', 'fr_amide', 'AROM', 'ALERTS', 'n_HM']
for desc in Desc_list_new:
  df[desc] = df[desc + '#' + '1'] * df['X#1 (molar fraction)'] + df[desc + '#' + '2'] * df['X#2 (molar fraction)']


df = df[['Tmelt, K', 'Smiles#1', 'Smiles#2', 'Component#1', 'Component#2',	'X#1 (molar fraction)',	'X#2 (molar fraction)',	'T#1', 'T#2', 'MW', 'HBD', 'fr_Al_COO', 'fr_Ar_COO','fr_Ar_N', 'fr_Ar_OH', 'fr_NH0', 'fr_NH1', 'fr_amide', 'AROM', 'ALERTS','n_HM']]
df.to_csv('Melt_temp_ML.csv')
