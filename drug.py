# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 12:50:12 2022

@author: VAGUE
"""

'''
"John S. Delaney" - ESOL: Estimating Aqueous Solubility Directly from Molecular Structure (research paper)
linear regression model for predicting molecular solubility

"Pat Walters" - Deep Learning for the Life Sciences (author)
'''

# DATA COLLECTION
#---------------------------------------------------------------------------------------------------------
# importing libraries
import pandas as pd
import numpy as np
import seaborn as sns
sns.set(style='ticks')
import matplotlib.pyplot as plt
# machine learning
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.feature_selection import VarianceThreshold
# lazypredict machine learning
import lazypredict
from lazypredict.Supervised import LazyRegressor
# rdkit for lipinski descriptors
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
# CHEMBL database
from chembl_webresource_client.new_client import new_client

# search for target protein
target = new_client.target
target_query = target.search('coronavirus')
targets = pd.DataFrame().from_dict(target_query)

# select and retreive bioactivity data for SARS coronavirus 3C-like proteinase
selected_target = targets.target_chembl_id[4]

# filtering by selected_target and standard_type
activity = new_client.activity
res = activity.filter(target_chembl_id = selected_target).filter(standard_type = 'IC50')

# creating dataframe
df = pd.DataFrame().from_dict(res)
df.head(3)
df.standard_type.unique()

# dataframe to csv
df.to_csv('bioactivity_data.csv', index = False)

# dropping missing value for the standard_value column
df = df[df.standard_value.notna()]

# data pre-processing
# labelling compound as either being active, inactive or intermediate
bioactivity_class = []
for i in df.standard_value:
    if float(i) >= 10000:
        bioactivity_class.append('inactive')
    elif float(i) <= 1000:
        bioactivity_class.append('active')
    else:
        bioactivity_class.append('intermediate')

# iterating through molecule_chembl_id
mol_cid = []
for i in df.molecule_chembl_id:
    mol_cid.append(i)
    
# iterating through canonical_smiles
canonical_smiles = []
for i in df.canonical_smiles:
    canonical_smiles.append(i)
    
# iterating through standard_value
standard_value = []
for i in df.standard_value:
    standard_value.append(i)

# combining the 4 lists into a dataframe
data_tuples = list(zip(mol_cid, canonical_smiles, bioactivity_class, standard_value))
df_bioclass = pd.DataFrame(data_tuples, columns=['molecule_chembl_id', 'canonical_smiles', 'bioactivity_class', 'standard_value'])
        
# =============================================================================
# # selecting different parameters to iterate through and make them unique
# selection = ['molecule_chembl_id', 'canonical_smiles', 'standard_value']
# df_filter = df[selection]
# # merging dataframes
# df_bioclass = pd.concat([df_filter, pd.Series(bioactivity_class)], axis=1)
# =============================================================================

# dataframe to csv
df_bioclass.to_csv('bioactivity_preprocessed_data.csv', index = False)
#---------------------------------------------------------------------------------------------------------


# DATA CLEANING & PROCESSING
#---------------------------------------------------------------------------------------------------------
# source - https://codeocean.com/explore/capsules?query=tag:data-curation
# calculating lipinski descriptors
def lipinski(smiles, verbose = False):
    moldata = []
    for elem in smiles:
        mol = Chem.MolFromSmiles(elem)
        moldata.append(mol)
        
    baseData = np.arange(1,1)
    i=0
    for mol in moldata:
        desc_MolWt = Descriptors.MolWt(mol)
        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_NumHDonors = Lipinski.NumHDonors(mol)
        desc_NumHAcceptors = Lipinski.NumHAcceptors(mol)
        
        row = np.array([desc_MolWt,
                        desc_MolLogP,
                        desc_NumHDonors,
                        desc_NumHAcceptors
                    ])
        
        if(i==0):
            baseData = row
        else:
            baseData = np.vstack([baseData, row])
        i=i+1
    
    columnNames = ['MW', 'LogP', 'NumHDonors', 'NumHAcceptors']
    descriptors = pd.DataFrame(data=baseData, columns=columnNames)
    
    return descriptors

# dataframe
df_lipinski = lipinski(df.canonical_smiles)
df_lipinski = pd.concat([df_bioclass, df_lipinski], axis=1)

# dataframe to csv
df_lipinski.to_csv('lipinski_descriptors.csv', index=False)

# NOTE
# values greater than 100,000,000 will be fixed at that value
# otherwise the negative logarithmic value will become negative
df_lipinski.standard_value.describe()
-np.log10((10**-9)*100000000)
-np.log10((10**-9)*10000000000)

# changing data type of 'standard_value'
print(df_lipinski.dtypes)
df_lipinski['standard_value'] = df_lipinski['standard_value'].astype(float)
print(df_lipinski.dtypes)

# capping the standard_value
def norm_value(input):
    norm = []
    
    for i in input['standard_value']:
        if i > 100000000:
            i = 100000000
        norm.append(i)
    
    input['standard_value_norm'] = norm
    #x = input.drop('standard_value', 1)
    x = input
    
    return x


# converting IC50 to pIC50
# source - https://github.com/chaninlab/estrogen-receptor-alpha-qsar/blob/master/02_ER_alpha_RO5.ipynb
def pIC50(input):
    pIC50 = []
    
    for i in input['standard_value_norm']:
        molar = float(i)*(10**-9)  # converts nM to M
        pIC50.append(-np.log10(molar))
        
    input['pIC50'] = pIC50
    x = input.drop('standard_value_norm', 1)
    
    return x

# applying normalization of values
df_final = norm_value(df_lipinski)
# applying conversion of IC50 to pIC50
df_final = pIC50(df_final)
# removing 'intermediate' bioactivity class
df_final = df_final[df_final.bioactivity_class != 'intermediate']

# dataframe to csv
df_final.to_csv('final_data.csv', index=False)
#---------------------------------------------------------------------------------------------------------


# EDA (CHEMICAL SPACE ANALYSIS)
#---------------------------------------------------------------------------------------------------------
'''
Jose Medina Franco (author)
each chemical compound could be thought of as stars, i.e. active molecules would be compared to as constellations
he developed an approach termed as "Constellation Plot" whereby one can perform chemical space analysis & create constellation plot where

[active molecule would be correspondingly have larger sizes compared to less active molecule]
'''
# frequency plot of bioactivity classes comparing inactive and active molecules
plt.figure(figsize = (5.5, 5.5))
sns.countplot(x = 'bioactivity_class', data = df_final, edgecolor = 'black')

plt.xlabel('Bioactivity class', fontsize = 14, fontweight = 'bold')
plt.ylabel('Frequency', fontsize = 14, fontweight = 'bold')

plt.savefig('plot_bioactivity_class.jpg')

# scatter plot of molecular weight(MW) v/s molecular solubility(logP)
plt.figure(figsize = (5.5, 5.5))
sns.scatterplot(x = 'MW', y = 'LogP', data = df_final, hue = 'bioactivity_class', size = 'pIC50', edgecolor = 'black', alpha = 0.7)

plt.xlabel('MW', fontsize = 14, fontweight = 'bold')
plt.ylabel('LogP', fontsize = 14, fontweight = 'bold')
plt.legend(bbox_to_anchor = (1.05, 1), loc = 2, borderaxespad = 0)

plt.savefig('plot_MW_vs_logP.jpg')

# mann-whitney U test
# source - https://machinelearningmastery.com/nonparametric-statistical-significance-tests-in-python
def mannwhitney(descriptor, verbose = False):
    from numpy.random import seed
    from scipy.stats import mannwhitneyu
    
    # seeding the random number generator
    seed(1)
    
    # actives and inactives
    selection = [descriptor, 'bioactivity_class']
    df = df_final[selection]
    active = df[df.bioactivity_class == 'active']
    active = active[descriptor]
    
    selection = [descriptor, 'bioactivity_class']
    df = df_final[selection]
    inactive = df[df.bioactivity_class == 'inactive']
    inactive = inactive[descriptor]
    
    # compare samples
    stat, p = mannwhitneyu(active, inactive)
    print('Statistics=%.3f, p=%.3f' % (stat, p))
    
    # interpret
    alpha = 0.05
    if p > alpha:
        interpretation = 'Same distribution (fail to reject H0)'
    else:
        interpretation = 'Different distribution (reject H0)'
        
    results = pd.DataFrame({'Descriptor': descriptor,
                            'Statistics': stat,
                            'p': p,
                            'alpha': alpha,
                            'Interpretation': interpretation}, index=[0])
    
    filename = 'mannwhitneyu_' + descriptor + '.csv'
    results.to_csv(filename)
    
    return results

# box plots for pIC50, MW, LogP, NumHDonors, NumHAcceptors
# performing mann-whitney analysis for each one of them

# pIC50
plt.figure(figsize = (5.5, 5.5))
sns.boxplot(x = 'bioactivity_class', y = 'pIC50', data = df_final)
plt.xlabel('Bioactivity_class', fontsize = 14, fontweight = 'bold')
plt.ylabel('pIC50 value', fontsize = 14, fontweight = 'bold')
plt.savefig('plot_ic50.jpg')
mannwhitney('pIC50')

# MW
plt.figure(figsize = (5.5, 5.5))
sns.boxplot(x = 'bioactivity_class', y = 'MW', data = df_final)
plt.xlabel('Bioactivity_class', fontsize = 14, fontweight = 'bold')
plt.ylabel('MW', fontsize = 14, fontweight = 'bold')
plt.savefig('plot_MW.jpg')
mannwhitney('MW')

# LogP
plt.figure(figsize = (5.5, 5.5))
sns.boxplot(x = 'bioactivity_class', y = 'LogP', data = df_final)
plt.xlabel('Bioactivity_class', fontsize = 14, fontweight = 'bold')
plt.ylabel('LogP', fontsize = 14, fontweight = 'bold')
plt.savefig('plot_LogP.jpg')
mannwhitney('LogP')

# NumHDonors
plt.figure(figsize = (5.5, 5.5))
sns.boxplot(x = 'bioactivity_class', y = 'NumHDonors', data = df_final)
plt.xlabel('Bioactivity_class', fontsize = 14, fontweight = 'bold')
plt.ylabel('NumHDonors', fontsize = 14, fontweight = 'bold')
plt.savefig('plot_NumHDonors.jpg')
mannwhitney('NumHDonors')

# NumHAcceptors
plt.figure(figsize = (5.5, 5.5))
sns.boxplot(x = 'bioactivity_class', y = 'NumHAcceptors', data = df_final)
plt.xlabel('Bioactivity_class', fontsize = 14, fontweight = 'bold')
plt.ylabel('NumHAcceptors', fontsize = 14, fontweight = 'bold')
plt.savefig('NumHAcceptors.jpg')
mannwhitney('NumHAcceptors')
#---------------------------------------------------------------------------------------------------------


# DATA PREPARATION
#---------------------------------------------------------------------------------------------------------
selection = ['canonical_smiles', 'molecule_chembl_id']
df_final_selection = df_final[selection]
df_final_selection.to_csv('molecule.smi', sep='\t', index=False, header=False)

# preparing X and Y matrices
# X
df_final_X = pd.read_csv('descriptors_output.csv')
df_final_X = df_final_X.drop(columns=['Name'])
# Y
df_final_Y = df_final['pIC50']
df_final_Y = df_final_Y.to_frame()

# combining X and Y variable
dataset = pd.concat([df_final_X, df_final_Y], axis=1)
# dataframe to csv
dataset.to_csv('dataset.csv', index=False)
#---------------------------------------------------------------------------------------------------------

# MACHINE LEARNING MODELS
#---------------------------------------------------------------------------------------------------------
# input features
X = dataset.drop('pIC50', axis=1)
# output features
Y = dataset.pIC50
# data dimension
print(X.shape)
print(Y.shape)

# =============================================================================
# # remove low variance features
# selection = VarianceThreshold(threshold=(.8*(1-.8)))
# X = selection.fit_transform(X)
# print(X.shape)
# =============================================================================

# data split
x_train, x_test, y_train, y_test = train_test_split(X, Y, test_size=0.2)
print(x_test.shape)
print(y_test.shape)

Y = Y.to_frame()
y_test = y_test.to_frame()
y_train = y_train.to_frame()

# building over 42 regression models
# using LazyPredict
# default parameters
clf = LazyRegressor(verbose=0, ignore_warnings=(True), custom_metric=(None))
train, test = clf.fit(x_train, x_test, y_train, y_test)
# performance table of training set
print(train)
# performance table of test set
print(test)
#---------------------------------------------------------------------------------------------------------

# DATA VISUALIZATION
#---------------------------------------------------------------------------------------------------------
# bar plots
# R-Squared values
# =============================================================================
# plt.figure(figsize = (5, 10))
# sns.set_theme(style='whitegrid')
# ax = sns.barplot(y=train.index, x='R-Squared', data=train)
# ax.set(xlim=(0,1))
# plt.savefig('R-Squared.jpg')
# 
# # RMSE values
# plt.figure(figsize = (5, 10))
# sns.set_theme(style='whitegrid')
# ax = sns.barplot(y=train.index, x='RMSE', data=train)
# ax.set(xlim=(0,10))
# plt.savefig('RMSE.jpg')
# 
# # calculation time
# plt.figure(figsize = (5, 10))
# sns.set_theme(style='whitegrid')
# ax = sns.barplot(y=train.index, x='Time Taken', data=train)
# ax.set(xlim=(0,10))
# plt.savefig('Calculation_Time.jpg')
# =============================================================================



