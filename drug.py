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

# importing libraries
import pandas as pd
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