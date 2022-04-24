## Description
Building a **QSAR** model to predict **pIC50** values of different biological compounds from their **pubChem fingerprints**\
Inspired by John S. Delaney & Pat Walters (article/research paper) dataset comprising of collection of compounds along with their molecular solubilities\
\
source - \
[original research paper](https://pubs.acs.org/doi/10.1021/ci034243x) \
[published article](https://www.tandfonline.com/doi/full/10.1517/17460441.2015.1016497)

## Bibliography
1. 
2. 
3. 

## Process
#### ChEMBL database
Retrieving and downloading biological activity data from ChEMBL database to construct machine learning models or QSARs. These data contains curated bioactivity data of more than 2 million compounds. It is compiled from more than 76,000 documents, 1.2 million essays and the data spans 13,000 targets, 1800 cells and 33,000 indications

#### Acetylcholinesterase
Cholinergic enzyme primarily found at postsynaptic neuromuscular junctions, especially in muscles and nerves \
Lack of acetylcholinesterase contributes to neuromuscular junction dysfunction in Type 1 diabetic neuropathy \
\
Search **"acetylcholinesterase"** in all targets :  referring to target proteins and target organisms \
Selected Target: CHEMBL220 \
standard_type: IC50 : inhibitory concentration @50% of target protein \
standard_value: potency of the drug \
Biologically, these compounds will come into contact with the protein/organism and induce a modulotory activity towards it {either to activate or to inhibit}

#### Activity Classification
* active compound : IC50 < 1000 nM
* inactive compound : IC50 > 10000 nM
* intermediate compound : 1000 < IC50 < 10000 

#### Lipinski Descriptors
**rdkit** allows us to compute molecular decriptors for the compounds in the dataset that we have compiled \
\
**Christopher Lipinski**, a scientist at Pfizer came up with a set of rule-of-five for evaluating the 'relative druglikeness of the compound' of compounds based on the key pharmaceutical kinetic properties:
* Absorption
* Distribution
* Metabolism
* Excretion
\
Lipinski analyzed all orally active FDA-approved drugs in the formulation whether it can be absorbed into the body, distributed to the proper tissue/organs and become metabolised. Following 4 descriptors that was used for his analysis has corresponding values in multiples of 5:
* Molecular Weight < 500 Dalton
* Octanol-Water partition coefficient (LogP) < 5
* Hydrogen bond donors < 5
* Hydrogen bond acceptors < 10
\
**Calculating Descriptors** \
Converting standard_value from IC50 to pIC50 to allow IC50 to be more uniformly distributed. We will convert the IC50 values to the negative logarithmic scale which is essentially -log10. Applying value normalization because values greater than negative logarithmic of 100,000,000 will become negative

#### Chemical Space Analysis
**Jose Medina Franco** (author) \
"Each chemical compound can be thought of as stars, i.e. active molecules would be compared to as constellations" \
He developed an approach termed as "constellation plot" whereby one can perform chemcial space analysis and create constellation plot where active molecule have correspondingly have larger size compared to less active molecule \

1. Frequency plot of 2 bioactivity classes comparing the inactive and active molecules
    ![freq_plot](https://github.com/subhashishansda4/Bio-Informatics/blob/main/assets/plots/plot_bioactivity_class.jpg)
2. Scatter plot of molecular weight (MW) v/s molecular solubilityt (LogP)
    ![scatter_plot](https://github.com/subhashishansda4/Bio-Informatics/blob/main/assets/plots/plot_MW_vs_logP.jpg) \
    It can be seen that the 2 bioactivity classes are spanning similar chemical spaces as evident by the scatter plot of MW v/s LogP \

#### Mann Whitney Analysis

    



