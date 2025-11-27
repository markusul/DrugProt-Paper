# drugResponse

#need
data
data/ProteinMatrix_sampleID_MapEC50_20240229.csv
Z
results
results/DrugEffects
results/ProteinEffects
results/anchorG
results/anchor_opt

R
figures

# data
R/data_preprocessing.R
needs
data/ProteinMatrix_sampleID_MapEC50_20240229.csv
generates 
data/drugLookup.RData
data/na_count.RData
data/prepData.RData

R/data_preparation.R
needs
data/prepData.RData
generates
data/aggData.RData
data/protNames.RData

R/lagged_time.R
needs
data/prepData.RData
data/protNames.RData
generates
data/order.RData
data/laggedData.RData

# Pvalue generation
R/getZ.R
needs
data/laggedData.RData
generates 
Z/6.RData
Z/24.RData
Z/48.RData


R/drugInteraction.R
needs
data/laggedData.RData
Z/6.RData
Z/24.RData
Z/48.RData

generates 
results/DrugEffects/...
results/ProteinEffects/...


R/pValOrganization.R
needs
data/order.RData
results/DrugEffects/...
results/ProteinEffects/...
generates
results/DrugEffects.RData
results/proteinNetworkPval.RData

R/pValVis.R
needs
data/drugLookup.RData
data/order.RData
results/DrugEffects.RData
results/proteinNetworkPval.RData
Z/6.RData
Z/24.RData
Z/48.RData
results/anchor_opt/proteinSelection.RData
generates
all P-value plots
results/P_Results.txt

# Anchor Forest
R/anchorG_CV.R
needs
R/utils.R
data/aggData.RData
data/protNames.RData
generates
results/anchorG/...

R/anchorG_opt.R
needs
R/utils.R
data/aggData.RData
data/protNames.RData
generates
results/anchorG_opt.RData

R/anchorG_opt_res.R
needs
results/anchorG_opt.RData
generates
results/anchor_opt/var_importance.RData
results/anchor_opt/regPath.RData
results/anchor_opt/stability_selection.RData
results/anchor_opt/partial_dependence.RData

R/anchorG_vis.R
needs
R/utils.R
data/aggData.RData
data/protNames.RData
data/order.RData
results/anchorG/...
results/anchor_opt/var_importance.RData
results/anchor_opt/regPath.RData
results/anchor_opt/stability_selection.RData
results/anchor_opt/partial_dependence.RData
generates
Anchor Forest plots
results/A_Results.txt
results/anchor_opt/proteinSelection.RData
results/most_important_proteins.txt



