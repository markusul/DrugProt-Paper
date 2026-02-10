# Drug-Prot: A query system for statistical inference of drug effects and interactions in dynamic proteomic networks

This repository contains all the code used for the paper "[*Drug-Prot: A query system for statistical inference of drug effects and interactions in dynamic proteomic networks*]()" by Markus Ulmer, Rui Sun, Liujia Qian, Ruedi Aebersold, Tiannan Guo, and Peter Bühlmann (2026).

All the statistical analysis was run on the [*Euler*](https://scicomp.ethz.ch/wiki/Euler) using the batch jobs in the 
`slurm` folder. See [*here*](https://scicomp.ethz.ch/wiki/Euler_applications_and_libraries_ubuntu) and sessionInfo.txt for details about the used libraries. The experiment will need the following folders and the proteomics measurement `ProteinMatrix_sampleID_MapEC50_20240229.csv`.

```
📂 DrugProt
├── 📂 R
├── 📂 Z
├── 📂 data
│   └── 📊 ProteinMatrix_sampleID_MapEC50_20240229.csv
├── 📂 figures
└── 📂 results
    ├── 📂 DrugEffects
    ├── 📂 ProteinEffects
    ├── 📂 anchorG
    └── 📂 anchor_opt
```

### 1. Analysis Pipeline
Run the following scripts in order to reproduce the model fits:

1.  `prepare.slurm`: **Preprocessing** of the data.
2.  `getZ.sh`: Calculates projections needed for the **de-sparsified Lasso regressions**.
3.  `fit.slurm`: Fits the models.
4.  `anchorG_CV.sh`: Performs out-of-distribution cross-validation for the **Anchor Forests**.
5.  `anchorG_opt.slurm`: Refits the Anchor Forests with the optimal gamma.
    > *Note: Determine the optimal gamma from `R/anchorG_vis.R` (Line 71) before running.*
6.  `anchorG_opt_res.slurm`: Calculates regularization paths and partial dependencies.

### 2. Visualization & Output
Once the pipeline is complete, run these R scripts to generate figures:

1.  `anchorG_vis.R`: **Visualizes** the results from Anchor Forests and saves findings to `results/A_Results.txt`.
2.  `pValOrganization.R`: **Organizes** the p-values from `DrugProt` for visualization.
3.  `pValVis.R`: **Visualizes** the p-values from `DrugProt` and saves findings to `results/P_Results.txt`.

## More details on the different scripts

<details>
<summary><strong>📂 Click to view detailed File Inputs & Outputs</strong></summary>

### 1. Data Processing
| Script | Needs (Input) | Generates (Output) |
| :--- | :--- | :--- |
| `R/data_preprocessing.R` | `data/ProteinMatrix_sampleID_MapEC50_20240229.csv` | `data/drugLookup.RData`<br>`data/na_count.RData`<br>`data/prepData.RData` |
| `R/data_preparation.R` | `data/prepData.RData` | `data/aggData.RData`<br>`data/protNames.RData` |
| `R/lagged_time.R` | `data/prepData.RData`<br>`data/protNames.RData` | `data/order.RData`<br>`data/laggedData.RData` |

### 2. P-Value Generation
| Script | Needs (Input) | Generates (Output) |
| :--- | :--- | :--- |
| `R/getZ.R` | `data/laggedData.RData` | `Z/6.RData`<br>`Z/24.RData`<br>`Z/48.RData` |
| `R/drugInteraction.R` | `data/laggedData.RData`<br>`Z/6.RData`, `Z/24.RData`, `Z/48.RData` | `results/DrugEffects/...`<br>`results/ProteinEffects/...` |
| `R/pValOrganization.R` | `data/order.RData`<br>`results/DrugEffects/...`<br>`results/ProteinEffects/...` | `results/DrugEffects.RData`<br>`results/proteinNetworkPval.RData` |
| `R/pValVis.R` | `data/drugLookup.RData`<br>`data/order.RData`<br>`Z/6.RData`, `Z/24.RData`, `Z/48.RData`<br>`results/proteinNetworkPval.RData`<br>`results/anchor_opt/proteinSelection.RData` | **All P-value Plots**<br>`results/P_Results.txt` |

### 3. Anchor Forest Analysis
| Script | Needs (Input) | Generates (Output) |
| :--- | :--- | :--- |
| `R/anchorG_CV.R` | `R/utils.R`<br>`data/aggData.RData`<br>`data/protNames.RData` | `results/anchorG/...` |
| `R/anchorG_opt.R` | `R/utils.R`<br>`data/aggData.RData`<br>`data/protNames.RData` | `results/anchorG_opt.RData` |
| `R/anchorG_opt_res.R` | `results/anchorG_opt.RData` | `results/anchor_opt/var_importance.RData`<br>`results/anchor_opt/regPath.RData`<br>`results/anchor_opt/stability_selection.RData`<br>`results/anchor_opt/partial_dependence.RData` |
| `R/anchorG_vis.R` | `R/utils.R`<br>`data/aggData.RData`<br>`data/protNames.RData`<br>`data/order.RData`<br>`results/anchorG/...`<br>`results/anchor_opt/var_importance.RData`<br>`results/anchor_opt/regPath.RData`<br>`results/anchor_opt/stability_selection.RData`<br>`results/anchor_opt/partial_dependence.RData` | **Anchor Forest Plots**<br>`results/A_Results.txt`<br>`results/anchor_opt/proteinSelection.RData`<br>`results/most_important_proteins.txt` |

</details>

```mermaid
graph TD
    %% Define Styles
    classDef script fill:#f9f,stroke:#333,stroke-width:2px;
    classDef file fill:#e1f5fe,stroke:#0277bd,stroke-width:2px;
    classDef output fill:#e8f5e9,stroke:#2e7d32,stroke-width:2px;

    subgraph Data_Prep ["**Data Preparation**"]
        S1(R/data_preprocessing.R):::script
        S2(R/data_preparation.R):::script
        S3(R/lagged_time.R):::script
        
        Raw[ProteinMatrix...20240229.csv]:::file --> S1
        S1 --> D1[data/drugLookup.RData]:::file
        S1 --> D2[data/na_count.RData]:::file
        S1 --> D3[data/prepData.RData]:::file
        
        D3 --> S2
        S2 --> D4[data/aggData.RData]:::file
        S2 --> D5[data/protNames.RData]:::file
        
        D3 & D5 --> S3
        S3 --> D6[data/order.RData]:::file
        S3 --> D7[data/laggedData.RData]:::file
    end

    subgraph P_Value ["**P-Value Generation**"]
        S4(R/getZ.R):::script
        S5(R/drugInteraction.R):::script
        S6(R/pValOrganization.R):::script
        S7(R/pValVis.R):::script
        
        D7 --> S4
        S4 --> Z1[Z/6.RData, Z/24.RData...]:::file
        
        D7 & Z1 --> S5
        S5 --> R1[results/DrugEffects/...]:::file
        S5 --> R2[results/ProteinEffects/...]:::file
        
        D6 & R1 & R2 --> S6
        S6 --> R3[results/DrugEffects.RData]:::file
        S6 --> R4[results/proteinNetworkPval.RData]:::file
        
        D1 & D6 & R3 & R4 & Z1 --> S7
        S7 --> Out1[All P-value Plots]:::output
        S7 --> Out2[results/P_Results.txt]:::output
    end

    subgraph Anchor_Forest ["**Anchor Forest**"]
        S8(R/anchorG_CV.R):::script
        S9(R/anchorG_opt.R):::script
        S10(R/anchorG_opt_res.R):::script
        S11(R/anchorG_vis.R):::script
        
        Utils[R/utils.R]:::file
        
        Utils & D4 & D5 --> S8
        S8 --> R5[results/anchorG/...]:::file
        
        Utils & D4 & D5 --> S9
        S9 --> R6[results/anchorG_opt.RData]:::file
        
        R6 --> S10
        S10 --> R7[results/anchor_opt/*]:::file
        
        Utils & D4 & D5 & D6 & R5 & R7 --> S11
        S11 --> Out3[Anchor Forest Plots]:::output
        S11 --> Out4[results/A_Results.txt]:::output
    end
    
    %% Connect visualization back to data if needed
    S11 -.->|Generates| Sel[results/anchor_opt/proteinSelection.RData]:::file
    Sel -.-> S7

    %% NEW CONNECTION: Optimal Gamma Feedback
    S11 -.->|Optimal Gamma| S9
```

