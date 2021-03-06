###### README for 'The multilayer community structure of Medulloblastoma' (Nuñez-Carpintero et al, 2021. DOI:https://doi.org/10.1016/j.isci.2021.102365) scripts. Please, feel free to contact the authors for any information and doubts at: iker.nunez@bsc.es, davide.cirillo@bsc.es .
###### To cite this article or its scripts, and CmmD package, use: Núñez-Carpintero, I. et al. (2021) ‘The multilayer community structure of medulloblastoma’, iScience, 24(4). doi: 10.1016/j.isci.2021.102365.

You can find the networks from the multilayer at: ***https://github.com/cirillodavide/gene_multilayer_network/tree/master/networks***

This Readme file contains the instructions to launch the scripts generated for the different multilayer network analysis

***DISCLAIMER: Accessing RNA-sequencing data from Forget et al. 2018 requires signing a Data Transfer Agreement (DPA) by the European Genome-phenome Archive (EGA). EGA accession number: EGAB00000001285. Therefore, this data is not available in the repository.***

#### Before launching any of the scripts, please take a moment to check if you have installed the library dependencies of these scripts:

      R packages:
            CmmD (available to install via devtools at GitHub repository ikernunezca/CmmD)
            BiocManager
            RCy3
            paxtoolsr
            igraph
            plyr
            biomaRt
            parallel
            rJava
            dendroextras
            dendextend
            circlize
            sigclust2
            pvclust
            fpc
            jaccard
            readr
            knitr
            data.table
            AnnotationDbi
            igraph
            stringr
            e1701
            mgcv
            MASS
            ggplot2
            ggrepel
            neat
      Python libraries:
            pandas
            mygene
            shutil
            urllib.parse
            urllib.request
            contextlib
            csv
            tqdm
            re
            os
            gzip
            io
            subprocess
            requests
            os.path
            numpy
            collections
            matplotlib
      
R version used: 3.6.1
Python version: 2.7
Molti version: https://github.com/gilles-didier/MolTi

### First steps in order to use this repository correctly:
    1. Download the full repository and access it via cd. 
    2. All the scripts should be run from the directory where we saved the repository:
          ~/Medulloblastoma $
          
### Preprocessing scripts @mpetrizzelli: 
***If you do not intend to produce this data again, you can jump to section 1***

These scripts generate the input data for the rest of the study, by preprocessing the data of both Medulloblastoma cohorts (Forget et al. 2018) (Archer et al. 2018).

You can use this scripts to reproduce the data, which then correponds to two files that are used by the rest of the scripts: ***https://raw.githubusercontent.com/iPC-project-H2020/ipcrg/master/scripts/CURIE2gr/multi.layer.net.gr*** (Forget et al. 2018 cohort, which is downloaded by the scripts themselves) and ***~/Medulloblastoma/data/multilayer_archer.gr*** (Archer et al. 2018). 

Files generated by the preprocessing scripts are saved at ***~/Medulloblastoma/studies_preprocessing/preprocessingCOHORT_2018/Networks/multi.layer.no.ggi.net***. Be sure to check the path dependencies on the script you are going to use if you want to run them based on this files, as they are by default configured to work with the files mentioned at the previous paragraph. 

### Reproducing preprocessing scripts :
    Forget et al. 2018:
            1. Run 'studies_preprocessing/preprocessing_Forget2018/Codes/multilayer.R:
                  ~/Medulloblastoma $: Rscript studies_preprocessing/preprocessing_Forget2018/Codes/multilayer.R
            2. This will produce the file '~/Medulloblastoma/studies_preprocessing/preprocessing_Forget2018/Networks/multi.layer.no.ggi.net'
            3. This file is a 3 column data frame, where the first is the patient ID, the second correspond to an associated gene (HGNC gene) and the third indicating the multi-omic data from where the relationship has been computed. 
            4. *IMPORTANT*: Convert the second column to Entrez gene ID before using this as input for the rest of the scripts. This can be done with a python script available at ~/Medulloblastoma/studies_preprocessing/preprocessing_Forget2018/Codes/CURIE2gr.py . An example of the final format you need for this file can be checked from https://raw.githubusercontent.com/iPC-project-H2020/ipcrg/master/scripts/CURIE2gr/multi.layer.net.gr. 
            5. With your table format ready in your file, you should now check the scripts you are going to run before doing so, in order to change one of the code lines. 
            THIS ONLY APPLIES TO SCRIPTS WORKING WITH Forget et al. 2018 COHORT. Most of this scripts include the following line that loads the cohort's data:
                     
                     tata <- as.data.frame(read.csv(file = "https://raw.githubusercontent.com/iPC-project-H2020/ipcrg/master/scripts/CURIE2gr/multi.layer.net.gr",sep = "\t",header = F, stringsAsFactors=F))
                     
            For example, in the script 'Get_best_accuracy_Ward.R' this corresponds to line 42. It may not be in the same position in other scripts (e.g., in the script Randomized_curie.R it is line 276) so be sure to check if this line exists before launching the script.
            If the line exists, change the URL directly to the directory where you saved the multiomics data file generated.
    
    Archer et al. 2018:
            1. Open a new R session: 
                        ~/Medulloblastoma $: R
            2. Run the content from the script 'studies_preprocessing/preprocessing_Archer2018/Codes/functions.R' within the session.
            3. Run the content from the script 'studies_preprocessing/preprocessing_Archer2018/Codes/packages.R' within the session if you may have some of the dependencies missing.
            4. Run the content from the script 'studies_preprocessing/preprocessing_Archer2018/Codes/MultilayerNetworkConstruction.R' within the session.  This script will generate the file '~/Medulloblastoma/studies_preprocessing/preprocessing_Archer2018/Networks/multi.layer.no.ggi.net'.
            5. This file is a 3 column data frame, where the first is the patient ID, the second correspond to an associated gene (HGNC gene) and the third indicating the multi-omic data from where the relationship has been computed. 
            6. *IMPORTANT*: Convert the second column to Entrez gene ID before using this as input for the rest of the scripts. This can be done with a python script available at ~/Medulloblastoma/studies_preprocessing/preprocessing_Archer2018/Codes/CURIE2gr_mod.py . An example of the final format you need for this file can be checked from ~/Medulloblastoma/data/multilayer_archer.gr.
            7. With your table format ready, you can now replace the file at ~/Medulloblastoma/data/multilayer_archer.gr .


### 1) Obtain full theta-lambda accuracy analysis (Figure 4, Supplementary Tables 1, 2 & 3) @ikernunezca:
    Run Get_best_accuracy_Ward.R script: 
          1. ~/Medulloblastoma $ Rscript Get_best_accuracy_Ward.R
          2. When finished, you will find 2 files at ~/Medulloblastoma/data/Output_Get_best_accuracy_Ward corresponding to supplementary table 1 (clusters_wardd2.csv), Supplementary table 2 (accuracy_wardd2.csv) and Matthew's Correlation Coefficient (matthews_wardd2.csv)
          3. Figure 4 is saved at ~/Medulloblastoma/data/Plots/Figure_4.pdf
          
### 2) Plot Figure 3 (Community trajectory dendrogram) @ikernunezca: 
    Run Figure_3.R script: 
          1. ~/Medulloblastoma $ Rscript Figure_3.R
          2. When finished, you will find the plot at ~/Medulloblastoma/Plots/Figure_3.pdf

### 3) Plot Figure 5, Supplementary 8 and 9 (Hierarchical clustering and bootstrap significance analysis) @ikernunezca: 
    Run Figure_5_Supplementary_8_&_9.R: 
          1. ~/Medulloblastoma $ Rscript Figure_5_Supplementary_8_&_9.R
          2. When finished, you will find a pdf file with the 4 plots (Figure 5, Supplementary 8, Supplementary 9A and 9B) at ~/Medulloblastoma/data/Plots/clustering_pv_shc_hclust.pdf

### 4) Plot Supplementary Figure 3 (Inflection point definition) @cirillodavide:
    Run Supplementary_Figure_3.ipynb script: 
          1. This file is a jupyter notebook written in python 2.7 version. Run a jupyter notebook server from ~/Medulloblastoma and access the script via web browser.

### 5) Plot Supplementary Figure 4 (Dynamic events analysis) @cirillodavide:
    Run dynamic_events.ipynb script: 
          1. This file is a jupyter notebook written in python 2.7 version. Run a jupyter notebook server from ~/Medulloblastoma and access the script via web browser:: ~/Medulloblastoma/dynamic_events_analysis/dynamic_events.ipynb
          
### 6) Plot Supplementary Figure 5 (Gene shuffling tests) @ikernunezca:
     Input randomization files used for the analysis represented in this plot can be accesed at ~/Medulloblastoma/data/Randomizations. In order to generate your own randomizations, run the Randomized_curie.R script. You can find the new randomizations at the very same directory.
                   ~/Medulloblastoma $ Rscript Randomized_curie.R 1 50 # In order to generate the first 50 randomizations
                   ~/Medulloblastoma $ Rscript Randomized_curie.R 51 100 # Randomizations 51 to 100
        Only 50 randomizations are allowed at each run of Randomized_curie as it is intented to be used within an HPC enviroment in a parallelized way.
        Once your randomizations are ready, Run Supplementary_Figure_5.R script. 
          1. ~/Medulloblastoma $ Rscript Supplementary_Figure_5.R
          2. When finished, you will find the pdf at: ~/Medulloblastoma/data/Plots/Supplementary_Figure_5.pdf

### 7) Plot Supplementary Figure 6 (Distributions of optimization accuracies) @ikernunezca:
      Run Supplementary_Figure_6.R script: 
          1. ~/Medulloblastoma $ Rscript Supplementary_Figure_6.R
          2. When finished, you will find the plot at ~/Medulloblastoma/Plots/Supplementary_Figure_6.png

### 8) Plot Supplementary Figure 7 (Recursive exclusion test) @ikernunezca:
     *NOTE* expected run time is approximately 5 hours.
        Part 1: Script for performing the analysis. Final output is a sequence of barplots:
            Run Supplementary_Figure_7.R script: 
                  1. ~/Medulloblastoma $ Rscript Supplementary_Figure_7.R
                  2. You will find your output plot at: ~/Medulloblastoma/data/Plots/Supp_Figure_7.png
        Part 2: For obtaining Supplementary figure as it is in the paper's final version
                  Run Supplementary_Figure_7_new.R:
                  1. ~/Medulloblastoma $ Rscript Supplementary_Figure_7_new.R
                  2. You will find your output plot at: ~/Medulloblastoma/data/Plots/Supp_Figure_7_paper_version.pdf
          
### 9) Supplementary Table 6 (Provenance Analysis) @cirillodavide
          1. This file is a jupyter notebook written in python 2.7 version. Run a jupyter notebook server from ~/Medulloblastoma and access the script via web browser: ~/Medulloblastoma/provenance_analysis/provenance_analysis.ipynb
          2. The output table file (Supplementary Table 6) is saved at ~/Medulloblastoma/provenance_analysis/provenance_analysis.tsv
          3. Supplementary Table 6 can also be found at ~/Medulloblastoma/Supplementary_Tables directory
          
### 10) Full Protein-Protein interaction Enrichment analysis @ikernunezca
      *NOTE* We generated a single script for the ppi enrichment analysis because of resource capabilities. This script is intended to be used within an HPC enviroment in a parallelized way. Run Ppi_enrichment.R script in the following way:
            1. ~/Medulloblastoma $ Run Ppi_enrichment.R 1 # This will perform the enrichment analysis for patient 1. Changing argument 1 to 2 will perform the analysis for patient 2 an so on. 
            2. This should be done for each one of the 38 patients.
            3. The output files are saved at data/Protein_Enrichments. All output files appended conform the full PPI enrichment analysis. 

### 11) Supplementary Table 4 and Full Network Enrichment analysis (Drugs, Diseases -Genetic variants-, Pathways and Metabolism) @ikernunezca
      *NOTE*: Excepted runtime is approximately 15-20 min.
      Run Full_Network_Enrichment_analysis.R script:
            1. ~/Medulloblastoma $ Rscript Full_Network_Enrichment_analysis.R
            2. Outputs are saved at: ~/Medulloblastoma/Supplementary_Tables/
            
### 12) Archer et al. 2018 cohort analysis: Supplementary Figures 10 and 11 + Corresponding accuracy/ Cluster length/ Matthew's Correlation Coefficients @ikernunezca
      #### For obtaining Supplementary Figure 10:
            Run Supplementary_Figure_10.R script:
            1. ~/Medulloblastoma $ Rscript Supplementary_Figure_7.R
            2. This script also generates 3 files, that respectively correspond to the accuracy table (saved as Supplementary Table 7), the number of suggested clusters (saved as Supplementary Table 8) and the values for MCC (saved as Supplementary Table 9). This files are not named in the manuscript to avoid redundancy.
      #### For obtaining Supplementary Figure 11, and Tables 10 to 12:
            Run Supplementary_Figure_11.R script:
            1. ~/Medulloblastoma $ Rscript Supplementary_Figure_11.R
            2. This script also generates 3 files, that respectively correspond to the accuracy table (saved as Supplementary Table 10), the number of suggested clusters (saved as Supplementary Table 11) and the values for MCC (saved as Supplementary Table 12). This files are not named in the manuscript to avoid redundancy.
      Supplementary table outputs are saved at: ~/Medulloblastoma/Supplementary_Tables/
      Supplementary figure outputs are saved at: data/Plots/
