
# MetDNA2InSilicoTool

<!-- badges: start -->
<!-- badges: end -->

[**Knowledge-guided multi-layer network (KGMN)**](http://metdna.zhulab.cn/) is a new approach leveraging knowledge-guided multi-layer networks to annotate known and unknown metabolites in untargeted metabolomics data. Although KGMN is an independent software tool, it can further integrate with other workflows to help users discover and validate metabolites. This tutorial aims to provide an easy instruction to integrated KGMN results with 3 common in-silico MS/MS tools (MetFrag, CFM-ID, MS-FINDER).<br>

The goal of MetDNA2InSilicoTool focus on providing ways to help users linking KGMN with other tools. It should be note that the parameters need to be adjusted according to their instrument settings and experimental designs. The detailed usage please refer their own tutorials. <br>

**Detailed tutorial** can be found [**here**](https://github.com/ZhuMetLab/MetDNA2_Web/blob/main/Tutorials/Tutorial_KGMN_and_insilico_ms2.pdf)


## Installation

You can install the development version of MetDNA2InSilicoTool like so:

``` r
# Install required packages
if(!require(devtools)){
install.packages("devtools")
}

if(!require(BiocManager)){
install.packages("BiocManager")
}

# Install CRAN/Bioconductor packages
required_pkgs <- c("dplyr","tidyr","readr","stringr","rcdk")
list_installed <- installed.packages()

new_pkgs <- required_pkgs[!(required_pkgs %in% list_installed[,'Package'])]
if (length(new_pkgs) > 0) {
  BiocManager::install(new_pkgs)
} else {
  cat('Required CRAN/Bioconductor packages installed\n')
}


# Install GitHub packages - call MetFrag
devtools::install_github("schymane/ReSOLUTION")

# Install GitHub packages
devtools::install_github("ZhuMetLab/MetDNA2InSilicoTool")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
# set working directory
setwd('G:/00_projects/03_MetDNA2/00_data/20220609_insilico_ms2_demo/NIST_urine_pos/')

# load packages
library(dplyr)
library(MetDNA2InSilicoTool)

# copy files
copyFiles4InsilicoTool(dir_path = '.')

# set working directory
setwd('G:/00_projects/03_MetDNA2/00_data/20220609_insilico_ms2_demo/NIST_urine_pos/07_insilico_msms/')

# reformat identification_table
reformatTable1(dir_path = '.')


# generate files for in-silico MS/MS match
# peak 'M196T420' as example
generateFiles4InsilicoMsMs(peak_id = 'M196T420')

# run MetFrag
runMetFragMatch(peak_id = 'M196T420', 
                metfrag_path = 'F:/software/metfrag/MetFrag2.4.5-CL.jar', 
                ppm = 25, 
                mzabs = 0.01, 
                frag_ppm = 25)


# run CFM-ID
runCfmIdMatch(peak_id = 'M196T420',
              cfmid_path = 'F:/software/cfm_id/cfm-id.exe',
              config_file = 'F:/software/cfm_id/metab_se_cfm/param_config.txt',
              param_file = 'F:/software/cfm_id/metab_se_cfm/param_output0.log',
              score_type = 'Jaccard',
              ppm = 25,
              mzabs = 0.01)

# run MS-FINDER
# note: the dir_path must be given 
runMsFinderMatch(peak_id = 'M196T420',
                 dir_path = 'G:/00_projects/03_MetDNA2/00_data/20220609_insilico_ms2_demo/NIST_urine_pos/07_insilico_msms',
                 msfinder_path = 'F:/software/MSFINDER/MSFINDER_ver_3.24/MsfinderConsoleApp.exe')

```

## Citation
This free open-source software implements academic research by the authors and co-workers. If you use it, please support the project by citing the appropriate journal articles. For other in-silico tools, please cite their origin papers.

Zhiwei Zhou†, Mingdu Luo†, Haosong Zhang, Yandong Yin, Yuping Cai, and Zheng-Jiang Zhu*, Metabolite annotation from knowns to unknowns through knowledge-guided multi-layer metabolic network, **Submitted**, 2022. [**bioRxiv**](https://doi.org/10.1101/2022.06.02.494523)

## License
<a rel="license" href="https://creativecommons.org/licenses/by-nc-nd/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-nd/4.0/88x31.png" /></a> 
This work is licensed under the Attribution-NonCommercial-NoDerivatives 4.0 International (CC BY-NC-ND 4.0)

