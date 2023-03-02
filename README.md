# SRRMLR

This repository contains all the necesaary codes for duplicating the results in the article ``Simultaneous Dimension Reduction and Variable Selection for Multinomial Logistic Regression'', authored by Canhong Wen, Zhenduo Li, Ruipeng Dong, Yijin Ni, and Wenliang Pan.


## Codes

- __simulation__ folder contains the code used for numerical experiments.

- __MNIST__ folder contains the code used for real data analysis in Section 5.

- __HAM__ folder contains the code used for HAM10000 data analysis in Appendix C of the supplementary material.

    
- __src_test__ folder contains the source files for implementing the SRRMLR algorithm. 

    * The source filefolder __SRRMLR__. It is only used when one want to modify the SRRMLR algorithm.
      
    * After downloading it, you need to run the following code in R to install it.

            install.packages("Your_download_path/SRRMLR_1.0.tar.gz", repos = NULL)
    * A synthetic test instance

- __warm_start__ folder provides some warm-start trials.


## Cite

To cite this data, please cite the [research article](https://doi.org/10.1287/????) and the data itself, using the following DOI.

[![DOI](https://zenodo.org/badge/524492449.svg)](https://zenodo.org/badge/latestdoi/524492449)


Below is the BibTex for citing the data.

```
@article{Data_SDO,
author =        {{\c{S}}eker, Oylum and Cevik, Mucahit and Bodur, Merve and Lee, Young and Ruschin, Mark},
publisher =     {INFORMS Journal on Computing},
title =         {Data for Sector Duration Optimization in Stereotactic Radiosurgery Treatment Planning},
year =          {2022},
doi =           {10.5281/zenodo.7048848},
url =           {https://github.com/INFORMSJoC/2021.0103},
}  
```



## Content

This repository includes

1. A synthetic test instance for sector duration optimization (actual patient data used in the paper is not shared due to privacy concerns)
1. A file describing instance data format  
1. A file containing sample results for the provided instance


Instance data files are located in the folder [data](data). The folder contains five text files.

* There are four different text files containing dose rate matrices for different structures, which are one tumor & ring, and two organs-at-risk (OARs) (dose rate data are generated using two predetermined isocenter locations).
    1. Dose rate file for the tumor: [doseRateMatrix_tumor.txt](data/doseRateMatrix_tumor.txt)
    1. Dose rate file for the ring: [doseRateMatrix_ring.txt](data/doseRateMatrix_ring.txt)
    1. Dose rate file for OAR1: [doseRateMatrix_OAR1.txt](data/doseRateMatrix_OAR1.txt)
    1. Dose rate file for OAR2: [doseRateMatrix_OAR2.txt](data/doseRateMatrix_OAR2.txt) 
