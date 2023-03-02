# SRRMLR

This repository contains all the necesaary codes for duplicating the results in the article ``Simultaneous Dimension Reduction and Variable Selection for Multinomial Logistic Regression'', authored by Canhong Wen, Zhenduo Li, Ruipeng Dong, Yijin Ni, and Wenliang Pan.


## Codes

- __simulation__ folder contains the code used for numerical experiments.

    * The filefolder __conv_plot__ and __aset_plot__ contain the replication code for convergence analysis:
        *  __conv.R__ and __conv.ipynb__ for Figures 1-2.
        *  __aset.R__ for Figure 3. 
        *  __aset-n.R__ for Figure 4.
    * The __Function.R__ file contails all the functions needed.
    * The __simulation.R__ file replicates the results in Sections 4.1-4.2.
    
 - __src_test__ folder contains the source files for implementing the SRRMLR algorithm. 

    * The source filefolder __SRRMLR__. It is only used when one want to modify the SRRMLR algorithm.    
    * The R package source file __SRRMLR_1.0.tar.gz__. After downloading it, you need to run the following code in R to install it.
    
            install.packages("Your_download_path/SRRMLR_1.0.tar.gz", repos = NULL)  
            
    * The __test.R__ file provides a synthetic test instance.

- __mnist__ folder contains the code used for real data analysis in Section 5.

- __HAM__ folder contains the code used for HAM10000 data analysis in Appendix C of the supplementary material.

    


- __warm_start__ folder provides some warm-start trials.


## Cite

Please cite the following publications if you make use of the material here.

- Canhong Wen, Zhenduo Li, Ruipeng Dong, Yijin Ni, and Wenliang Pan. Simultaneous Dimension Reduction and Variable Selection for Multinomial Logistic Regression. 

Below is the BibTex for citing the data.

```
@article{srrmlr,
author =        {Wen, Canhong and Li, Zhenduo and Dong, Ruipeng and Ni, Yijin and Pan, Wenliang},
title =         {Simultaneous Dimension Reduction and Variable Selection for Multinomial Logistic Regression},
year =          {2023}
}  
```


## Contact
Please direct questions and comments to the [issues page](https://github.com/C2S2-HF/SRRMLR/issues).
