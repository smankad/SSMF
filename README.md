# SSMF
Data and Code for the paper "Single Stage Prediction with Embedded Topic Modeling of Online Reviews for Mobile App Management"

The "Data" subfolder contains the raw data used in the paper. 

The "Code" subfolder contains 3 implementations of the proposed matrix factorization. 
  + ssmf.R  performs the factorization given in equation (1) of the paper; 
  + ssmf_glmnet.R performs the factorization given in equation (1) of the paper PLUS a lasso penalty on the beta; 
  + ssmf_continuation_ratio.R performs the full continuation ratio model with potentially multiple apps and time points; 
  + tester.R runs a simple check that the functions compile and run using the "crude" dataset contained in the "tm" package. 

If you find the code useful or in need of improvement, please let me know at spm263@cornell.edu.

Referenced Paper: https://arxiv.org/abs/1607.07515

@ARTICLE{2016arXiv160707515M,
   author = {{Mankad}, S. and {Hu}, S. and {Gopal}, A.},
    title = "{Single Stage Prediction with Embedded Topic Modeling of Online Reviews for Mobile App Management}",
  journal = {ArXiv e-prints},
archivePrefix = "arXiv",
   eprint = {1607.07515},
 primaryClass = "stat.AP",
 keywords = {Statistics - Applications, Computer Science - Information Retrieval, Computer Science - Software Engineering},
     year = 2016,
    month = jul,
   adsurl = {http://adsabs.harvard.edu/abs/2016arXiv160707515M},
  adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

