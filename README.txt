!*************************************************************************
!*                                                                       *
!*     IviaCLR - Imputation method for incomplete data preprocessing     *
!*                 using Nonsmooth Clusterwise Linear Regression         *                                                                   
!*                                                                       *
!*     by Napsu Karmitsa 2017 (last modified 19.11.2019)                 *
!*                                                                       *
!*     The software is free for academic teaching and research           *
!*     purposes but I ask you to refer the appropriate references        *
!*     given below, if you use it.                                       *
!*                                                                       *
!*     *******************************************************************
!*     The IviaCLR method proposed in the following paper:               *
!*     "Clusterwise Linear Regression based Missing Value Imputation     * 
!*     for Data Preprocessing" by authours:                              *
!*     "N. Karmitsa, S. Taheri, A. Bagirov, P. Mäkinen"                  *
!*     is submitted to a journal and its author version is:              *
!*     Iviaclr_tucs.pdf which is attached.                               *
!*************************************************************************
!*
!*
!*     Codes included:
!*
!*
!*     iviaclr.f03           - Mainprogram for imputation via clusterwise
!*                             linear regression (this file).
!*     initiviaclr.f03       - initialization of parameters for imputation,
!*                             clusterwise linear regression, and LMBM.
!*                             Includes modules:
!*                               - initimp   - Initialization of imputation procedure.
!*                               - initclr   - Initialization of parameters for
!*                                             clusterwise linear regeression.
!*                               - initlmbm  - Initialization of LMBM.
!*     parameters.f03        - Parameters. Inludes modules:
!*                               - r_precision - Precision for reals,
!*                               - param - Parameters
!*                               - exe_time - Execution time.
!*     imput.f03             - initial imputations and final predictions.
!*     clrsubpro.f03         - Subroutines for clusterwise linear regression.
!*     clrfunctionmod.f03    - Computation of function and (sub)gradients values.
!*     lmbm.f03              - LMBM - limited memory bundle method.
!*     lmbmsubpro.f03        - subprograms for LMBM.
!*
!*     Makefile              - makefile.
!*
!*
!*     Calling sequence:
!*
!*       ./iviaclr infile outfile nrecord nft maxclust [init_imp predi misvalue]
!*
!*     with the following arguments:
!*
!*     infile                - the data set with missing values (-9999 by default).
!*     outfile               - name of the imputed data set.
!*     nrecord               - number of records in the data set.
!*     nft                   - number of features in the data set.
!*     maxclust              - maximum number of linear regression functions.
!*     init_imp              - optional, initial imputation method, mean by default.
!*                             Possible values:
!*                               0 - mean imputation (default),
!*                               1 - regression on complete objects,
!*                               2 - recursive regression.
!*     predi                 - optional, prediction method.
!*                             Possible values:
!*                               0 - k-NN with popularity based weights on clusters.
!*                               1 - k-NN with rmsq based weights on clusters (default).
!*                               2 - the weight based on popularity of clusters (don't use).
!*     misvalue              - optional, value for missing value, must be a number,
!*                             -9999 by default.
!*
!*
!*     An example >> ./iviaclr iris_missing.txt iris_imputed.txt 150 4 3 0 1 -9999
!*
!*
!*     References:
!*     for IviaCLR:
!*       N. Karmitsa, S. Taheri,  A. Bagirov, P. Mäkinen, "Missing Value Imputation 
!*       via Nonsmooth Clusterwise Linear Regression", submitted.
!*
!*       N. Karmitsa, S. Taheri, A. Bagirov, P. Mäkinen, "Clusterwise Linear Regression 
!*       based Missing Value Imputation for Data Preprocessing", TUCS Technical Report, 
!*       No. 1193, Turku Centre for Computer Science, Turku, 2018.
!*
!*     for LMBM-CLR:
!*       N. Karmitsa, A. Bagirov, S. Taheri, "Clusterwise Linear Regression using LMBM"
!*
!*     for LMBM:
!*       N. Haarala, K. Miettinen, M.M. Mäkelä, "Globally Convergent Limited Memory Bundle Method  
!*       for Large-Scale Nonsmooth Optimization", Mathematical Programming, Vol. 109, No. 1,
!*       pp. 181-205, 2007. DOI 10.1007/s10107-006-0728-2.
!*
!*       M. Haarala, K. Miettinen, M.M. Mäkelä, "New Limited Memory Bundle Method for Large-Scale 
!*       Nonsmooth Optimization", Optimization Methods and Software, Vol. 19, No. 6, pp. 673-692, 2004. 
!*       DOI 10.1080/10556780410001689225.
!*
!*       N. Karmitsa, "Numerical Methods for Large-Scale Nonsmooth Optimization" in Big Data
!*       Optimization. A. Emrouznejad (eds.), Springer, 2016.
!*
!*     for NSO clusterwise linear regression:
!*       A. Bagirov, N. Karmitsa, M.M. Mäkelä, "Introduction to nonsmooth optimization: theory, 
!*       practice and software", Springer, 2014.
!*
!*
!*     Acknowledgements:
!*
!*     The work was financially supported by the Academy of Finland (Project No. 289500).
!*
