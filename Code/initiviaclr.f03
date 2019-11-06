!*************************************************************************
!*                                                                       *
!*     Initialization of parameters for IviaCLR (imputation via          *
!*     nonsmooth clusterwise linear regression) software                 *
!*     (last modified 12.02.2019)                                        *
!*                                                                       *
!*************************************************************************
!*
!*     Modules included:
!*
!*     initimp          ! Initialization of parameters for imputation.
!*     initclr          ! Initialization of parameters for regerssion.
!*     initlmbm         ! Initialization of LMBM -solver.
!*


MODULE initimp  ! Initialization of parameters for imputation.

    USE r_precision, ONLY : prec            ! Precision for reals.
    IMPLICIT NONE

    ! Names of input and output files:
    CHARACTER(LEN=80), SAVE :: &
        infile, &                           ! Name of the dataset file. Given as commandline argument.
        outfile1 = 'hyperplanes.txt', &     ! Regression results with hyperplanes.
        outfile2 = 'predictions.txt', &     ! Regression results with predictions.
        outfile3 = 'imputed_values.txt', &  ! Result file with imputed values.
        outfile4                            ! Complete data set with imputed values. Given as commandline argument.

    ! Input real parameters. Give values here or as on optional commandline argument.
    REAL(KIND=prec), SAVE :: &
        misval = -9999_prec                 ! Missing value: every value in infile equal to misval
                                            !   is considered to be missing (needs to be numeric).

    ! Input integer parameters.
    INTEGER, SAVE :: &
        inimp = 0, &                        ! Initial imputation:
                                            !   0 - mean imputation (default),
                                            !   1 - regression on complete objects,
                                            !   2 - recursive regression.
        ipredi = 1, &                       ! Prediction method:
                                            !   0 - k-NN with popularity based weights on clusters,
                                            !   1 - k-NN with rmsq based weights on clusters (default),
                                            !   2 - the weight based on popularity of clusters (don't use).
        knn = 5, &                          ! Number of nearest neighbours.
        maxrounds = 5, &                    ! Maximum number of outerloops maxrounds >=1.
        nft, &                              ! Number of features in data, from user.
        nrecord, &                          ! Number of records in data, from user.
        nclust, &                           ! Maximum number of clusters, from user.
        nsample                             ! Size of the subset taken to compute k-NN: nsample <= nrecord.
                                            ! Default nsample = MIN(150,nrecord). May be changed in def_par.


    ! Other Real parameters.
    REAL(KIND=prec), SAVE, DIMENSION(:,:), allocatable :: &
        a, &                                ! Data set, from input file.
        a2, &                               ! Original data set.
        a3, &                               ! "Complete" deleted data set.
        a4                                  ! Auxiliary data matrix.

    ! Other integer parameters.
    INTEGER, SAVE :: & !
        noutcom, &                          ! Index of the output column.
        nrecordtmp, &                       ! Number of complete objects.
        nmis                                ! Number of missing values in data set.

    ! Allocatable tables
    INTEGER, SAVE, DIMENSION(:), allocatable :: &
        nalk, &                             ! nalk(i) = number of data points in cluster i.
        clalk, &                            ! clalk(i) = cluster where data point i belongs to.
        mis_ins, &                          ! Missing instances.
        mis_ft, &                           ! Missing features.
        nmis_ins, &                         ! Number of missing features in each data point
        nmis_ft                             ! Number of missing values in each feature

CONTAINS


    SUBROUTINE init_imppar()  ! Initialization of parameters (when needed). May be left empty.

        IMPLICIT NONE


    END SUBROUTINE init_imppar

    SUBROUTINE def_imppar()  ! Default values for parameters.

        USE param, ONLY : zero,one          ! Parameters.

        IMPLICIT NONE

        nsample = MIN(150,nrecord)          ! Size of the subset taken to compute k-NN. nsample <= nrecord.
        IF (nclust == 0) maxrounds = 1

    END SUBROUTINE def_imppar

END MODULE initimp


MODULE initclr  ! Initialization of parameters for clusterwise linear regression codes.

    USE r_precision, ONLY : prec            ! Precision for reals.
    USE initimp
    IMPLICIT NONE

    ! Input real parameters, give value here.
    REAL, SAVE :: &
        tlimit    = 7000.0                  ! Maximum CPU-time in seconds.

    ! Input integer parameters.
    INTEGER, SAVE :: &
        iscale = 1                          ! Choose:
                                            !   1: if scaling is used (default),
                                            !   2: if no scaling is used.

    ! Other Real parameters.
    REAL(KIND=prec), SAVE :: & !
        gamma1, &                           ! Parameter affecting the size of the set of starting point,
        gamma2, &                           ! Parameter affecting the size of the set of starting point,
        gamma3, &                           ! Parameter affecting the size of the set of starting point,
        gamma4, &                           ! Parameter affecting the size of the set of starting point,
        gam1, &                             ! Multiplier for gamma1
        gam2, &                             ! Multiplier for gamma2
        gam3, &                             ! Multiplier for gamma3
        gam4, &                             ! Multiplier for gamma4
        toler

    REAL(KIND=prec), SAVE, DIMENSION(:), allocatable :: & !
        dminim, &       !
        clin, &
        dlin


    ! Other integer parameters.
    INTEGER, SAVE :: & !
        maxdim, &                           ! Maximum number of variables in optimization,
                                            !   maxdim = MAX(nft,nclust*nft).
        maxinit, &                          ! Maximum number of candidate points of data set,
                                            !   maxinit >> maxdim.
                                            ! Default 10000*maxdim. May be changed in def_par.
        nupdate, &                          ! Number of updates.
        naux2, &                            ! Number of first problem evaluations.
        naux, &                             ! Number of auxiliary function evaluations.
        nfit, &                             ! Number of function evaluations.
        nc, &                               ! Current number of clusters, loops from 1 to nclust
        ns, &                               ! Switch for auxiliary and real clustering problem.
        m, &                                ! Number of variables in optimization:
                                            !   m = nft    if ns = 1,
                                            !   m = nft*nc if ns = 2.
        nel1                                !

    INTEGER, SAVE, DIMENSION(:,:), allocatable :: & !
        nob

    INTEGER, SAVE, DIMENSION(:), allocatable :: & !
        nob1, &                             !
        nel                                 !



CONTAINS


    SUBROUTINE init_clrpar()  ! User supplied subroutine for further initialization of parameters (when needed).
                              ! May be left empty.
        IMPLICIT NONE

    END SUBROUTINE init_clrpar

    SUBROUTINE def_clrpar()   ! Default values for parameters.

        USE param, ONLY : zero,one          ! Parameters.

        IMPLICIT NONE

        maxinit  = 10000*maxdim             ! Maximum number of candidate points of data set,
                                            !   maxinit >> maxdim. Default 10000*maxdim.

        gam1 = one
        gam2 = one
        gam3 = one
        gam4 = one

        IF (nrecord <  200) THEN
            gamma1 = 0.00_prec  ! gamma1 in [0,1], smaller value -> more data points
            gamma2 = 20.0_prec  ! gamma2 in [1, infty) -> smaller value less minimizers are chosen
            gamma3 = 20.0_prec  ! gamma3 in [1, infty) -> smaller value less minimizers are chosen
            gamma4 = 0.10_prec

        ELSE IF (nrecord < 1000) THEN !
            gamma1 = 0.50_prec
            gamma2 = 2.5_prec
            gamma3 = 2.25_prec
            gamma4 = 0.50_prec

        ELSE IF (nrecord < 4000) THEN
            gamma1=0.6_prec
            gamma2=2.50_prec
            gamma3=1.250_prec
            gamma4=1.0_prec

        ELSE IF (nrecord < 15000) THEN
            gamma1=0.95_prec
            gamma2=2.0_prec
            gamma3=1.150_prec
            gamma4=1.0_prec

        ELSE IF (nrecord < 50000) THEN
            gamma1=0.9750_prec
            gamma2=2.0_prec
            gamma3=1.050_prec
            gamma4=1.0_prec

        ELSE
            gamma1=0.99_prec
            gamma2=1.01_prec
            gamma3=1.01_prec
            gamma4=1.0_prec

        END IF

        RETURN

    END SUBROUTINE def_clrpar

END MODULE initclr


MODULE initlmbm  ! Initialization of parameters for LMBM.

    USE r_precision, ONLY : prec           ! Precision for reals.
    USE param, ONLY : zero, one            ! Parameters.
    USE initimp

    IMPLICIT NONE

    ! Parameters
    INTEGER, PARAMETER :: &
        na      = 2, &                     ! Maximum bundle dimension, na >= 2.
        mcu     = 15, &                    ! Maximum number of stored corrections, mcu >=1.
        mcinit  = 7, &                     ! Initial maximum number of stored corrections, mcu >= mcinit >= 3.
                                           ! If mcinit <= 0, the default value mcinit = 3 will be used.
                                           ! However, the value mcinit = 7 is recommented.
        inma    = 3, &                     ! Selection of line search method:
                                           !   inma = 0, Armijo line search,
                                           !   inma = 1, nonmonotone Armijo line search.
                                           !   inma = 2, weak Wolfe line search.
                                           !   inma = 3, nonmonotone  weak Wolfe line search.
        mnma    = 10, &                    ! Maximum number of function values used in nonmonotone line search.
        maxnin  = 20                       ! Maximum number of interpolations, maxnin >= 0.
                                           ! The value maxnin = 2 - 20 is recommented with inma=0,
                                           ! maxnin >= 20 with inma=1 and 3, and maxnin =200 with inma=2.
                                           ! For example:
                                           !   inma = 0, maxin = 20.
                                           !   inma = 1, mnma = 20, maxin = 30.
                                           !   inma = 2, maxnin = 200.
                                           !   inma = 3, mnma=10, maxnin = 20.


    ! Real parameters (if parameter value <= 0.0 the default value of the parameter will be used).
    REAL(KIND=prec), SAVE :: &
        tolb    = zero, &                  ! Tolerance for the function value (default = -large).
        tolf    = zero, &                  ! Tolerance for change of function values (default = 1.0E-8).
        tolf2   = -1.0_prec, &             ! Second tolerance for change of function values.
                                           !   - If tolf2 < 0 the the parameter and the corresponding termination
                                           !   criterion will be ignored (recommended with inma=1,3).
                                           !   - If tolf2 = 0 the default value 1.0E+4 will be used.
        tolg    = 1.0E-5_prec, &           ! Tolerance for the termination criterion (default = 1.0E-5).
        tolg2   = 1.0E-3_prec, &           ! Tolerance for the second termination criterion (default = 1.0E-3).
        eta     = 1.0E-3_prec, &           ! Distance measure parameter, eta > 0.
                                           !   - If eta < 0  the default value 0.0001 will be used.
        epsl    = 0.22E+00, &              ! Line search parameter, 0 < epsl < 0.25 (default = 0.24).
                                           ! use small value with inma=2 and larger value with inma=3
        xmax    = 1000.0_prec              ! Maximum stepsize, 1 < XMAX (default = 1000).

    ! Integer parameters (if value <= 0 the default value of the parameter will be used).
    INTEGER, SAVE :: &
        n,  &                              ! Number of variables.
        mit     = 5000, &                  ! Maximun number of iterations.
        mfe     = 5000, &                  ! Maximun number of function evaluations.
        mtesf   = 20, &                    ! Maximum number of iterations with changes of
                                           ! function values smaller than tolf (default = 20).
        iprint  = 0, &                     ! Printout specification:
                                           !    -1  - No printout.
                                           !     0  - Only the error messages.
                                           !     1  - The final values of the objective function
                                           !          (default used if iprint < -1).
                                           !     2  - The final values of the objective function and the
                                           !          most serious warning messages.
                                           !     3  - The whole final solution.
                                           !     4  - At each iteration values of the objective function.
                                           !     5  - At each iteration the whole solution
        iscale  = 0                        ! Selection of the scaling with LMBM:
                                           !     0  - Scaling at every iteration with STU/UTU (default).
                                           !     1  - Scaling at every iteration with STS/STU.
                                           !     2  - Interval scaling with STU/UTU.
                                           !     3  - Interval scaling with STS/STU.
                                           !     4  - Preliminary scaling with STU/UTU.
                                           !     5  - Preliminary scaling with STS/STU.
                                           !     6  - No scaling.

    REAL(KIND=prec), SAVE, DIMENSION(:), allocatable :: & !
         x                                 ! Vector of variables



CONTAINS

    SUBROUTINE defaults()  ! Default values for parameters.

        USE param, ONLY: small, large, zero, one, half
        IMPLICIT NONE

        IF (iprint < -1) iprint  = 1               ! Printout specification.
        IF (mit   <= 0) mit      = 5000            ! Maximum number of iterations.
        IF (mfe   <= 0) mfe      = n*mit           ! Maximum number of function evaluations.
        IF (tolf  <= zero) tolf  = 1.0E-08_prec    ! Tolerance for change of function values.
        IF (tolf2 == zero) tolf2 = 1.0E+04_prec    ! Second tolerance for change of function values.
        IF (tolb  == zero) tolb  = -large + small  ! Tolerance for the function value.
        IF (tolg  <= zero) tolg  = 1.0E-05_prec    ! Tolerance for the termination criterion.
        IF (tolg2 <= zero) tolg = 1.0E-03_prec     ! Tolerance for the second termination criterion.
        IF (xmax  <= zero) xmax  = 1000.0_prec     ! Maximum stepsize.
        IF (eta   <  zero) eta   = 1.0E-03_prec    ! Distance measure parameter
        IF (epsl  <= zero) epsl  = 0.22_prec       ! Line search parameter,
        IF (mtesf <= 0) mtesf    = 20              ! Maximum number of iterations with changes
                                                   ! of function values smaller than tolf.
        IF (iscale > 6 .OR. iscale < 0) iscale = 0 ! Selection of the scaling.


    END SUBROUTINE defaults

    SUBROUTINE init_lmbmpar()  ! User supplied subroutine for further initialization of parameters
                               ! (when needed) for LMBM. May be left empty.

        USE initclr, ONLY : ns,nrecord
        IMPLICIT NONE

        ! For big data sets use larger values
        IF (nrecord < 500) THEN ! small data sets

            IF (ns == 3) THEN !
                mit=5000
                tolg2= 1.0E-12
                tolg=1.0E-7
            ELSE IF (ns == 1) THEN
                mit=5000
                tolg2= 1.0E-12
                tolg = 1.0E-5
            ELSE
                mit=5000
                tolg2 = 1.0E-8
                tolg = 1.0E-5
            END IF

        ELSE IF (nrecord < 4000) THEN ! medium data sets
            IF (ns == 3) THEN !
                mit=5000
                tolg2= 1.0E-12
                tolg=1.0E-7
            ELSE IF (ns == 1) THEN
                mit=5000
                tolg2 = 1.0E-6
                tolg = 1.0E-5 !
            ELSE
                mit=5000
                tolg2 = 1.0E-6
                tolg = 1.0E-5 !
            END IF

        ELSE IF (nrecord < 50000) THEN ! large datasets
            IF (ns == 3) THEN !
                mit=5000
                tolg2= 1.0E-12 !32
                tolg=1.0E-7
            ELSE IF (ns == 1) THEN
                mit=5000
                tolg2 = 1.0E-5 !64
                tolg = 1.0E-6 !
            ELSE
                mit=5000
                tolg2 = 1.0E+1
                tolg = 1.0E-5 !
            END IF

        ELSE ! very large datasets
            IF (ns == 3) THEN !
                mit=5000
                tolg2= 1.0E-8 !
                tolg=1.0E-7
            ELSE IF (ns == 1) THEN
                mit=5000
                tolg2 = 1.0E-3 !
                tolg = 1.0E-6 !
            ELSE
                mit=5000
                tolg2 = 1.0E+2
                tolg = 1.0E-4 !
            END IF

        END IF


    END SUBROUTINE init_lmbmpar

END MODULE initlmbm



