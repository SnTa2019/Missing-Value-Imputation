!*************************************************************************
!*                                                                       *
!*     IviaCLR - Imputation method for incomplete datasets using         *
!*               Nonsmooth Clusterwise Linear Regression                 *
!*                                                                       *
!*     by Napsu Karmitsa 2017 (last modified 12.02.2019).                *
!*                                                                       *
!*     The software is free for academic teaching and research           *
!*     purposes but I ask you to refer the appropriate references        *
!*     given below, if you use it.                                       *
!*                                                                       *
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
!*       ./iviaclr infile outfile nrecord nft maxclust [init_imp predi misvalue]
!*
!*     with the following arguments:
!*
!*     infile                - the data set with missing values (-9999 by default).
!*     outfile               - name of the imputed data set.
!*     nrecord               - number of records in the data set.
!*     nft                   - number of features in the data set.
!*     maxclust              - maximum number of linear regression functions.
!*     init_imp              - initial imputation method, optional, mean by default.
!*                             Possible values:
!*                               0 - mean imputation (default),
!*                               1 - regression on complete objects,
!*                               2 - recursive regression.
!*     predi                 - prediction method, optional.
!*                             Possible values:
!*                               0 - k-NN with popularity based weights on clusters.
!*                               1 - k-NN with rmsq based weights on clusters (default).
!*                               2 - the weight based on popularity of clusters (don't use).
!*     misvalue              - value for missing value, must be a number,
!*                             optional, -9999 by default.
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
!*****************************************************************************************
!*
!*     * PROGRAM iviaclr *
!*
!*     Main program for imputation via nonsmooth clusterwise linear regression.
!*
!*****************************************************************************************

PROGRAM iviaclr

    USE r_precision, ONLY : prec        ! Precision for reals.
    USE param, ONLY : zero, one, large  ! Parameters.

    USE initimp, ONLY : &               ! Initialization of parameters.
        infile, &                       ! Name of the dataset file with missing values.
        outfile1, &                     ! Regression results with hyperplanes.
        outfile2, &                     ! Regression results with predictions.
        outfile3, &                     ! Result file with imputed values.
        outfile4, &                     ! Complete imputed data set.
        inimp, &                        ! Initial imputation:
                                        !   0 - mean imputation,
                                        !   1 - regression on complete objects.
                                        !   2 - recursive regression.
        ipredi, &                       ! Prediction method:
                                        !   0 - k-NN with popularity based weights on clusters.
                                        !   1 - k-NN with rmsq based weights on clusters.
                                        !   2 - weight based on popularity of clusters.
        misval, &                       ! Missing value: every value in infile equal to misval
                                        !   is considered to be missing (needs to be numeric).
        a, &                            ! Data matrix a(nrecord,nft), from input file.
        a2,a3,a4, &                     ! Auxiliary data files.
        nclust, &                       ! Maximum number of clusters, from user.
        nft, &                          ! Number of features in data, from user.
        nrecord, &                      ! Number of instances in data, from user.
        noutcom, &                      ! Output column.
        nmis_ft, &                      ! Number of missing values in features.
        maxrounds, &                    ! Maximum number of loops.
        init_imppar, &                  ! Further initilaization of imputation.
        def_imppar, &                   ! Default parameters for imputation.
        nalk, &                         ! nalk(i) = number of data points in cluster i
        clalk, &                        ! clalk(i) = cluster where data point i belongs to.
        mis_ins, &                      ! Missing instances.
        mis_ft, &                       ! Missing features.
        nmis_ins, &                     ! Number of missing features in each data point.
        nmis_ft                         ! Number of missing values in each feature.

    USE firstimput, ONLY : &            ! Initial imputation.
        read_data, &                    ! Reading input data.
        mean_imput, &                   ! Mean imputation.
        reg_imput, &                    ! Regression on complete objects.
        rreg_imput                      ! Recursive regression.

    USE prediction, ONLY : &            ! Predictions.
        pred_weight, &                  ! Form complete data set by weighted imputations.
        pred_knn, &                     ! Form complete data set by imputing missing values
                                        !   using k-NN with popularity based weights on clusters.
        pred_knn2, &                    ! (not in use)
        pred_knnr                       ! Form complete data set by imputing missing values
                                        !   using k-NN with rmsq based weights on clusters.

    USE initclr, ONLY : &               ! Initialization of parameters.
        maxdim, &                       ! Maximum number of variables for optimization.
        maxinit, &                      ! Maximum number of candidate points of data set.
        gamma1, &                       ! Parameter affecting the size of the set of starting point.
        gamma2, &                       ! Parameter affecting the size of the set of starting point.
        gamma3, &                       ! Parameter affecting the size of the set of starting point.
        gamma4, &                       ! Parameter affecting the size of the set of starting point.
        gam1, &                         ! Multiplier for gamma1.
        gam2, &                         ! Multiplier for gamma2.
        gam3, &                         ! Multiplier for gamma3.
        gam4, &                         ! Multiplier for gamma4.
        dminim, &                       !
        toler, &                        !
        tlimit, &                       ! Time limit, from user.
        iscale, &                       ! Selection of scaling:
                                        !   1: scaling (default),
                                        !   2: no scaling.
        nupdate, &                      ! Number of updates.
        naux, &                         ! Number of auxiliary function evaluations.
        naux2, &                        ! Number of first auxiliary function evaluations.
        nfit, &                         ! Number of function evaluations.
        nc, &                           ! Current number of clusters, loops from 1 to nclust.
        ns, &                           ! Switch for auxiliary and real clustering problem.
        nel1, &                         !
        nob1, &                         !
        nob, &                          !
        nel, &                          !
        clin, &                         !
        dlin, &                         !
        init_clrpar, &                  ! Furher initialization of parameters.
        def_clrpar                      ! Default values of clustering parameters.

    USE clrsubpro, ONLY : &             ! Subprograms for clusterwise linear regression.
        scaling, &                      ! Scaling for CLR procedure.
        rescaling, &                    ! Rescaling.
        distribution, &                 !
        step1, &                        !
        cleaning, &                     !
        checkfinal                      ! Final result from CLR.

    USE clrfunctionmod, ONLY : &
        clrfunc, &                      ! Function calculus.
        auxfunc                         ! Function calculus.

    USE initlmbm, ONLY : xlmbm => x     ! x for optimization purposis.

    USE lmbm_mod, ONLY : optim2         ! LMBM method for optimization.

    USE exe_time, ONLY : getime         ! Execution time.

    IMPLICIT NONE


    REAL(KIND=prec), DIMENSION(:), allocatable :: &
        x, &
        x1, &
        xbest, &
        x4, &
        z, &
        d8, &
        fv, &
        fv2, &
        xc

    REAL(KIND=prec) :: &
        atmp, &
        tol2, &
        p2,p3,p4, &
        d7,d9, &
        f,f2,f3, &
        fv1, &
        fbest,fb, &
        fbest1, &
        fmin,fmin1, &
        ff,rmse1,ce1
    REAL :: &
        time1, &                        ! Starting time
        time4, &                        ! Finishing time
        timef                           ! CPU time used
    INTEGER :: &
        count, &                        ! Number of commandline arguments
        ifea, &
        i,j,k, &
        jpoints, &
        ncount, &
        neltotal, &
        ngood0, &
        naux2tot, &
        nauxtot, &
        nfittot, &
        iround

    INTEGER, DIMENSION(:), allocatable :: &
        nmis_ft2

    CHARACTER(len=50) :: arg

    ! Intrinsic Functions
    INTRINSIC MAX,MAXLOC

    count = command_argument_count()
    IF (count < 1) THEN
    PRINT*
    PRINT*,'IviaCLR - Imputation method for incomplete datasets using Nonsmooth Clusterwise Linear Regression'
    PRINT*
    PRINT*,'Calling sequence:'
    PRINT*
    PRINT*,'   ./iviaclr infile outfile nrecord nft maxclust [init_imp predi misvalue]'
    PRINT*
    PRINT*,'with the following arguments:'
    PRINT*
    PRINT*,'   infile    - the data set with missing values'
    PRINT*,'               (missing value = -9999 by default).'
    PRINT*,'   outfile   - name of the imputed data set.'
    PRINT*,'   nrecord   - number of records in the data set.'
    PRINT*,'   nft       - number of features in the data set.'
    PRINT*,'   nclust    - maximum number of linear regression functions.'
    PRINT*,'   initimp   - optional, initial imputation method, mean by default,'
    PRINT*,'               possible values:'
    PRINT*,'                 0 - mean imputation (default),'
    PRINT*,'                 1 - regression on complete objects,'
    PRINT*,'                 2 - recursive regression.'
    PRINT*,'   predi     - optional, prediction method, k-NN with rmsq based weights on clusters by default,'
    PRINT*,'               possible values:'
    PRINT*,'                 0 - k-NN with popularity based weights on clusters,'
    PRINT*,'                 1 - k-NN with rmsq based weights on clusters (default),'
    PRINT*,'                 2 - the weight based on popularity of clusters (do not use).'
    PRINT*,'   misvalue  - optional, value for missing value, must be a number, -9999 by default.'
    PRINT*
    STOP

    END IF

    IF (count < 5) STOP '  Too few commandline arguments. Type ./iviaclr for help.'

    CALL get_command_argument(3, arg)
    READ (arg,*) nrecord
    CALL get_command_argument(4, arg)
    READ (arg,*) nft
    CALL get_command_argument(5, arg)
    READ (arg,*) nclust

    IF (count >= 6) THEN
        CALL get_command_argument(6, arg)
        READ (arg,*) inimp
    END IF
    IF (count >= 7) THEN
        CALL get_command_argument(7, arg)
        READ (arg,*) ipredi
    END IF
    IF (count >= 8) THEN
        CALL get_command_argument(8, arg)
        READ (arg,*) misval
    END IF


    CALL get_command_argument(1, infile)
    CALL get_command_argument(2, outfile4)

    PRINT*,' Input file: ',infile
    PRINT*,' Output file: ',outfile4
    PRINT*,' Number of records = ',nrecord
    PRINT*,' Number of features = ',nft
    PRINT*,' Number of regression functions = ',nclust

    maxdim=MAX(nft,nclust*nft)
    CALL init_imppar()
    CALL def_imppar()


    CALL init_clrpar()
    CALL def_clrpar()

    allocate(x(maxdim),xc(maxinit),x1(maxdim),x4(maxdim),xbest(maxdim),z(maxdim), &
        d8(nrecord),fv(nrecord),fv2(nrecord))
    allocate(a(nrecord,nft),a2(nrecord,nft),a3(nrecord,nft),a4(nrecord,nft))
    allocate(nalk(nclust),clalk(nrecord),mis_ins((nrecord-1)*nft),mis_ft((nrecord-1)*nft), &
        nmis_ins(nrecord),nmis_ft(nft))
    allocate(clin(nft),dlin(nft),dminim(nrecord),nob1(nrecord),nob(nclust,nrecord),nel(nclust))
    allocate(xlmbm(maxdim),nmis_ft2(nft))

    tol2=1.0E-7_prec*REAL(nrecord,prec)
    IF (nrecord > 15000) tol2 = 1.0E-4_prec*REAL(nrecord,prec)

    dminim=large

        OPEN(40,file=outfile1)
        OPEN(41,file=outfile2)
        OPEN(43,file=outfile3)
        OPEN(44,file=outfile4)

        SELECT CASE (inimp)

            CASE(0)
                WRITE(40, *) 'IviaCLR with mean as initial imputation and'
                WRITE(41, *) 'IviaCLR with mean as initial imputation and'
                WRITE(43, *) 'IviaCLR with mean as initial imputation and'
            CASE(1)
                WRITE(40, *) 'IviaCLR: initial imputation by regression of complete data and'
                WRITE(41, *) 'IviaCLR: initial imputation by regression of complete data and'
                WRITE(43, *) 'IviaCLR: initial imputation by regression of complete data and'
            CASE DEFAULT
                WRITE(40, *) 'IviaCLR: initial imputation by recursive regression and'
                WRITE(41, *) 'IviaCLR: initial imputation by recursive regression and'
                WRITE(43, *) 'IviaCLR: initial imputation by recursive regression and'

        END SELECT

        SELECT CASE (ipredi)

            CASE(0)
                WRITE(40, *) '  and k-NN with popularity based weights on clusters as predictions.'
                WRITE(41, *) '  and k-NN with popularity based weights on clusters as predictions.'
                WRITE(43, *) '  and k-NN with popularity based weights on clusters as predictions.'
            CASE(1)
                WRITE(40, *) '  and k-NN with rmsq based weights on clusters as predictions.'
                WRITE(41, *) '  and k-NN with rmsq based weights on clusters as predictions.'
                WRITE(43, *) '  and k-NN with rmsq based weights on clusters as predictions.'
            CASE DEFAULT
                WRITE(40, *) '  the weight based on popularity of clusters as predictions.'
                WRITE(41, *) '  the weight based on popularity of clusters as predictions.'
                WRITE(43, *) '  the weight based on popularity of clusters as predictions.'

        END SELECT

        CALL getime(time1)

        CALL read_data()

        WRITE(40,*)
        WRITE(40,*) 'Number of records: ',nrecord
        WRITE(40,*) 'Number of features: ',nft
        WRITE(41,*)
        WRITE(41,*) 'Number of records: ',nrecord
        WRITE(41,*) 'Number of features: ',nft
        WRITE(43,*)
        WRITE(43,*) 'Number of records: ',nrecord
        WRITE(43,*) 'Number of features: ',nft

        IF (inimp == 0) THEN
            CALL mean_imput()
        ELSE IF (inimp == 1) THEN
            CALL reg_imput()
        ELSE
            CALL rreg_imput()
        END IF

        iround = 1
        DO
            IF (iround > maxrounds) EXIT
            PRINT*,' Outer iteration: ',iround
            WRITE(40,*)
            WRITE(40,*) 'Outer iteration: ',iround
            WRITE(41,*)
            WRITE(41,*) 'Outer iteration: ',iround
            WRITE(43,*)
            WRITE(43,*) 'Outer iteration: ',iround

            iround = iround + 1
            IF (timef > tlimit) EXIT


            nmis_ft2 = nmis_ft

            featureloop: DO ifea = 1,nft ! Go through all features

                IF (nclust == 0) EXIT featureloop ! we need not to do this, if nclust=0.
                noutcom = MAXLOC(nmis_ft2,dim=1)

                IF(nmis_ft2(noutcom) == 0) EXIT featureloop ! no more missing values.
                nmis_ft2(noutcom) = 0

                WRITE(40,*)
                WRITE(40,*) 'Feature: ',noutcom
                PRINT*,'   Feature: ',noutcom


                IF (noutcom < nft) THEN ! swap features

                    DO i = 1,nrecord
                        atmp = a(i,noutcom)
                        a(i,noutcom) = a(i,nft)
                        a(i,nft) = atmp
                    END DO
                END IF


                DO k=1,nft
                    DO i=1,nrecord
                        a4(i,k)=a(i,k)
                    END DO
                END DO

                !================================================
                IF (iscale == 1) CALL scaling
                !===============================================

                p2 = zero
                p3 = zero
                p4 = zero
                ncount = 0
                nfittot=0
                nauxtot=0
                naux2tot=0

                innerloop2: DO nc=1,nclust

                    naux = 0
                    naux2 = 0
                    nfit = 0
                    neltotal = 0
                    PRINT*,'     The number of hyperplanes:',nc

                    IF (nc > 1) THEN

                        gamma1=gam1*gamma1
                        gamma2=gam2*gamma2
                        gamma3=gam3*gamma3
                        gamma4=gam4*gamma4
                        IF (gamma1 > one) gamma1=one
                        IF (gamma2 < 1.005) gamma2=1.005_prec
                        IF (gamma3 < 1.005_prec) gamma3=1.005_prec
                        IF (gamma4 > one) gamma4=one

                        CALL step1(x,xc,ngood0)
                        fmin=large
                        DO i=1,ngood0
                            DO j=1,nft
                                z(j)=xc(j+(i-1)*nft)
                            END DO
                            d7=large
                            d9=zero
                            d8(:)=zero

                            DO k=1,nrecord
                                d8(k)=zero
                                DO j=1,nft-1
                                    d8(k)=d8(k)+(a4(k,j)-z(j))**2
                                END DO
                                IF (d8(k) <= d7) d7=d8(k)
                                IF (d8(k) > d9) d9=d8(k)
                            END DO



                            d7=d7+gamma4*(d9-d7)

                            nel1=0
                            DO k=1,nrecord
                                f2=z(nft)
                                DO j=1,nft-1
                                    f2=f2+a4(k,j)*z(j)
                                END DO
                                f3=(a4(k,nft)-f2)**2

                                IF (f3 < dminim(k)) THEN
                                    IF (d8(k) <= d7) THEN
                                        nel1=nel1+1
                                        nob1(nel1)=k
                                    END IF
                                END IF
                            END DO
                            IF (nel1 > 0) THEN
                                ! ns=1
                                nupdate=0
                                ! Solving auxiliaty problem
                                ns=3
                                CALL optim2(z,xc((i-1)*nft+1:),f)
                                CALL auxfunc(xc((i-1)*nft+1:),f)
                                neltotal=neltotal+nel1*nupdate+nrecord
                                ncount=ncount+1
                            ELSE ! nel1 == 0
                                f=large
                            END IF
                            fv2(i) = f
                            IF (f < fmin) fmin=f
                        END DO


                        fmin1=gamma2*fmin
                        jpoints=0
                        DO i=1,ngood0
                            fv1=fv2(i)
                            IF (fv1 <=fmin1) THEN
                                jpoints=jpoints+1
                                DO j=1,nft
                                    xc(j+(jpoints-1)*nft)=xc(j+(i-1)*nft)
                                END DO
                            END IF
                        END DO
                        CALL cleaning(jpoints,xc)
                        ngood0=jpoints

                        fbest=large
                        DO i=1,ngood0
                            ns=1
                            DO j=1,nft
                                z(j)=xc(j+(i-1)*nft)
                            END DO
                            CALL optim2(z,xc((i-1)*nft+1:),f)

                            ncount=ncount+1
                            fv(i)=f
                            IF (f < fbest) fbest=f
                        END DO
                        fbest1=gamma3*fbest
                        jpoints=0
                        DO i=1,ngood0
                            fv1=fv(i)
                            IF (fv1 <= fbest1) THEN
                                jpoints=jpoints+1
                                DO j=1,nft
                                    xc(j+(jpoints-1)*nft)=xc(j+(i-1)*nft)
                                END DO
                            END IF
                        END DO
                        CALL cleaning(jpoints,xc)
                        ngood0=jpoints


                        fb=large
                        DO i=1,ngood0
                            ns=2
                            DO j=1,nft
                                x(j+(nc-1)*nft)=xc(j+(i-1)*nft)
                            END DO
                            CALL optim2(x,x1,f)
                            ncount=ncount+nc
                            IF (f < fb) THEN
                                fb=f
                                DO k=1,nc
                                    DO j=1,nft
                                        xbest(j+(k-1)*nft)=x1(j+(k-1)*nft)
                                    END DO
                                END DO
                            END IF
                        END DO

                        DO k=1,nc
                            DO j=1,nft
                                x(j+(k-1)*nft)=xbest(j+(k-1)*nft)
                            END DO
                        END DO

                    ELSE ! nc == 1
                        DO j=1,nft
                            z(j)=one
                        END DO
                        nel1=nrecord
                        DO k=1,nrecord
                            nob1(k)=k
                        END DO
                        nupdate=0
                        ns=3
                        CALL optim2(z,x,f)
                        CALL clrfunc(x,f)

                        ncount=ncount+1
                        neltotal = nel1*nupdate + nrecord
                    END IF

                    CALL distribution(x,f)
                    toler=MAX(tol2,tol2*f/REAL(nrecord,prec))
                    IF (nrecord > 50000) toler = MAX(tol2,100.0_prec*tol2*f/REAL(nrecord,prec))
                    IF (iscale == 1) THEN
                        call rescaling(x,x4)
                    ELSE
                        DO k=1,nc
                            DO j=1,nft
                                x4(j+(k-1)*nft)=x(j+(k-1)*nft)
                            END DO
                        END DO
                    END IF
                    CALL checkfinal(x4,ff,rmse1,ce1)
                    IF (nc == nclust) THEN
                        IF (ipredi == 0) THEN
                            CALL pred_knn(x4)
                        ELSE IF (ipredi == 1) THEN
                            CALL pred_knnr(x4)
                        ELSE
                            CALL pred_knn2(x4)
                        END IF
                    END IF
                    naux2tot=naux2tot+naux2
                    nauxtot=nauxtot+naux
                    nfittot=nfittot+nfit
                    WRITE(40,*)
                    WRITE(40,573) nc
573                 FORMAT('No of hyperplanes:',I4)
                    WRITE(40,*)
                    WRITE(41,*)
                    WRITE(41,528) noutcom,nc
528                 FORMAT(' Regression coefficiences for feature ', i5,' with ',i3,' regression functions are')

                    WRITE(41,*)
                    DO k=1,nc
                        WRITE(41,722) (x4(j+(k-1)*nft),j=1,nft)
                    END DO

722                 FORMAT(20f12.5)
                    WRITE(40,543) ff
543                 FORMAT('Fit function value:',f24.8)
                    WRITE(40,*)
                    WRITE(40,541) ncount
541                 FORMAT('The number of linear regression problems solved:',i10)
                    WRITE(40,542) naux2tot,nauxtot,nfittot
542                 FORMAT('The numbers of evaluations:',i10,i10,i10)
                    CALL getime(time4)
                    timef=time4-time1
                    WRITE(40,*)
                    WRITE(40,141) timef
                    WRITE(40,*)
                    WRITE(40,*)

                    p2=p2+REAL(neltotal,prec)
                    p3=p3+REAL(naux*nrecord,prec)
                    IF(nc == 1) THEN
                        p4=p4+REAL(neltotal,prec)
                    ELSE
                        p4=p4+REAL(nfit*nrecord*nc,prec)
                    END IF
                    IF (timef > tlimit) EXIT innerloop2
                END DO innerloop2
141             FORMAT('CPU time:',f12.3)


                ! Changing the order back
                IF (noutcom < nft) THEN

                    DO i=1,nrecord
                        atmp = a(i,noutcom)
                        a(i,noutcom)=a(i,nft)
                        a(i,nft)=atmp
                    END DO
                END IF

            END DO featureloop
        END DO

        DO i=1,nrecord
            WRITE(44,230) (a(i,k),k=1,nft)
        END DO
230     FORMAT(10f24.8)


        CLOSE(40)
        CLOSE(41)
        CLOSE(43)
        CLOSE(44)

    !END DO
    !END DO
    !END DO
    !END DO
    !END DO
    !CLOSE(78)
    deallocate(x,x1,x4,xbest,z,d8,fv,fv2)
    !deallocate(x,x1,x4,xbest,z,d8,fv,fv2,ftr,rmsetr,cetr)
    !deallocate(a1,x,x1,x4,xbest,z,z1,d8,fv,fv2)
    deallocate(a,a2,a3,a4) ! a täällä uusi
    deallocate(nalk,clalk,mis_ins,mis_ft,nmis_ins,nmis_ft) !uusi
    deallocate(clin,dlin,dminim,nob1,nob,nel)
    deallocate(xlmbm,nmis_ft2)


    STOP


END PROGRAM iviaclr
