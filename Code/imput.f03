!*************************************************************************
!*                                                                       *
!*     Initial imputations and predictions for IviaCLR                   *
!*     (last modified 14.02.2019)                                        *
!*                                                                       *
!*************************************************************************
!*
!*     Modules included:
!*
!*     firstimput       ! Initial imputation for IviaCLR.
!*     prediction       ! Subprograms for prediction for IviaCLR.
!*

MODULE firstimput

    USE r_precision, ONLY : prec  ! Precision for reals.
    USE initimp, ONLY : &
        nft, &                    ! Number of features in data.
        nrecord, &                ! Number of data points in data.
        misval, &                 ! Value of missing value.
        a                         ! Data set, a(nrecord,nft).
                                  !   Input: incomplete data set.
                                  !   Output: imputed data set.



    IMPLICIT NONE

    ! MODULE firstimput includes the following subroutines (S).
    PUBLIC :: &
        read_data, &              ! S Reading of incomplete data.
        mean_imput, &             ! S Initial imputation of mean to all missing values.
        reg_imput, &              ! S Initial imputation by regressing complete data.
        rreg_imput                ! S Initial imputation by recursively regressing data.
    PRIVATE :: &
        deldata                   ! S Form complete data set by deletion.

CONTAINS

    SUBROUTINE read_data()        ! Reading input file.

        USE initimp, ONLY : &
            infile, &             ! Input data file.
            a2                    ! Data file with missing values


        IMPLICIT NONE

        INTEGER :: &
            i,j

        OPEN(78,file=infile,status='old',form='formatted')

        DO i=1,nrecord
            READ(78,*) (a(i,j),j=1,nft)
        END DO

        CLOSE(78)

        a2 = a

    END SUBROUTINE read_data



    SUBROUTINE mean_imput()       ! Initial imputation of mean to all missing values.

        USE param, ONLY : zero    ! Parameters.
        USE initimp, ONLY : &
            nmis, &               ! number of missing values in data set.
            nmis_ins, &           ! number of missing features in each data point, nmis_ins(nrecord)
                                  ! nmis_ins(i) = 0 --- complete object
            nmis_ft, &            ! number of missing values in each feature, nmis_ft(nft)
                                  ! nmis_ft(i) = 0 --- the feature has no missing values
                                  !   -> no need for regression with this feature as an output.
            mis_ins, &            ! missing instances, mis_ins((nrecord-1)*nft)
            mis_ft                ! missing features, mis_ft((nrecord-1)*nft)


        IMPLICIT NONE

        REAL(KIND=prec) :: &
            newvalue
        INTEGER :: &
            i,j,i1,i2



        i1       = 0
        i2       = 0
        nmis     = 0
        nmis_ins = 0
        nmis_ft  = 0

        DO j = 1,nft
            newvalue = zero

            DO i = 1,nrecord
                IF (a(i,j) == misval) THEN

                    i1 = i1 + 1
                    mis_ins(i1) = i             ! missing instances
                    mis_ft(i1)  = j             ! missing features
                    nmis_ins(i) = nmis_ins(i)+1 ! number of missing features in each data point
                    nmis_ft(j)  = nmis_ft(j)+1  ! number of missing values in each feature
                    nmis        = nmis+1        ! number of missing values in data set
                ELSE
                    newvalue = newvalue + a(i,j)
                END IF
            END DO
            DO i = 1,nmis_ft(j)
                i2=i2+1
                a(mis_ins(i2),j)= newvalue / (nrecord - nmis_ft(j))
            END DO
        END DO

        WRITE(40,*) 'Number of missing values: ',nmis,'(',nmis*100 / (nrecord*nft),'%)'
        WRITE(41,*) 'Number of missing values: ',nmis,'(',nmis*100 / (nrecord*nft),'%)'
        WRITE(43,*) 'Number of missing values: ',nmis,'(',nmis*100 / (nrecord*nft),'%)'

        WRITE(41,*)
        i2=1
        loop1: DO j=1,nft
            IF (nmis_ft(j) == 0) CYCLE loop1
            WRITE(41,*) 'Mean value for feature ',j,' is ',a(mis_ins(i2),j)
            i2=i2+nmis_ft(j)
        END DO loop1

        PRINT*,' There are ',nmis,' missing values. That is',nmis*100 / (nrecord*nft),'%.'

    END SUBROUTINE mean_imput


    SUBROUTINE reg_imput()          ! Initial imputation by regressing complete data.

        USE param, ONLY : zero,one  ! Parameters.
        USE initimp, ONLY : &
            a2, &                   ! Data file with missing values
            a3, &                   ! "Complete" deleted data set, a3(nrecord,nft).
            a4, &                   ! Auxiliary data matrix, a4(nrecord,nft).
            nmis, &                 ! number of missing values in data set.
            nmis_ins, &             ! number of missing features in each data point, nmis_ins(nrecord)
                                    ! nmis_ins(i) = 0 --- complete object
            nmis_ft, &              ! number of missing values in each feature, nmis_ft(nft)
                                    ! nmis_ft(i) = 0 --- the feature has no missing values
                                    !   -> no need for regression with this feature as an output.
            mis_ins, &              ! missing instances, mis_ins((nrecord-1)*nft)
            mis_ft, &               ! missing features, mis_ft((nrecord-1)*nft)
            misval                  ! Missing value: every value in infile equal to misval
                                    !   is considered to be missing.
        USE initclr, ONLY : &
            ns, &
            nel1, &
            nob1
        USE lmbm_mod, ONLY : optim2 ! LMBM method for optimization.


        IMPLICIT NONE

        REAL(KIND=prec) :: &
            atmp, &
            f, &
            z(nft), &
            x(nft), &
            x4(nft)
        INTEGER :: &
            nrecordtmp, &
            noutcom, &
            i,j,i1,k,k1

        CALL deldata(a3,nmis,nmis_ins,nmis_ft,mis_ins,mis_ft,nrecordtmp)

        WRITE(40,*) 'Number of missing values: ',nmis,'(',nmis*100 / (nrecord*nft),'%)'
        WRITE(41,*) 'Number of missing values: ',nmis,'(',nmis*100 / (nrecord*nft),'%)'
        WRITE(43,*) 'Number of missing values: ',nmis,'(',nmis*100 / (nrecord*nft),'%)'

        IF (nrecordtmp == 0) THEN
            PRINT*,' All data points have missing values: no complete data with deletion. '
            PRINT*,' Try mean as initial imputation.'
            STOP
        END IF



        ! Regression needs to be done to all features with missing values
        ! but we can skip the feature if nmis_ft(j)=0.
        k1       = 0
        loop1: DO noutcom = 1,nft
            IF(nmis_ft(noutcom) == 0) CYCLE loop1

            DO k=1,nft
                DO i=1,nrecordtmp
                    a4(i,k)=a3(i,k)
                END DO
            END DO

            IF (noutcom < nft) THEN

                DO i=1,nrecordtmp
                    atmp = a4(i,noutcom)
                    a4(i,noutcom) = a4(i,nft)
                    a4(i,nft) = atmp
                END DO
            END IF

            z=one

            nel1=nrecordtmp
            DO k=1,nrecordtmp
                nob1(k)=k
            END DO
            ns=3
            CALL optim2(z,x,f)
            DO j=1,nft
                x4(j)=x(j)
            END DO
            WRITE(40,*)
            WRITE(41,*)
            WRITE(41,*) ' Initial regression coefficiences for feature ', noutcom
            WRITE(41,722) (x4(j),j=1,nft)
722         FORMAT(20f12.5)

            DO k=1,nmis_ft(noutcom) !

                k1=k1+1
                a(mis_ins(k1),noutcom)=x4(nft)

                i1=0
                loop2: DO i=1,nft
                    IF (i == noutcom) CYCLE loop2
                    IF (a2(mis_ins(k1),i)==misval) CYCLE loop2
                    i1=i1+1
                    a(mis_ins(k1),noutcom)=a(mis_ins(k1),noutcom) + x4(i1)*a(mis_ins(k1),i)
                END DO loop2
            END DO

        END DO loop1

        nrecordtmp = nrecord

    END SUBROUTINE reg_imput


    SUBROUTINE rreg_imput()    ! Initial imputation by recursive regression.

        ! Algortihm:
        ! 1. Regress the variable with the fewest number of missing values
        !     using complete data set (i.e. data set produced with deletion).
        ! 2. Impute values
        ! 3. Repeat the process for the variable with the next fewest missing
        !    values using updated data set with previously imputed values.

        USE param, ONLY : zero,one  ! Parameters.
        USE initimp, ONLY : &
            a2, &               ! Data file with missing values.
            a3, &               ! "Complete" deleted data set, a3(nrecord,nft).
            a4, &               ! Auxiliary data matrix, a4(nrecord,nft).
            nmis, &             ! number of missing values in data set.
            nmis_ins, &         ! number of missing features in each data point, nmis_ins(nrecord)
                                ! nmis_ins(i) = 0 --- complete object
            nmis_ft, &          ! number of missing values in each feature, nmis_ft(nft)
                                ! nmis_ft(i) = 0 --- the feature has no missing values
                                !   -> no need for regression with this feature as an output.
            mis_ins, &          ! missing instances, mis_ins((nrecord-1)*nft)
            mis_ft, &           ! missing features, mis_ft((nrecord-1)*nft)
            misval              ! Missing value: every value in infile equal to misval
                                !   is considered to be missing.

        USE initclr, ONLY : &
            ns, &
            nel1, &
            nob1
        USE lmbm_mod, ONLY : optim2         ! LMBM method for optimization.


        IMPLICIT NONE

        REAL(KIND=prec) :: &
            atmp, &
            f, &
            z(nft), &
            x(nft), &
            x4(nft)
        INTEGER :: &
            nmis2, &
            nmis_ins2(nrecord), &
            nmis_ft2(nft), &
            nmis_ft3(nft), &
            mis_ins2((nrecord-1)*nft), &
            mis_ft2((nrecord-1)*nft), &
            nrecordtmp, &
            ifea, &
            noutcom, &
            i,j,i1,k,k1


        CALL deldata(a3,nmis,nmis_ins,nmis_ft,mis_ins,mis_ft,nrecordtmp)

        nmis2=nmis
        nmis_ins2 = nmis_ins
        nmis_ft2 = nmis_ft
        mis_ins2=mis_ins
        mis_ft2=mis_ft
        WRITE(40,*) 'Number of missing values: ',nmis2,'(',nmis2*100 / (nrecord*nft),'%)'
        WRITE(41,*) 'Number of missing values: ',nmis2,'(',nmis2*100 / (nrecord*nft),'%)'
        WRITE(43,*) 'Number of missing values: ',nmis2,'(',nmis2*100 / (nrecord*nft),'%)'

        IF (nrecordtmp == 0) THEN
            PRINT*,' All data points have missing values: no complete data with deletion. '
            PRINT*,' Try mean as initial imputation.'
            STOP
        END IF

        nmis_ft3=nmis_ft

        loop1: DO ifea = 1,nft ! Go through all features

            noutcom = MINLOC(nmis_ft3,dim=1)

            IF(nmis_ft3(noutcom) == 0) THEN
                nmis_ft3(noutcom) = nrecord+1
                CYCLE loop1
            END IF

            DO k=1,nft
                DO i=1,nrecordtmp
                    a4(i,k)=a3(i,k)
                END DO
            END DO

            IF (noutcom < nft) THEN

                DO i=1,nrecordtmp
                    atmp = a4(i,noutcom)
                    a4(i,noutcom) = a4(i,nft)
                    a4(i,nft) = atmp

                END DO
            END IF

            z=one

            nel1=nrecordtmp
            DO k=1,nrecordtmp
                nob1(k)=k
            END DO
            ns=3
            CALL optim2(z,x,f)
            DO j=1,nft
                x4(j)=x(j)
            END DO

            WRITE(40,*)
            WRITE(41,*)
            WRITE(41,*) ' Initial regression coefficiences for feature ', noutcom
            WRITE(41,722) (x4(j),j=1,nft)
722         FORMAT(20f12.5)

            k1=0
            DO i1=1,noutcom-1
                k1=k1+nmis_ft2(i1)
            END DO
            DO k=1,nmis_ft3(noutcom)

                k1=k1+1
                a(mis_ins2(k1),noutcom)=x4(nft)

                i1=0

                loop2: DO i=1,nft
                    IF (i == noutcom) CYCLE loop2
                    IF (a2(mis_ins2(k1),i)==misval) CYCLE loop2
                    i1=i1+1
                    a(mis_ins2(k1),noutcom)=a(mis_ins2(k1),noutcom) + x4(i1)*a(mis_ins2(k1),i)

                END DO loop2

            END DO

            nmis_ft3(noutcom) = nrecord+1

            ! Update data set with deletion
            ! The easiest (but not efficient) way to do this is to do new deleted data set from a.
            CALL deldata(a3,nmis2,nmis_ins2,nmis_ft2,mis_ins2,mis_ft2,nrecordtmp)

        END DO loop1

        nrecordtmp = nrecord

    END SUBROUTINE rreg_imput


    SUBROUTINE deldata(a3,nmis,nmis_ins,nmis_ft,mis_ins,mis_ft,nrecordtmp)    ! Form complete data set by deletion.

        USE param, ONLY : zero,one  ! Parameters.

        IMPLICIT NONE

        ! Array Arguments
        REAL(KIND=prec), DIMENSION(:,:), INTENT(OUT) :: &
            a3                  ! Complete data set.
        INTEGER, DIMENSION(:), INTENT(OUT) :: &
            nmis_ins, &         ! number of missing features in each data point, nmis_ins(nrecord)
                                ! nmis_ins(i) = 0 --- complete object
            nmis_ft, &          ! number of missing values in each feature, nmis_ft(nft)
                                ! nmis_ft(i) = 0 --- the feature has no missing values
                                !   -> no need for regression with this feature as an output.
            mis_ins, &          ! missing instances, mis_ins((nrecord-1)*nft)
            mis_ft              ! missing features, mis_ft((nrecord-1)*nft)

        ! Scalar Arguments
        INTEGER, INTENT(OUT) :: &
            nrecordtmp, &       ! Number of complete objects.
            nmis                ! Number of missing values in data set.

        INTEGER :: &
            i,j,i1

        i1       = 0
        nmis     = 0
        nmis_ins = 0
        nmis_ft  = 0

        DO j = 1,nft

            DO i = 1,nrecord
                IF (a(i,j) == misval) THEN

                    i1 = i1 + 1
                    mis_ins(i1) = i   ! missing instances
                    mis_ft(i1)  = j   ! missing features
                    nmis_ins(i) = nmis_ins(i)+1 ! number of missing features in each data point
                    nmis_ft(j)  = nmis_ft(j)+1  ! number of missing values in each feature
                    nmis       = nmis+1 ! number of missing values in data set
                END IF
            END DO
        END DO

        ! complete data set with deletion
        DO j = 1,nft
            i1 = 0
            DO i = 1,nrecord
                IF (nmis_ins(i) == 0) THEN ! no missing values
                    i1 = i1+1
                    a3(i1,j) = a(i,j)
                END IF
            END DO
        END DO

        PRINT*,' There are ',nmis,' missing values. That is',nmis*100 / (nrecord*nft),'%.'

        nrecordtmp = i1
        !IF (nrecordtmp == 0) THEN
        !    PRINT*,' All data points have missing values: no complete data with deletion. '
        !    PRINT*,' The resulting imputed data may be strongly biased. Try mean as initial imputation. '
        !END IF

    END SUBROUTINE deldata

END MODULE firstimput


MODULE prediction

    USE r_precision, ONLY : prec   ! Precision for reals.
    USE initimp, ONLY : &
        nft, &                     ! Number of features in data.
        nrecord, &                 ! Number of data points in data.
        noutcom, &                 ! Index of original output column.
        misval, &                  ! Missing value:  every value equal to misval
                                   !   is considered to be missing.
        a                          ! Data set, a(nrecord,nft).
                                   !   Input: incomplete data set.
                                   !   Output: imputed data set.
                                   ! Note, this data set has been rearranged such that
                                   ! the last column is an output column.


    IMPLICIT NONE

    ! MODULE prediction includes the following subroutines (S).
    PUBLIC :: &
        pred_weight, &             ! S Form complete data set by imputing missing values
                                   !   with weighs based on popularity of clusters.
        pred_knn, &                ! S Form complete data set by imputing missing values
                                   !   using k-NN with popularity based weights on clusters.
        pred_knn2, &               ! S Form complete data set by imputing missing values
                                   !   using k-NN without weights on clusters.
        pred_knnr                  ! S Form complete data set by imputing missing values
                                   !   using k-NN with rmsq based weights on clusters.

        !PRIVATE :: &

CONTAINS


    SUBROUTINE pred_weight(x)      ! Impute data set using weighs based on popularity of clusters.

        USE param, ONLY : zero,one ! Parameters.
        USE initimp, ONLY : &
            mis_ins, &             ! missing instances.
            mis_ft, &              ! missing features.
            nmis, &                ! number of missing values in data set.
            nalk                   ! nalk(i) = number of data points in cluster i.
        USE initclr, ONLY : &
            nc                     ! Current number of clusters, loops from 1 to nclust.

        IMPLICIT NONE

        REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
            x                      ! regression coeffiences.

        REAL(KIND=prec) :: &
            w, &                   ! weight.
            z(nc)                  ! clusterwise predictions.

        INTEGER :: &
            i,j,k


        ! Computing imputation for missing value
        DO i=1,nmis ! Go trough every missing value. Do nothing if there is no missing value at noutcom.
            IF (mis_ft(i)==noutcom) THEN ! misins(i) is the index of the data point that has missing value at noutcom.
                a(mis_ins(i),nft) = zero ! Missing values are at column nft
                DO k=1,nc !
                    w = nalk(k)/REAL(nrecord,prec) !
                    z(k) = x(k*nft) !
                    DO j=1,nft-1
                        z(k)=z(k)+a(mis_ins(i),j)*x(j+(k-1)*nft)
                    END DO
                    a(mis_ins(i),nft) = a(mis_ins(i),nft) + w*z(k)
                END DO

                WRITE(43,*)
                WRITE(43,*) ' Missing value indices ', mis_ins(i),noutcom
                WRITE(43,*) ' New value ',a(mis_ins(i),nft)

            END IF
        END DO

    END SUBROUTINE pred_weight


    SUBROUTINE pred_knn(x)         ! Imputation using random sample k-nn
                                   ! with popularity based weights on clusters.

        USE param, ONLY : zero,one ! Parameters.
        USE initimp, ONLY : &
            a2, &                  ! data set with missing values.
            mis_ins, &             ! missing instances.
            mis_ft, &              ! missing features.
            nmis_ft, &             ! number of missing values in each feature.
            nmis, &                ! number of missing values in data set.
            clalk, &               ! clalk(i) = cluster where data point i belongs to.
            knn, &                 ! number of nearest neighbours.
            nsample                ! size of the subset taken to compute k-NN. nsample <= nrecord.
        USE initclr, ONLY : &
            nc                     ! Current number of clusters, loops from 1 to nclust.


        IMPLICIT NONE

        REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
            x


        INTEGER :: &
            j_rand, &              ! random index.
            j_rand_used(nsample+1), &   ! used random indices.
            i,j,k,itmp, &
            isamples, &            ! current number of samples.
            nn_index(knn), &       ! indices of data points close to target.
            nalk(nc)               ! nalk(i) = number of data points (close to target) in cluster i.

        REAL(KIND=prec) :: &
            u, &                   ! random number.
            w, &                   ! auxiliary weight.
            z(nc), &               ! predictions for different clusters.
            rmsq_tmp, &            ! rmsq of data points and target.
            rmsq(knn)              ! knn smallest rmsq of data points and target.



        ! Compute knn nearest neighbours from nsample random data points of data set
        ! exclude data points with missing value in noutcom.

        loop1: DO i=1,nmis ! Loop trough all missing values.
                           ! Do nothing if there is no missing value at noutcom.
            j = 0
            isamples = 0
            j_rand   = 0
            j_rand_used = 0
            IF (mis_ft(i)==noutcom) THEN ! data point a(mis_ins(i),noutcom) is missing.

                loop2: DO

                    IF (nsample >= nrecord - nmis_ft(noutcom)) THEN ! Go trough all indices

                        j_rand = j_rand+1
                        IF (j_rand > nrecord) EXIT loop2 ! No more samples to search.

                    ELSE ! Generate a discrete uniform distribution on the integers {1, 2, ... , nrecord}

                        call RANDOM_NUMBER(u)
                        j_rand = 1 + FLOOR(nrecord*u)  ! Choose one from nrecord integers
                        IF (ANY(j_rand_used == j_rand)) CYCLE loop2 ! Do not use the same index multiple times.

                    END IF

                    IF (j_rand==mis_ins(i)) CYCLE loop2 ! Do not compare data point to itself or to
                    IF (a2(j_rand,noutcom) == misval) CYCLE loop2  ! data points that have a missing value in noutcom.

                    isamples = isamples + 1
                    IF (isamples > nsample) EXIT loop2 ! Enough random samples.
                    IF (j_rand > nrecord+1) EXIT loop2 ! No more samples to search. Miksi täällä on +1???
                    j_rand_used(isamples)=j_rand


                    rmsq_tmp = zero
                    DO k = 1,nft-1
                        rmsq_tmp = rmsq_tmp + (a(mis_ins(i),k)-a(j_rand,k))**2
                    END DO
                    !   rmse_tmp=SQRT(f/REAL(nrecord,prec)) ! We do not need to get real result.

                    IF (isamples > knn) THEN
                        IF (rmsq_tmp < MAXVAL(rmsq)) THEN
                            itmp=MAXLOC(rmsq,DIM=1)
                            rmsq(itmp)=rmsq_tmp
                            nn_index(itmp)=j_rand
                        END IF

                    ELSE
                        rmsq(isamples) = rmsq_tmp
                        nn_index(isamples) = j_rand
                    END IF


                END DO loop2


                nalk=0
                DO j=1,knn

                    DO k=1,nc
                        IF (clalk(nn_index(j))==k) nalk(k)=nalk(k)+1
                    END DO
                END DO

                a(mis_ins(i),nft) = zero ! The missing values are at nft column
                DO k=1,nc !
                    w = nalk(k)/REAL(knn,prec)
                    z(k) = x(k*nft)
                    DO j=1,nft-1
                        z(k)=z(k)+a(mis_ins(i),j)*x(j+(k-1)*nft)
                    END DO
                    a(mis_ins(i),nft) = a(mis_ins(i),nft) + w*z(k)
                END DO
                WRITE(43,*)
                WRITE(43,*) ' Missing value indices ', mis_ins(i),noutcom
                WRITE(43,*) ' New value ',a(mis_ins(i),nft)


            END IF

        END DO loop1

    END SUBROUTINE pred_knn



    SUBROUTINE pred_knn2(x)     ! Form complete data set using random sample k-nn
                                ! without weighting.

        USE param, ONLY : zero,one  ! Parameters.
        USE initimp, ONLY : &
            a2, &              ! data set with missing values.
            misval, &          ! Missing value:  every value equal to misval
                               !   is considered to be missing.
            mis_ins, &         ! missing instances.
            mis_ft, &          ! missing features.
            nmis_ft, &         ! number of missing values in each feature.
            nmis, &            ! number of missing values in data set.
            clalk, &           ! clalk(i) = cluster where data point i belongs to.
            knn, &             ! number of nearest neighbours.
            nsample            ! size of the subset taken to compute k-NN. nsample <= nrecord.
        USE initclr, ONLY : &
            nc                 ! Current number of clusters, loops from 1 to nclust.


        IMPLICIT NONE

        REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
            x


        INTEGER :: &
            i,j,k, &
            j_rand, &          ! random index.
            j_rand_used(nsample+1), &   ! used random indices.
            itmp, &
            isamples, &        ! current number of samples.
            nn_index(knn), &   ! indices of data points close to target.
            k_array(1), &
            nalk(nc)           ! nalk(i) = number of data points (close to target) in cluster i.

        REAL(KIND=prec) :: &
            u, &               ! random number.
            z(nc), &           ! predictions for different clusters.
            rmsq_tmp, &        ! rmsq of data points and target.
            rmsq(knn)          ! knn smallest rmsq of data points and target.



        ! Compute knn nearest neighbours from nsample random data points of data set
        ! exclude data points with missing value in noutcom.

        loop1: DO i=1,nmis ! loop trough all missing values
            j = 0
            isamples = 0
            j_rand   = 0
            j_rand_used = 0
            IF (mis_ft(i)==noutcom) THEN ! data point a(mis_ins(i),noutcom) is missing.

                loop2: DO

                    IF (nsample >= nrecord - nmis_ft(noutcom)) THEN ! Go trough all indices

                        j_rand = j_rand+1
                        IF (j_rand > nrecord) EXIT loop2 ! No more samples to search.

                    ELSE ! Generate a discrete uniform distribution on the integers {1, 2, ... , nrecord}

                        call RANDOM_NUMBER(u)
                        j_rand = 1 + FLOOR(nrecord*u)  ! Choose one from nrecord integers
                        IF (ANY(j_rand_used == j_rand)) CYCLE loop2 ! Do not use the same index multiple times.

                    END IF

                    IF (j_rand==mis_ins(i)) CYCLE loop2 ! Do not compare data point to itself or to
                    IF (a2(j_rand,noutcom) == misval) CYCLE loop2  ! data points that have a missing value in noutcom.

                    isamples = isamples + 1
                    IF (isamples > nsample) EXIT loop2 ! Enough random samples.
                    IF (j_rand > nrecord+1) EXIT loop2 ! No more samples to search.
                    j_rand_used(isamples)=j_rand


                    rmsq_tmp = zero
                    DO k = 1,nft-1
                        rmsq_tmp = rmsq_tmp + (a(mis_ins(i),k)-a(j_rand,k))**2
                    END DO
                    !   rmse_tmp=SQRT(f/REAL(nrecord,prec)) ! We do not need to get real result.

                    IF (isamples > knn) THEN
                        IF (rmsq_tmp < MAXVAL(rmsq)) THEN
                            itmp=MAXLOC(rmsq,DIM=1)
                            rmsq(itmp)=rmsq_tmp
                            nn_index(itmp)=j_rand
                        END IF

                    ELSE
                        rmsq(isamples) = rmsq_tmp
                        nn_index(isamples) = j_rand
                    END IF


                END DO loop2


                nalk=0
                DO j=1,knn

                    DO k=1,nc
                        IF (clalk(nn_index(j))==k) nalk(k)=nalk(k)+1
                    END DO
                END DO
                k_array = MAXLOC(nalk)
                k=k_array(1)

                z(k) = x(k*nft)
                DO j=1,nft-1
                   z(k)=z(k)+a(mis_ins(i),j)*x(j+(k-1)*nft)
                END DO
                a(mis_ins(i),nft) = z(k)

                WRITE(43,*)
                WRITE(43,*) ' Missing value indices ', mis_ins(i),noutcom
                WRITE(43,*) ' New value ',a(mis_ins(i),nft)


            END IF

        END DO loop1

    END SUBROUTINE pred_knn2


    SUBROUTINE pred_knnr(x)    ! Form complete data set using random sample k-nn
                               ! with rmsq based weights on clusters.

        USE param, ONLY : zero,one,small  ! Parameters.
        USE initimp, ONLY : &
            a2, &              ! Data set with missing values.
            misval, &          ! Missing value:  every value equal to misval
                               !   is considered to be missing.
            mis_ins, &         ! Missing instances.
            mis_ft, &          ! Missing features.
            nmis_ft, &         ! Number of missing values in each feature.
            nmis, &            ! Number of missing values in data set.
            clalk, &           ! clalk(i) = cluster where data point i belongs to.
            knn, &             ! Number of nearest neighbours.
            nsample            ! Size of the subset taken to compute k-NN. nsample <= nrecord.
        USE initclr, ONLY : &
            nc                 ! Current number of clusters, loops from 1 to nclust.


        IMPLICIT NONE

        REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: &
            x


        INTEGER :: &
            i,j,k, &
            j_rand, &          ! Random index.
            j_rand_used(nsample+1), &   ! Used random indices.
            itmp, &
            isamples, &        ! Current number of samples.
            nn_index(knn), &   ! Indises of data points close to target.
            nalk(nc)           ! nalk(i) = number of data points in cluster i.



        REAL(KIND=prec) :: &
            u, &               ! Random number.
            w, &               ! Auxiliary weight.
            z(nc), &           ! Predictions for different clusters.
            ralk(nc), &        ! ralk(i) = rmsq based weight for cluster i.
            rmsq_tmp, &        ! RMSE of data points and target.
            rmsq(knn)          ! knn smallest RMSE of data points and target.



        ! Compute knn nearest neighbours from nsample random data points of data set
        ! exclude data points with missing value in noutcom.


        loop1: DO i=1,nmis ! loop trough all missing values in the data
            isamples = 0
            j_rand   = 0
            j_rand_used = 0
            IF (mis_ft(i)==noutcom) THEN ! data point a(mis_ins(i),noutcom) is missing.

                loop2: DO
                    IF (nsample >= nrecord - nmis_ft(noutcom)) THEN ! Go trough all indices
                        j_rand = j_rand+1
                        IF (j_rand > nrecord) EXIT loop2 ! No more samples to search.

                    ELSE ! Generate a discrete uniform distribution on the integers {1, 2, ... , nrecord}
                        call RANDOM_NUMBER(u)
                        j_rand = 1 + FLOOR(nrecord*u)  ! Choose one from nrecord integers
                        IF (ANY(j_rand_used == j_rand)) CYCLE loop2 ! Do not use the same index multiple times.

                    END IF


                    IF (j_rand==mis_ins(i)) CYCLE loop2 ! Do not compare data point to itself or to
                    IF (a2(j_rand,noutcom) == misval) CYCLE loop2  ! data points that have a missing value in noutcom.


                    isamples = isamples + 1

                    IF (isamples > nsample) EXIT loop2 ! Enough random samples.
                    IF (j_rand > nrecord+1) EXIT loop2 ! No more samples to search.
                    j_rand_used(isamples)=j_rand

                    rmsq_tmp = zero
                    DO k = 1,nft-1
                        rmsq_tmp = rmsq_tmp + (a(mis_ins(i),k)-a(j_rand,k))**2
                    END DO
                    !   rmse_tmp=SQRT(f/REAL(nrecord,prec)) ! We do not need to get real result.

                    IF (isamples > knn) THEN
                        IF (rmsq_tmp < MAXVAL(rmsq)) THEN
                            itmp=MAXLOC(rmsq,DIM=1)
                            rmsq(itmp)=rmsq_tmp
                            nn_index(itmp)=j_rand
                        END IF

                    ELSE
                        rmsq(isamples) = rmsq_tmp
                        nn_index(isamples) = j_rand
                    END IF


                END DO loop2


                nalk=0
                rmsq_tmp = SUM(rmsq)
                ralk=zero

                DO j=1,knn
                    IF (rmsq_tmp > small .AND. knn /= 1) THEN
                        DO k=1,nc
                            IF (clalk(nn_index(j))==k) THEN
                                ralk(k)=ralk(k)+(rmsq_tmp-rmsq(j))/((knn-1)*rmsq_tmp)
                                nalk(k)=nalk(k)+1
                            END IF
                        END DO

                    ELSE
                        DO k=1,nc
                            IF (clalk(nn_index(j))==k) THEN
                                ralk(k) = ralk(k)+one/REAL(knn,prec)
                                nalk(k)=nalk(k)+1
                            END IF
                        END DO
                    END IF
                END DO


                a(mis_ins(i),nft) = zero ! The missing values are at nft column
                DO k=1,nc !
                    w = ralk(k)
                    z(k) = x(k*nft)
                    DO j=1,nft-1
                        z(k)=z(k)+a(mis_ins(i),j)*x(j+(k-1)*nft)
                    END DO
                    a(mis_ins(i),nft) = a(mis_ins(i),nft) + w*z(k)
                END DO
                WRITE(43,*)
                WRITE(43,*) ' Missing value indices ', mis_ins(i),noutcom
                WRITE(43,*) ' New value ',a(mis_ins(i),nft)


            END IF

        END DO loop1

    END SUBROUTINE pred_knnr


END MODULE prediction




