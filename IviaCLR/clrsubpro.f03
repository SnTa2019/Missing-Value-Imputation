!*************************************************************************
!*                                                                      *
!*     Subroutines for nonsmooth optimization based clusterwise linear   *
!*     regression software (last modified 18.06.2016).                   *
!*                                                                       *
!*************************************************************************
!*
!*     Modules included:
!*
!*     clrsubpro       !
!*

MODULE clrsubpro       ! Subroutines for clusterwise linear regression software

    USE r_precision, ONLY : prec      ! Precision for reals.
    IMPLICIT NONE

    ! MODULE clusteringmod includes the following subroutines.

    PUBLIC :: &
        distribution, &    ! Distribution over clusters.
        checkfinal, &      ! Final distribution over clusters.
        step1, &           ! Initialization of clusters.
        cleaning, &        ! Cleaning arrays.
        scaling, &
        rescaling
!    PRIVATE :: &

CONTAINS

    !==================================================================
    !     distribution over clusters
    !===================================================================

    SUBROUTINE distribution(x,f)

        USE param, ONLY : zero, large
        USE initimp, ONLY : &
            nft, &         ! Number of used features.
            nrecord        ! Number of records in data.
        USE initclr, ONLY : &
            maxdim, &
            a4, &                           ! Auxiliary data matrix a4(nrecord,nft).
            dminim, &                       !
            nc, &                           ! Current number of clusters, loops from 1 to nclust.
            nel, &
            nob

        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec) :: &
            f, &
            x(maxdim)

        ! Local variables
        REAL(KIND=prec) :: &
            d1,d2,d3
        INTEGER ::  &
            i,i1,j,k

        i1 = 1
        f = zero
        DO i=1,nc
            nel(i)=0
        END DO

        DO k=1,nrecord
            d1 = large
            DO i=1,nc
                d2=x(i*nft)
                DO j=1,nft-1
                    d2=d2+a4(k,j)*x(j+(i-1)*nft)
                END DO
                d3=(d2-a4(k,nft))**2
                IF (d3 < d1) THEN
                    d1=d3
                    i1=i
                END IF
            END DO
            f=f+d1
            nel(i1)=nel(i1)+1
            nob(i1,nel(i1))=k
            dminim(k)=d1
            END DO

        RETURN

    END SUBROUTINE distribution


    !===================================================================
    !  final distribution over clusters
    !===================================================================
    SUBROUTINE checkfinal(x,f,rmse1,ce1)

        USE param, ONLY : zero, one, large
        USE initimp, ONLY : &
            a, &                 ! Data matrix a(nrecord,nft), from input file.
            nft, &               ! Number of features in data, from user.
            nrecord, &           ! Number of instances in data, from user.
            nalk, &              ! nalk(i) = number of data points in cluster i.
            clalk                ! clalk(i) = cluster where data point i belongs to.
        USE initclr, ONLY : &
            maxdim, &            ! Maximum number of variables for optimization.
            nc                   ! Current number of clusters, loops from 1 to nclust.

        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec) :: &
            f,             &     !
            x(maxdim),     &     ! Regression coefficiences.
            rmse1,         &     ! Root-mean-square-error of data points versus regression line.
            ce1

        ! Local variables
        REAL(KIND=prec) :: &
            r1,r2, &
            fe1,fe2, &
            d2,d3,d4
        INTEGER ::  &
            kalk, &             ! kalk index of cluster (temporary variable)
            i,j,k

        ! Intrinsic Functions
        INTRINSIC SQRT,REAL

        kalk = 0
        nalk = 0
        rmse1 = zero
        ce1 = zero

        r1 = zero
        loop1: DO i=1,nrecord
            r1=r1+a(i,nft)
        END DO loop1

        r2= zero

        r1=r1/REAL(nrecord,prec)

        fe1 = zero
        loop2: DO i=1,nrecord
            fe1=fe1+(a(i,nft)-r1)**2
        END DO loop2

        fe2 = zero

        f = zero
        loop3: DO i=1,nrecord
            d4 = large

            DO k=1,nc
                d2=x(k*nft)
                DO j=1,nft-1
                    d2=d2+a(i,j)*x(j+(k-1)*nft)
                END DO
                d3 = (d2-a(i,nft))**2 ! error to last element (i.e. missing value), need some checking
                IF (d3 < d4) THEN
                   d4=d3
                   kalk=k
                   clalk(i)=k         !
                END IF
            END DO
            nalk(kalk)=nalk(kalk)+1
            f=f+d4
        END DO loop3

        rmse1=SQRT(f/REAL(nrecord,prec))
        ce1=one-f/fe1

        RETURN

    END SUBROUTINE checkfinal

    !===================================================================
    !  Initialize clusters
    !===================================================================

    SUBROUTINE step1(x,z,jpoints)

        USE param, ONLY : zero
        USE initimp, ONLY : &
            nft, &                          ! Number of used features.
            nrecord                         ! Number of records in data.
        USE initclr, ONLY : &
            maxdim, &                       ! Maximum number of variables for optimization.
            maxinit, &                      ! Maximum number of candidate points of data set.
            a4, &                           ! Auxiliary data matrix a4(nrecord,nft).
            gamma1, &                       ! Parameter affecting the size of the set of starting point,
            dminim, &                       !
            nclust, &                       ! Maximum number of clusters, from user.
            nc, &                           ! Current number of clusters, loops from 1 to nclust.
            nel, &
            nob

        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec) ::      &
            x(maxdim), &
            z(maxinit)
        INTEGER :: &
            jpoints

        ! Local variables
        REAL(KIND=prec) ::      &
            y, &
            d0,d2,d3,d4,d6,d7, &
            d1(nclust,nrecord), &
            d5(nrecord)
        INTEGER ::  &
            l1(nrecord), &
            i,i1,j,k,k1


        DO i=1,nc-1
            DO k=1,nel(i)
                k1=nob(i,k)
                l1(k1)=i
            END DO
        END DO

        DO i=1,nc-1
            DO k=1,nrecord
                d0=zero
                DO j=1,nft-1
                    d0=d0+x(j+(i-1)*nft)*a4(k,j)
                END DO
                d1(i,k)=d0
            END DO
        END DO

        d6=zero
        DO i=1,nrecord
            i1=l1(i)
            y=a4(i,nft)-d1(i1,i)
            d5(i)=zero
            DO k=1,nrecord
                d2=d1(i1,k)+y
                d3=(d2-a4(k,nft))**2
                d4=dmin1(zero,d3-dminim(k))
                d5(i)=d5(i)+d4
            END DO
            IF (d5(i) <= d6) THEN
                d6=d5(i)
            END IF
        END DO

        jpoints=0
        d7=gamma1*d6
        DO k=1,nrecord
            i1=l1(k)

            IF (d5(k) <= d7) THEN
                jpoints=jpoints+1
                DO j=1,nft-1
                    z(j+(jpoints-1)*nft)=x(j+(i1-1)*nft)
                END DO
                z(jpoints*nft)=a4(k,nft)-d1(i1,k)

            END IF

        END DO

        CALL cleaning(jpoints,z)

        RETURN

    END SUBROUTINE step1


    !===================================================================
    !  Cleaning arrays
    !===================================================================

    SUBROUTINE cleaning(jpoints,z)

        USE param, ONLY : zero
        USE initclr, ONLY : &
            maxinit, &                   ! Maximum number of candidate points of data set.
            toler, &                     !
            nft                          ! Number of features in data, from user.

        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec) ::      &
            z(maxinit)
        INTEGER :: &
            jpoints

        ! Local variables
        REAL(KIND=prec) ::      &
            z1(nft), &
            z2(maxinit), &
            d2
        INTEGER ::  &
            i,j,jpoints0,k

        ! Intrinsic Functions
        INTRINSIC ABS

        jpoints0=0
        loop1: DO i=1,jpoints
            DO j=1,nft
                z1(j)=z(j+(i-1)*nft)
            END DO
            IF (i == 1) THEN
                jpoints0=1
                DO j=1,nft
                    z2(j)=z1(j)
                END DO

            ELSE
                DO k=1,jpoints0
                    d2=zero
                    DO j=1,nft
                        d2=d2+ABS(z1(j)-z2(j+(k-1)*nft))
                    END DO
                    IF(d2 <= toler) CYCLE loop1
                END DO
                jpoints0=jpoints0+1
                DO j=1,nft
                    z2(j+(jpoints0-1)*nft)=z1(j)
                END DO
            END IF
        END DO loop1

        DO k=1,jpoints0
            DO j=1,nft
                z(j+(k-1)*nft)=z2(j+(k-1)*nft)
            END DO
        END DO
        jpoints=jpoints0

        RETURN

    END SUBROUTINE cleaning

    !===================================================================
    !  Scaling
    !===================================================================

    SUBROUTINE scaling

        USE param, ONLY : zero, one
        USE initclr, ONLY : &
            a4, &                           ! Auxiliary data matrix a4(nrecord,nft).
            nft, &                          ! Number of features in data, from user.
            nrecord, &                      ! Number of instances in data, from user.
            clin, &
            dlin

        IMPLICIT NONE

        ! Local variables
        REAL(KIND=prec) ::      &
            dm,var
        INTEGER ::  &
            i,j

        ! Intrinsic Functions
        INTRINSIC ABS,SQRT

        DO i=1,nft
            dm=zero
            DO j=1,nrecord
                dm=dm+a4(j,i)
            END DO
            dm=dm/REAL(nrecord,prec)
            var=zero
            DO j=1,nrecord
                var=var+(a4(j,i)-dm)**2
            END DO
            var=var/REAL(nrecord-1,prec)
            var=SQRT(var)

            IF (var >= 1.0E-05_prec) THEN
                clin(i)=one/var
                dlin(i)=-dm/var
                DO j=1,nrecord
                    a4(j,i)=clin(i)*a4(j,i)+dlin(i)
                END DO
            ELSE
                IF (ABS(dm) <= one) THEN
                    clin(i)=one
                    dlin(i)=zero
                ELSE
                    clin(i)=one/dm
                    dlin(i)=zero
                END IF
                DO j=1,nrecord
                    a4(j,i)=clin(i)*a4(j,i)+dlin(i)
                END DO
            END IF
        END DO

        RETURN

    END SUBROUTINE scaling


    !===================================================================
    !  Rescaling
    !===================================================================

    SUBROUTINE rescaling(x,x4)

        USE param, ONLY : zero, one, large
        USE initclr, ONLY : &
            maxdim, &                 ! Maximum number of variables for optimization.
            nft, &                    ! Number of features in data, from user.
            clin, &
            dlin, &
            nc                        ! Current number of clusters, loops from 1 to nclust.

        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec) ::      &
            x(maxdim),x4(maxdim)

        ! Local variables
        REAL(KIND=prec) ::      &
             d0,rabs
        INTEGER ::  &
            i,j

        ! Intrinsic Functions
        INTRINSIC ABS

        rabs=ABS(clin(nft))
        if (rabs > 1.0E-08_prec) THEN
            DO i=1,nc
                d0=zero
                DO j=1,nft-1
                    x4(j+(i-1)*nft)=x(j+(i-1)*nft)*clin(j)/rabs
                    d0=d0+x(j+(i-1)*nft)*dlin(j)
                END DO
                x4(i*nft)=(x(i*nft)+d0-dlin(nft))/rabs
            END DO
        END IF

        RETURN

    END SUBROUTINE rescaling



END MODULE clrsubpro
