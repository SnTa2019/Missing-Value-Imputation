!*************************************************************************
!*                                                                       *
!*     Computation of function and subgradients values of clusterwise    *
!*     linear regression problem. (last modified 28.02.2017).            *
!*                                                                       *
!*************************************************************************
!*
!*     Modules included:
!*
!*     clrfunctionmod         !
!*     clrobjfun              !
!*


MODULE clrfunctionmod ! Computation of function and (sub)gradients values for
                      ! clusterwise linear regression problem.

    USE r_precision, ONLY : prec      ! Precision for reals.
    IMPLICIT NONE

    ! MODULE clrfunctionmod includes the following subroutines.

    PUBLIC :: &
        auxfunc, &         ! Computation of aux clusterwise linear regression problem.
        clrfunc, &         ! Computation of clusterwise linear regression problem.
        fgrad, &           ! Computation of subgradient.
        fgrad1, &          ! Computation of gradient of the first DC component.
        fgrad2, &          ! Computation of subgradient of the second DC component.
        fvb, &
        gradient

!    PRIVATE :: &


CONTAINS

    !==========================================================
    ! Computation of aux clusterwise linear regression problem
    !==========================================================
    SUBROUTINE auxfunc(x,fval)

        USE param, ONLY : zero
        USE initimp, ONLY : &
            nft, &         ! Number of used features.
            nrecord        ! Number of records in data.
        USE initclr, ONLY : &
            maxdim, &      ! Maximum number of variables in optimization,
                           !   maxdim >= nclust * nft.
            a4, &          ! Data matrix a4(maxrec,maxnft)
            dminim, &      ! dminim(nrecord)
            naux           ! The total number of auxiliary problems computed.

        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec) ::  &
            x(maxdim), &
            fval

        ! Local variables
        REAL(KIND=prec) ::  &
            f1,f2,f3
        INTEGER :: &
            i,j

        INTRINSIC ::  DMIN1

        fval=zero

        DO i=1,nrecord
            f1=dminim(i)
            f2=x(nft)

            DO j=1,nft-1
                f2=f2+a4(i,j)*x(j)
            END DO

            f3=(a4(i,nft)-f2)**2
            fval=fval+DMIN1(f1,f3)
        END DO

        naux=naux+1
        RETURN

    END SUBROUTINE auxfunc


    !======================================================
    ! Computation of clusterwise linear regression problem
    !======================================================
    SUBROUTINE clrfunc(x,f)

        USE param, ONLY : zero, large
        USE initimp, ONLY : &
            nft, &         ! Number of used features.
            nrecord        ! Number of records in data.
        USE initclr, ONLY : &
            maxdim, &      ! Maximum number of variables in optimization,
                           !   maxdim >= nclust * nft.
            a4, &          ! Data matrix a4(maxrec,maxnft).
            nc, &          ! Current number of clusters, loops from 1 to nclust.
            nfit           ! The number of function evaluations.

        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec) ::  &
            x(maxdim), &
            f

        ! Local variables
        REAL(KIND=prec) ::  &
            f2,d1,d2
        INTEGER :: &
            i,j,k

        INTRINSIC ::  DMIN1

        f=zero
        DO i=1,nrecord
            f2=large
            DO k=1,nc
                d1=x(k*nft)
                DO j=1,nft-1
                    d1=d1+a4(i,j)*x(j+(k-1)*nft)
                END DO
                d2=(d1-a4(i,nft))**2
                f2=DMIN1(f2,d2)
            END DO
            f=f+f2
        END DO
        nfit=nfit+1
        RETURN

    END SUBROUTINE clrfunc

    !================================
    ! Computation of the subgradient
    !================================

    SUBROUTINE fgrad(x,grad)

        USE param, ONLY : zero, two, large
        USE initimp, ONLY : &
            nft, &         ! Number of used features.
            nrecord        ! Number of records in data.
        USE initclr, ONLY : &
            maxdim, &      ! Maximum number of variables in optimization,
                           !   maxdim >= nclust * nft.
            a4, &          ! Data matrix a4(maxrec,maxnft).
            dminim, &      !
            nc, &          ! Current number of clusters, loops from 1 to nclust.
            ns, &          ! Switch for auxiliary and real clustering problem.
            m              ! Number of variables in optimization:
                           !   m = nft    if ns = 1,
                           !   m = nft*nc if ns = 2.

        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec) ::  &
            x(maxdim), &
            grad(maxdim)

        ! Local variables
        REAL(KIND=prec) ::  &
            d1,d2,d4,f12
        INTEGER :: &
            i,j,k,jmin

        jmin = 1
        d4 = zero

        DO i=1,m
            grad(i)=zero
        END DO

        IF (ns == 1) THEN
            DO i=1,nrecord
                d1 = x(nft)
                DO j=1,nft-1
                    d1=d1+a4(i,j)*x(j)
                END DO
                d1 = d1-a4(i,nft)
                d2=d1*d1
                IF (d2 <= dminim(i)) THEN
                    DO j=1,nft-1
                        grad(j)=grad(j)+two*d1*a4(i,j)
                    END DO
                    grad(nft)=grad(nft)+two*d1
                END IF
            END DO
        ELSE

            DO i=1,nrecord
                f12 = large
                DO k=1,nc
                    d1=x(k*nft)
                    DO j=1,nft-1
                        d1=d1+a4(i,j)*x(j+(k-1)*nft)
                    END DO
                    d2=(d1-a4(i,nft))**2
                    IF (f12 > d2) THEN
                        f12=d2
                        jmin=k
                        d4=d1-a4(i,nft)
                    END IF
                END DO
                DO j=1,nft-1
                    grad(j+(jmin-1)*nft)=grad(j+(jmin-1)*nft) + two*d4*a4(i,j)
                END DO
                grad(jmin*nft)=grad(jmin*nft)+two*d4

            END DO
        END IF
        !grad=grad/REAL(nrecord,prec)

        RETURN

    END SUBROUTINE fgrad

    !=====================================================
    ! Computation of the subgradient (first DC component)
    !=====================================================
    SUBROUTINE fgrad1(x,grad)

        USE param, ONLY : zero, two
        USE initimp, ONLY : &
            nft, &         ! Number of used features.
            nrecord        ! Number of records in data.
        USE initclr, ONLY : &
            maxdim, &      ! Maximum number of variables in optimization,
                           !   maxdim >= nclust * nft.
            a4, &          ! Data matrix a4(maxrec,maxnft).
            nc, &          ! Current number of clusters, loops from 1 to nclust.
            ns, &          ! Switch for auxiliary and real clustering problem.
            m              ! Number of variables in optimization:
                           !   m = nft    if ns = 1,
                           !   m = nft*nc if ns = 2.

        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec) ::  &
            x(maxdim), &
            grad(maxdim)

        ! Local variables
        REAL(KIND=prec) ::  &
            d1
        INTEGER :: &
            i,j,k

        DO i=1,m
            grad(i)=zero
        END DO

        IF (ns == 1) THEN
            DO i=1,nrecord
                d1=x(nft)
                DO j=1,nft-1
                    d1=d1+a4(i,j)*x(j)
                END DO
                DO j=1,nft-1
                    grad(j)=grad(j)+two*(d1-a4(i,nft))*a4(i,j)
                END DO
                grad(nft)=grad(nft)+two*(d1-a4(i,nft))
            END DO
        ELSE
            DO i=1,nrecord
                DO k=1,nc
                    d1=x(k*nft)
                    DO j=1,nft-1
                        d1=d1+a4(i,j)*x(j+(k-1)*nft)
                    END DO
                    d1=d1-a4(i,nft)
                    DO j=1,nft-1
                        grad(j+(k-1)*nft)=grad(j+(k-1)*nft)+two*d1*a4(i,j)
                    END DO
                    grad(k*nft)=grad(k*nft)+two*d1
                END DO
            END DO
        END IF

        RETURN

    END SUBROUTINE fgrad1

    !======================================================
    ! Computation of the subgradient (second DC component)
    !======================================================
    SUBROUTINE fgrad2(x,grad)

        USE param, ONLY : zero, two
        USE initimp, ONLY : &
            nft, &         ! Number of used features.
            nrecord        ! Number of records in data.
        USE initclr, ONLY : &
            maxdim, &      ! Maximum number of variables in optimization,
                           !   maxdim >= nclust * nft.
            a4, &          ! Data matrix a4(maxrec,maxnft).
            nc, &          ! Current number of clusters, loops from 1 to nclust.
            dminim, &      ! dminim(nrecord).
            ns, &          ! Switch for auxiliary and real clustering problem.
            m              ! Number of variables in optimization:
                           !   m = mf    if ns = 1,
                           !   m = mf*nc if ns = 2.


        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec) ::  &
            x(maxdim), &
            grad(maxdim)

        ! Local variables
        REAL(KIND=prec), DIMENSION(nc) ::  &
            d2,d3          ! Auxiliary arrays.

        REAL(KIND=prec) ::  &
            d0,d1,d4,d5
        INTEGER :: &
            i,j,k,jindex

        jindex = 0
        DO i=1,m
            grad(i)=zero
        END DO

        IF (ns == 1) THEN
            DO i=1,nrecord
                d1=x(nft)
                DO j=1,nft-1
                    d1=d1+a4(i,j)*x(j)
                END DO
                d0=(d1-a4(i,nft))**2
                IF (d0 >= dminim(i)) THEN
                    DO j=1,nft-1
                        grad(j)=grad(j)+two*(d1-a4(i,nft))*a4(i,j)
                    END DO
                    grad(nft)=grad(nft)+two*(d1-a4(i,nft))
                END IF
            END DO
        ELSE
            DO i=1,nrecord
                DO k=1,nc
                    d3(k)=x(k*nft)
                    DO j=1,nft-1
                        d3(k)=d3(k)+a4(i,j)*x(j+(k-1)*nft)
                    END DO
                    d2(k)=(d3(k)-a4(i,nft))**2
                END DO
                d4=zero
                DO j=1,nc
                    d5=zero
                    DO k=1,nc
                        IF (k /= j) THEN
                            d5=d5+d2(k)
                        END IF
                    END DO
                    IF (d4 < d5) THEN
                        d4=d5
                        jindex=j
                    END IF
                END DO

                DO j=1,nc
                    IF (j /= jindex) THEN
                        DO k=1,nft-1
                            grad(k+(j-1)*nft)=grad(k+(j-1)*nft)+two*(d3(j)-a4(i,nft))*a4(i,k)
                        END DO
                        grad(j*nft)=grad(j*nft)+two*(d3(j)-a4(i,nft))
                    END IF
                END DO
            END DO
        END IF

        RETURN

    END SUBROUTINE fgrad2


    SUBROUTINE fvb(x,f)

        USE param, ONLY : zero, one
        USE initclr, ONLY : &
            maxdim, &      ! Maximum number of variables.
            nft, &
            a4, &
            nel1,nob1, &
            nupdate,naux2

        ! Subroutine arguments
        REAL(KIND=prec) ::  &
            x(maxdim), &
            f

        ! Local variables
        REAL(KIND=prec) ::      &
            d1,d2
        INTEGER :: &
            j,k

        d2=zero
        do k=1,nel1
            d1=x(nft)
            do j=1,nft-1
                d1=d1+a4(nob1(k),j)*x(j)
            end do
            d2=d2+(d1-a4(nob1(k),nft))**2
        end do
        f=d2
        nupdate=nupdate+1
        naux2=naux2+1

        RETURN
    END SUBROUTINE fvb

    SUBROUTINE gradient(x,grad)

        USE param, ONLY : zero, two
        USE initclr, ONLY : &
            maxdim, &      ! Maximum number of variables.
            nft, &
            a4, &
            nel1,nob1

        ! Subroutine arguments
        REAL(KIND=prec) ::  &
            x(maxdim), &
            grad(maxdim)

        ! Local variables
        REAL(KIND=prec) ::      &
            d1
        INTEGER :: &
            j,k

        DO j=1,nft
            grad(j)=zero
        END DO
        DO k=1,nel1
            d1=x(nft)
            DO j=1,nft-1
                d1=d1+a4(nob1(k),j)*x(j)
            END DO
            d1=d1-a4(nob1(k),nft)
            DO j=1,nft-1
                grad(j)=grad(j)+two*d1*a4(nob1(k),j)
            END DO
            grad(nft)=grad(nft)+two*d1
        END DO

        RETURN
    END SUBROUTINE gradient


END MODULE clrfunctionmod

MODULE clrobjfun  ! Computation of the value and the subgradient of the
                  ! objective function.

    USE r_precision, ONLY : prec  ! Precision for reals.
    IMPLICIT NONE

    PUBLIC :: &
        myf,    &   ! Computation of the value of the objective.
        myg         ! Computation of the subgradient of the objective.

CONTAINS
    !************************************************************************
    !*                                                                      *
    !*     * SUBROUTINE myf *                                               *
    !*                                                                      *
    !*     Computation of the value of the objective.                       *
    !*                                                                      *
    !************************************************************************

    SUBROUTINE myf(n,x,f,iterm)

        USE param, ONLY : large    ! Parameters.
        USE initclr, ONLY : ns     ! Switch for auxiliary and real clustering problem.

        USE clrfunctionmod, ONLY : &
            auxfunc, &             ! S Computation of auxiliary clustering problem
            clrfunc, &                ! S Computation of clustering problem
            fvb

        IMPLICIT NONE

        ! Array Arguments
        REAL(KIND=prec), DIMENSION(n), INTENT(IN) :: &
            x  ! Vector of variables.

        ! Scalar Arguments
        REAL(KIND=prec), INTENT(OUT) :: f  ! Value of the function.
        INTEGER, INTENT(IN) :: n           ! Number of variables.
        INTEGER, INTENT(OUT) :: iterm      ! Cause of termination:
                                           !   0  - Everything is ok.
                                           !  -3  - Failure in function calculations

        iterm = 0

        ! Function evaluation
        IF (ns == 1) THEN  ! Auxiliary problem
            CALL auxfunc(x,f)
        ELSE IF (ns==3) THEN
            CALL fvb(x,f)
        ELSE  ! Clustering problem
            CALL clrfunc(x,f)
        END IF

        ! Error checking.
        !    IF (f > large) iterm = -3  !
        !    IF (f < -large) iterm = -3 !
        RETURN

    END SUBROUTINE myf

    !************************************************************************
    !*                                                                      *
    !*     * SUBROUTINE myg *                                               *
    !*                                                                      *
    !*     Computation of the subgradient of the objective function.        *
    !*                                                                      *
    !************************************************************************

    SUBROUTINE myg(n,x,g,iterm)

        USE initclr, ONLY : ns       ! Switch for auxiliary and real clustering problem.
        USE clrfunctionmod, ONLY : &
            fgrad, &
            fgrad1, &                ! S Computation of the subgradient of the
            fgrad2, &                !   clusterwise linear regression function
            gradient

        IMPLICIT NONE

        ! Array Arguments
        REAL(KIND=prec), DIMENSION(n), INTENT(IN) :: x  ! Vector of variables.
        REAL(KIND=prec), DIMENSION(n), INTENT(OUT) :: g ! Subgradient.

        ! Scalar Arguments
        INTEGER, INTENT(IN) :: n                        ! Number of variables.
        INTEGER, INTENT(OUT) :: iterm                   ! Cause of termination:
                                                        !   0  - Everything is ok.
                                                        !  -3  - Failure in subgradient calculations
                                                        !        (assigned by the user).


        iterm = 0

        ! Gradient evaluation.
        IF (ns == 3) THEN
            CALL gradient(x,g)
        ELSE
            CALL fgrad(x,g)
        END IF
        RETURN

    END SUBROUTINE myg

END MODULE clrobjfun

