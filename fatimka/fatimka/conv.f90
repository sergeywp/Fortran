MODULE VAR
IMPLICIT NONE
INTEGER, PARAMETER :: NI = 20, NJ = 20, NK = 20, NIJ = 20
REAL(8), PARAMETER :: PI = 3.1415926535897932384626433832795
INTEGER NF,NP,L1,L2,L3,M1,M2,M3,N1,N2,N3,ITER,LAST, &
        NTIMES_UVW,NTIMES_P,NTIMES_T,IWRIT,ISHOW,IFLAG,CLOCK1,CLOCK2,II,LK
REAL(8) :: TIME, DT, XL, YL, ZL, thc, visk, dens, SMAX, SMAX_B
REAL(8), DIMENSION (:),       ALLOCATABLE :: X, Y, Z, XU, YV, ZW, XDIF, YDIF, ZDIF, XCV, YCV, ZCV, fx, fy, fz
REAL(8), DIMENSION (:, :),    ALLOCATABLE :: ARI, ARJ, ARK
REAL(8), DIMENSION (:, :, :), ALLOCATABLE :: U, V, W, P, T, &
                                             U0, V0, W0, T0, UHAT, VHAT, WHAT, U1, V1, DU, DV, DW, &
                                             UA, VA, WA, &
                                             BU, BV, BW, PHI, BU0, BV0, BW0, &
                                             VCV, VCVU, VCVV, VCVW, GAM, RHO, &
                                             AIP, AIM, AJP, AJM, AKP, AKM, CON, AP
!----------------------------------------------------------------------
REAL(8) :: Ra, Pr, E, VOLD, Ekin, viskm, Re, B0, Pm, Am, R0

END MODULE VAR
!*********************************************************************

MODULE USER
IMPLICIT NONE
CONTAINS


FUNCTION ACOF1 (X1,Y1)
REAL(8)  ACOF, ACOF1, X1, Y1

  IF(X1.EQ.0.) THEN
    ACOF = 0.0_8
  ELSE
    ACOF = DMAX1(0.0_8,((1.0_8 - 0.1_8 * DABS(Y1/X1))**5))
  END IF

  ACOF1 = X1*ACOF + DMAX1(0.0_8,Y1)

END  FUNCTION ACOF1

FUNCTION visk1 (X1,Y1,Z1)
USE VAR
REAL(8)  visk1, X1, Y1, Z1, X2, Y2, Z2, R

  X2 = 0.5_8 * XL;
  Y2 = 0.5_8 * YL;
  Z2 = 0.0_8

  R = DSQRT( (X1 - X2)**2 + (Y1 - Y2)**2 + (Z1 - Z2)**2 )

  IF( R >= R0 .AND. R <= (R0 + 1.0_8) ) THEN
    visk1 = 1.0_8
  ELSE
    visk1 = 1.D+30
  END IF
END  FUNCTION visk1
!*************************************************
             
SUBROUTINE GRID
USE VAR

  LAST = 1000

  IWRIT = 100
  ISHOW = 1
  NTIMES_UVW = 2
  NTIMES_P = 10
  NTIMES_T = 2

  !R0 = 7.0_8 / 13.0_8
  XL = 2.0_8 * DSQRT(3.0_8) * PI
  YL = 2.0_8 * DSQRT(3.0_8) * PI    
  ZL = 2.0_8 * DSQRT(3.0_8) * PI

  L1 = NI
  M1 = NJ
  N1 = NK

  DT = 1.D-5

  !---------------------------
  Re  = 50.0_8
  B0 = 1.0_8  ! Ha
  Pm = 1.0_8

  Pr = 1.0_8 	
  Ra = 1.D+5
  E = 5.D-4

  !---------------------------

  Am = Re / (DCOSH(0.5_8 * B0) - 1.0_8)
  visk = 1.0_8
  dens = 1.0_8
  viskm = 1.0_8 / Pm   
  thc = 1.0_8 / Pr     

  CALL UGRID1(XL, L1, L2, L3, XU, X, XDIF, XCV, fx)
  CALL UGRID1(YL, M1, M2, M3, YV, Y, YDIF, YCV, fy)
  CALL UGRID1(ZL, N1, N2, N3, ZW, Z, ZDIF, ZCV, fz)

END SUBROUTINE GRID

!-------------------------------------------------

SUBROUTINE START
USE VAR
INTEGER I, J, K
REAL(8) C1, C2, C3

  C1 = 1.0_8 / DSQRT(3.0_8)
  C2 = 1.0_8 / 3.0_8
  C3 = 2.0_8 / 3.0_8

! !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J)
!  DO K = 1, N1
!    DO J = 1, M1
!      DO I = 1, L1
!        U(I,J,K) = C * DSIN(C * XU(I)) * DCOS(C * Y(J)) * DSIN(C * Z(K)) &
!                   + (1.0_8 / 3.0_8) * DCOS(C * XU(I)) * DSIN(C * Y(J)) * DCOS(C * Z(K))
!        V(I,J,K) = - C * DCOS(C * X(I)) * DSIN(C * YV(J)) * DSIN(C * Z(K)) &
!                   + (1.0_8 / 3.0_8) * DSIN(C * X(I)) * DCOS(C * YV(J)) * DCOS(C * Z(K))
!        W(I,J,K) = (2.0_8 / 3.0_8) * DSIN(C * X(I)) * DSIN(C * Y(J)) * DSIN(C * ZW(K))           
!      ENDDO
!    ENDDO
!  ENDDO


  DO K = 1, N1
    DO J = 1, M1
      DO I = 2, L1
        U(I,J,K) = C1 * DSIN(C1 * XU(I)) * DCOS(C1 * Y(J)) * DSIN(C1 * Z(K)) + &
                   C2 * DCOS(C1 * XU(I)) * DSIN(C1 * Y(J)) * DCOS(C1 * Z(K))
      ENDDO
    ENDDO
  ENDDO

  DO K = 1, N1
    DO J = 2, M1
      DO I = 1, L1
        V(I,J,K) = - C1 * DCOS(C1 * X(I)) * DSIN(C1 * YV(J)) * DSIN(C1 * Z(K)) + &
                     C2 * DSIN(C1 * X(I)) * DCOS(C1 * YV(J)) * DCOS(C1 * Z(K))
      ENDDO
    ENDDO
  ENDDO

  DO K = 2, N1
    DO J = 1, M1
      DO I = 1, L1
        W(I,J,K) = C3 * DSIN(C1 * X(I)) * DSIN(C1 * Y(J)) * DSIN(C1 * ZW(K))           
      ENDDO
    ENDDO
  ENDDO


END SUBROUTINE START

!===================================================

SUBROUTINE BOUND
USE VAR
INTEGER I, J, K
REAL(8) C1, C2, C3

  C1 = 1.0_8 / DSQRT(3.0_8)
  C2 = 1.0_8 / 3.0_8
  C3 = 2.0_8 / 3.0_8


  DO K = 1, N1
    DO J = 1, M1
    
        U(2,J,K) = ( C1 * DSIN(C1 * XU(2)) * DCOS(C1 * Y(J)) * DSIN(C1 * Z(K)) + &
                     C2 * DCOS(C1 * XU(2)) * DSIN(C1 * Y(J)) * DCOS(C1 * Z(K)) ) * DEXP(-TIME)

        V(1,J,K) = ( - C1 * DCOS(C1 * X(1)) * DSIN(C1 * YV(J)) * DSIN(C1 * Z(K)) + &
                       C2 * DSIN(C1 * X(1)) * DCOS(C1 * YV(J)) * DCOS(C1 * Z(K)) ) * DEXP(-TIME)

        W(1,J,K) = ( C3 * DSIN(C1 * X(1)) * DSIN(C1 * Y(J)) * DSIN(C1 * ZW(K)) ) * DEXP(-TIME)
           
        U(L1,J,K) = ( C1 * DSIN(C1 * XU(L1)) * DCOS(C1 * Y(J)) * DSIN(C1 * Z(K)) + &
                      C2 * DCOS(C1 * XU(L1)) * DSIN(C1 * Y(J)) * DCOS(C1 * Z(K)) ) * DEXP(-TIME)

        V(L1,J,K) = ( - C1 * DCOS(C1 * X(L1)) * DSIN(C1 * YV(J)) * DSIN(C1 * Z(K)) + &
                        C2 * DSIN(C1 * X(L1)) * DCOS(C1 * YV(J)) * DCOS(C1 * Z(K)) ) * DEXP(-TIME)

        W(L1,J,K) = ( C3 * DSIN(C1 * X(L1)) * DSIN(C1 * Y(J)) * DSIN(C1 * ZW(K)) ) * DEXP(-TIME)         

    ENDDO
  ENDDO

  DO K = 1, N1
    DO I = 1, L1
    
        U(I,1,K)  = ( C1 * DSIN(C1 * XU(I)) * DCOS(C1 * Y(1)) * DSIN(C1 * Z(K)) + &
                     C2 * DCOS(C1 * XU(I)) * DSIN(C1 * Y(1)) * DCOS(C1 * Z(K)) ) * DEXP(-TIME)

        V(I,2,K)  = ( - C1 * DCOS(C1 * X(I)) * DSIN(C1 * YV(2)) * DSIN(C1 * Z(K)) + &
                       C2 * DSIN(C1 * X(I)) * DCOS(C1 * YV(2)) * DCOS(C1 * Z(K)) ) * DEXP(-TIME)

        W(I,1,K)  = ( C3 * DSIN(C1 * X(I)) * DSIN(C1 * Y(1)) * DSIN(C1 * ZW(K)) ) * DEXP(-TIME)   
        
        U(I,M1,K) = ( C1 * DSIN(C1 * XU(I)) * DCOS(C1 * Y(M1)) * DSIN(C1 * Z(K)) + &
                      C2 * DCOS(C1 * XU(I)) * DSIN(C1 * Y(M1)) * DCOS(C1 * Z(K)) ) * DEXP(-TIME)

        V(I,M1,K) = ( - C1 * DCOS(C1 * X(I)) * DSIN(C1 * YV(M1)) * DSIN(C1 * Z(K)) + &
                        C2 * DSIN(C1 * X(I)) * DCOS(C1 * YV(M1)) * DCOS(C1 * Z(K)) ) * DEXP(-TIME)

        W(I,M1,K) = ( C3 * DSIN(C1 * X(I)) * DSIN(C1 * Y(M1)) * DSIN(C1 * ZW(K)) ) * DEXP(-TIME)      

 
    ENDDO
  ENDDO

  DO J = 1, M1
   DO I = 1, L1
    
        U(I,J,1) = ( C1 * DSIN(C1 * XU(I)) * DCOS(C1 * Y(J)) * DSIN(C1 * Z(1)) + &
                     C2 * DCOS(C1 * XU(I)) * DSIN(C1 * Y(J)) * DCOS(C1 * Z(1)) ) * DEXP(-TIME)

        V(I,J,1) = ( - C1 * DCOS(C1 * X(I)) * DSIN(C1 * YV(J)) * DSIN(C1 * Z(1)) + &
                       C2 * DSIN(C1 * X(I)) * DCOS(C1 * YV(J)) * DCOS(C1 * Z(1)) ) * DEXP(-TIME)

        W(I,J,2) = ( C3 * DSIN(C1 * X(I)) * DSIN(C1 * Y(J)) * DSIN(C1 * ZW(2)) ) * DEXP(-TIME)   
        
        U(I,J,N1) = ( C1 * DSIN(C1 * XU(I)) * DCOS(C1 * Y(J)) * DSIN(C1 * Z(N1)) + &
                      C2 * DCOS(C1 * XU(I)) * DSIN(C1 * Y(J)) * DCOS(C1 * Z(N1)) ) * DEXP(-TIME)

        V(I,J,N1) = ( - C1 * DCOS(C1 * X(I)) * DSIN(C1 * YV(J)) * DSIN(C1 * Z(N1)) + &
                        C2 * DSIN(C1 * X(I)) * DCOS(C1 * YV(J)) * DCOS(C1 * Z(N1)) ) * DEXP(-TIME)

        W(I,J,N1) = ( C3 * DSIN(C1 * X(I)) * DSIN(C1 * Y(J)) * DSIN(C1 * ZW(N1)) ) * DEXP(-TIME)          
 
    ENDDO
  ENDDO

END SUBROUTINE BOUND

SUBROUTINE GAMSOR_U
USE VAR
INTEGER I, J, K


  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J)
  DO K = 1, N1
    DO J = 1, M1
      DO I = 2, L1
        GAM(I,J,K) = visk             
      ENDDO
    ENDDO
  ENDDO

  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J)
  DO K = 2, N2
    DO J = 2, M2
      DO I = 3, L2

        CON(I,J,K) = 0.0_8

        AP(I,J,K) = 0.0_8

      ENDDO
    ENDDO
  ENDDO
                   
END SUBROUTINE GAMSOR_U


SUBROUTINE GAMSOR_V
USE VAR
INTEGER I, J, K
REAL(8) Bz_t, Bz_b, BVz_t, BVz_b, Bx_t, Bx_b, BVx_t, BVx_b, S1, S2, FcorY, FbY

  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J)
  DO K = 1, N1
    DO J = 2, M1
      DO I = 1, L1
        GAM(I,J,K) = visk          
      ENDDO
    ENDDO
  ENDDO

  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J)
  DO K = 2, N2
    DO J = 3, M2
      DO I= 2, L2

        CON(I,J,K) = 0.0_8

        AP(I,J,K) = 0.0_8
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE GAMSOR_V


SUBROUTINE GAMSOR_W
USE VAR
INTEGER I, J, K

  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J)
  DO K = 2, N1
    DO J = 1, M1
      DO I = 1, L1
        GAM(I,J,K) = visk          
      ENDDO
    ENDDO
  ENDDO

  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J)
  DO K = 3, N2
    DO J = 2, M2
      DO I = 2, L2

        CON(I,J,K) = 0.0_8

        AP(I,J,K) = 0.0_8

      ENDDO
    ENDDO
  ENDDO
                   
END SUBROUTINE GAMSOR_W

SUBROUTINE GAMSOR_T
USE VAR
INTEGER I, J, K

  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J)
  DO K = 1, N1
    DO J = 1, M1
      DO I = 1, L1
        RHO(I,J,K) = 1.0_8
        GAM(I,J,K) = thc
        CON(I,J,K) = 0.0_8
        AP(I,J,K)  = 0.0_8  
      ENDDO
    ENDDO
  ENDDO
                                                             
END SUBROUTINE GAMSOR_T
!-------------------------------------------------

SUBROUTINE OUTPUT
USE VAR
REAL(8) Ux, Uy, Uz, b
INTEGER I, J, K
REAL(8) C1, C2, C3, delU, delV, delW


  IF(MOD(ITER,ISHOW).EQ.0) THEN



  C1 = 1.0_8 / DSQRT(3.0_8)
  C2 = 1.0_8 / 3.0_8
  C3 = 2.0_8 / 3.0_8

  delU = 0.0_8; delV = 0.0_8; delW = 0.0_8;
  DO K = 1, N1
    DO J = 1, M1
      DO I = 2, L1
        UA(I,J,K) = (C1 * DSIN(C1 * XU(I)) * DCOS(C1 * Y(J)) * DSIN(C1 * Z(K)) + &
                   C2 * DCOS(C1 * XU(I)) * DSIN(C1 * Y(J)) * DCOS(C1 * Z(K)) ) * DEXP(-TIME)
        delU = DMAX1(DABS(U(I, J, K) - UA(I, J, K)), delU)

      ENDDO
    ENDDO
  ENDDO

  DO K = 1, N1
    DO J = 2, M1
      DO I = 1, L1
        VA(I,J,K) = (- C1 * DCOS(C1 * X(I)) * DSIN(C1 * YV(J)) * DSIN(C1 * Z(K)) + &
                     C2 * DSIN(C1 * X(I)) * DCOS(C1 * YV(J)) * DCOS(C1 * Z(K)) )  * DEXP(-TIME)
		delV = DMAX1(DABS(V(I, J, K) - VA(I, J, K)), delV)
      ENDDO
    ENDDO
  ENDDO

  DO K = 2, N1
    DO J = 1, M1
      DO I = 1, L1
        WA(I,J,K) = (C3 * DSIN(C1 * X(I)) * DSIN(C1 * Y(J)) * DSIN(C1 * ZW(K)) ) * DEXP(-TIME)
        delW = DMAX1(DABS(W(I, J, K) - WA(I, J, K)), delW)       
      ENDDO
    ENDDO
  ENDDO



  WRITE(*,501)
  WRITE(1,501)
  501	FORMAT(5X,'ITER',8X,'SMAX',12X,'TIME')

  WRITE(*,503)ITER,SMAX,TIME
  WRITE(1,503)ITER,SMAX,TIME

  503	FORMAT (I8,1P6E16.4)

  WRITE(*,502)
  WRITE(1,502)
  502	FORMAT(7X,'delU',8X,'delV',8X,'delW')

  WRITE(*,504)delU, delV,delW
  WRITE(1,504)delU, delV,delW
  WRITE(*,*)
  WRITE(1,*)

504	FORMAT (1P7E16.4)

      WRITE(*,9)
      WRITE(1,9)
9     FORMAT(1X,'--------------------------------------------------------------')

      	END IF      

      	IF (ITER.NE.LAST) RETURN


	OPEN(UNIT=40,FILE='ALL.DAT',STATUS='UNKNOWN')
	WRITE(40,*)'VARIABLES = "X", "Y", "Z","P", "U", "V", "W", "b"'
	WRITE(40,*)'ZONE I=',L3,', J=',M3,', K=',N3,', F=POINT'
                
	DO K = 2, N2
		DO J = 2, M2
			DO I = 2, L2
             
			  Ux = 0.5_8*(U(I,J,K) + U(I+1,J,K))
			  Uy = 0.5_8*(V(I,J,K) + V(I,J+1,K))
			  Uz = 0.5_8*(W(I,J,K) + W(I,J,K+1))	

			  b = ARI(J,K)*(U(I,J,K) - U(I+1,J,K)) + &
			      ARJ(I,K)*(V(I,J,K) - V(I,J+1,K)) + &
                              ARK(I,J)*(W(I,J,K) - W(I,J,K+1))

			  b = b/VCV(I,J,K)

			  WRITE(40,'(1P7E15.5)') X(I),Y(J),Z(K),P(I,J,K),Ux,Uy,Uz,b

      ENDDO           
    ENDDO
  ENDDO

	ENDFILE 40
	CLOSE(40)


	OPEN(UNIT=41,FILE='ALL1.DAT',STATUS='UNKNOWN')
	WRITE(41,*)'VARIABLES = "X", "Y", "Z","U","V","W","UA","VA","WA","delU","delV","delW"'    
	WRITE(41,*)'ZONE I=',L1,', J=',M1,', K=',N1,', F=POINT'
                
	DO K = 1, N1
		DO J = 1, M1
			DO I = 1, L1
             
			  WRITE(41,'(1P12E15.5)') X(I),Y(J),Z(K), U(I, J, K), V(I, J, K), W(I, J, K), UA(I, J, K), VA(I, J, K), WA(I, J, K), &
                        DABS(U(I, J, K) - UA(I, J, K)), DABS(V(I, J, K) - VA(I, J, K)), DABS(W(I, J, K) - WA(I, J, K))

      ENDDO           
    ENDDO
  ENDDO

	ENDFILE 41
	CLOSE(41)

END SUBROUTINE OUTPUT
!==============================================
END MODULE USER
