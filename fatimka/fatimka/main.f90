PROGRAM MAIN
USE VAR
USE USER
IMPLICIT NONE

  INTEGER :: clock_start, clock_end, clock_rate, IERR
  REAL(8) :: elapsed_time

  CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate)
  CALL SYSTEM_CLOCK(COUNT=clock_start)

  ALLOCATE( U(NI, NJ, NK), V(NI, NJ, NK), W(NI, NJ, NK), P(NI, NJ, NK), T(NI, NJ, NK), &
            U0(NI, NJ, NK), V0(NI, NJ, NK), W0(NI, NJ, NK), T0(NI, NJ, NK), UHAT(NI, NJ, NK), &
            VHAT(NI, NJ, NK), WHAT(NI, NJ, NK), U1(NI, NJ, NK), V1(NI, NJ, NK), DU(NI, NJ, NK), &
            DV(NI, NJ, NK), DW(NI, NJ, NK), BU(NI, NJ, NK), BV(NI, NJ, NK), BW(NI, NJ, NK), &
            PHI(NI, NJ, NK), BU0(NI, NJ, NK), BV0(NI, NJ, NK), BW0(NI, NJ, NK), VCV(NI, NJ, NK), &
            VCVU(NI, NJ, NK), VCVV(NI, NJ, NK), VCVW(NI, NJ, NK), GAM(NI, NJ, NK), RHO(NI, NJ, NK), &
            AIP(NI, NJ, NK), AIM(NI, NJ, NK), AJP(NI, NJ, NK), AJM(NI, NJ, NK), AKP(NI, NJ, NK), &
            AKM(NI, NJ, NK), CON(NI, NJ, NK), AP(NI, NJ, NK), UA(NI, NJ, NK), VA(NI, NJ, NK), WA(NI, NJ, NK), STAT = IERR )


  ALLOCATE( X(NI), XU(NI), XDIF(NI), XCV(NI), fx(NI), Y(NJ), YV(NJ), YDIF(NJ), YCV(NJ), fy(NJ), &
            Z(NK), ZW(NK), ZDIF(NK), ZCV(NK), fz(NK), ARI(NJ, NK), ARJ(NI, NK), ARK(NI, NJ), STAT = IERR )

	OPEN(UNIT=1,FILE='Q.OUT',STATUS='UNKNOWN')

	CALL GRID; CALL GEOM; CALL START; CALL OUTPUT; CALL F0_F

  DO ITER = 1, LAST
    TIME = TIME+DT
    SMAX = 1.0_8
    II = 0
    DO  WHILE ( SMAX > 1.D-4 )
    CALL BOUND
    CALL PISO
    !CALL GAMSOR_T;  CALL COF_T; CALL SOLVE(2, 2, 2, NTIMES_T, T)
    II = II + 1            
  ENDDO !WHILE

  CALL F0_F; CALL OUTPUT

  ENDDO



  DEALLOCATE( X, Y, Z, XU, YV, ZW, XDIF, YDIF, ZDIF, XCV, YCV, ZCV, fx, fy, fz, &
              ARI, ARJ, ARK, U, V, W, P, T, U0, V0, W0, T0, UHAT, VHAT, U1, V1, DU, DV, DW, &  
              BU, BV, BW, PHI, BU0, BV0, BW0, VCV, VCVU, VCVV, VCVW, GAM, RHO, AIP, AIM, AJP, &
              AJM, AKP, AKM, CON, AP, UA, VA, WA, STAT = IERR )




  CALL SYSTEM_CLOCK(COUNT=clock_end)

  elapsed_time=DBLE((clock_end-clock_start)/clock_rate)
  WRITE(*, *)elapsed_time

END PROGRAM MAIN

 
!*********************************************************************
SUBROUTINE UGRID1(XL, L1, L2, L3, XU, X, XDIF, XCV, fx)
IMPLICIT NONE
INTEGER  L1, L2, L3, I
REAL(8) XL, DX, XU(L1), X(L1), XDIF(L1), XCV(L1), fx(L1)

  L2 = L1-1
  L3 = L1-2

  DX=XL/DBLE(L3)
  XU(2) = 0.0_8

  DO I = 3, L1
    XU(I) = XU(I-1) + DX
  ENDDO

  X(1) = XU(2)
  
  DO I= 2, L2
    X(I) = 0.5_8*( XU(I) + XU(I+1) )
  ENDDO

  X(L1) = XU(L1)
  DO I = 2, L1
    XDIF(I) = X(I) - X(I-1)
  ENDDO

  DO I = 2, L2
    XCV(I) = XU(I+1) - XU(I)
  ENDDO

  fx(3) = (XU(3) - X(1))/(X(3) - X(1))                            
  
  DO I=4,L3
    fx(I) = (XU(I) - X(I-1))/(X(I) - X(I-1))
  ENDDO

  fx(L2) = (XU(L2) - X(L3))/(X(L1) - X(L3))                            

END SUBROUTINE UGRID1
!*********************************************************************
SUBROUTINE GEOM
USE VAR
IMPLICIT NONE
INTEGER I, J, K
                                     	
  DO K = 2, N2
    DO J = 2, M2
      ARI(J,K) = YCV(J)*ZCV(K)
    ENDDO
  ENDDO

  DO K = 2, N2
    DO I = 2, L2
      ARJ(I,K) = XCV(I)*ZCV(K)
    ENDDO
  ENDDO

  DO J = 2, M2
    DO I = 2, L2
      ARK(I,J) = XCV(I)*YCV(J)                                                          
    ENDDO
  ENDDO

  DO K = 2, N2
    DO J = 2, M2
      DO I = 2, L2
        VCV(I,J,K) = XCV(I)*YCV(J)*ZCV(K)
      ENDDO
    ENDDO
  ENDDO

  DO K = 2, N2
    DO J = 2, M2
      VCVU(3,J,K) = (X(3)-X(1))*YCV(J)*ZCV(K)
      DO I = 4, L3
        VCVU(I,J,K) = XDIF(I)*YCV(J)*ZCV(K)
      ENDDO
        VCVU(L2,J,K) = (X(L1)-X(L1-2))*YCV(J)*ZCV(K)
    ENDDO
  ENDDO

  DO K = 2, N2
    DO J = 4, M3
      DO I = 2, L2
        VCVV(I,J,K) = XCV(I)*YDIF(J)*ZCV(K)
      ENDDO
    ENDDO

    DO I = 2, L2
      VCVV(I,3,K) = XCV(I)*ZCV(K)*(Y(3)-Y(1))
      VCVV(I,M2,K) = XCV(I)*ZCV(K)*(Y(M1)-Y(M3))
    ENDDO
  ENDDO

  DO K = 4, N3
    DO J = 2, M2
      DO I = 2, L2
        VCVW(I,J,K) = XCV(I)*YCV(J)*ZDIF(K)
      ENDDO
    ENDDO
  ENDDO

  DO J = 2, M2
    DO I = 2, L2
      VCVW(I,J,3) = XCV(I)*YCV(J)*(Z(3)-Z(1))
      VCVW(I,J,N2) = XCV(I)*YCV(J)*(Z(N1)-Z(N3))
    ENDDO
  ENDDO

END SUBROUTINE GEOM
!*********************************************************************
SUBROUTINE COF_T
USE VAR
USE USER
IMPLICIT NONE
INTEGER I, J, K
REAL(8) DIFF, FLOW, R, AP0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J,DIFF,FLOW)
  DO K = 2, N2
    DO J = 2, M2
      DO I = 2, L2
        !----------------- AIM, AIP -----------------
        DIFF = GAM(I-1,J,K)*GAM(I,J,K)*ARI(J,K)/( GAM(I-1,J,K)*(X(I)-XU(I)) + GAM(I,J,K)*(XU(I)-X(I-1)) )
        FLOW = ( RHO(I-1,J,K)*(XU(I)-X(I-1)) + RHO(I,J,K)*(X(I)-XU(I)) )*U(I,J,K)*ARI(J,K)/XDIF(I)
        AIM(I,J,K) = ACOF1(DIFF,FLOW)                                        
        DIFF = GAM(I,J,K)*GAM(I+1,J,K)*ARI(J,K)/( GAM(I,J,K)*(X(I+1)-XU(I+1)) + GAM(I+1,J,K)*(XU(I+1)-X(I)) )         	    
        FLOW = ( RHO(I,J,K)*(XU(I+1)-X(I)) + RHO(I+1,J,K)*(X(I+1)-XU(I+1)) )*U(I+1,J,K)*ARI(J,K)/XDIF(I+1)	          			    
        AIP(I,J,K) = ACOF1(DIFF,-FLOW)  
        !----------------- AJM, AJP -----------------
        DIFF = GAM(I,J-1,K)*GAM(I,J,K)*ARJ(I,K)/( GAM(I,J-1,K)*(Y(J)-YV(J)) + GAM(I,J,K)*(YV(J)-Y(J-1)) )
        FLOW = ( RHO(I,J-1,K)*(YV(J)-Y(J-1)) + RHO(I,J,K)*(Y(J)-YV(J)) )*V(I,J,K)*ARJ(I,K)/YDIF(J)	               
        AJM(I,J,K) = ACOF1(DIFF,FLOW)                                          
        DIFF = GAM(I,J,K)*GAM(I,J+1,K)*ARJ(I,K)/( GAM(I,J,K)*(Y(J+1)-YV(J+1)) + GAM(I,J+1,K)*(YV(J+1)-Y(J)) )	
        FLOW = ( RHO(I,J,K)*(YV(J+1)-Y(J)) + RHO(I,J+1,K)*(Y(J+1)-YV(J+1)) )*V(I,J+1,K)*ARJ(I,K)/YDIF(J+1)	
        AJP(I,J,K) = ACOF1(DIFF,-FLOW)
        !----------------- AKM, AKP -----------------
        DIFF = GAM(I,J,K-1)*GAM(I,J,K)*ARK(I,J)/( GAM(I,J,K-1)*(Z(K)-ZW(K)) + GAM(I,J,K)*(ZW(K)-Z(K-1)) )
        FLOW = ( RHO(I,J,K-1)*(ZW(K)-Z(K-1)) + RHO(I,J,K)*(Z(K)-ZW(K)) )*W(I,J,K)*ARK(I,J)/ZDIF(K)          
        AKM(I,J,K) = ACOF1(DIFF,FLOW)
        DIFF = GAM(I,J,K)*GAM(I,J,K+1)*ARK(I,J)/( GAM(I,J,K)*(Z(K+1)-ZW(K+1)) + GAM(I,J,K+1)*(ZW(K+1)-Z(K)) )     
        FLOW = ( RHO(I,J,K)*(ZW(K+1)-Z(K)) + RHO(I,J,K+1)*(Z(K+1)-ZW(K+1)) )*W(I,J,K+1)*ARK(I,J)/ZDIF(K+1)          
        AKP(I,J,K) = ACOF1(DIFF,-FLOW)
      ENDDO
    ENDDO
  ENDDO
!--------------------------- Boundary conditions of the second kind on z = 0 ------------------
  
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I)
  DO J = 2, M2
    DO I = 2, L2                        
      R = DSQRT( (X(I) - 0.5_8 * XL )**2 + ( Y(J) - 0.5_8 * YL)**2 )
      IF( R >= R0 .AND. R <= R0 + 1.0_8) THEN	                 
        AKM(I,J,2) = 0.0_8                                          
      END IF               
    ENDDO
  ENDDO
  
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J,AP0)
  DO K = 2, N2
    DO J = 2, M2
      DO I = 2, L2
        !----------------- AP0, CON, AP -----------------
        AP0 = RHO(I,J,K)*VCV(I,J,K)/DT
        CON(I,J,K) = CON(I,J,K)+AP0*T0(I,J,K)
        AP(I,J,K) = -AP(I,J,K)*VCV(I,J,K)+AP0+ &
        AIM(I,J,K)+AIP(I,J,K)+ &
        AJM(I,J,K)+AJP(I,J,K)+ &
        AKM(I,J,K)+AKP(I,J,K)    
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE COF_T             
!*********************************************************************

SUBROUTINE SOLVE(IST, JST, KST, NTIMES, F)
USE VAR
IMPLICIT NONE
INTEGER I, J, K, IST, JST, KST, NTIMES, MT
REAL(8) :: DENOM
REAL(8), DIMENSION (NIJ) :: PT, QT
REAL(8), DIMENSION (NI,NJ,NK) :: F, TEMPS

  DO MT = 1, NTIMES
  !---------------- X_dir -------------
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J)
    DO K = KST, NK-1
      DO J = JST, NJ-1
        DO I = IST, NI-1
          TEMPS(I,J,K) = CON(I,J,K) &
                         + AJP(I,J,K)*F(I,J+1,K) &
                         + AJM(I,J,K)*F(I,J-1,K) &
                         + AKP(I,J,K)*F(I,J,K+1) &
                         + AKM(I,J,K)*F(I,J,K-1)
        ENDDO
      ENDDO
    ENDDO                      

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J,PT,QT,DENOM)
    DO K = KST, NK-1
      DO J = JST, NJ-1
        PT(IST-1) = 0.0_8
        QT(IST-1) = F(IST-1,J,K)
        DO I = IST, NI-1
          DENOM = AP(I,J,K)-PT(I-1)*AIM(I,J,K)
          PT(I) = AIP(I,J,K)/(DENOM+1.D-30)
          QT(I) = (TEMPS(I,J,K)+AIM(I,J,K)*QT(I-1))/(DENOM+1.D-30)
        ENDDO
        
        DO I = NI-1, IST, -1  
          F(I,J,K) = F(I+1,J,K)*PT(I)+QT(I)
        ENDDO

      ENDDO
    ENDDO

    !------------ Y_dir ---------------
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J)
    DO K = KST, NK-1
      DO J = JST, NJ-1
        DO I = IST, NI-1
          TEMPS(I,J,K) = CON(I,J,K) &
                         +AIP(I,J,K)*F(I+1,J,K) &
                         +AIM(I,J,K)*F(I-1,J,K) &
                         +AKP(I,J,K)*F(I,J,K+1) &
                         +AKM(I,J,K)*F(I,J,K-1)
        ENDDO
      ENDDO
    ENDDO

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J,PT,QT,DENOM)
    DO K = KST, NK-1
      DO I = IST, NI-1
        PT(JST-1) = 0.0_8
        QT(JST-1) = F(I,JST-1,K)
        DO J = JST, NJ-1
          DENOM = AP(I,J,K)-PT(J-1)*AJM(I,J,K)
          PT(J) = AJP(I,J,K)/(DENOM+1.D-30)
          QT(J) = (TEMPS(I,J,K)+AJM(I,J,K)*QT(J-1))/(DENOM+1.D-30)
        ENDDO

        DO J = NJ-1, JST, -1
          F(I,J,K) = F(I,J+1,K)*PT(J)+QT(J)
        ENDDO
      ENDDO
    ENDDO

    !------------- Z_dir ---------------
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J)
    DO K = KST, NK-1
      DO J = JST, NJ-1
        DO I = IST, NI-1
          TEMPS(I,J,K) = CON(I,J,K) &
                         + AIP(I,J,K)*F(I+1,J,K) &
                         + AIM(I,J,K)*F(I-1,J,K) &
                         + AJP(I,J,K)*F(I,J+1,K) &
                         + AJM(I,J,K)*F(I,J-1,K)
        ENDDO
      ENDDO
    ENDDO

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,K,PT,QT,DENOM)
    DO J = JST, NJ-1
      DO I = IST, NI-1
        PT(KST-1) = 0.0_8
        QT(KST-1) = F(I,J,KST-1)  
        DO K = KST, NK-1
          DENOM = AP(I,J,K)-PT(K-1)*AKM(I,J,K)
          PT(K) = AKP(I,J,K)/(DENOM+1.D-30)
          QT(K) = (TEMPS(I,J,K)+AKM(I,J,K)*QT(K-1))/(DENOM+1.D-30)
        ENDDO
    
        DO K = NK-1, KST, -1
          F(I,J,K) = F(I,J,K+1)*PT(K)+QT(K)
        ENDDO
      ENDDO
    ENDDO

  ENDDO ! NTIMES   		
END SUBROUTINE SOLVE
!**********************************************************************
SUBROUTINE F0_F
USE VAR
USE USER
IMPLICIT NONE
  U0 = U; V0 = V; W0 = W; T0 = T; BU0 = BU; BV0 = BV; BW0 = BW
END SUBROUTINE F0_F
