SUBROUTINE PISO
USE VAR
USE USER
IMPLICIT NONE

  U1 = U; V1 = V
  CALL GAMSOR_U; CALL COF_U;  CALL SOLVE(3, 2, 2, NTIMES_UVW, U1)
  CALL GAMSOR_V; CALL COF_V;  CALL SOLVE(2, 3, 2, NTIMES_UVW, V1)
  CALL GAMSOR_W; CALL COF_W;  CALL SOLVE(2, 2, 3, NTIMES_UVW, W)
  U = U1; V = V1
  !--------------------------------------------------------------------------------------------------
  CALL PRESS; CALL SOLVE(2, 2, 2, NTIMES_P, P)
  !--------------------------------------------------------------------------------------------------
  CALL CORRECT_UVW
  CALL GAMSOR_U; CALL COF_U
  CALL GAMSOR_V; CALL COF_V
  CALL GAMSOR_W; CALL COF_W        
  !--------------------------------------------------------------------------------------------------
  CALL PRESS; CALL SOLVE(2, 2, 2, NTIMES_P, P)
  !--------------------------------------------------------------------------------------------------
  CALL CORRECT_UVW
  CALL COMPUTE_SMAX

END SUBROUTINE PISO

!***********************************************************************************
SUBROUTINE COF_U
USE VAR
USE USER
IMPLICIT NONE

INTEGER I, J, K
REAL(8) DIFF, FLOW, AP0

  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J,DIFF,FLOW)
  DO K = 2, N2
    DO J = 2, M2	
      !---------------- AIM, AIP,..., APK  FOR THE U EQUATION (I = 4,...,L3) ------------		
      DO I = 4, L3
        DIFF = GAM(I-1,J,K)*GAM(I,J,K)*ARI(J,K)/( GAM(I-1,J,K)*(XU(I)-X(I-1)) + GAM(I,J,K)*(X(I-1)-XU(I-1)) )
        FLOW = dens*0.5_8*(U(I-1,J,K) + U(I,J,K))*ARI(J,K)
        AIM(I,J,K) = ACOF1(DIFF,FLOW)
        DIFF = GAM(I,J,K)*GAM(I+1,J,K)*ARI(J,K)/( GAM(I,J,K)*(XU(I+1)-X(I)) + GAM(I+1,J,K)*(X(I)-XU(I)) )                
        FLOW = dens*0.5_8*(U(I,J,K) + U(I+1,J,K))*ARI(J,K)
        AIP(I,J,K) = ACOF1(DIFF,-FLOW)

        DIFF = GAM(I,J-1,K)*GAM(I,J,K)*ZCV(K)*XDIF(I)/( GAM(I,J-1,K)*(Y(J)-YV(J)) + GAM(I,J,K)*(YV(J)-Y(J-1)) )
        FLOW = dens*(V(I-1,J,K)*(XU(I) - X(I-1)) + V(I,J,K)*(X(I) - XU(I)))*ZCV(K)
        AJM(I,J,K) = ACOF1(DIFF,FLOW)
        DIFF = GAM(I,J,K)*GAM(I,J+1,K)*ZCV(K)*XDIF(I)/( GAM(I,J,K)*(Y(J+1)-YV(J+1)) + GAM(I,J+1,K)*(YV(J+1)-Y(J)) )
        FLOW = dens*(V(I-1,J+1,K)*(XU(I) - X(I-1)) + V(I,J+1,K)*(X(I) - XU(I)))*ZCV(K)
        AJP(I,J,K) = ACOF1(DIFF,-FLOW)

        DIFF = GAM(I,J,K-1)*GAM(I,J,K)*YCV(J)*XDIF(I)/( GAM(I,J,K-1)*(Z(K)-ZW(K)) + GAM(I,J,K)*(ZW(K)-Z(K-1)) )
        FLOW = dens*(W(I-1,J,K)*(XU(I) - X(I-1)) + W(I,J,K)*(X(I) - XU(I)))*YCV(J)
        AKM(I,J,K) = ACOF1(DIFF,FLOW)
        DIFF = GAM(I,J,K)*GAM(I,J,K+1)*YCV(J)*XDIF(I)/( GAM(I,J,K)*(Z(K+1)-ZW(K+1)) + GAM(I,J,K+1)*(ZW(K+1)-Z(K)) )
        FLOW = dens*(W(I-1,J,K+1)*(XU(I) - X(I-1)) + W(I,J,K+1)*(X(I) - XU(I)))*YCV(J)
        AKP(I,J,K) = ACOF1(DIFF,-FLOW)
      ENDDO
      !---------------- AIM, AIP,..., APK  FOR THE U EQUATION (I = 3) ----------------
      DIFF = GAM(2,J,K)*GAM(3,J,K)*ARI(J,K)/( GAM(2,J,K)*(XU(3)-X(2)) + GAM(3,J,K)*(X(2)-XU(2)) )
      FLOW = dens*U(2,J,K)*ARI(J,K)
      AIM(3,J,K) = ACOF1(DIFF,FLOW)
      DIFF = GAM(3,J,K)*GAM(4,J,K)*ARI(J,K)/( GAM(3,J,K)*(XU(4)-X(3)) + GAM(4,J,K)*(X(3)-XU(3)) )
      FLOW = dens*0.5_8*(U(3,J,K) + U(4,J,K))*ARI(J,K)
      AIP(3,J,K) = ACOF1(DIFF,-FLOW)

      DIFF = GAM(3,J-1,K)*GAM(3,J,K)*ZCV(K)*(X(3)-X(1))/( GAM(3,J-1,K)*(Y(J)-YV(J)) + GAM(3,J,K)*(YV(J)-Y(J-1)) )
      FLOW = dens*(V(2,J,K)*XCV(2) + V(3,J,K)*(X(3) - XU(3)))*ZCV(K)
      AJM(3,J,K) = ACOF1(DIFF,FLOW)
      DIFF = GAM(3,J,K)*GAM(3,J+1,K)*ZCV(K)*(X(3)-X(1))/( GAM(3,J,K)*(Y(J+1)-YV(J+1)) + GAM(3,J+1,K)*(YV(J+1)-Y(J)) )
      FLOW = dens*(V(2,J+1,K)*XCV(2) + V(3,J+1,K)*(X(3) - XU(3)))*ZCV(K)
      AJP(3,J,K) = ACOF1(DIFF,-FLOW)

      DIFF = GAM(3,J,K-1)*GAM(3,J,K)*YCV(J)*(X(3)-X(1))/( GAM(3,J,K-1)*(Z(K)-ZW(K)) + GAM(3,J,K)*(ZW(K)-Z(K-1)) )
      FLOW = dens*(W(2,J,K)*XCV(2) + W(3,J,K)*(X(3) - XU(3)))*YCV(J)
      AKM(3,J,K) = ACOF1(DIFF,FLOW)
      DIFF = GAM(3,J,K)*GAM(3,J,K+1)*YCV(J)*(X(3)-X(1))/( GAM(3,J,K)*(Z(K+1)-ZW(K+1)) + GAM(3,J,K+1)*(ZW(K+1)-Z(K)) )
      FLOW = dens*(W(2,J,K+1)*XCV(2) + W(3,J,K+1)*(X(3) - XU(3)))*YCV(J)
      AKP(3,J,K) = ACOF1(DIFF,-FLOW)
      !---------------- AIM, AIP,..., APK  FOR THE U EQUATION (I = L2) ----------------
      DIFF = GAM(L2-1,J,K)*GAM(L2,J,K)*ARI(J,K)/( GAM(L2-1,J,K)*(XU(L2)-X(L2-1)) + GAM(L2,J,K)*(X(L2-1)-XU(L2-1)) )
      FLOW = dens*0.5_8*(U(L3,J,K) + U(L2,J,K))*ARI(J,K)
      AIM(L2,J,K) = ACOF1(DIFF,FLOW)
      DIFF = GAM(L2,J,K)*GAM(L1,J,K)*ARI(J,K)/( GAM(L2,J,K)*(XU(L1)-X(L2)) + GAM(L1,J,K)*(X(L2)-XU(L2)) )
      FLOW = dens*U(L1,J,K)*ARI(J,K)
      AIP(L2,J,K) = ACOF1(DIFF,-FLOW)

      DIFF = GAM(L2,J-1,K)*GAM(L2,J,K)*ZCV(K)*(X(L1)-X(L1-2))/( GAM(L2,J-1,K)*(Y(J)-YV(J)) + GAM(L2,J,K)*(YV(J)-Y(J-1)) )
      FLOW = dens*(V(L3,J,K)*(XU(L2) - X(L3)) + V(L2,J,K)*XCV(L2))*ZCV(K)
      AJM(L2,J,K) = ACOF1(DIFF,FLOW)
      DIFF = GAM(L2,J,K)*GAM(L2,J+1,K)*ZCV(K)*(X(L1)-X(L1-2))/( GAM(L2,J,K)*(Y(J+1)-YV(J+1)) + GAM(L2,J+1,K)*(YV(J+1)-Y(J)) )
      FLOW = dens*(V(L3,J+1,K)*(XU(L2) - X(L3)) + V(L2,J+1,K)*XCV(L2))*ZCV(K)
      AJP(L2,J,K) = ACOF1(DIFF,-FLOW)

      DIFF = GAM(L2,J,K-1)*GAM(L2,J,K)*YCV(J)*(X(L1)-X(L1-2))/( GAM(L2,J,K-1)*(Z(K)-ZW(K)) + GAM(L2,J,K)*(ZW(K)-Z(K-1)) )
      FLOW = dens*(W(L3,J,K)*(XU(L2) - X(L3)) + W(L2,J,K)*XCV(L2))*YCV(J)
      AKM(L2,J,K) = ACOF1(DIFF,FLOW)
      DIFF = GAM(L2,J,K)*GAM(L2,J,K+1)*YCV(J)*(X(L1)-X(L1-2))/( GAM(L2,J,K)*(Z(K+1)-ZW(K+1)) + GAM(L2,J,K+1)*(ZW(K+1)-Z(K)) )
      FLOW = dens*(W(L3,J,K+1)*(XU(L2) - X(L3)) + W(L2,J,K+1)*XCV(L2))*YCV(J)
      AKP(L2,J,K) = ACOF1(DIFF,-FLOW)	
    ENDDO !J
  ENDDO !K

  !--------------------------- Boundary conditions of the second kind on z = 0 ------------------   	
  !!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I)
  !DO J = 2, M2
  !  DO I = 3, L2	   
  !    AKM(I,J,2) = 0.0_8               
  !  ENDDO
  !ENDDO    
  !-------------------- Sc, Sp -----------------------------------------------------
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J,AP0)
  DO K = 2, N2
    DO J = 2, M2
      DO I = 3, L2
        AP0 = dens*VCVU(I,J,K)/DT
        CON(I,J,K) = CON(I,J,K) + AP0*U0(I,J,K)
        AP(I,J,K) = AP0 +  &
        AIM(I,J,K) + AIP(I,J,K) + &
        AJM(I,J,K) + AJP(I,J,K) + &
        AKM(I,J,K) + AKP(I,J,K)
        DU(I,J,K) = ARI(J,K)/AP(I,J,K)
      ENDDO

      DO I = 2, L2
        UHAT(I,J,K) = &
        (AIP(I,J,K)*U(I+1,J,K)+AIM(I,J,K)*U(I-1,J,K)+ &
        AJP(I,J,K)*U(I,J+1,K)+AJM(I,J,K)*U(I,J-1,K)+ &
        AKP(I,J,K)*U(I,J,K+1)+AKM(I,J,K)*U(I,J,K-1)+ &
        CON(I,J,K))/AP(I,J,K)
      ENDDO

      DO I = 3, L2
        CON(I,J,K) = CON(I,J,K) + ARI(J,K)*(P(I-1,J,K)-P(I,J,K))
      ENDDO
    ENDDO !J
  ENDDO !K
END SUBROUTINE COF_U
!***********************************************************************************

SUBROUTINE COF_V
USE VAR
USE USER
IMPLICIT NONE
INTEGER I, J, K
REAL(8) DIFF, FLOW, AP0
           
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J,DIFF,FLOW)
  DO K = 2, N2
  !---------------- AIM, AIP,..., APK  FOR THE V EQUATION (J = 4,...,M3) ------------
    DO J = 4, M3
      DO I = 2, L2
        DIFF = GAM(I-1,J,K)*GAM(I,J,K)*ZCV(K)*YDIF(J)/( GAM(I-1,J,K)*(X(I)-XU(I)) + GAM(I,J,K)*(XU(I)-X(I-1)) )	
        FLOW = dens*(U(I,J-1,K)*(YV(J) - Y(J-1)) + U(I,J,K)*(Y(J) - YV(J)))*ZCV(K)
        AIM(I,J,K) = ACOF1(DIFF,FLOW)
        DIFF = GAM(I,J,K)*GAM(I+1,J,K)*ZCV(K)*YDIF(J)/( GAM(I,J,K)*(X(I+1)-XU(I+1)) + GAM(I+1,J,K)*(XU(I+1)-X(I)) )		
        FLOW = dens*(U(I+1,J-1,K)*(YV(J) - Y(J-1)) + U(I+1,J,K)*(Y(J) - YV(J)))*ZCV(K)
        AIP(I,J,K) = ACOF1(DIFF,-FLOW)

        DIFF = GAM(I,J-1,K)*GAM(I,J,K)*ARJ(I,K)/( GAM(I,J-1,K)*(YV(J)-Y(J-1)) + GAM(I,J,K)*(Y(J-1)-YV(J-1)) )	
        FLOW = dens*0.5_8*(V(I,J,K) + V(I,J-1,K))*ARJ(I,K)
        AJM(I,J,K) = ACOF1(DIFF,FLOW)
        DIFF = GAM(I,J,K)*GAM(I,J+1,K)*ARJ(I,K)/( GAM(I,J,K)*(YV(J+1)-Y(J)) + GAM(I,J+1,K)*(Y(J)-YV(J)) )	
        FLOW = dens*0.5_8*(V(I,J,K) + V(I,J+1,K))*ARJ(I,K)
        AJP(I,J,K) = ACOF1(DIFF,-FLOW)

        DIFF = GAM(I,J,K-1)*GAM(I,J,K)*XCV(I)*YDIF(J)/( GAM(I,J,K-1)*(Z(K)-ZW(K)) + GAM(I,J,K)*(ZW(K)-Z(K-1)) )
        FLOW = dens*(W(I,J-1,K)*(YV(J) - Y(J-1)) + W(I,J,K)*(Y(J) - YV(J)))*XCV(I)
        AKM(I,J,K) = ACOF1(DIFF,FLOW)
        DIFF = GAM(I,J,K)*GAM(I,J,K+1)*XCV(I)*YDIF(J)/( GAM(I,J,K)*(Z(K+1)-ZW(K+1)) + GAM(I,J,K+1)*(ZW(K+1)-Z(K)) )
        FLOW = dens*(W(I,J-1,K+1)*(YV(J) - Y(J-1)) + W(I,J,K+1)*(Y(J) - YV(J)))*XCV(I)
        AKP(I,J,K) = ACOF1(DIFF,-FLOW)
      ENDDO !I
    ENDDO !J
 
    DO I = 2, L2
      !---------------- AIM, AIP,..., APK  FOR THE V EQUATION (J = 3 ) ------------
      DIFF = GAM(I-1,3,K)*GAM(I,3,K)*ZCV(K)*(Y(3)-Y(1))/( GAM(I-1,3,K)*(X(I)-XU(I)) + GAM(I,3,K)*(XU(I)-X(I-1)) )
      FLOW = dens*(U(I,2,K)*YCV(2) + U(I,3,K)*(Y(3) - YV(3)))*ZCV(K)
      AIM(I,3,K) = ACOF1(DIFF,FLOW)
      DIFF = GAM(I,3,K)*GAM(I+1,3,K)*ZCV(K)*(Y(3)-Y(1))/( GAM(I,3,K)*(X(I+1)-XU(I+1)) + GAM(I+1,3,K)*(XU(I+1)-X(I)) )
      FLOW = dens*(U(I+1,2,K)*YCV(2) + U(I+1,3,K)*(Y(3) - YV(3)))*ZCV(K)
      AIP(I,3,K) = ACOF1(DIFF,-FLOW)

      DIFF = GAM(I,2,K)*GAM(I,3,K)*ARJ(I,K)/( GAM(I,2,K)*(YV(3)-Y(2)) + GAM(I,3,K)*(Y(2)-YV(2)) )
      FLOW = dens*V(I,2,K)*ARJ(I,K)
      AJM(I,3,K) = ACOF1(DIFF,FLOW)
      DIFF = GAM(I,3,K)*GAM(I,4,K)*ARJ(I,K)/( GAM(I,3,K)*(YV(4)-Y(3)) + GAM(I,4,K)*(Y(3)-YV(3)) )			
      FLOW = dens*0.5_8*(V(I,3,K) + V(I,4,K))*ARJ(I,K)
      AJP(I,3,K) = ACOF1(DIFF,-FLOW)

      DIFF = GAM(I,3,K-1)*GAM(I,3,K)*XCV(I)*(Y(3)-Y(1))/( GAM(I,3,K-1)*(Z(K)-ZW(K)) + GAM(I,3,K)*(ZW(K)-Z(K-1)) )  
      FLOW = dens*(W(I,2,K)*YCV(2) + W(I,3,K)*(Y(3) - YV(3)))*XCV(I)
      AKM(I,3,K) = ACOF1(DIFF,FLOW)
      DIFF = GAM(I,3,K)*GAM(I,3,K+1)*XCV(I)*(Y(3)-Y(1))/( GAM(I,3,K)*(Z(K+1)-ZW(K+1)) + GAM(I,3,K+1)*(ZW(K+1)-Z(K)) )
      FLOW = dens*(W(I,2,K+1)*YCV(2) + W(I,3,K+1)*(Y(3) - YV(3)))*XCV(I)
      AKP(I,3,K) = ACOF1(DIFF,-FLOW)
      !---------------- AIM, AIP,..., APK  FOR THE V EQUATION (J = M3 ) ------------	
      DIFF = GAM(I-1,M2,K)*GAM(I,M2,K)*ZCV(K)*(Y(M1)-Y(M3))/( GAM(I-1,M2,K)*(X(I)-XU(I)) + GAM(I,M2,K)*(XU(I)-X(I-1)) )			
      FLOW = dens*(U(I,M3,K)*(YV(M2) - Y(M3)) + U(I,M2,K)*YCV(M2))*ZCV(K)
      AIM(I,M2,K) = ACOF1(DIFF,FLOW)
      DIFF = GAM(I,M2,K)*GAM(I+1,M2,K)*ZCV(K)*(Y(M1)-Y(M3))/( GAM(I,M2,K)*(X(I+1)-XU(I+1)) + GAM(I+1,M2,K)*(XU(I+1)-X(I)) )	
      FLOW = dens*(U(I+1,M3,K)*(YV(M2) - Y(M3)) + U(I+1,M2,K)*YCV(M2))*ZCV(K)
      AIP(I,M2,K) = ACOF1(DIFF,-FLOW)                                                                           

      DIFF = GAM(I,M2-1,K)*GAM(I,M2,K)*ARJ(I,K)/( GAM(I,M2-1,K)*(YV(M2)-Y(M2-1)) + GAM(I,M2,K)*(Y(M2-1)-YV(M2-1)) )
      FLOW = dens*0.5_8*(V(I,M2,K) + V(I,M3,K))*ARJ(I,K)
      AJM(I,M2,K) = ACOF1(DIFF,FLOW)
      DIFF = GAM(I,M2,K)*GAM(I,M1,K)*ARJ(I,K)/( GAM(I,M2,K)*(YV(M1)-Y(M2)) + GAM(I,M1,K)*(Y(M2)-YV(M2)) )	
      FLOW = dens*V(I,M1,K)*ARJ(I,K)
      AJP(I,M2,K) = ACOF1(DIFF,-FLOW)

      DIFF = GAM(I,M2,K-1)*GAM(I,M2,K)*XCV(I)*(Y(M1)-Y(M3))/( GAM(I,M2,K-1)*(Z(K)-ZW(K)) + GAM(I,M2,K)*(ZW(K)-Z(K-1)) )  
      FLOW = dens*(W(I,M3,K)*(YV(M2) - Y(M3)) + W(I,M2,K)*YCV(M2))*XCV(I)
      AKM(I,M2,K) = ACOF1(DIFF,FLOW)
      DIFF = GAM(I,M2,K)*GAM(I,M2,K+1)*XCV(I)*(Y(M1)-Y(M3))/( GAM(I,M2,K)*(Z(K+1)-ZW(K+1)) + GAM(I,M2,K+1)*(ZW(K+1)-Z(K)) )
      FLOW = dens*(W(I,M3,K+1)*(YV(M2) - Y(M3)) + W(I,M2,K+1)*YCV(M2))*XCV(I)
      AKP(I,M2,K) = ACOF1(DIFF,-FLOW)
    ENDDO !I
  ENDDO !K

  !--------------------------- Boundary conditions of the second kind on z = 0 ------------------   	
  !!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I)
  !DO J = 3, M2
  !  DO I = 2, L2	   
  !    AKM(I,J,2) = 0.0_8               
  !  ENDDO
  !ENDDO
    
  !-------------------- Sc, Sp -----------------------------------------------------
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J,AP0)
  DO K = 2, N2	
    DO J = 3, M2
      DO I = 2, L2
        AP0 = dens*VCVV(I,J,K)/DT
        CON(I,J,K) = CON(I,J,K) + AP0*V0(I,J,K)
        AP(I,J,K) = AP0 + &
        AIM(I,J,K) + AIP(I,J,K) + &
        AJM(I,J,K) + AJP(I,J,K) + &
        AKM(I,J,K) + AKP(I,J,K)			
        DV(I,J,K) = ARJ(I,K)/AP(I,J,K)
      ENDDO
    ENDDO

    DO J = 2, M2
      DO I = 2, L2
        VHAT(I,J,K) = &
        (AIP(I,J,K)*V(I+1,J,K)+AIM(I,J,K)*V(I-1,J,K)+ &
        AJP(I,J,K)*V(I,J+1,K)+AJM(I,J,K)*V(I,J-1,K)+ &
        AKP(I,J,K)*V(I,J,K+1)+AKM(I,J,K)*V(I,J,K-1)+ &
        CON(I,J,K))/AP(I,J,K)
      ENDDO
    ENDDO

    DO J = 3, M2
      DO I = 2, L2
        CON(I,J,K) = CON(I,J,K) + ARJ(I,K)*(P(I,J-1,K)-P(I,J,K))
      ENDDO
    ENDDO

  ENDDO !K

END SUBROUTINE COF_V

!***********************************************************************************
SUBROUTINE COF_W
USE VAR
USE USER
IMPLICIT NONE
INTEGER I, J, K
REAL(8) DIFF, FLOW, AP0

  !---------------- AIM, AIP,..., APK  FOR THE W EQUATION (K = 4,...,N3) ------------
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,K,DIFF,FLOW)
  DO K = 4, N3
    DO J = 2, M2
      DO I = 2, L2
        DIFF = GAM(I-1,J,K)*GAM(I,J,K)*YCV(J)*ZDIF(K)/( GAM(I-1,J,K)*(X(I)-XU(I)) + GAM(I,J,K)*(XU(I)-X(I-1)) )                        
        FLOW = dens*(U(I,J,K-1)*(ZW(K) - Z(K-1)) + U(I,J,K)*(Z(K) - ZW(K)))*YCV(J)
        AIM(I,J,K) = ACOF1(DIFF,FLOW)
        DIFF = GAM(I,J,K)*GAM(I+1,J,K)*YCV(J)*ZDIF(K)/( GAM(I,J,K)*(X(I+1)-XU(I+1)) + GAM(I+1,J,K)*(XU(I+1)-X(I)) )    
        FLOW = dens*(U(I+1,J,K-1)*(ZW(K) - Z(K-1)) + U(I+1,J,K)*(Z(K) - ZW(K)))*YCV(J)
        AIP(I,J,K) = ACOF1(DIFF,-FLOW)

        DIFF = GAM(I,J-1,K)*GAM(I,J,K)*XCV(I)*ZDIF(K)/( GAM(I,J-1,K)*(Y(J)-YV(J)) + GAM(I,J,K)*(YV(J)-Y(J-1)) )	    
        FLOW = dens*(V(I,J,K-1)*(ZW(K) - Z(K-1)) + V(I,J,K)*(Z(K) - ZW(K)))*XCV(J)
        AJM(I,J,K) = ACOF1(DIFF,FLOW)
        DIFF = GAM(I,J,K)*GAM(I,J+1,K)*XCV(I)*ZDIF(K)/( GAM(I,J,K)*(Y(J+1)-YV(J+1)) + GAM(I,J+1,K)*(YV(J+1)-Y(J)) )
        FLOW = dens*(V(I,J+1,K-1)*(ZW(K) - Z(K-1)) + V(I,J+1,K)*(Z(K) - ZW(K)))*XCV(J)
        AJP(I,J,K) = ACOF1(DIFF,-FLOW)

        DIFF = GAM(I,J,K-1)*GAM(I,J,K)*ARK(I,J)/( GAM(I,J,K-1)*(ZW(K)-Z(K-1)) + GAM(I,J,K)*(Z(K-1)-ZW(K-1)) )
        FLOW = dens*0.5_8*(W(I,J,K) + W(I,J,K-1))*ARK(I,J)
        AKM(I,J,K) = ACOF1(DIFF,FLOW)
        DIFF = GAM(I,J,K)*GAM(I,J,K+1)*ARK(I,J)/( GAM(I,J,K)*(ZW(K+1)-Z(K)) + GAM(I,J,K+1)*(Z(K)-ZW(K)) )
        FLOW = dens*0.5_8*(W(I,J,K) + W(I,J,K+1))*ARK(I,J)
        AKP(I,J,K) = ACOF1(DIFF,-FLOW)
      ENDDO
    ENDDO
  ENDDO

  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,DIFF,FLOW)
  DO J = 2, M2
    DO I = 2, L2
      !---------------- AIM, AIP,..., APK  FOR THE W EQUATION (K = 3) ------------
      DIFF = GAM(I-1,J,3)*GAM(I,J,3)*YCV(J)*(Z(3)-Z(1))/( GAM(I-1,J,3)*(X(I)-XU(I)) + GAM(I,J,3)*(XU(I)-X(I-1)) )                                   
      FLOW = dens*(U(I,J,2)*ZCV(2) + U(I,J,3)*(Z(3) - ZW(3)))*YCV(J)
      AIM(I,J,3) = ACOF1(DIFF,FLOW)
      DIFF = GAM(I,J,3)*GAM(I+1,J,3)*YCV(J)*(Z(3)-Z(1))/( GAM(I,J,3)*(X(I+1)-XU(I+1)) + GAM(I+1,J,3)*(XU(I+1)-X(I)) )
      FLOW = dens*(U(I+1,J,2)*ZCV(2) + U(I+1,J,3)*(Z(3) - ZW(3)))*YCV(J)
      AIP(I,J,3) = ACOF1(DIFF,-FLOW)

      DIFF = GAM(I,J-1,3)*GAM(I,J,3)*XCV(I)*(Z(3)-Z(1))/( GAM(I,J-1,3)*(Y(J)-YV(J)) + GAM(I,J,3)*(YV(J)-Y(J-1)) )	
      FLOW = dens*(V(I,J,2)*ZCV(2) + V(I,J,3)*(Z(3) - ZW(3)))*XCV(J)
      AJM(I,J,3) = ACOF1(DIFF,FLOW)
      DIFF = GAM(I,J,3)*GAM(I,J+1,3)*XCV(I)*(Z(3)-Z(1))/( GAM(I,J,3)*(Y(J+1)-YV(J+1)) + GAM(I,J+1,3)*(YV(J+1)-Y(J)) )
      FLOW = dens*(V(I,J+1,2)*ZCV(2) + V(I,J+1,3)*(Z(3) - ZW(3)))*XCV(J)
      AJP(I,J,3) = ACOF1(DIFF,-FLOW)

      DIFF = GAM(I,J,2)*GAM(I,J,3)*ARK(I,J)/( GAM(I,J,2)*(ZW(3)-Z(2)) + GAM(I,J,3)*(Z(2)-ZW(2)) )		    
      FLOW = dens*W(I,J,2)*ARK(I,J)
      AKM(I,J,3) = ACOF1(DIFF,FLOW)
      DIFF = GAM(I,J,3)*GAM(I,J,4)*ARK(I,J)/( GAM(I,J,3)*(ZW(4)-Z(3)) + GAM(I,J,4)*(Z(3)-ZW(3)) )
      FLOW = dens*0.5_8*(W(I,J,3) + W(I,J,4))*ARK(I,J)
      AKP(I,J,3) = ACOF1(DIFF,-FLOW)
      !---------------- AIM, AIP,..., APK  FOR THE W EQUATION (K = N2) ------------
      DIFF = GAM(I-1,J,N2)*GAM(I,J,N2)*YCV(J)*(Z(N1)-Z(N3))/( GAM(I-1,J,N2)*(X(I)-XU(I)) + GAM(I,J,N2)*(XU(I)-X(I-1)) )
      FLOW = dens*(U(I,J,N3)*(ZW(N2) - Z(N3)) + U(I,J,N2)*ZCV(N2))*YCV(J)
      AIM(I,J,N2) = ACOF1(DIFF,FLOW)			    
      DIFF = GAM(I,J,N2)*GAM(I+1,J,N2)*YCV(J)*(Z(N1)-Z(N3))/( GAM(I,J,N2)*(X(I+1)-XU(I+1)) + GAM(I+1,J,N2)*(XU(I+1)-X(I)) )
      FLOW = dens*(U(I+1,J,N3)*(ZW(N2) - Z(N3)) + U(I+1,J,N2)*ZCV(N2))*YCV(J)
      AIP(I,J,N2) = ACOF1(DIFF,-FLOW)


      DIFF = GAM(I,J-1,N2)*GAM(I,J,N2)*XCV(I)*(Z(N1)-Z(N3))/( GAM(I,J-1,N2)*(Y(J)-YV(J)) + GAM(I,J,N2)*(YV(J)-Y(J-1)) )
      FLOW = dens*(V(I,J,N3)*(ZW(N2) - Z(N3)) + V(I,J,N2)*ZCV(N2))*XCV(J)
      AJM(I,J,N2) = ACOF1(DIFF,FLOW)
      DIFF = GAM(I,J,N2)*GAM(I,J+1,N2)*XCV(I)*(Z(N1)-Z(N3))/( GAM(I,J,N2)*(Y(J+1)-YV(J+1)) + GAM(I,J+1,N2)*(YV(J+1)-Y(J)) )
      FLOW = dens*(V(I,J+1,N3)*(ZW(N2) - Z(N3)) + V(I,J+1,N2)*ZCV(N2))*XCV(J)
      AJP(I,J,N2) = ACOF1(DIFF,-FLOW)

      DIFF = GAM(I,J,N2-1)*GAM(I,J,N2)*ARK(I,J)/( GAM(I,J,N2-1)*(ZW(N2)-Z(N2-1)) + GAM(I,J,N2)*(Z(N2-1)-ZW(N2-1)) )
      FLOW = dens*0.5_8*(W(I,J,N2) + W(I,J,N3))*ARK(I,J)
      AKM(I,J,N2) = ACOF1(DIFF,FLOW)
      DIFF = GAM(I,J,N2)*GAM(I,J,N1)*ARK(I,J)/( GAM(I,J,N2)*(ZW(N1)-Z(N2)) + GAM(I,J,N1)*(Z(N2)-ZW(N2)) )                     
      FLOW = dens*W(I,J,N1)*XCV(I)*ARK(I,J)
      AKP(I,J,N2) = ACOF1(DIFF,-FLOW)
    ENDDO
  ENDDO
  !---------------- Sc, Sp ------------------

  !-----AP0, CON, AP FOR THE W EQUATION -----------------------------------------------
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J,AP0)
  DO K = 3, N2
    DO J = 2, M2
      DO I = 2, L2

        AP0 = dens*VCVW(I,J,K)/DT

        CON(I,J,K) = CON(I,J,K) + AP0*W0(I,J,K)

        AP(I,J,K) = AP0 + &
        AIM(I,J,K) + AIP(I,J,K) + &
        AJM(I,J,K) + AJP(I,J,K) + &
        AKM(I,J,K) + AKP(I,J,K)

        DW(I,J,K) = ARK(I,J)/AP(I,J,K)

      ENDDO
    ENDDO
  ENDDO


  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J)
  DO K = 2, N2
    DO J = 2, M2
      DO I = 2, L2
        WHAT(I,J,K)=(AIP(I,J,K)*W(I+1,J,K)+AIM(I,J,K)*W(I-1,J,K)+ &
        AJP(I,J,K)*W(I,J+1,K)+AJM(I,J,K)*W(I,J-1,K)+ &
        AKP(I,J,K)*W(I,J,K+1)+AKM(I,J,K)*W(I,J,K-1)+ &
        CON(I,J,K))/AP(I,J,K)
      ENDDO
    ENDDO
  ENDDO

  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J)
  DO K = 3, N2
    DO J = 2, M2
      DO I = 2, L2
        CON(I,J,K) = CON(I,J,K)+ARK(I,J)*(P(I,J,K-1)-P(I,J,K))
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE COF_W
!***********************************************************************************

SUBROUTINE PRESS
USE VAR
USE USER
IMPLICIT NONE
INTEGER I, J, K

  !--------------- AIM, AIP ------------------
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J)
  DO K = 2, N2
    DO J = 2, M2
      DO I = 3, L3

        AIM(I,J,K) = ARI(J,K)*DU(I,J,K)
        AIP(I,J,K) = ARI(J,K)*DU(I+1,J,K)
        CON(I,J,K) = ARI(J,K)*(UHAT(I,J,K) - UHAT(I+1,J,K))

      ENDDO          
    ENDDO
  ENDDO

  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(J)
  DO K = 2, N2
    DO J = 2, M2

      AIM(2,J,K) = 0.0_8
      AIP(2,J,K) = ARI(J,K)*DU(3,J,K)
      CON(2,J,K) = ARI(J,K)*(U(2,J,K) - UHAT(3,J,K))
      AIM(L2,J,K) = ARI(J,K)*DU(L2,J,K)
      AIP(L2,J,K) = 0.0_8
      CON(L2,J,K)=ARI(J,K)*(UHAT(L2,J,K) - U(L1,J,K))

    ENDDO
  ENDDO

  !-------------- AJM, AJP -------------------
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J)
  DO K = 2, N2
    DO J = 3, M3
      DO I = 2, L2

        AJM(I,J,K) = ARJ(I,K)*DV(I,J,K)
        AJP(I,J,K) = ARJ(I,K)*DV(I,J+1,K)
        CON(I,J,K) = CON(I,J,K) + ARJ(I,K)*(VHAT(I,J,K) - VHAT(I,J+1,K))

      ENDDO          
    ENDDO
  ENDDO

  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I)
  DO K = 2, N2
    DO I = 2, L2	

      AJM(I,2,K) = 0.0_8
      AJP(I,2,K) = ARJ(I,K)*DV(I,3,K)
      CON(I,2,K) = CON(I,2,K) + ARJ(I,K)*(V(I,2,K) - VHAT(I,3,K))
      AJM(I,M2,K) = ARJ(I,K)*DV(I,M2,K)
      AJP(I,M2,K) = 0.0_8
      CON(I,M2,K) = CON(I,M2,K) + ARJ(I,K)*(VHAT(I,M2,K) - V(I,M1,K))

    ENDDO
  ENDDO

  !--------------- AKM, AKP --------------------
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I)
  DO J = 2, M2
    DO I = 2, L2
      AKM(I,J,2) = 0.0_8
      AKP(I,J,2) = ARK(I,J)*DW(I,J,3)
      CON(I,J,2) = CON(I,J,2) + ARK(I,J)*(W(I,J,2) - WHAT(I,J,3))

      AKM(I,J,N2) = ARK(I,J)*DW(I,J,N2)
      AKP(I,J,N2) = 0.0_8
      CON(I,J,N2) = CON(I,J,N2) + ARK(I,J)*(WHAT(I,J,N2) - W(I,J,N1))
    ENDDO
  ENDDO

  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J)
  DO K = 3, N3	
    DO J = 2, M2
      DO I = 2, L2
        AKM(I,J,K) = ARK(I,J)*DW(I,J,K)
        AKP(I,J,K) = ARK(I,J)*DW(I,J,K+1)		
        CON(I,J,K) = CON(I,J,K) + ARK(I,J)*(WHAT(I,J,K) - WHAT(I,J,K+1))
      ENDDO          
    ENDDO
  ENDDO

  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J)
  DO K = 2, N2
    DO J = 2, M2
      DO I = 2, L2

        AP(I,J,K) = AIM(I,J,K)+AIP(I,J,K)+AJM(I,J,K)+AJP(I,J,K)+AKM(I,J,K)+AKP(I,J,K)

      ENDDO          
    ENDDO
  ENDDO

END SUBROUTINE PRESS

!***********************************************************************************
SUBROUTINE COMPUTE_SMAX
USE VAR
USE USER
IMPLICIT NONE
INTEGER I, J, K

  SMAX = 0.0_8
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J)
  DO K = 2, N2
    DO J = 2, M2
      DO I = 2, L2	

        CON(I,J,K)= DABS( ARI(J,K)*(U(I,J,K) - U(I+1,J,K))+ &
                          ARJ(I,K)*(V(I,J,K) - V(I,J+1,K))+ &
                          ARK(I,J)*(W(I,J,K) - W(I,J,K+1)) )/VCV(I,J,K)

      ENDDO          
    ENDDO
  ENDDO

  SMAX = MAXVAL(CON)

END SUBROUTINE COMPUTE_SMAX
!***********************************************************************************

SUBROUTINE CORRECT_UVW
USE VAR
USE USER
IMPLICIT NONE
INTEGER I, J, K
!------------ COME HERE TO CORRECT THE VELOCITIES ------------

  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J)
  DO K = 2, N2
    DO J = 2, M2
      DO I = 3, L2

        U(I,J,K) =  UHAT(I,J,K) + DU(I,J,K)*(P(I-1,J,K)-P(I,J,K))

      ENDDO          
    ENDDO
  ENDDO

  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J)
  DO K = 2, N2
    DO J = 3, M2
      DO I = 2, L2

        V(I,J,K) = VHAT(I,J,K) + DV(I,J,K)*(P(I,J-1,K)-P(I,J,K))

      ENDDO          
    ENDDO
  ENDDO

  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J)
  DO K = 3, N2
    DO J = 2, M2
      DO I = 2, L2

        W(I,J,K) = WHAT(I,J,K) + DW(I,J,K)*(P(I,J,K-1)-P(I,J,K))         

      ENDDO          
    ENDDO
  ENDDO

END SUBROUTINE CORRECT_UVW
