PROGRAM MAIN
USE VAR
USE USER

	CALL SYSTEM_CLOCK(CLOCK1)

	OPEN(UNIT=1,FILE='Q.OUT',STATUS='UNKNOWN')
	CALL GRID
	CALL SETUP1
	CALL START
10	CALL DENSE
	CALL BOUND
	CALL OUTPUT
	
	IF (LSTOP) THEN
	    CALL SYSTEM_CLOCK(CLOCK2)
	    WRITE (*,*),(CLOCK2 - CLOCK1)/10000.
	    WRITE (1,*),(CLOCK2 - CLOCK1)/10000.
	    STOP
	ENDIF

	CALL SETUP2
	GO TO 10
END

!=====================================================

SUBROUTINE UGRID
USE VAR

   XU(2) = 0.
   DX = XL/DBLE(L1-2)
   DO I = 3,L1
     XU(I) = XU(I-1)+DX
   END DO

   YV(2) = 0.
   DY = YL/DBLE(M1-2)
   DO J = 3,M1
     YV(J) = YV(J-1)+DY
   ENDDO

   ZW(2) = 0.
   DZ = ZL/DBLE(N1-2)
   DO K = 3,N1
     ZW(K) = ZW(K-1)+DZ
   ENDDO

END SUBROUTINE UGRID
!=====================================================

SUBROUTINE SETUP1
USE VAR

      	L2 = L1-1
      	L3 = L2-1
      	M2 = M1-1
      	M3 = M2-1
		N2 = N1-1
      	N3 = N2-1

      	X(1) = XU(2)
      	 DO I = 2,L2
       		X(I) = 0.5*(XU(I+1)+XU(I))
      	 ENDDO
      	X(L1) = XU(L1)

      	Y(1) = YV(2)
      	 DO J = 2,M2
       		Y(J) = 0.5*(YV(J+1)+YV(J))
      	 ENDDO 
      	Y(M1) = YV(M1)

		Z(1) = ZW(2)
      	 DO K = 2,N2
       		Z(K) = 0.5*(ZW(K+1)+ZW(K))
      	 ENDDO 
      	Z(N1) = ZW(N1)		

      	DO I = 2,L1
       		XDIF(I) = X(I)-X(I-1)
      	ENDDO
      	DO I = 2,L2
       		XCV(I) = XU(I+1)-XU(I)
      	ENDDO

      	DO J = 2,M1
       		YDIF(J) = Y(J)-Y(J-1)
     	ENDDO
      	DO J = 2,M2
       		YCV(J) = YV(J+1)-YV(J)
      	ENDDO

	DO K = 2,N1
       		ZDIF(K) = Z(K)-Z(K-1)
     	ENDDO
      	DO K = 2,N2
       		ZCV(K) = ZW(K+1)-ZW(K)
      	ENDDO
     
       	DO K = 2,N2
	 DO J = 2,M2
          DO I = 2,L2
             VCV(I,J,K) = XCV(I)*YCV(J)*ZCV(K) 
          ENDDO
	 ENDDO
       	ENDDO

!-----  SC,SP,RHO ARRAY ARE INITIALIZED HERE --------

      	DO K = 1,N1
           DO J = 1, M1
       	      DO I = 1, L1
          		SC(I,J,K) = 0.
          		SP(I,J,K) = 0.
          		RHO(I,J,K) = 1.
       	      ENDDO
      	   ENDDO
	ENDDO

END SUBROUTINE SETUP1

!=====================================================

SUBROUTINE SETUP2
USE VAR
USE USER

    IST = 2
    JST = 2
    KST = 2

      	CALL RESET
      	CALL GAMSOR

!---------- AIM, AIP, AJM, AJP -------------

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(J,I)

	DO K = 2, N2 
	 DO J = 2, M2
	  DO I = 2, L2

	      AIM(I,J,K) = GAMI(I,J,K)*YCV(J)*ZCV(K)/XDIF(I)            
	      AIP(I,J,K) = GAMI(I+1,J,K)*YCV(J)*ZCV(K)/XDIF(I+1)  
              AJM(I,J,K) = GAMJ(I,J,K)*XCV(I)*ZCV(K)/YDIF(J)
              AJP(I,J,K) = GAMJ(I,J+1,K)*XCV(I)*ZCV(K)/YDIF(J+1)
              AKM(I,J,K) = GAMK(I,J,K)*XCV(I)*YCV(J)/ZDIF(K)
              AKP(I,J,K) = GAMK(I,J,K+1)*XCV(I)*YCV(J)/ZDIF(K+1)
	      AP0(I,J,K) = RHO(I,J,K)*VCV(I,J,K)/DT
   
	  ENDDO
	 ENDDO
	ENDDO

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(J,I)   

	  DO K = 2, N2  
	    DO J = 2, M2
              DO I = 2, L2

                 CON(I,J,K) = SC(I,J,K)+AP0(I,J,K)*F(I,J,K)           		
                 AP(I,J,K) = -SP(I,J,K)+AP0(I,J,K)+AIM(I,J,K)+AIP(I,J,K)+AJM(I,J,K)+AJP(I,J,K)+AKM(I,J,K)+AKP(I,J,K)
     
              ENDDO
	    ENDDO
	  ENDDO


!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J)
      DO K = 2, N2
        DO J = 2, M2
          DO I = 2, L2


           AP(I,J,K) = AP(I,J,K)/RELAX
           CON(I,J,K) = CON(I,J,K)+(1.0-RELAX)*AP(I,J,K)*F(I,J,K)

     
          ENDDO
	ENDDO
      ENDDO

CALL SOLVE

 TIME = TIME+DT
 ITER = ITER+1
 IF (ITER.GE.LAST) LSTOP = .TRUE.

END SUBROUTINE SETUP2

!=====================================================

SUBROUTINE SOLVE   !SOR
USE VAR

                       
      DO K = 2, N2
        DO J = 2, M2
          DO I = 2, L2

            F(I,J,K) =  (AIP(I,J,K)*F(I+1,J,K)+AIM(I,J,K)*F(I-1,J,K)+ &
                         AJP(I,J,K)*F(I,J+1,K)+AJM(I,J,K)*F(I,J-1,K)+ &
                         AKP(I,J,K)*F(I,J,K+1)+AKM(I,J,K)*F(I,J,K-1)+ &
                         CON(I,J,K))/AP(I,J,K)
                        
	  ENDDO
        ENDDO
      ENDDO


END SUBROUTINE SOLVE

!-------------------------------

SUBROUTINE RESET
USE VAR
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J)

	DO K = 2, N2
       	    DO J = 2, M2
        	DO I = 2, L2
         		CON(I,J,K) = 0.
         		AP(I,J,K) = 0.
        	ENDDO
	    ENDDO
        ENDDO
END SUBROUTINE RESET
