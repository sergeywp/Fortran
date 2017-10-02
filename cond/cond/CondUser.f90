MODULE VAR

INTEGER, PARAMETER :: NI=99,NJ=99,NK=99,NIJ=99
REAL*8, PARAMETER :: pi=3.1415926535897932384626433832795

LOGICAL LSTOP

INTEGER L1,L2,L3,M1,M2,M3,N1,N2,N3,LAST,ITER,ISHOW,IST,JST,KST,ISTF,JSTF,KSTF
INTEGER :: ITR=1, IT=1,CLOCK1,CLOCK2

REAL*8 PT(NIJ),QT(NIJ),TEMP(NI,NJ,NK)

REAL*8 XL,YL,ZL,SMAX,SSUM,TIME,DT,DIFF,RELAX

REAL*8 X(NI),XU(NI),XDIF(NI),XCV(NI), &
       Y(NJ),YV(NJ),YDIF(NJ),YCV(NJ), &    
       Z(NK),ZW(NK),ZDIF(NK),ZCV(NK), &      
       VCV(NI,NJ,NK),RHO(NI,NJ,NK),GAMI(NI,NJ,NK),GAMJ(NI,NJ,NK),GAMK(NI,NJ,NK)

REAL*8 F(NI,NJ,NK),CON(NI,NJ,NK),AIP(NI,NJ,NK),AIM(NI,NJ,NK), &
       AJP(NI,NJ,NK),AJM(NI,NJ,NK),AKP(NI,NJ,NK),AKM(NI,NJ,NK),	&
       AP(NI,NJ,NK),AP0(NI,NJ,NK), SC(NI,NJ,NK),SP(NI,NJ,NK), &
       F1(NI,NJ,NK)

CHARACTER(10) DAT,BEG_TIM,END_TIM,ZON

INTEGER T_BEG,T_END, T_PARALLEL,T_UNPARALLEL

CHARACTER(LEN=40) :: NAME_FILE

DATA LSTOP /.FALSE./
DATA LAST,TIME,ITER/5,0.,0/
                                                                           
END MODULE VAR

!=====================================================

MODULE USER
CONTAINS
             
SUBROUTINE GRID
USE VAR

  LAST=800
  ISHOW=50

  XL=2.*pi
  YL=2.*pi
  ZL=2.*pi

  L1=99
  M1=99
  N1=99

  RELAX = 1.8

  DT=1.E+30


  CALL UGRID

END SUBROUTINE GRID
!=====================================================

SUBROUTINE START

USE VAR

   DO I=1,L1
    DO J=1,M1
     DO K=1,N1
       F(I,J,K)=0.
     ENDDO
    ENDDO
   ENDDO

END SUBROUTINE START
!=====================================================

SUBROUTINE DENSE
USE VAR

END SUBROUTINE DENSE
!=====================================================

SUBROUTINE BOUND
USE VAR

END SUBROUTINE BOUND
!=====================================================

SUBROUTINE OUTPUT
USE VAR

  IF(MOD(ITER,ISHOW).EQ.0) THEN

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(J,I)
   DO K=1,N1
    DO J=1,M1
     DO I=1,L1
       F1(I,J,K)=DSIN(X(I))*DSIN(Y(J))*DSIN(Z(K))
     ENDDO
    ENDDO
   ENDDO

DELTA=0.

 DO K=1,N1
  DO J=1,M1
   DO I=1,L1

     DELTA =DMAX1(DELTA, DABS(F(I,J,K)-F1(I,J,K))) 
   
   ENDDO
  ENDDO
 ENDDO


      WRITE(*,502)
      WRITE(1,502)
 502  FORMAT(5X,'ITER',8X,'TIME',8X,'F(5,5,5)',10X,'F1(5,5,5)',8X,'F(10,10,10)',&
             8X,'F1(10,10,10)',8X,'DELTA')

      WRITE(*,504)ITER,TIME,F(5,5,5),F1(5,5,5),F(10,10,10),F1(10,10,10),DELTA
      WRITE(1,504)ITER,TIME,F(5,5,5),F1(5,5,5),F(10,10,10),F1(10,10,10),DELTA

 504  FORMAT (I8,1P6E16.4)

  END IF

 
      	IF (ITER.NE.LAST) RETURN


    	OPEN(UNIT=3,FILE='ALL99_800.DAT',STATUS='UNKNOWN')
      	WRITE(3,*)'VARIABLES = "X", "Y","Z", "T","Tan","DT"'
	WRITE(3,*)'ZONE I=99, J=99, K=99, F=POINT'
                
        DO J=1,M1
         DO I=1,L1
	  DO K=1,N1 

            WRITE(3,'(1P6E15.5)') X(I),Y(J),Z(K),F(I,J,K),F1(I,J,K), &
                                  DABS(F(I,J,K)-F1(I,J,K))

	  ENDDO
         ENDDO           
        ENDDO
      
      	ENDFILE 3
      	CLOSE(3)

END SUBROUTINE OUTPUT

!=====================================================

SUBROUTINE GAMSOR
USE VAR

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(J,I)
  DO K=1,N1
    DO J=1,M1
      DO I=1,L1
              GAMI(I,J,K)=1.
              GAMJ(I,J,K)=1.
              GAMK(I,J,K)=1.
       ENDDO
    ENDDO 
  ENDDO

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(J,I)
  DO K=2,N2
    DO J=2,M2
      DO I=2,L2     
              SC(I,J,K)=3.*(DCOS(XU(I))-DCOS(XU(I+1)))* &
                           (DCOS(YV(J))-DCOS(YV(J+1)))* &   
                           (DCOS(ZW(K))-DCOS(ZW(K+1)))

       ENDDO
    ENDDO 
  ENDDO

                   
END SUBROUTINE GAMSOR

END MODULE USER