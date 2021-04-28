C
C-----------------------------------------------------------------------
C
C   Model for the plasmaspheric electron temperature based on Akebono
C   TED data:
C     I. Kutiev, K. Oyama, T. Abe, and P. Marinov, Plasamsphere 
C         electron temperature model based on Akebono data, Advances in 
C         Space Research 33, 975-979, 2004
C
C   This file includes the model subroutine TE3ARG and a short driver
C   program 
C
C-----------------------------------------------------------------------
C   Subroutine TE3ARG and main program
C    
C   Examples:
C      IN: TLOC, GMLAT, HEIGHT     OUT: TLOC, GMLAT, Height, Te 
C               14,-50,1000        14.00   -50.00  1000.00   5283.16
C               14,-70,5000        14.00   -70.00  5000.00   3935.66
C               14,-71,5000        14.00   -71.00  5000.00    3935.7 warning 1
C               14,70,5000         14.00    70.00  5000.00   3935.66
C               14,-80,0           14.00   -80.00     0.00    2764.0 warning 2
C               warning 1 indicates one input parameter out of range
C               warning 2 indicates two input parameters out of range
C-----------------------------------------------------------------------
C
C     2007/01/10
C     2005/07/17
C     2002/04/28
C   		Pencho G. Marinov                                                    
C   		Institute for Parallel Processing                           
C   		Bulgarian Academy of Science                                         
C   		E-mail: pencho@bas.bg                                           
C
C     1) The Polynomial Approximation
C
C     2) 3-dimensioned X_space
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
      PROGRAM PGMIK
C
      IMPLICIT REAL*4 (A-H,O-Z)
      IMPLICIT INTEGER*2 (I-N)
C
      COMMON /BLBAS/ XL(3),XR(3),SL(3),SR(3),NBAS(2,3),MBAS(3)
      COMMON /BLDAT/ BUFF(50)
      COMMON /BLWRK/ WA(216),WB(216),WC(6,6,6)
      COMMON /BLCFN/ Z(216),A(216),CFN(6,6,6),NCFN(20)
      COMMON /BLEPS/ ERS(20), NITER(20), O, E
      COMMON /BLIOT/ NINP1, NINP2, NOUT1, NOUT2, NTAB1, NTAB2
C
      LOGICAL IOCHECK
      CHARACTER*64 NAME(10), NAMP(10)
      integer*4 n,n0
c
      data NR, NW, npgm /1,2,4/
C
      DATA  NAME(1)/'AKEBONO.INP'/
      DATA  NAME(2)/'AKEBONO.OUT'/
      DATA  NAME(3)/'AKEBONO.INF'/
      DATA  NAME(4)/'AKEBONO.CFN'/
      DATA  NAMP(7)/'(1X,A)                  '/
      DATA  NAMP(8)/'(1X,3F9.2, F10.2  )     '/
      DATA  NAMP(9)/'(1X,3F9.2,F10.1,8H warning,I2)'/
C
         n=0
         n0=0
         write(*,*) 'LT [0.0 - 24.0], GMLAT [-70 to +70], HEIGHT ',
     &   '[1000-10,000 km] (# to end)'
C
  200 CONTINUE
         READ (* , *, ERR=201, END=201) TLOC, GMLAT, HEIGHT
         n=n+1
         call TE3ARG(TLOC, GMLAT, HEIGHT, TE, ierr)
         if (ierr.GT.0) then
            n0=n0+1
            WRITE( *,NAMP(9)) TLOC, GMLAT, HEIGHT, TE, ierr
         else
            WRITE(* ,NAMP(8)) TLOC, GMLAT, HEIGHT, TE
         endif
         GOTO 200
  201 CONTINUE
C
      write(*,*) n , ' - rows processed  '
      write(*,*) n0, ' - warnings occured '
C
 301  continue
      STOP '  '
C
      END
C
C
       SUBROUTINE TE3ARG(TLOC,GMLAT,HEIGHT,TE,ierr)
C--------------------------------------------------------------------------
C      Input parameters are: 
C          TLOC solar local time in hours (0-24) 
C          GMLAT magnetic latitude in degrees (-70 - +70) 
C          HEIGHT altitude in km (1000 - 10,000)
C          If input parameters is outside valid range it will be set to the
C          nearest range boundary, e.g., +80 to +70
C      Output parameters are: 
C          TE  electron temperature in K
C          ierr  # of cases of input parameter ouside valid range
C--------------------------------------------------------------------------
C
C
      IMPLICIT REAL*4 (A-H,O-Z)
      IMPLICIT INTEGER*2  (I-N)
C
      COMMON /BLBAS/ XL(3),XR(3),SL(3),SR(3),NBAS(2,3),MBAS(3)
      COMMON /BLDAT/ BUFF(50)
      COMMON /BLCFN/ Z(216),A(216),CFN(6,6,6),NCFN(20)
      COMMON /BLEPS/ ERS(20), NITER(20), O, E
C
      CHARACTER*64 FORMD,FORME
C
      FFL(V)=(2.0D0-V-3*V*V+2*V**3)*0.5D0
      FFR(V)=(      V+3*V*V-2*V**3)*0.5D0
C
      A(11) = TLOC
      A(12) = GMLAT
      A(13) = HEIGHT
C
C     Test if input data are in correct range.
C          if NOT then ierr>0
C
         ierr=0
         A(1)=A(11)
         A(2)=A(12)
         A(3)=A(13)
         IF (A(11).LT.(0.0)) THEN
            A(1)=0.0
            ierr=ierr+1
         ENDIF
         IF (A(11).GT.(24.)) THEN
            A(1)=24.0
            ierr=ierr+1
         ENDIF
         IF ( (A(12).LT.XL(2)).AND.(A(12).LT.O) ) A(2)=-A(12)
         IF (A( 2).GT.XR(2)) THEN
            A(2)=XR(2)
            ierr=ierr+1
         ENDIF
         IF (A(13).LT.XL(3)) THEN
            A(3)=XL(3)
            ierr=ierr+1
         ENDIF
         IF (A(13).GT.XR(3)) THEN
            A(3)=XR(3)
            ierr=ierr+1
         ENDIF
C
         DO J=2,3
            BUFF(J)=SL(J)+(A(J)-XL(J))/(XR(J)-XL(J))*(SR(J)-SL(J))
         ENDDO
         BUFF(1)=SL(1)
         R1 = R(BUFF)
         BUFF(1)=SR(1)
         R2 = R(BUFF)
         IF ((A(1)-9.0)*(16.0-A(1)).GE.O) THEN
            BUFF(4) = R2
            GOTO 202
         ENDIF
         IF ( (A(1).GE.(22.0)).OR.(A(1).LE.(4.0)) ) THEN
            BUFF(4) = R1
            GOTO 202
         ENDIF
         IF ((A(1)-4.0)*(9.0-A(1)).GE.O) THEN
            SS = (A(1)-4.0)/5.0
            BUFF(4) = R1*FFL(SS)+R2*FFR(SS)
            GOTO 202
         ENDIF
         IF ((A(1)-16.0)*(22.0-A(1)).GE.O) THEN
            SS = (A(1)-16.0)/6.0
            BUFF(4) = R2*FFL(SS)+R1*FFR(SS)
            GOTO 202
         ENDIF
C
  202 CONTINUE
         A   (4)=BUFF(4)
         A  (14)=   A(4)*FLOAT(MBAS(1))
         TE = A(14) 
C
C
      RETURN
      END
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
      REAL*4 FUNCTION R(S)
C
      IMPLICIT REAL*4 (A-H,O-Z)
      IMPLICIT INTEGER*2  (I-N)
C
      COMMON /BLCFN/ Z(6,6,6),A(216),CFN(6,6,6),NCFN(20)
      COMMON /BLEPS/ ERS(20), NITER(20), O, E
      DIMENSION S(5)
C
      R=O
      CALL PVGB(S)
              DO  I3 = 1, NCFN(3)
                  DO  I2 = 1, NCFN(2)
                      DO  I1 = 1, NCFN(1)
                          R=R+Z(I1,I2,I3)*CFN(I1,I2,I3)
                      END DO
                  END DO
              END DO
C
      RETURN
      END
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C     NBAS  1 2 3  4 5 6  7 8 9
C     NB1   1 1 1  2 2 2  3 3 3    1= S**   2=Tn  3= sin cos
C     NB2   1 2 3  1 2 3  1 2 3    1=0,1,2,3... 2=0,2,4... 3=0,1,3,5...
C
C
      SUBROUTINE PVGB(S)
C
      IMPLICIT REAL*4 (A-H,O-Z)
      IMPLICIT INTEGER*2  (I-N)
C
      COMMON /BLBAS/ XL(3),XR(3),SL(3),SR(3),NBAS(2,3),MBAS(3)
      COMMON /BLWRK/ WA(36,6),WB(216),WC(6,6,6)
      COMMON /BLCFN/ Z(6,6,6),A(216),CFN(6,6,6),NCFN(20)
      COMMON /BLEPS/ ERS(20), NITER(20), O, E
C
      DIMENSION S(5)
C
      DO  I=1,3
          CALL PVGB1( S(I), WA(1,I), NCFN(I), NBAS(1,I), NBAS(2,I) )
      END DO
C
      WB(4)=E
              DO  I3 = 1, NCFN(3)
                  WB(3)=WB(4)*WA(I3,3)
                  DO  I2 = 1, NCFN(2)
                      WB(2)=WB(3)*WA(I2,2)
                      DO  I1 = 1, NCFN(1)
                          WB(1)=WB(2)*WA(I1,1)
                          Z(I1,I2,I3)=WB(1)
                      END DO
                  END DO
              END DO
C
      RETURN
      END
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
      SUBROUTINE PVGB1(S,T,N,NBAS1,NBAS2)
C
C
C     NBAS  1 2 3  4 5 6  7 8 9
C     NB1   1 1 1  2 2 2  3 3 3    1= S**   2=Tn  3= sin cos
C     NB2   1 2 3  1 2 3  1 2 3    1=0,1,2,3... 2=0,2,4... 3=0,1,3,5...
C
      IMPLICIT REAL*4 (A-H,O-Z)
      IMPLICIT INTEGER*2  (I-N)
C
      COMMON /BLEPS/ ERS(20), NITER(20), O, ONE
      DIMENSION T(5)
C
      GOTO (10,20,30),NBAS1
  10  T(1)=ONE
      GOTO (11,13,14),NBAS2
  11  DO 12 I=2,N
  12  T(I)=T(I-1)*S
      RETURN
  13  T(2)=S*S
      GOTO 15
  14  T(2)=S
  15  SS=S*S
      DO 16 I=3,N
  16  T(I)=T(I-1)*SS
      RETURN
  20  T(1)=ONE
      T(2)=S
      SS=S+S
      NN=N
      IF(NBAS2.GT.1) NN=N+N
      DO 21 I=3,NN
  21  T(I)=T(I-1)*SS-T(I-2)
      GOTO (26,22,24),NBAS2
  22  DO 23 I=2,N
      J=2*I-1
  23  T(I)=T(J)
      GOTO 26
  24  DO 25 I=2,N
      J=2*I-2
  25  T(I)=T(J)
  26  RETURN
  30  T(1)=ONE
      SS= SIN(S)
      SC= COS(S)
      NN=N
      IF(NBAS2.GT.1) NN=N+N
      T(2)=SS
      T(3)=SC
      DO 31 I=4,NN,2
      T(I)  =T(I-2)*SC+T(I-1)*SS
  31  T(I+1)=T(I-1)*SC-T(I-2)*SS
      GOTO (36,32,34),NBAS2
  32  DO 33 I=2,N
      J=2*I-1
  33  T(I)=T(J)
      GOTO 36
  34  DO 35 I=2,N
      J=2*I-2
  35  T(I)=T(J)
  36  RETURN
      END
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
      BLOCK DATA PMIK
C
      IMPLICIT REAL*4 (A-H,O-Z)
      IMPLICIT INTEGER*2  (I-N)
C
      COMMON /BLBAS/ XL(3),XR(3),SL(3),SR(3),NBAS(2,3),MBAS(3)
      COMMON /BLDAT/ BUFF(50)
      COMMON /BLWRK/ WA(216),WB(216),WC(6,6,6)
      COMMON /BLCFN/ Z(216),A(216),CFN(6,6,6),NCFN(20)
      COMMON /BLEPS/ ERS(20), NITER(20), O, E
      COMMON /BLIOT/ NINP1, NINP2, NOUT1, NOUT2, NTAB1, NTAB2
C
      DATA NITER/20*0/,ERS/20*.0E0/,O/0.0E0/,E/1.000E0/
      DATA NINP1, NINP2 /1,2/
      DATA NOUT1, NOUT2 /5,6/
      DATA NTAB1, NTAB2 /3,4/
C
      DATA XL/1.0, -70.0, 1000. /
      DATA XR/13.,  70.0, 10000./
      DATA SL/0.0, -1.00, -1.00 /
      DATA SR/1.0,  1.00,  1.00 /
      DATA NBAS/ 1, 1, 2,  2, 2, 1 /
      DATA MBAS/1000, 1, 0 /
      DATA NCFN/ 2, 3, 4, 0, 0, 1, 2, 6, 24, 0, 24, 9*0 /
      DATA CFN/ 0.312742E+01, 0.267892E+01, 4*0.0,
     *         -0.142170E+00,-0.935806E+00, 4*0.0,
     *         -0.719833E+00,-0.372147E+00,22*0.0,
     *          0.323231E+00, 0.683911E+00, 4*0.0,
     *         -0.797365E+00,-0.605986E+00, 4*0.0,
     *         -0.100220E+00, 0.527238E+00,22*0.0, 
     *         -0.127842E+00,-0.321196E+00, 4*0.0,
     *         -0.676206E-01, 0.212179E+00, 4*0.0,
     *          0.136602E+00, 0.126485E+00,22*0.0, 
     *          0.149113E+00, 0.217611E+00, 4*0.0,
     *          0.166502E+00, 0.223654E+00, 4*0.0,
     *          0.834719E-01,-0.402023E-01,22*0.0,72*0.0/
      END
