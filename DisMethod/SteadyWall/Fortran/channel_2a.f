c
c     FULLY-IMPLICIT METHOD OF  
c     THE INCOMPRESSIBLE NAVIER-STOKES EQUATIONS
c     FOR THREE DIMENSIONAL PERIODIC TURBULET CHANNEL FLOW
c
c     1.No specific BCs for intermediate velocity and pressure
c
c     2.Crank-Nicolson mehtod for all terms
c     
c     3.LU decomposition with 2nd order accuracy 
c       for velocity and pressure decoupling
c
c     4.Keeping constant mass flow rate 
c
c     5.Convection terms are linearized with 2nd order accuracy
c
c     6.2D-FFT Solver for Poisson equation JIChoi 99/10/30
c
c     7.Streamwise and Spanwise averaging is added. JIChoi 99/11/2
c
c                                                      
c                                                   Kyoungyoun KIM
c                                          Flow Control Laboratory
c                             Department of Mechanical Engineering
c                                                            KAIST 
C-----------------------------------------------------------------------
C      Run test.                                         2010.11.16
C-----------------------------------------------------------------------
C      Replace the fft subroutines                       2015.01.19
C-----------------------------------------------------------------------

      Program MAIN
c	USE MSIMSL

      Include 'channel.h'
      Common/NSTEP/NTST
      COMMON/TSTEP/DT,CFLMAX
      COMMON/FINOUT/INCODE,IDTOPT,NWRITE,NREAD,IAVG,NPRN,INSF,NINS

      Real U(3,0:M1,0:M2,0:M3)
      Real P(M1,M2,M3)


      OPEN(31,FILE='Wallss.plt',STATUS='UNKNOWN')
      OPEN(32,FILE='FLOW.plt',STATUS='UNKNOWN')
      OPEN(33,FILE='CFL.plt',STATUS='UNKNOWN')

      Call SETUP
      IF(NREAD.EQ.0) Call INIUP(U,P,PRESG,PRESG3)
      IF(NREAD.EQ.1) Call READUP(U,P,PRESG,PRESG3)

      TIME=0.0
      Call DIVCHECK(U,TIME,DIVMAX)
      PRINT*, 'DIVMAX = ', DIVMAX
      Call CHKMF(U,TIME)

      IMORE=0                ! Index of written file

      NAVG=0
      DTR=DT

      Do 10 NTIME=1,NTST

          Call CFL(U,CFLM)
          IF (CFLM*DTR.GE.CFLMAX.AND.IDTOPT.EQ.1) DT=CFLMAX/CFLM
          IF (CFLM*DTR.LE.CFLMAX.OR.IDTOPT.NE.1) DT=DTR

          TIME=TIME+DT

          NAVG=NAVG+1

          Call GETUP(U,P,TIME,PRESG,PRESG3)
          Call DIVCHECK(U,TIME,DIVMAX)
          PRINT*, 'DIVMAX = ', DIVMAX
          Call CHKMF(U,TIME)
          CALL WALLSS(U,P,TIME)
          Call ENERGY(U,P)

          IF(IAVG.EQ.1) THEN
              IF(NTIME.EQ.1) Call SAVER
              IF(NTIME.NE.1) Call TAVER
          ENDIF

          IF(NWRITE.EQ.1.AND.MOD(NTIME,NINS).EQ.0)
     >    CALL WRITEUP(U,P,PRESG,PRESG3,IMORE)
          IF(IAVG.EQ.1.AND.MOD(NTIME,NPRN).EQ.0) Call OUTPUT(TIME,NAVG)

          IF(INSF.EQ.1.AND.MOD(NTIME,NINS).EQ.0) Call INSFIELD(U,P)

          Call CFL(U,CFLM)
          IF(MOD(NTIME,NPRN).EQ.0) Write(*,100) NTIME,DT,TIME,DIVMAX,CFLM*DT
          IF(MOD(NTIME,NPRN).EQ.0) PRINT*, PRESG, PRESG3
          !PAUSE
          WRITE(33,*) TIME,DT,CFLM*DT

 10   Continue

      CALL WRITEUP(U,P,PRESG,PRESG3,IMORE)
      IF (IAVG.EQ.1) CALL OUTPUT(TIME,NAVG)
      CALL PROFILE(U,P)

 100  Format(/'STEP=',I5,2X,'DT=',E11.4,2X,'TIME=',E11.4,
     >2X,'DIVMAX=',E11.4,2X,'CFL=',E11.4)
      Stop
      End

c***************** SETUP ***********************     
      Subroutine SETUP
      Common/NSTEP/NTST
      COMMON/TSTEP/DT,CFLMAX
      Common/DIM/N1,N2,N3,N1M,N2M,N3M
      Common/PARA/RE
      COMMON/VPERIN/VPER
      Common/SIZE/ALX,ALY,ALZ,VOL
      COMMON/FINOUT/INCODE,IDTOPT,NWRITE,NREAD,IAVG,NPRN,INSF,NINS
      COMMON/FILENAME/fileini,filegrd,fileout,fileavg
      CHARACTER*11 fileini,filegrd,fileout,fileavg

      CHARACTER*65 DUMMY

      OPEN(1,FILE='parame.ter',STATUS='OLD')
      READ (1,300) DUMMY
      WRITE(*,300) DUMMY
      READ (1,300) DUMMY
      WRITE(*,300) DUMMY
      READ (1,300) DUMMY
      WRITE(*,300) DUMMY
      READ (1,301) DUMMY,N1
      WRITE(*,301) DUMMY,N1
      READ (1,301) DUMMY,N2
      WRITE(*,301) DUMMY,N2
      READ (1,301) DUMMY,N3
      WRITE(*,301) DUMMY,N3
      READ (1,302) DUMMY,RE
      WRITE(*,302) DUMMY,RE
      READ (1,302) DUMMY,ALX
      WRITE(*,302) DUMMY,ALX
      READ (1,302) DUMMY,ALZ
      WRITE(*,302) DUMMY,ALZ
      READ (1,301) DUMMY,INCODE
      WRITE(*,301) DUMMY,INCODE
      READ (1,301) DUMMY,NTST
      WRITE(*,301) DUMMY,NTST
      READ (1,302) DUMMY,VPER
      WRITE(*,302) DUMMY,VPER
      READ (1,302) DUMMY,DT
      WRITE(*,302) DUMMY,DT
      READ (1,301) DUMMY,IDTOPT
      WRITE(*,301) DUMMY,IDTOPT
      READ (1,302) DUMMY,CFLMAX
      WRITE(*,302) DUMMY,CFLMAX
      READ (1,301) DUMMY,NWRITE
      WRITE(*,301) DUMMY,NWRITE
      READ (1,301) DUMMY,NREAD
      WRITE(*,301) DUMMY,NREAD
      READ (1,301) DUMMY,IAVG
      WRITE(*,301) DUMMY,IAVG
      READ (1,301) DUMMY,NPRN
      WRITE(*,301) DUMMY,NPRN
      READ (1,301) DUMMY,INSF
      WRITE(*,301) DUMMY,INSF
      READ (1,301) DUMMY,NINS
      WRITE(*,301) DUMMY,NINS
      READ (1,300) DUMMY
      WRITE(*,300) DUMMY
      READ (1,300) DUMMY
      WRITE(*,300) DUMMY
      READ (1,303) DUMMY,fileini
      WRITE(*,303) DUMMY,fileini
      READ (1,303) DUMMY,filegrd
      WRITE(*,303) DUMMY,filegrd
      READ (1,300) DUMMY
      WRITE(*,300) DUMMY
      READ (1,300) DUMMY
      WRITE(*,300) DUMMY
      CLOSE(1)

  300 FORMAT(A65)
  301 FORMAT(A45,I15)
  302 FORMAT(A45,E15.7)
  303 FORMAT(A45,A11)

C------------------------------------
C     PHYSICAL LENGTH
      IF(INCODE.EQ.1) THEN
          PI=ACOS(-1.0)
          ALX=PI
          ALZ=0.289*PI
      ENDIF
C------------------------------------

      N1M=N1-1
      N2M=N2-1
      N3M=N3-1

      Call MESH
      Call INDICES
      Call INIWAVE
      Return
      End

c***************** MESH ***********************     
      Subroutine MESH
      Include 'channel.h'
      Common/DIM/N1,N2,N3,N1M,N2M,N3M
      Common/MESH1/DX1,DX1Q,DX3,DX3Q
      Common/MESH2/Y(0:M2)
      Common/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      Common/MESH4/DYM(M2),DYC(M2),DYP(M2)
      Common/SIZE/ALX,ALY,ALZ,VOL
      COMMON/FILENAME/fileini,filegrd,fileout,fileavg
      CHARACTER*11 fileini,filegrd,fileout,fileavg
      PI=ACOS(-1.0)

c      CREATE THE UNIFORM GRID IN X2 DIRECTION      
c      ALY=2.0
c      Y(0)=0.0 
c      Do 10 J=1,N2 
c      Y(J)=ALY*DBLE(J-1)/DBLE(N2M) 
c 10   Continue 
      OPEN(2,FILE=filegrd,STATUS='OLD')
      READ(2,*) (Y(J),J=0,N2)
      CLOSE(2)

      ALY=Y(N2)

      VOL=ALX*ALY*ALZ

      DX1=DBLE(N1M)/ALX
      DX3=DBLE(N3M)/ALZ

      DX1Q=DX1**2.0
      DX3Q=DX3**2.0

      DY(1)=Y(2)
      Do 20 J=2,N2M
          DY(J)=Y(J+1)-Y(J)
          H(J)=0.5*(DY(J)+DY(J-1))
 20   Continue
      H(1)=0.5*DY(1)
      H(N2)=DY(N2M)*0.5

      Do 30 J=2,N2M-1
          HP(J)=1.0/H(J+1)/DY(J)
          HC(J)=(H(J+1)+H(J))/H(J+1)/H(J)/DY(J)
          HM(J)=1.0/H(J)/DY(J)
 30   Continue
C     J=1
      HP(1)=1.0/H(2)/DY(1)
      HC(1)=(1./H(2)+2./DY(1))/DY(1)
      HM(1)=2.0/DY(1)/DY(1)
C     J=N2M
      HP(N2M)=2.0/DY(N2M)/DY(N2M)
      HC(N2M)=(1./H(N2M)+2./DY(N2M))/DY(N2M)
      HM(N2M)=1.0/DY(N2M)/H(N2M)

      Do 40 J=2,N2M
          DYP(J)=1.0/H(J)/DY(J)
          DYC(J)=1.0/H(J)*(1.0/DY(J)+1.0/DY(J-1))
          DYM(J)=1.0/H(J)/DY(J-1)
 40   Continue

      Return
      End

c***************** INDICES ***********************     
      Subroutine INDICES
      Include 'channel.h'
      Common/DIM/N1,N2,N3,N1M,N2M,N3M
      Common/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      Common/INDX2/JPA(M2),JMU(M2),JMV(M2)

      Do 10 IC=1,N1M
          IPA(IC)=IC+1
 10   IMA(IC)=IC-1
      IPA(N1M)=1
      IMA(1)=N1M

      Do 20 KC=1,N3M
          KPA(KC)=KC+1
 20   KMA(KC)=KC-1
      KPA(N3M)=1
      KMA(1)=N3M

      Do 30 JC=1,N2M
          JPA(JC)=JC+1
          JMU(JC)=JC-1
 30   JMV(JC)=JC-1
      JPA(N2M)=N2M
      JMU(1)=1
      JMV(2)=2

      Return
      End


C  ****************************** INIUP **********************
C     THIS ROUTINE IS TO GIVE INITIAL FLOW FIELDS 
C     ON A PARABOLIC PROFILES
C     THIS SUBROUTINE IS MODIFIEDY BY JICHOI 99/07/27
C     TO MAKE CONSTANT FLOW RATE AT EACH PLANE
C     U IN YZ-PLANE, W IN XY-PLANE

      Subroutine INIUP(U,P,PRESG,PRESG3)
      Include 'channel.h'
      Common/DIM/N1,N2,N3,N1M,N2M,N3M
      Common/MESH1/DX1,DX1Q,DX3,DX3Q
      Common/MESH2/Y(0:M2)
      Common/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      Common/PARA/RE
      COMMON/VPERIN/VPER
      Common/SIZE/ALX,ALY,ALZ,VOL

      Real U(3,0:M1,0:M2,0:M3)
      Real P(M1,M2,M3)
      INTEGER*4 SEED
      REAL NXZ
      PI=ACOS(-1.0)
      DO 10 K=1,N3M
          DO 10 J=0,N2
              DO 10 I=1,N1M
                  U(1,I,J,K)=VPER*(RAN(SEED)-RAN(SEED))
                  U(2,I,J,K)=VPER*(RAN(SEED)-RAN(SEED))
                  U(3,I,J,K)=VPER*(RAN(SEED)-RAN(SEED))
c      U(1,I,J,K)=SIN(REAL(J)/REAL(N2)*PI)*U(1,I,J,K)
c      U(2,I,J,K)=SIN(REAL(J)/REAL(N2)*PI)*U(2,I,J,K)
c      U(3,I,J,K)=SIN(REAL(J)/REAL(N2)*PI)*U(3,I,J,K)
10    CONTINUE

C     IMPOSE ZERO VELOCITY AT BOUNDARY
      Do 15 I=1,N1M
          Do 15 K=1,N3M
              U(1,I,0,K)=0.0
              U(1,I,N2,K)=0.0
              U(2,I,1,K)=0.0
              U(2,I,N2,K)=0.0
              U(3,I,0,K)=0.0
              U(3,I,N2,K)=0.0
15    Continue

C     ELIMINATE MEAN QUANTITIES OF RANDOM FLUCTUATIONS
      DO 21 I=1,N1M
          V1M=0.
          DO 22 J=1,N2M
              DO 22 K=1,N3M
                  V1M=V1M+U(1,I,J,K)*DY(J)/DX3
22        CONTINUE
          V1M=V1M/ALY/ALZ
          DO 23 J=1,N2M
              DO 23 K=1,N3M
                  U(1,I,J,K)=U(1,I,J,K)-V1M
23        CONTINUE
21    CONTINUE

      DO 31 J=2,N2M
          V2M=0.
          DO 32 I=1,N1M
              DO 32 K=1,N3M
                  V2M=V2M+U(2,I,J,K)/DX1/DX3
32        CONTINUE
          V2M=V2M/ALX/ALZ
          DO 33 I=1,N1M
              DO 33 K=1,N3M
                  U(2,I,J,K)=U(2,I,J,K)-V2M
33        CONTINUE
31    CONTINUE

      DO 41 K=1,N3M
          V3M=0.
          DO 42 J=1,N2M
              DO 42 I=1,N1M
                  V3M=V3M+U(3,I,J,K)*DY(J)/DX1
42        CONTINUE
          V3M=V3M/ALX/ALY
          DO 43 I=1,N1M
              DO 43 J=1,N2M
                  U(3,I,J,K)=U(3,I,J,K)-V3M
43        CONTINUE
41    CONTINUE

C     IMPOSE LAMINAR VELOCITY PROFIELS IN U VELOCITIES
      DO 30 J=1,N2M
          JP=J+1
          DO 30 K=1,N3M
              DO 30 I=1,N1M
                  YH=0.5*(Y(J)+Y(JP))
                  REALU=YH*(ALY-YH)
                  U(1,I,J,K)=U(1,I,J,K)+REALU
30    CONTINUE

C     CHECK FLOW RATE
      FLOW1=0.0
      FLOW3=0.0
      DO 40 I=1,N1M
          DO 40 J=1,N2M
              DO 40 K=1,N3M
                  FLOW1=FLOW1+U(1,I,J,K)*DY(J)/DX1/DX3
                  FLOW3=FLOW3+U(3,I,J,K)*DY(J)/DX1/DX3
40    CONTINUE
      FLOW1=FLOW1/ALX/ALZ
      FLOW3=FLOW3/ALX/ALZ/ALY

      FLOWR=ALY**3.0/6.0       ! LAMINAR MASS FLOW RATE IN X1 DIRECTION
      DO 45 I=1,N1M
          DO 45 J=1,N2M
              DO 45 K=1,N3M
                  U(1,I,J,K)=FLOWR/FLOW1*U(1,I,J,K)
                  U(3,I,J,K)=U(3,I,J,K)-FLOW3 ! ZERO MASS FLOW RATE IN X3 DIRECTION
   45 CONTINUE

C     IMPOSE ZERO-PRESSURE FLUCTUATIONS
      DO 60 I=1,N1M
          DO 60 J=1,N2M
              DO 60 K=1,N3M
                  P(I,J,K)=0.0
   60 CONTINUE

C     INITIAL MEAN PRESSURE GRADIENT AT LAMINAR FLOW FIELD
      PRESG=0.0
      PRESG3=0.0      ! Z-direction pressure gradient

      RETURN
      END

C  ************************  READUP **********************
C     READ FLOW FIELD AND BOUNDARY CONDITIONS
C     AND MEAN PRESSURE GRADIENT

      Subroutine READUP(U,P,PRESG1,PRESG3)
      Include 'channel.h'
      Common/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/FILENAME/fileini,filegrd,fileout,fileavg
      CHARACTER*11 fileini,filegrd,fileout,fileavg

      Real U(3,0:M1,0:M2,0:M3)
      Real P(M1,M2,M3)

      OPEN(3,FILE=fileini,FORM='UNFORMATTED',STATUS='OLD')
      READ(3) PRESG1,PRESG3
      READ(3) (((U(1,I,J,K),U(2,I,J,K),U(3,I,J,K)
     >,K=1,N3M),J=0,N2),I=1,N1M)
      READ(3) (((P(I,J,K),K=1,N3M),J=1,N2M),I=1,N1M)
      CLOSE(3)

      RETURN
      END



C  ************************  WRITEUP **********************
C     WRITE FLOW FIELD AND BOUNDARY CONDITIONS
C     AND MEAN PRESSURE GRADIENT

      Subroutine WRITEUP(U,P,PRESG,PRESG3,IMORE)
      Include 'channel.h'
      Common/DIM/N1,N2,N3,N1M,N2M,N3M

      Real U(3,0:M1,0:M2,0:M3)
      Real P(M1,M2,M3)

      NFILE=100+IMORE
      WRITE(NFILE) PRESG,PRESG3
      WRITE(NFILE) (((U(1,I,J,K),U(2,I,J,K),U(3,I,J,K)
     >,K=1,N3M),J=0,N2),I=1,N1M)
      WRITE(NFILE) (((P(I,J,K),K=1,N3M),J=1,N2M),I=1,N1M)
      IMORE=IMORE+1
      RETURN
      END


c***************** GETUP ***********************     
      Subroutine GETUP(U,P,TIME,PRESG,PRESG3)
      Include 'channel.h'
      Common/DIM/N1,N2,N3,N1M,N2M,N3M

      Real U(3,0:M1,0:M2,0:M3)
      Real UH(3,0:M1,0:M2,0:M3)
      Real P(M1,M2,M3)
      Real DP(M1,M2,M3)

C INITIALIZE THE INTERMEDIATE VELOCITY AND PRESSURE
      DO 10 NV=1,3
          Do 10 I=1,N1
              Do 10 J=1,N2
                  Do 10 K=1,N3
 10   UH(NV,I,J,K)=0.0
      Do 20 I=1,N1
          Do 20 J=1,N2
              Do 20 K=1,N3
 20   DP(I,J,K)=0.0
      Call BCOND(U,TIME)

C CALCULATE THE INTERMEDIATE VELOCITY
      Call UHCALC(U,UH,P,PRESG,PRESG3)

C CALCULATE DP
      Call DPCALC(U,UH,DP)

C UPDATE THE N+1 TIME STEP VELOCITY AND PRESSURE
      Call UPCALC(U,P,UH,DP,PRESG,PRESG3)

      Return
      End

c***************** BCOND ***********************     
      Subroutine BCOND(U,TIME)
      Include 'channel.h'
      Common/DIM/N1,N2,N3,N1M,N2M,N3M
      Common/BCON/UBC(2,3,0:M1,0:M3)
      Real U(3,0:M1,0:M2,0:M3)

      Do 10 I=1,N1M
          Do 10 K=1,N3M
              UBC(1,1,I,K)=0.0  ! lower wall : no-slip condition
              UBC(1,2,I,K)=0.0
              UBC(1,3,I,K)=0.0
              UBC(2,1,I,K)=0.0  ! upper wall : no-slip condition
              UBC(2,2,I,K)=0.0
              UBC(2,3,I,K)=0.0
 10   Continue

      Return
      End

c***************** UHCALC ***********************     
      Subroutine UHCALC(U,UH,P,PRESG,PRESG3)
      Include 'channel.h'
      Common/DIM/N1,N2,N3,N1M,N2M,N3M

      Real U(3,0:M1,0:M2,0:M3)
      Real P(M1,M2,M3)

      Real UH(3,0:M1,0:M2,0:M3)
      Real RUH1(0:M1,0:M2,0:M3)
      Real RUH2(0:M1,0:M2,0:M3)
      Real RUH3(0:M1,0:M2,0:M3)

      Call RHS1(U,P,RUH1,PRESG)
      Call RHS2(U,P,RUH2)
      Call RHS3(U,P,RUH3,PRESG3)

      Call GETUH1(U,UH,RUH1)
      Call GETUH2(U,UH,RUH2)
      Call GETUH3(U,UH,RUH3)

      Return
      End

c***************** DPCALC ***********************     
      Subroutine DPCALC(U,UH,DP)
      Include 'channel.h'

      Real U(3,0:M1,0:M2,0:M3)
      Real UH(3,0:M1,0:M2,0:M3)
      Real DP(M1,M2,M3)
      Real RDP(M1,M2,M3)

      Call RHSDP(RDP,UH,U)
      Call GETDP(DP,RDP)

      Return
      End

c***************** UPCALC ***********************     
      Subroutine UPCALC(U,P,UH,DP,PRESG,PRESG3)
      Include 'channel.h'
      Common/DIM/N1,N2,N3,N1M,N2M,N3M
      Common/BCON/UBC(2,3,0:M1,0:M3)
      Common/MESH1/DX1,DX1Q,DX3,DX3Q
      Common/MESH2/Y(0:M2)
      Common/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      Common/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      Common/INDX2/JPA(M2),JMU(M2),JMV(M2)
      Common/NSTEP/NTST
      COMMON/TSTEP/DT,CFLMAX
      Common/SIZE/ALX,ALY,ALZ,VOL

      Real U(3,0:M1,0:M2,0:M3)
      Real P(M1,M2,M3)
      Real UH(3,0:M1,0:M2,0:M3)
      Real DP(M1,M2,M3)
      Real DMPRESG
      Real DMPRESG3

C     U1 VELOCITY UPDATE
C     To keep a constant mass flow rate,
C     calculate the mean pressure gradient difference
C     between n+1/2 time step And n-1/2 time step 

      DMPRESG=0.0
      FLOW1=0.0
      Do 1 I=1,N1M
          IM=IMA(I)
          Do 1 K=1,N3M
              Do 1 J=1,N2M
                  FLOW1=FLOW1+U(1,I,J,K)*DY(J)/DX3/DX1
                  DMPRESG=DMPRESG
     >            +UH(1,I,J,K)*DY(J)/DX3/DX1
     >            -DT*(DP(I,J,K)-DP(IM,J,K))*DY(J)/DX3
  1   Continue

      DMPRESG=(DMPRESG-FLOW1)/VOL/DT

      Do 10 I=1,N1M
          IM=IMA(I)
          Do 10 K=1,N3M
              Do 10 J=1,N2M
                  U(1,I,J,K)=UH(1,I,J,K)
     >            -DT*(DP(I,J,K)-DP(IM,J,K))*DX1
     >            -DT*DMPRESG
  10  Continue

      Do 11 I=1,N1M
          Do 11 K=1,N3M
              U(1,I,0,K)=UBC(1,1,I,K)
              U(1,I,N2,K)=UBC(2,1,I,K)
  11  Continue

c     U2 VELOCITY UPDATE

      Do 20 I=1,N1M
          Do 20 K=1,N3M
              Do 20 J=2,N2M
                  U(2,I,J,K)=UH(2,I,J,K)
     >            -DT*(DP(I,J,K)-DP(I,J-1,K))/H(J)
  20  Continue

      Do 21 I=1,N1M
          Do 21 K=1,N3M
              U(2,I,1,K)=UBC(1,2,I,K)
              U(2,I,N2,K)=UBC(2,2,I,K)
  21  Continue

C     U3 VELOCITY UPDATE
C     To make the mass flow rate in x3 direction zero
C     calculate the mean pressure gradient in x3 direction
C     this term may be zero at quasi-steady state

      DMPRESG3=0.0
      FLOW3=0.0
      Do 3 K=1,N3M
          KM=KMA(K)
          Do 3 I=1,N1M
              Do 3 J=1,N2M
                  FLOW3=FLOW3+U(3,I,J,K)*DY(J)/DX1/DX3
                  DMPRESG3=DMPRESG3
     >            +UH(3,I,J,K)*DY(J)/DX3/DX1
     >            -DT*(DP(I,J,K)-DP(I,J,KM))*DY(J)/DX1
  3   Continue
      DMPRESG3=(DMPRESG3-FLOW3)/VOL/DT

      Do 30 I=1,N1M
          Do 30 K=1,N3M
              KM=KMA(K)
              Do 30 J=1,N2M
                  U(3,I,J,K)=UH(3,I,J,K)
     >            -DT*(DP(I,J,K)-DP(I,J,KM))*DX3
     >            -DT*DMPRESG3
  30  Continue

      Do 31 I=1,N1M
          Do 31 K=1,N3M
              U(3,I,0,K)=UBC(1,3,I,K)
              U(3,I,N2,K)=UBC(2,3,I,K)
  31  Continue

C     PRESSURE UPDATE     

      Do 40 I=1,N1M
          Do 40 J=1,N2M
              Do 40 K=1,N3M
                  P(I,J,K)=P(I,J,K)+DP(I,J,K)
  40  Continue

      PRESG=PRESG+DMPRESG
      PRESG3=PRESG3+DMPRESG3

      Return
      End

c***************** RHS1 ***********************     
      Subroutine RHS1(U,P,RUH1,PRESG)
      Include 'channel.h'
      Common/DIM/N1,N2,N3,N1M,N2M,N3M
      Common/PARA/RE
      Common/BCON/UBC(2,3,0:M1,0:M3)
      Common/MESH1/DX1,DX1Q,DX3,DX3Q
      Common/MESH2/Y(0:M2)
      Common/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      Common/NSTEP/NTST
      COMMON/TSTEP/DT,CFLMAX
      Common/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      Common/INDX2/JPA(M2),JMU(M2),JMV(M2)

      Real U(3,0:M1,0:M2,0:M3)
      Real P(M1,M2,M3)
      Real RUH1(0:M1,0:M2,0:M3)

      Do 10 J=1,N2M
          JP=J+1
          JM=J-1
          JUM=J-JMU(J)
          JUP=JPA(J)-J
          Do 10 I=1,N1M
              IP=IPA(I)
              IM=IMA(I)
              Do 10 K=1,N3M
                  KP=KPA(K)
                  KM=KMA(K)

                  VISCOS=0.5*DX1Q/RE*(U(1,IP,J,K )-2.0*U(1,I,J,K)+U(1,IM,J,K ))
     >            +0.5*DX3Q/RE*(U(1,I ,J,KP)-2.0*U(1,I,J,K)+U(1,I ,J,KM))
     >            +0.5/RE*(HP(J)*U(1,I,JP,K)
     >            -HC(J)*U(1,I,J ,K)
     >            +HM(J)*U(1,I,JM,K))

                  PRESSG1=DX1*(P(I,J,K)-P(IM,J,K))
     >            +PRESG

                  BC_DOWN=0.5/RE*HM(1)*UBC(1,1,I,K)
     >            +0.5/DY(1)*0.5*(U(2,I,1,K)+U(2,IM,1,K))*UBC(1,1,I,K)
     >            +0.5/DY(1)*U(1,I,0,K)*0.5*(UBC(1,2,I,K)+UBC(1,2,IM,K))

                  BC_UP  =0.5/RE*HP(N2M)*UBC(2,1,I,K)
     >            -0.5/DY(N2M)*0.5*(U(2,I,N2,K)+U(2,IM,N2,K))*UBC(2,1,I,K)
     >            -0.5/DY(N2M)*U(1,I,N2,K)*0.5*(UBC(2,2,I,K)+UBC(2,2,IM,K))

                  BC=(1-JUM)*BC_DOWN
     >            +(1-JUP)*BC_UP

                  RUH1(I,J,K)=1./DT*U(1,I,J,K)
     >            -PRESSG1+VISCOS
     >            +BC

C     R1=r1-AU^n

C     M11U^N      
                  V2=0.5*(U(2,I,JP,K)+U(2,IM,JP,K))
                  V1=0.5*(U(2,I,J ,K)+U(2,IM,J ,K))
                  APJ=JUP*(
     >            -0.5*HP(J)/RE
     >            +0.5/DY(J)*V2/H(JP)*DY(J)/2.0
     >            )
                  ACJ=   0.5*HC(J)/RE
     >            +0.5/DY(J)*(JUP*V2/H(JP)*DY(JP)/2.0
     >            -JUM*V1/H(J )*DY(JM)/2.0)
                  AMJ=JUM*(
     >            -0.5*HM(J)/RE
     >            -0.5/DY(J)*V1/H(J)*DY(J)/2.0
     >            )
                  U2=0.5*(U(1,IP,J,K)+U(1,I ,J,K))
                  U1=0.5*(U(1,I ,J,K)+U(1,IM,J,K))
                  API=  -0.5*DX1Q/RE
     >            +DX1*U2*0.5
                  ACI=   DX1Q/RE
     >            +DX1*(U2*0.5-U1*0.5)
                  AMI=  -0.5*DX1Q/RE
     >            -DX1*U1*0.5
                  W2=0.5*(U(3,IM,J,KP)+U(3,I,J,KP))
                  W1=0.5*(U(3,IM,J,K )+U(3,I,J,K ))
                  APK=   -0.5*DX3Q/RE
     >            +0.5*DX3*W2*0.5
                  ACK=    DX3Q/RE
     >            +0.5*DX3*(W2*0.5-W1*0.5)
                  AMK=   -0.5*DX3Q/RE
     >            -0.5*DX3*W1*0.5
                  RM11U_N=APJ*U(1,I,JP,K)
     >            +ACJ*U(1,I,J ,K)
     >            +AMJ*U(1,I,JM,K)
     >            +API*U(1,IP,J,K)
     >            +ACI*U(1,I ,J,K)
     >            +AMI*U(1,IM,J,K)
     >            +APK*U(1,I,J,KP)
     >            +ACK*U(1,I,J,K )
     >            +AMK*U(1,I,J,KM)
C     M12V^N
                  U2=1.0/H(JP)*(DY(J )/2.*U(1,I,JP,K)+DY(JP)/2.*U(1,I,J ,K))
                  U1=1.0/H(J )*(DY(JM)/2.*U(1,I,J ,K)+DY(J )/2.*U(1,I,JM,K))
                  V2=0.5*(U(2,I,JP,K)+U(2,IM,JP,K))
                  V1=0.5*(U(2,I,J ,K)+U(2,IM,J ,K))
                  RM12V_N=JUP*0.5/DY(J)*U2*V2
     >            -JUM*0.5/DY(J)*U1*V1
C     M13W^N
                  U2=0.5*(U(1,I,J,KP)+U(1,I,J,K ))
                  U1=0.5*(U(1,I,J,K )+U(1,I,J,KM))
                  W2=0.5*(U(3,IM,J,KP)+U(3,I,J,KP))
                  W1=0.5*(U(3,IM,J,K )+U(3,I,J,K ))
                  RM13W_N=0.5*DX3*(U2*W2-U1*W1)

                  RUH1(I,J,K)=RUH1(I,J,K)
     >            -1./DT*U(1,I,J,K)
     >            -RM11U_N-RM12V_N-RM13W_N

 10   Continue

      Return
      End

c***************** RHS2 ***********************     
      Subroutine RHS2(U,P,RUH2)
      Include 'channel.h'
      Common/DIM/N1,N2,N3,N1M,N2M,N3M
      Common/PARA/RE
      Common/BCON/UBC(2,3,0:M1,0:M3)
      Common/MESH1/DX1,DX1Q,DX3,DX3Q
      Common/MESH2/Y(0:M2)
      Common/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      Common/MESH4/DYM(M2),DYC(M2),DYP(M2)
      Common/NSTEP/NTST
      COMMON/TSTEP/DT,CFLMAX
      Common/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      Common/INDX2/JPA(M2),JMU(M2),JMV(M2)

      Real U(3,0:M1,0:M2,0:M3)
      Real P(M1,M2,M3)
      Real RUH2(0:M1,0:M2,0:M3)

      Do 10 J=2,N2M
          JP=J+1
          JM=J-1
          JUM=J-JMV(J)
          JUP=JPA(J)-J

          Do 10 I=1,N1M
              IP=IPA(I)
              IM=IMA(I)
              Do 10 K=1,N3M
                  KP=KPA(K)
                  KM=KMA(K)

                  VISCOS= 0.5*DX1Q/RE*(U(2,IP,J,K)-2.0*U(2,I,J,K)+U(2,IM,J,K))
     >            +0.5*DX3Q/RE*(U(2,I,J,KP)-2.0*U(2,I,J,K)+U(2,I,J,KM))
     >            +0.5/RE*(DYP(J)*U(2,I,JP,K)
     >            -DYC(J)*U(2,I,J ,K)
     >            +DYM(J)*U(2,I,JM,K))

                  PRESSG2=(P(I,J,K)-P(I,JM,K))/H(J)

                  BC_DOWN=0.5/RE*DYM(2)*UBC(1,2,I,K)
     >            +1./H(2)*0.5*(U(2,I,2,K)+U(2,I,1,K))*0.5*UBC(1,2,I,K)

                  BC_UP  =0.5/RE*DYP(N2M)*UBC(2,2,I,K)
     >            -1./H(N2M)*0.5*(U(2,I,N2,K)+U(2,I,N2M,K))*0.5*UBC(2,2,I,K)

                  BC=(1-JUM)*BC_DOWN
     >            +(1-JUP)*BC_UP

                  RUH2(I,J,K)=1./DT*U(2,I,J,K)
     >            -PRESSG2+VISCOS
     >            +BC

C     R2=r2-AU^n      

C     M22V^N

                  V2=0.5*(U(2,I,JP,K)+U(2,I,J ,K))
                  V1=0.5*(U(2,I,J ,K)+U(2,I,JM,K))
                  APJ=JUP*(
     >            -0.5*DYP(J)/RE
     >            +1.0/H(J)*V2*0.5
     >            )
                  ACJ=
     >            +0.5*DYC(J)/RE
     >            +1.0/H(J)*(V2*0.5-V1*0.5)
                  AMJ=JUM*(
     >            -0.5*DYM(J)/RE
     >            -1.0/H(J)*V1*0.5
     >            )
                  U2=1.0/H(J)*(DY(J)/2.0*U(1,IP,JM,K)+DY(JM)/2.0*U(1,IP,J,K))
                  U1=1.0/H(J)*(DY(J)/2.0*U(1,I ,JM,K)+DY(JM)/2.0*U(1,I ,J,K))
                  API=
     >            -0.5*DX1Q/RE
     >            +0.5*DX1*U2*0.5
                  ACI=
     >            +DX1Q/RE
     >            +0.5*DX1*(U2*0.5-U1*0.5)
                  AMI=
     >            -0.5*DX1Q/RE
     >            -0.5*DX1*U1*0.5
                  W2=1.0/H(J)*(DY(J)/2.0*U(3,I,JM,KP)+DY(JM)/2.0*U(3,I,J,KP))
                  W1=1.0/H(J)*(DY(J)/2.0*U(3,I,JM,K )+DY(JM)/2.0*U(3,I,J,K ))
                  APK=
     >            -0.5*DX3Q/RE
     >            +0.5*DX3*W2/2.0
                  ACK=
     >            +DX3Q/RE
     >            +0.5*DX3*(W2/2.0-W1/2.0)
                  AMK=
     >            -0.5*DX3Q/RE
     >            -0.5*DX3*W1/2.0
                  RM22V_N=APJ*U(2,I,JP,K)
     >            +ACJ*U(2,I,J ,K)
     >            +AMJ*U(2,I,JM,K)
     >            +API*U(2,IP,J,K)
     >            +ACI*U(2,I ,J,K)
     >            +AMI*U(2,IM,J,K)
     >            +APK*U(2,I,J,KP)
     >            +ACK*U(2,I,J,K )
     >            +AMK*U(2,I,J,KM)
C     M21U^N
                  V2=0.5*(U(2,IP,J,K)+U(2,I ,J,K))
                  V1=0.5*(U(2,I ,J,K)+U(2,IM,J,K))
                  U2=1.0/H(J)*(DY(J)/2.*U(1,IP,JM,K)+DY(JM)/2.*U(1,IP,J,K))
                  U1=1.0/H(J)*(DY(J)/2.*U(1,I ,JM,K)+DY(JM)/2.*U(1,I ,J,K))
                  RM21U_N=0.5*DX1*(V2*U2-V1*U1)
C     M23W^N
                  V2=0.5*(U(2,I,J,KP)+U(2,I,J,K ))
                  V1=0.5*(U(2,I,J,K )+U(2,I,J,KM))
                  W2=1.0/H(J)*(DY(J)/2.*U(3,I,JM,KP)+DY(JM)/2.*U(3,I,J,KP))
                  W1=1.0/H(J)*(DY(J)/2.*U(3,I,JM,K )+DY(JM)/2.*U(3,I,J,K ))
                  RM23W_N=0.5*DX3*(V2*W2-V1*W1)

                  RUH2(I,J,K)=RUH2(I,J,K)
     >            -1./DT*U(2,I,J,K)
     >            -RM21U_N-RM22V_N-RM23W_N

 10   Continue


      Return
      End

c***************** RHS3 ***********************     
      Subroutine RHS3(U,P,RUH3,PRESG3)
      Include 'channel.h'
      Common/DIM/N1,N2,N3,N1M,N2M,N3M
      Common/PARA/RE
      Common/BCON/UBC(2,3,0:M1,0:M3)
      Common/MESH1/DX1,DX1Q,DX3,DX3Q
      Common/MESH2/Y(0:M2)
      Common/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      Common/NSTEP/NTST
      COMMON/TSTEP/DT,CFLMAX
      Common/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      Common/INDX2/JPA(M2),JMU(M2),JMV(M2)

      Real U(3,0:M1,0:M2,0:M3)
      Real P(M1,M2,M3)
      Real RUH3(0:M1,0:M2,0:M3)

      Do 10 J=1,N2M
          JP=J+1
          JM=J-1
          JUM=J-JMU(J)
          JUP=JPA(J)-J

          Do 10 I=1,N1M
              IP=IPA(I)
              IM=IMA(I)

              Do 10 K=1,N3M
                  KP=KPA(K)
                  KM=KMA(K)

                  VISCOS=0.5*DX1Q/RE*(U(3,IP,J,K )-2.0*U(3,I,J,K)+U(3,IM,J,K ))
     >            +0.5*DX3Q/RE*(U(3,I ,J,KP)-2.0*U(3,I,J,K)+U(3,I ,J,KM))
     >            +0.5/RE*(HP(J)*U(3,I,JP,K)
     >            -HC(J)*U(3,I,J ,K)
     >            +HM(J)*U(3,I,JM,K))

                  PRESSG3=DX3*(P(I,J,K)-P(I,J,KM))
     >            +PRESG3

                  BC_DOWN=0.5/RE*HM(1)*UBC(1,3,I,K)
     >            +0.5/DY(1)*0.5*(U(2,I,1,KM)+U(2,I,1,K))*UBC(1,3,I,K)
     >            +0.5/DY(1)*U(3,I,0,K)*0.5*(UBC(1,2,I,K)+UBC(1,2,I,KM))

                  BC_UP  =0.5/RE*HP(N2M)*UBC(2,3,I,K)
     >            -0.5/DY(N2M)*0.5*(U(2,I,N2,KM)+U(2,I,N2,K))*UBC(2,3,I,K)
     >            -0.5/DY(N2M)*U(3,I,N2,K)*0.5*(UBC(2,2,I,K)+UBC(2,2,I,KM))

                  BC=(1-JUM)*BC_DOWN
     >            +(1-JUP)*BC_UP

                  RUH3(I,J,K)=1./DT*U(3,I,J,K)
     >            -PRESSG3+VISCOS
     >            +BC

C     R3=r3-AU^n      

C     M33W^N

                  V2=0.5*(U(2,I,JP,K)+U(2,I,JP,KM))
                  V1=0.5*(U(2,I,J ,K)+U(2,I,J ,KM))
                  APJ=JUP*(
     >            -0.5*HP(J)/RE
     >            +0.5/DY(J)*V2/H(JP)*DY(J)/2.0
     >            )
                  ACJ=
     >            +0.5*HC(J)/RE
     >            +0.5/DY(J)*(JUP*V2/H(JP)*DY(JP)/2.0
     >            -JUM*V1/H(J)*DY(JM)/2.0)
                  AMJ=JUM*(
     >            -0.5*HM(J)/RE
     >            -0.5/DY(J)*V1/H(J)*DY(J)/2.0
     >            )
                  W2=0.5*(U(3,I,J,KP)+U(3,I,J,K ))
                  W1=0.5*(U(3,I,J,K )+U(3,I,J,KM))
                  APK=
     >            -0.5*DX3Q/RE
     >            +DX3*W2/2.0
                  ACK=
     >            +DX3Q/RE
     >            +DX3*(W2/2.0-W1/2.0)
                  AMK=
     >            -0.5*DX3Q/RE
     >            -DX3*W1/2.0
                  U2=0.5*(U(1,IP,J,K)+U(1,IP,J,KM))
                  U1=0.5*(U(1,I ,J,K)+U(1,I ,J,KM))
                  API=
     >            -0.5*DX1Q/RE
     >            +0.5*DX1*U2/2.0
                  ACI=
     >            +DX1Q/RE
     >            +0.5*DX1*(U2/2.0-U1/2.0)
                  AMI=
     >            -0.5*DX1Q/RE
     >            -0.5*DX1*U1/2.0
                  RM33W_N=APJ*U(3,I,JP,K)
     >            +ACJ*U(3,I,J ,K)
     >            +AMJ*U(3,I,JM,K)
     >            +API*U(3,IP,J,K)
     >            +ACI*U(3,I ,J,K)
     >            +AMI*U(3,IM,J,K)
     >            +APK*U(3,I,J,KP)
     >            +ACK*U(3,I,J,K )
     >            +AMK*U(3,I,J,KM)
C     M31U^N
                  W2=0.5*(U(3,IP,J,K)+U(3,I ,J,K))
                  W1=0.5*(U(3,I ,J,K)+U(3,IM,J,K))
                  U2=0.5*(U(1,IP,J,K)+U(1,IP,J,KM))
                  U1=0.5*(U(1,I ,J,K)+U(1,I ,J,KM))
                  RM31U_N=0.5*DX1*(W2*U2-W1*U1)
C     M32V^N
                  W2=1.0/H(JP)*(DY(J )/2.*U(3,I,JP,K)+DY(JP)/2.*U(3,I,J ,K))
                  W1=1.0/H(J )*(DY(JM)/2.*U(3,I,J ,K)+DY(J )/2.*U(3,I,JM,K))
                  V2=0.5*(U(2,I,JP,K)+U(2,I,JP,KM))
                  V1=0.5*(U(2,I,J ,K)+U(2,I,J ,KM))
                  RM32V_N=JUP*0.5/DY(J)*W2*V2
     >            -JUM*0.5/DY(J)*W1*V1

                  RUH3(I,J,K)=RUH3(I,J,K)
     >            -1./DT*U(3,I,J,K)
     >            -RM31U_N-RM32V_N-RM33W_N

 10   Continue


      Return
      End

c***************** RHSDP ***********************     
      Subroutine RHSDP(RDP,UH,U)
      Include 'channel.h'
      Common/DIM/N1,N2,N3,N1M,N2M,N3M
      Common/BCON/UBC(2,3,0:M1,0:M3)
      Common/MESH1/DX1,DX1Q,DX3,DX3Q
      Common/MESH2/Y(0:M2)
      Common/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      Common/NSTEP/NTST
      COMMON/TSTEP/DT,CFLMAX
      Common/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      Common/INDX2/JPA(M2),JMU(M2),JMV(M2)

      Real U(3,0:M1,0:M2,0:M3)
      Real UH(3,0:M1,0:M2,0:M3)
      Real RDP(M1,M2,M3)

      Do 10 J=1,N2M
          JP=J+1
          JM=J-1
          JUM=J-JMU(J)
          JUP=JPA(J)-J

          Do 10 I=1,N1M
              IP=IPA(I)
              IM=IMA(I)

              Do 10 K=1,N3M
                  KP=KPA(K)
                  KM=KMA(K)

                  DIVUH=(UH(1,IP,J,K)-UH(1,I,J,K))*DX1
     >            +(JUP*UH(2,I,JP,K)-JUM*UH(2,I,J,K))/DY(J)
     >            +(UH(3,I,J,KP)-UH(3,I,J,K))*DX3

                  CBC=(1-JUM)*UBC(1,2,I,K)/DY(J)
     >            -(1-JUP)*UBC(2,2,I,K)/DY(J)

                  RDP(I,J,K)=(DIVUH-CBC)/DT
  10  Continue

      Return
      End


c***************** GETUH1 ***********************     
      Subroutine GETUH1(U,UH,RUH1)
      Include 'channel.h'
      Common/DIM/N1,N2,N3,N1M,N2M,N3M
      Common/PARA/RE
      Common/MESH1/DX1,DX1Q,DX3,DX3Q
      Common/MESH2/Y(0:M2)
      Common/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      Common/NSTEP/NTST
      COMMON/TSTEP/DT,CFLMAX
      Common/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      Common/INDX2/JPA(M2),JMU(M2),JMV(M2)

      Real U(3,0:M1,0:M2,0:M3)
      Real UH(3,0:M1,0:M2,0:M3)
      Real RUH1(0:M1,0:M2,0:M3)
      Real API(M1),ACI(M1),AMI(M1)
      Real APJ(M2),ACJ(M2),AMJ(M2)
      Real APK(M3),ACK(M3),AMK(M3)
      Real R1(M1),R2(M2),R3(M3)


      DO 2 I=1,N1M
          IP=IPA(I)
          IM=IMA(I)
          Do 2 K=1,N3M
              Do 20 J=1,N2M
                  JP=J+1
                  JM=J-1
                  JUM=J-JMU(J)
                  JUP=JPA(J)-J

                  V2=0.5*(U(2,I,JP,K)+U(2,IM,JP,K))
                  V1=0.5*(U(2,I,J ,K)+U(2,IM,J ,K))

                  APJ(J)=JUP*(
     >            -0.5*HP(J)/RE
     >            +0.5/DY(J)*V2/H(JP)*DY(J)/2.0
     >            )*DT
                  ACJ(J)=1.0+(
     >            0.5*HC(J)/RE
     >            +0.5/DY(J)*(JUP*V2/H(JP)*DY(JP)/2.0
     >            -JUM*V1/H(J )*DY(JM)/2.0)
     >            )*DT
                  AMJ(J)=JUM*(
     >            -0.5*HM(J)/RE
     >            -0.5/DY(J)*V1/H(J)*DY(J)/2.0
     >            )*DT
                  R2(J)=RUH1(I,J,K)*DT
  20          Continue
              Call TDMA1(AMJ,ACJ,APJ,R2,UH,I,K,1)
  2   Continue

      Do 1 J=1,N2M
          Do 1 K=1,N3M
              Do 10 I=1,N1M
                  IP=IPA(I)
                  IM=IMA(I)

                  U2=0.5*(U(1,IP,J,K)+U(1,I ,J,K))
                  U1=0.5*(U(1,I ,J,K)+U(1,IM,J,K))

                  API(I)=(
     >            -0.5*DX1Q/RE
     >            +DX1*U2*0.5
     >            )*DT
                  ACI(I)=1.0+(
     >            +DX1Q/RE
     >            +DX1*(U2*0.5-U1*0.5)
     >            )*DT
                  AMI(I)=(
     >            -0.5*DX1Q/RE
     >            -DX1*U1*0.5
     >            )*DT
                  R1(I)=UH(1,I,J,K)
  10          Continue
              Call CTDMA1(AMI,ACI,API,R1,UH,J,K,1,N1M)
  1   Continue

      DO 3 I=1,N1M
          IP=IPA(I)
          IM=IMA(I)
          Do 3 J=1,N2M
              Do 30 K=1,N3M
                  KP=KPA(K)
                  KM=KMA(K)

                  W2=0.5*(U(3,IM,J,KP)+U(3,I,J,KP))
                  W1=0.5*(U(3,IM,J,K )+U(3,I,J,K ))

                  APK(K)=(
     >            -0.5*DX3Q/RE
     >            +0.5*DX3*W2*0.5
     >            )*DT
                  ACK(K)=1.+(
     >            +DX3Q/RE
     >            +0.5*DX3*(W2*0.5-W1*0.5)
     >            )*DT
                  AMK(K)=(
     >            -0.5*DX3Q/RE
     >            -0.5*DX3*W1*0.5
     >            )*DT
                  R3(K)=UH(1,I,J,K)
  30          Continue
              Call CTDMA3(AMK,ACK,APK,R3,UH,I,J,1,N3M)
  3   Continue

      Return
      End

c***************** GETUH2 ***********************     
      Subroutine GETUH2(U,UH,RUH2)
      Include 'channel.h'
      Common/DIM/N1,N2,N3,N1M,N2M,N3M
      Common/PARA/RE
      Common/MESH1/DX1,DX1Q,DX3,DX3Q
      Common/MESH2/Y(0:M2)
      Common/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      Common/MESH4/DYM(M2),DYC(M2),DYP(M2)
      Common/NSTEP/NTST
      COMMON/TSTEP/DT,CFLMAX
      Common/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      Common/INDX2/JPA(M2),JMU(M2),JMV(M2)

      Real U(3,0:M1,0:M2,0:M3)
      Real UH(3,0:M1,0:M2,0:M3)
      Real RUH2(0:M1,0:M2,0:M3)
      Real API(M1),ACI(M1),AMI(M1)
      Real APJ(M2),ACJ(M2),AMJ(M2)
      Real APK(M3),ACK(M3),AMK(M3)
      Real R1(M1),R2(M2),R3(M3)

      Do 2 I=1,N1M
          IP=IPA(I)
          IM=IMA(I)
          Do 2 K=1,N3M
              KP=KPA(K)
              KM=KMA(K)
              Do 20 J=2,N2M
                  JP=J+1
                  JM=J-1
                  JUM=J-JMV(J)
                  JUP=JPA(J)-J

                  V2=0.5*(U(2,I,JP,K)+U(2,I,J ,K))
                  V1=0.5*(U(2,I,J ,K)+U(2,I,JM,K))

                  APJ(J)=JUP*(
     >            -0.5*DYP(J)/RE
     >            +1.0/H(J)*V2*0.5
     >            )*DT
                  ACJ(J)=1.0+(
     >            +0.5*DYC(J)/RE
     >            +1.0/H(J)*(V2*0.5-V1*0.5)
     >            )*DT
                  AMJ(J)=JUM*(
     >            -0.5*DYM(J)/RE
     >            -1.0/H(J)*V1*0.5
     >            )*DT

C     M21UH
                  V2=0.5*(U(2,IP,J,K)+U(2,I ,J,K))
                  V1=0.5*(U(2,I ,J,K)+U(2,IM,J,K))
                  UH2=1.0/H(J)*(DY(J)/2.*UH(1,IP,JM,K)+DY(JM)/2.*UH(1,IP,J,K))
                  UH1=1.0/H(J)*(DY(J)/2.*UH(1,I ,JM,K)+DY(JM)/2.*UH(1,I ,J,K))
                  RM21UH=0.5*DX1*(V2*UH2-V1*UH1)

                  R2(J)=DT*(RUH2(I,J,K)-RM21UH)

  20          Continue
              Call TDMA2(AMJ,ACJ,APJ,R2,UH,I,K)
  2   Continue

      DO 1 J=2,N2M
          JP=J+1
          JM=J-1
          Do 1 K=1,N3M
              KP=KPA(K)
              KM=KMA(K)
              Do 10 I=1,N1M
                  IP=IPA(I)
                  IM=IMA(I)

                  U2=1.0/H(J)*(DY(J)/2.0*U(1,IP,JM,K)+DY(JM)/2.0*U(1,IP,J,K))
                  U1=1.0/H(J)*(DY(J)/2.0*U(1,I ,JM,K)+DY(JM)/2.0*U(1,I ,J,K))

                  API(I)=(
     >            -0.5*DX1Q/RE
     >            +0.5*DX1*U2*0.5
     >            )*DT
                  ACI(I)=1.+(
     >            +DX1Q/RE
     >            +0.5*DX1*(U2*0.5-U1*0.5)
     >            )*DT
                  AMI(I)=(
     >            -0.5*DX1Q/RE
     >            -0.5*DX1*U1*0.5
     >            )*DT
                  R1(I)=UH(2,I,J,K)
  10          CONTINUE
              Call CTDMA1(AMI,ACI,API,R1,UH,J,K,2,N1M)
  1   Continue

      DO 3 I=1,N1M
          IP=IPA(I)
          IM=IMA(I)
          Do 3 J=2,N2M
              JP=J+1
              JM=J-1
              Do 30 K=1,N3M
                  KP=KPA(K)
                  KM=KMA(K)

                  W2=1.0/H(J)*(DY(J)/2.0*U(3,I,JM,KP)+DY(JM)/2.0*U(3,I,J,KP))
                  W1=1.0/H(J)*(DY(J)/2.0*U(3,I,JM,K )+DY(JM)/2.0*U(3,I,J,K ))

                  APK(K)=(
     >            -0.5*DX3Q/RE
     >            +0.5*DX3*W2/2.0
     >            )*DT
                  ACK(K)=1.+(
     >            +DX3Q/RE
     >            +0.5*DX3*(W2/2.0-W1/2.0)
     >            )*DT
                  AMK(K)=(
     >            -0.5*DX3Q/RE
     >            -0.5*DX3*W1/2.0
     >            )*DT
                  R3(K)=UH(2,I,J,K)
  30          Continue
              Call CTDMA3(AMK,ACK,APK,R3,UH,I,J,2,N3M)
  3   Continue

      Return
      End

c***************** GETUH3 ***********************     
      Subroutine GETUH3(U,UH,RUH3)
      Include 'channel.h'
      Common/DIM/N1,N2,N3,N1M,N2M,N3M
      Common/PARA/RE
      Common/MESH1/DX1,DX1Q,DX3,DX3Q
      Common/MESH2/Y(0:M2)
      Common/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      Common/NSTEP/NTST
      COMMON/TSTEP/DT,CFLMAX
      Common/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      Common/INDX2/JPA(M2),JMU(M2),JMV(M2)

      Real U(3,0:M1,0:M2,0:M3)
      Real UH(3,0:M1,0:M2,0:M3)
      Real RUH3(0:M1,0:M2,0:M3)
      Real API(M1),ACI(M1),AMI(M1)
      Real APJ(M2),ACJ(M2),AMJ(M2)
      Real APK(M3),ACK(M3),AMK(M3)
      Real R1(M1),R2(M2),R3(M3)

      DO 2 I=1,N1M
          IP=IPA(I)
          IM=IMA(I)
          Do 2 K=1,N3M
              KP=KPA(K)
              KM=KMA(K)
              Do 20 J=1,N2M
                  JP=J+1
                  JM=J-1
                  JUM=J-JMU(J)
                  JUP=JPA(J)-J

                  V2=0.5*(U(2,I,JP,K)+U(2,I,JP,KM))
                  V1=0.5*(U(2,I,J ,K)+U(2,I,J ,KM))

                  APJ(J)=JUP*(
     >            -0.5*HP(J)/RE
     >            +0.5/DY(J)*V2/H(JP)*DY(J)/2.0
     >            )*DT
                  ACJ(J)=1.+(
     >            +0.5*HC(J)/RE
     >            +0.5/DY(J)*(JUP*V2/H(JP)*DY(JP)/2.0
     >            -JUM*V1/H(J)*DY(JM)/2.0)
     >            )*DT
                  AMJ(J)=JUM*(
     >            -0.5*HM(J)/RE
     >            -0.5/DY(J)*V1/H(J)*DY(J)/2.0
     >            )*DT

C     M31UH
                  W2=0.5*(U(3,IP,J,K)+U(3,I ,J,K))
                  W1=0.5*(U(3,I ,J,K)+U(3,IM,J,K))
                  UH2=0.5*(UH(1,IP,J,K)+UH(1,IP,J,KM))
                  UH1=0.5*(UH(1,I ,J,K)+UH(1,I ,J,KM))
                  RM31UH=0.5*DX1*(W2*UH2-W1*UH1)
C     M32VH
                  W2=1.0/H(JP)*(DY(J )/2.*U(3,I,JP,K)+DY(JP)/2.*U(3,I,J ,K))
                  W1=1.0/H(J )*(DY(JM)/2.*U(3,I,J ,K)+DY(J )/2.*U(3,I,JM,K))
                  VH2=0.5*(UH(2,I,JP,K)+UH(2,I,JP,KM))
                  VH1=0.5*(UH(2,I,J ,K)+UH(2,I,J ,KM))
                  RM32VH=JUP*0.5/DY(J)*W2*VH2
     >            -JUM*0.5/DY(J)*W1*VH1

                  R2(J)=DT*(RUH3(I,J,K)-RM31UH-RM32VH)

  20          Continue
              Call TDMA1(AMJ,ACJ,APJ,R2,UH,I,K,3)
  2   CONTINUE

      Do 3 I=1,N1M
          IP=IPA(I)
          IM=IMA(I)
          Do 3 J=1,N2M
              JP=J+1
              JM=J-1
              JUM=J-JMU(J)
              JUP=JPA(J)-J
              Do 30 K=1,N3M
                  KP=KPA(K)
                  KM=KMA(K)

                  W2=0.5*(U(3,I,J,KP)+U(3,I,J,K ))
                  W1=0.5*(U(3,I,J,K )+U(3,I,J,KM))

                  APK(K)=(
     >            -0.5*DX3Q/RE
     >            +DX3*W2/2.0
     >            )*DT
                  ACK(K)=1.+(
     >            +DX3Q/RE
     >            +DX3*(W2/2.0-W1/2.0)
     >            )*DT
                  AMK(K)=(
     >            -0.5*DX3Q/RE
     >            -DX3*W1/2.0
     >            )*DT
                  R3(K)=UH(3,I,J,K)
  30          Continue
              Call CTDMA3(AMK,ACK,APK,R3,UH,I,J,3,N3M)
  3   Continue

      DO 1 J=1,N2M
          Do 1 K=1,N3M
              KP=KPA(K)
              KM=KMA(K)
              Do 10 I=1,N1M
                  IP=IPA(I)
                  IM=IMA(I)

                  U2=0.5*(U(1,IP,J,K)+U(1,IP,J,KM))
                  U1=0.5*(U(1,I ,J,K)+U(1,I ,J,KM))

                  API(I)=(
     >            -0.5*DX1Q/RE
     >            +0.5*DX1*U2/2.0
     >            )*DT
                  ACI(I)=1.+(
     >            +DX1Q/RE
     >            +0.5*DX1*(U2/2.0-U1/2.0)
     >            )*DT
                  AMI(I)=(
     >            -0.5*DX1Q/RE
     >            -0.5*DX1*U1/2.0
     >            )*DT
                  R1(I)=UH(3,I,J,K)
  10          Continue
              Call CTDMA1(AMI,ACI,API,R1,UH,J,K,3,N1M)
  1   Continue

C     DVH UPDATE        
      Do 110 J=2,N2M
          JP=J+1
          JM=J-1
          JUM=J-JMV(J)
          JUP=JPA(J)-J

          Do 110 I=1,N1M
              IP=IPA(I)
              IM=IMA(I)
              Do 110 K=1,N3M
                  KP=KPA(K)
                  KM=KMA(K)

C     M23WH
                  V2=0.5*(U(2,I,J,KP)+U(2,I,J,K ))
                  V1=0.5*(U(2,I,J,K )+U(2,I,J,KM))
                  WH2=1.0/H(J)*(DY(J)/2.*UH(3,I,JM,KP)+DY(JM)/2.*UH(3,I,J,KP))
                  WH1=1.0/H(J)*(DY(J)/2.*UH(3,I,JM,K )+DY(JM)/2.*UH(3,I,J,K ))
                  RM23WH=0.5*DX3*(V2*WH2-V1*WH1)

                  UH(2,I,J,K)=UH(2,I,J,K)-DT*RM23WH
110   CONTINUE

C     DUH UPDATE
      Do 210 J=1,N2M
          JP=J+1
          JM=J-1
          JUM=J-JMU(J)
          JUP=JPA(J)-J
          Do 210 I=1,N1M
              IP=IPA(I)
              IM=IMA(I)
              Do 210 K=1,N3M
                  KP=KPA(K)
                  KM=KMA(K)

C     M12VH
                  U2=1.0/H(JP)*(DY(J )/2.*U(1,I,JP,K)+DY(JP)/2.*U(1,I,J ,K))
                  U1=1.0/H(J )*(DY(JM)/2.*U(1,I,J ,K)+DY(J )/2.*U(1,I,JM,K))
                  VH2=0.5*(UH(2,I,JP,K)+UH(2,IM,JP,K))
                  VH1=0.5*(UH(2,I,J ,K)+UH(2,IM,J ,K))
                  RM12VH=JUP*0.5/DY(J)*U2*VH2
     >            -JUM*0.5/DY(J)*U1*VH1
C     M13WH
                  U2=0.5*(U(1,I,J,KP)+U(1,I,J,K ))
                  U1=0.5*(U(1,I,J,K )+U(1,I,J,KM))
                  WH2=0.5*(UH(3,IM,J,KP)+UH(3,I,J,KP))
                  WH1=0.5*(UH(3,IM,J,K )+UH(3,I,J,K ))
                  RM13WH=0.5*DX3*(U2*WH2-U1*WH1)

                  UH(1,I,J,K)=UH(1,I,J,K)-DT*(RM12VH+RM13WH)
210   CONTINUE

C     INTERMEDIATE VELOCITY UPDATE
      DO 300 I=1,N1M
          DO 300 J=1,N2M
              DO 300 K=1,N3M
                  UH(1,I,J,K)=U(1,I,J,K)+UH(1,I,J,K)
300   UH(3,I,J,K)=U(3,I,J,K)+UH(3,I,J,K)

      DO 310 I=1,N1M
          DO 310 J=2,N2M
              DO 310 K=1,N3M
310   UH(2,I,J,K)=U(2,I,J,K)+UH(2,I,J,K)


      Return
      End

c***************** TDMA1 ***********************     
      Subroutine TDMA1(A,B,C,R,X,I,K,NV)
      Include 'channel.h'
      Common/DIM/N1,N2,N3,N1M,N2M,N3M
      Real A(M2),B(M2),C(M2),R(M2)
      Real X(3,0:M1,0:M2,0:M3)
      Real BET(M2),GAM(M2)

      BET(1)=B(1)
      GAM(1)=C(1)/BET(1)
      X(NV,I,1,K)=R(1)/BET(1)

      Do 10 J=2,N2M
          BET(J)=B(J)-A(J)*GAM(J-1)
          GAM(J)=C(J)/BET(J)
 10   X(NV,I,J,K)=(R(J)-A(J)*X(NV,I,J-1,K))/BET(J)

      Do 20 J=N2M,1,-1
 20   X(NV,I,J,K)=X(NV,I,J,K)-X(NV,I,J+1,K)*GAM(J)

      Return
      End

c***************** TDMA2 ***********************     
      Subroutine TDMA2(A,B,C,R,X,I,K)
      Include 'channel.h'
      Common/DIM/N1,N2,N3,N1M,N2M,N3M
      Real A(M2),B(M2),C(M2),R(M2)
      Real X(3,0:M1,0:M2,0:M3)
      Real BET(M2),GAM(M2)

      GAM(1)=0.
      Do 10 J=2,N2M
          BET(J)=B(J)-A(J)*GAM(J-1)
          GAM(J)=C(J)/BET(J)
 10   X(2,I,J,K)=(R(J)-A(J)*X(2,I,J-1,K))/BET(J)

      Do 20 J=N2M,2,-1
 20   X(2,I,J,K)=X(2,I,J,K)-X(2,I,J+1,K)*GAM(J)

      Return
      End

c***************** CTDMA1 ***********************    
      Subroutine CTDMA1(A,B,C,R,X,J,K,NV,N)
      Include 'channel.h'
      Common/DIM/N1,N2,N3,N1M,N2M,N3M
      Real A(M1),B(M1),C(M1),R(M1)
      Real X(3,0:M1,0:M2,0:M3)
      Real BET(M1),GAM(M1)
      Real P(M1),Q(M1)

      BET(1)=B(1)
      GAM(1)=C(1)/BET(1)
      X(NV,1,J,K)=R(1)/BET(1)
      P(1)=C(N)
      Q(1)=A(1)/BET(1)

      Do 10 I=2,N-1
          BET(I)=B(I)-A(I)*GAM(I-1)
          GAM(I)=C(I)/BET(I)
          P(I)=-P(I-1)*GAM(I-1)
          Q(I)=-A(I)/BET(I)*Q(I-1)
 10   X(NV,I,J,K)=(R(I)-A(I)*X(NV,I-1,J,K))/BET(I)

      P(N-1)=A(N)-P(N-2)*GAM(N-2)
      Q(N-1)=(C(N-1)-A(N-1)*Q(N-2))/BET(N-1)

      X(NV,N,J,K)=R(N)
      P(N)=B(N)
      Do 20 I=1,N-1
          X(NV,N,J,K)=X(NV,N,J,K)-P(I)*X(NV,I,J,K)
 20   P(N)=P(N)-P(I)*Q(I)
      X(NV,N,J,K)=X(NV,N,J,K)/P(N)

      GAM(N-1)=0.0
      Do 30 I=N-1,1,-1
 30   X(NV,I,J,K)=X(NV,I,J,K)-GAM(I)*X(NV,I+1,J,K)-Q(I)*X(NV,N,J,K)

      Return
      End

c***************** CTDMA3 ***********************    
      Subroutine CTDMA3(A,B,C,R,X,I,J,NV,N)
      Include 'channel.h'
      Common/DIM/N1,N2,N3,N1M,N2M,N3M
      Real A(M3),B(M3),C(M3),R(M3)
      Real X(3,0:M1,0:M2,0:M3)
      Real BET(M3),GAM(M3)
      Real P(M3),Q(M3)

      BET(1)=B(1)
      GAM(1)=C(1)/BET(1)
      X(NV,I,J,1)=R(1)/BET(1)
      P(1)=C(N)
      Q(1)=A(1)/BET(1)

      Do 10 K=2,N-1
          BET(K)=B(K)-A(K)*GAM(K-1)
          GAM(K)=C(K)/BET(K)
          P(K)=-P(K-1)*GAM(K-1)
          Q(K)=-A(K)/BET(K)*Q(K-1)
 10   X(NV,I,J,K)=(R(K)-A(K)*X(NV,I,J,K-1))/BET(K)

      P(N-1)=A(N)-P(N-2)*GAM(N-2)
      Q(N-1)=(C(N-1)-A(N-1)*Q(N-2))/BET(N-1)

      X(NV,I,J,N)=R(N)
      P(N)=B(N)
      Do 20 K=1,N-1
          X(NV,I,J,N)=X(NV,I,J,N)-P(K)*X(NV,I,J,K)
 20   P(N)=P(N)-P(K)*Q(K)
      X(NV,I,J,N)=X(NV,I,J,N)/P(N)

      GAM(N-1)=0.0
      Do 30 K=N-1,1,-1
 30   X(NV,I,J,K)=X(NV,I,J,K)-GAM(K)*X(NV,I,J,K+1)-Q(K)*X(NV,I,J,N)
      Return
      End


c***************** DIVCHECK ***********************    
      Subroutine DIVCHECK(U,TIME,DIVMAX)
      Include 'channel.h'
      Common/DIM/N1,N2,N3,N1M,N2M,N3M
      Common/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      Common/INDX2/JPA(M2),JMU(M2),JMV(M2)
      Common/MESH1/DX1,DX1Q,DX3,DX3Q
      Common/MESH2/Y(0:M2)
      Common/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      Common/NSTEP/NTST
      COMMON/TSTEP/DT,CFLMAX

      Real U(3,0:M1,0:M2,0:M3)

      DIVMAX=0.0

      Do 20 I=1,N1M
          IP=IPA(I)
          IM=IMA(I)
          Do 20 J=1,N2M
              JP=J+1
              JM=J-1
              JUM=J-JMU(J)
              JUP=JPA(J)-J
              Do 20 K=1,N3M
                  KP=KPA(K)
                  KM=KMA(K)
                  DIV=(U(1,IP,J,K)-U(1,I,J,K))*DX1
     >            +(U(2,I,JP,K)-U(2,I,J,K))/DY(J)
     >            +(U(3,I,J,KP)-U(3,I,J,K))*DX3
                  DIVMAX=AMAX1(DIV,DIVMAX)
  20  Continue

      Return
      End

c*********************** CFL ***********************
c     This subroutine calculate the maximum local CFL number
c     devided by DT
c     AT THE CELL CENTER
c
      Subroutine CFL(U,CFLM)
      Include 'channel.h'
      Common/DIM/N1,N2,N3,N1M,N2M,N3M
      Common/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      Common/INDX2/JPA(M2),JMU(M2),JMV(M2)
      Common/MESH1/DX1,DX1Q,DX3,DX3Q
      Common/MESH2/Y(0:M2)
      Common/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      Common/NSTEP/NTST
      COMMON/TSTEP/DT,CFLMAX

      Real U(3,0:M1,0:M2,0:M3)
      CFLM=0.0
      DO 10 K=1,N3M
          KP=KPA(K)
          DO 10 I=1,N1M
              IP=IPA(I)
              DO 10 J=1,N2M
                  JP=J+1
                  CFLL=ABS(U(1,I,J,K)+U(1,IP,J,K))*0.5*DX1
     >            +ABS(U(2,I,J,K)+U(2,I,JP,K))*0.5/DY(J)
     >            +ABS(U(3,I,J,K)+U(3,I,J,KP))*0.5*DX3
                  CFLM=AMAX1(CFLM,CFLL)
 10   CONTINUE

      RETURN
      END

C-JIChoi 99/10/30 X1/X3 FOURIER COMPUTATIONS FOR PERIODIC CHANNEL FLOW
C**********************************************************************
C
C     POISSON EQUATION IS SOLVED BY
C     1. FFT IN X1 & X3-DIRECTION
C     3. TRIDIAGONAL SOLVER IN X2-DIRECTION
C
C
C   ********************* INIWAVE ******************************
C
C     GET THE MODIFIED WAVENUMBERS AND THE LOCATION
      SUBROUTINE INIWAVE

      Include 'channel.h'
      PARAMETER (M3M=M3-1,M3MH=M3M/2+1)

      Common/DIM/N1,N2,N3,N1M,N2M,N3M
      Common/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      Common/INDX2/JPA(M2),JMU(M2),JMV(M2)
      Common/MESH1/DX1,DX1Q,DX3,DX3Q
      Common/MESH2/Y(0:M2)
      Common/DAK/N3MH,N3MD,N3MDU,N3MU
      Common/WAVK13/AK1(M1),AK3(M3MH)

      REAL ANI(M1),ANK(M3)

      N3MH=N3M/2+1
      PI=ACOS(-1.0)

C     MODIFIED WAVE NUMBER DEFINITION NECESSARY FOR THE X1 DIRECTION
C     IT MUST BE CONSIDERED ABOUT FREQUENCY SHIFTING IN X1 DIRECTION
      DO 10 I=1,N1M
          II=I
          IF(I.GT.N1M/2+1) II=I-N1M
   10 ANI(I)=(II-1)*2.0*PI
      WRITE(6,764) (ANI(I),I=1,N1M)
  764 FORMAT(1X,'ANI',2X,11F10.3)

      DO 11 II=1,N1M
  11  AK1(II)=2.*(1.-COS(ANI(II)/N1M))*DX1Q
      WRITE(6,764) (AK1(I),I=1,N1M)

C     MODIFIED WAVE NUMBER DEFINITION NECESSARY FOR THE X3 DIRECTION
      DO 20 K=1,N3MH
   20 ANK(K)=(K-1)*2.*PI
      WRITE(6,765) (ANK(K),K=1,N3MH)
  765 FORMAT(1X,'ANK',2X,11F10.3)

      DO 21 KK=1,N3MH
  21  AK3(KK)=2.*(1.-COS(ANK(KK)/N3M))*DX3Q
      WRITE(6,765) (AK3(K),K=1,N3MH)

      CALL METRICPOISSON

      RETURN
      END

C  ************************ METRICPOISSON  **********************
C
C     CALCULATE THE COEFFICIENTS OF THE POISSON EQUATION FOR DPH.
C     COEFFICIENTS FOR THE THREE POINTS STENCIL OF THE POISSON EQ.
C

      SUBROUTINE METRICPOISSON

      Include 'channel.h'
      Common/DIM/N1,N2,N3,N1M,N2M,N3M
      Common/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      Common/INDX2/JPA(M2),JMU(M2),JMV(M2)
      Common/MESH1/DX1,DX1Q,DX3,DX3Q
      Common/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      Common/METPOI/PMJ(M2),PCJ(M2),PPJ(M2)

      DO 100 JC=1,N2M
          JP=JC+1
          JM=JC-1
          JUM=JC-JMU(JC)
          JUP=JPA(JC)-JC
          PMJ(JC)=JUM*JUP*(1.0/DY(JC)/H(JC))+(1-JUP)*(1.0/DY(N2M)/H(N2M))
          PCJ(JC)=JUM*JUP*(-1.0/DY(JC)/H(JC)-1.0/DY(JC)/H(JP))
     >    +(1-JUP)*(-1.0/DY(N2M)/H(N2M))
     >    +(1-JUM)*(-1.0/DY(1)/H(2))
          PPJ(JC)=JUM*JUP*(1.0/DY(JC)/H(JP))+(1-JUM)*(1.0/DY(1)/H(2))
  100 CONTINUE

      RETURN
      END

      Include 'fft4f2d.f'   ! For fast transform
C   ********************* GETDP ********************************
C     CALCULATE DP.
C     PERIODIC DIRECTION ALONG X1 & X3 ALLOWS
C     TO USE THE COMPLEX FOURIER TRANSFORM.

      Subroutine GETDP(DP,RDP)

      Include 'channel.h'
      PARAMETER (M1M=M1-1,M3M=M3-1,M3MH=M3M/2+1,M3MD=M3M+2)

      Common/DIM/N1,N2,N3,N1M,N2M,N3M
      Common/WAVK13/AK1(M1),AK3(M3MH)
      Common/METPOI/PMJ(M2),PCJ(M2),PPJ(M2)

C---- FFT4F2D CDFT2D PARAMETERS --------------------------------------
      COMPLEX COEF(N1M,N3M)
      REAL  PF(2*M1M,M3M)
      REAL*8 PF_TEMP(0:2*M1M-1,0:M3M-1)
      INTEGER IPF(M1M*2), IPF_TEMP(0:M1M*2-1)
      REAL WF(M1M*2)
      REAL*8 WF_TEMP(0:M1M*2-1)
C---------------------------------------------------------------------
      REAL RDP(M1,M2,M3),DP(M1,M2,M3)
      REAL FRDP(M1,M2,M3MD),FDP(M1,M2,M3MD)
      REAL CMJ(M2),CCJ(M2),CPJ(M2),CFJ(M2)

      N1MH=N1M/2+1
      N3MH=N3M/2+1

      DO 10 J=1,N2M
          DO 1 I=1,N1M
              DO 1 K=1,N3M
                  PF(2*I-1,K)=RDP(I,J,K)
                  PF(2*I  ,K)=0.
    1     CONTINUE

          IPF(1)=0    ! Initialize
          IPF_TEMP = IPF
          PF_TEMP  = PF
          WF_TEMP  = WF
          CALL CDFT2D(2*M1M,2*M1M,M3M,1,PF_TEMP,IPF_TEMP,WF_TEMP)
          PF = PF_TEMP
          WF = WF_TEMP
          IPF= IPF_TEMP

          RATIO=DBLE(N1M*N3M)
          DO I=1,N1M
              DO K=1,N3M
                  COEF(I,K)=CMPLX(PF(2*I-1,K),PF(2*I,K))/RATIO
              ENDDO
          ENDDO

          DO 11 I=1,N1M
              DO 11 K=1,N3MH
                  KR=2*K-1
                  KI=2*K
                  FRDP(I,J,KR)=DBLE(COEF(I,K))
                  FRDP(I,J,KI)=IMAG(COEF(I,K))
   11     CONTINUE

   10 CONTINUE

C     SOLVE FDP BY USING TDMA IN X2-DIRECTION
C     SPLIT REAL AND IMAGINARY PART
C---- REAL PART CALCULATION
      DO 110 I=1,N1M
          DO 110 K=1,N3MH
              KR=2*K-1
              DO 120 J=1,N2M
                  CMJ(J)=PMJ(J)
                  CCJ(J)=PCJ(J)-AK3(K)-AK1(I)
                  CPJ(J)=PPJ(J)
                  CFJ(J)=FRDP(I,J,KR)
                  IF(I.EQ.1.AND.K.EQ.1) THEN
                      CMJ(1)=0.0
                      CCJ(1)=1.0
                      CPJ(1)=0.0
                      CFJ(1)=0.0
                  ENDIF
  120         CONTINUE
              CALL TDMAP(CMJ,CCJ,CPJ,CFJ,N2M,CFJ)

              DO 130 J=1,N2M
                  FDP(I,J,KR)=CFJ(J)
  130         CONTINUE

  110 CONTINUE

C---- IMAG PART CALCULATION
      DO 140 I=1,N1M
          DO 140 K=1,N3MH
              KI=2*K
              DO 150 J=1,N2M
                  CMJ(J)=PMJ(J)
                  CCJ(J)=PCJ(J)-AK3(K)-AK1(I)
                  CPJ(J)=PPJ(J)
                  CFJ(J)=FRDP(I,J,KI)
                  IF(((I.EQ.1).AND.(K.EQ.1)).OR.((I.EQ.1).AND.(K.EQ.N3MH)).OR.
     >            ((I.EQ.N1MH).AND.(K.EQ.1)).OR.((I.EQ.N1MH).AND.(K.EQ.N3MH))) THEN
                      CMJ(J)=0.0
                      CCJ(J)=1.0
                      CPJ(J)=0.0
                      CFJ(J)=0.0
                  ENDIF
  150         CONTINUE
              CALL TDMAP(CMJ,CCJ,CPJ,CFJ,N2M,CFJ)

              DO 160 J=1,N2M
                  FDP(I,J,KI)=CFJ(J)
  160         CONTINUE

  140 CONTINUE

C     DO THE INVERSE FFT.
      DO 20 J=1,N2M
          DO 21 I=1,N1M
              DO 21 K=1,N3MH
                  KR=2*K-1
                  KI=2*K
                  COEF(I,K)=CMPLX(FDP(I,J,KR),FDP(I,J,KI))
   21     CONTINUE

          DO 22 I=1,N1M
              DO 22 K=N3MH+1,N3M
                  II=N1M-I+2
                  KK=N3M-K+2
                  IF(I.EQ.1) II=1
                  COEF(I,K)=CONJG(COEF(II,KK))
   22     CONTINUE

          DO 23 I=1,N1M
              DO 23 K=1,N3M
                  PF(2*I-1,K)=DBLE(COEF(I,K))
                  PF(2*I  ,K)=IMAG(COEF(I,K))
   23     CONTINUE

          IPF(1)=0    ! Initialize
          IPF_TEMP = IPF
          PF_TEMP  = PF
          WF_TEMP  = WF
          CALL CDFT2D(2*M1M,2*M1M,M3M,-1,PF_TEMP,IPF_TEMP,WF_TEMP)
          IPF = IPF_TEMP
          PF  = PF_TEMP
          WF  = WF_TEMP

          DO 24 I=1,N1M
              DO 24 K=1,N3M
                  DP(I,J,K)=PF(2*I-1,K)
   24     CONTINUE

   20 CONTINUE

      RETURN
      END

C  ****************************** TDMAP **********************
      SUBROUTINE TDMAP(A,B,C,R,N,X)
      Include 'channel.h'
      Real GAM(M2),A(M2),B(M2),C(M2),R(M2),X(M2)

      BET=B(1)
      X(1)=R(1)/BET
      DO 11 J=2,N
          GAM(J)=C(J-1)/BET
          BET=B(J)-A(J)*GAM(J)
          X(J)=(R(J)-A(J)*X(J-1))/BET
   11 CONTINUE
      DO 12 J=N-1,1,-1
          X(J)=X(J)-GAM(J+1)*X(J+1)
   12 CONTINUE
      RETURN
      END


c***************** CHKMF ***********************    
      Subroutine CHKMF(U,TIME)
      Include 'channel.h'
      Common/DIM/N1,N2,N3,N1M,N2M,N3M
      Common/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      Common/INDX2/JPA(M2),JMU(M2),JMV(M2)
      Common/MESH1/DX1,DX1Q,DX3,DX3Q
      Common/MESH2/Y(0:M2)
      Common/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      Common/SIZE/ALX,ALY,ALZ,VOL
      Common/NSTEP/NTST
      COMMON/TSTEP/DT,CFLMAX

      Real U(3,0:M1,0:M2,0:M3)

C     Calculate the mass flow rate in x1 direction      
      FLOW1=0.0
      DO 1 K=1,N3M
          DO 1 J=1,N2M
              DO 1 I=1,N1M
                  FLOW1=FLOW1
     >            +U(1,I,J,K)*DY(J)/DX3/DX1
    1 CONTINUE
      FLOW1=FLOW1/ALX/ALZ    ! NORMALIZED LENGTH IS Y

C     Calculate the mass flow rate in x3 direction      
      FLOW3=0.0
      DO 2 I=1,N1M
          DO 2 J=1,N2M
              DO 2 K=1,N3M
                  FLOW3=FLOW3
     >            +U(3,I,J,1)*DY(J)/DX1/DX3
    2 CONTINUE
      FLOW3=FLOW3/ALX/ALZ

      WRITE(32,100) TIME,FLOW1,FLOW3
 100  Format(3(E12.5,2X))

      Return
      End

C***************** PROFILE ***********************    
C     SPATIALLY AVERAGED FLOW FILED
C     AT A CERTAIN INSTANTANEOUS TIME

      Subroutine PROFILE(U,P)
      Include 'channel.h'
      Common/DIM/N1,N2,N3,N1M,N2M,N3M
      Common/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      Common/INDX2/JPA(M2),JMU(M2),JMV(M2)
      Common/MESH1/DX1,DX1Q,DX3,DX3Q
      Common/MESH2/Y(0:M2)
      Common/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      Common/SIZE/ALX,ALY,ALZ,VOL
      Common/NSTEP/NTST
      COMMON/TSTEP/DT,CFLMAX

      Real U(3,0:M1,0:M2,0:M3)
      Real P(M1,M2,M3)
      REAL NXZ

      Real UM(0:M2),VM(M2),WM(0:M2),PM(0:M2)

      NXZ=1.0/DBLE(N1M*N3M)


      DO 1 J=1,N2M
          UM(J)=0.0
          WM(J)=0.0
          PM(J)=0.0
          DO 10 I=1,N1M
              DO 10 K=1,N3M
                  UM(J)=UM(J)+U(1,I,J,K)
                  WM(J)=WM(J)+U(3,I,J,K)
                  PM(J)=PM(J)+P(I,J,K)
   10     CONTINUE
          UM(J)=UM(J)*NXZ
          WM(J)=WM(J)*NXZ
          PM(J)=PM(J)*NXZ
    1 CONTINUE

      DO 2 J=1,N2
          VM(J)=0.0
          DO 20 I=1,N1M
              DO 20 K=1,N3M
                  VM(J)=VM(J)+U(2,I,J,K)
   20     CONTINUE
          VM(J)=VM(J)*NXZ
    2 CONTINUE

      UM(0)=0.0
      WM(0)=0.0
      PM(0)=PM(1)
      PM(N2)=PM(N2M)

      OPEN(36,FILE='U.plt',STATUS='UNKNOWN')
      OPEN(37,FILE='V.plt',STATUS='UNKNOWN')
      OPEN(38,FILE='W.plt',STATUS='UNKNOWN')
      OPEN(39,FILE='P.plt',STATUS='UNKNOWN')
      X2=0.0
      WRITE(36,100) X2,UM(0),U(1,1,0,1)
      WRITE(38,100) X2,WM(0),U(3,1,0,1)
      DO 30 J=1,N2
          X2=X2+H(J)
          WRITE(36,100) X2,UM(J),U(1,1,J,1)
          WRITE(38,100) X2,WM(J),U(3,1,J,1)
          IF (J.NE.N2) WRITE(39,100) X2,PM(J),P(1,J,1)
   30 CONTINUE

      DO 40 J=1,N2
          WRITE(37,100) Y(J),VM(J),U(2,1,J,1)
   40 CONTINUE

  100 FORMAT(3(E12.5,2X))

      CLOSE(36)
      CLOSE(37)
      CLOSE(38)
      CLOSE(39)

      RETURN
      END


C  ***************************** INSFIELD **********************

      SUBROUTINE INSFIELD(U,P)

      Include 'channel.h'
      Common/DIM/N1,N2,N3,N1M,N2M,N3M
      Common/PARA/RE
      Common/MESH1/DX1,DX1Q,DX3,DX3Q
      Common/MESH2/Y(0:M2)
      Common/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      Common/NSTEP/NTST
      COMMON/TSTEP/DT,CFLMAX
      Common/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      Common/INDX2/JPA(M2),JMU(M2),JMV(M2)
      COMMON/VMEAN/VM(3,M2),VRMS(6,M2)
      COMMON/WMEAN/WM(3,M2),WRMS(3,M2)
      COMMON/PMEAN/PM(M2),PRMS(M2)
      COMMON/WSMEM/WSM(3),WSMO(3),ENEJ,ENEJO
      Common/SIZE/ALX,ALY,ALZ,VOL

      Real U(3,0:M1,0:M2,0:M3),P(M1,M2,M3)
      Real VOR(3,M1,M2,M3)
      Real WS(M1,2,M3)


      OPEN(20,FILE='INSYZ.plt',STATUS='UNKNOWN')
      OPEN(21,FILE='INSXZ.plt',STATUS='UNKNOWN')
      OPEN(22,FILE='INSXY.plt',STATUS='UNKNOWN')
      OPEN(23,FILE='INSXZW.plt',STATUS='UNKNOWN')

C     CALCULATE VORTICITY at Cell Center Points
      DO 10 K=1,N3M
          KP=KPA(K)
          KM=KMA(K)
          DO 10 J=1,N2M
              JP=J+1
              JM=J-1
              JUM=J-JMU(J)
              JUP=JPA(J)-J
              DO 10 I=1,N1M
                  IP=IPA(I)
                  IM=IMA(I)

                  U1KP=0.5*(U(1,I,J,KP)+U(1,IP,J,KP))
                  U1KM=0.5*(U(1,I,J,KM)+U(1,IP,J,KM))
                  U2KP=0.5*(U(2,I,J,KP)+U(2,I,JP,KP))
                  U2KM=0.5*(U(2,I,J,KM)+U(2,I,JP,KM))
                  U2IP=0.5*(U(2,IP,J,K)+U(2,IP,JP,K))
                  U2IM=0.5*(U(2,IM,J,K)+U(2,IM,JP,K))
                  U3IP=0.5*(U(3,IP,J,K)+U(3,IP,J,KP))
                  U3IM=0.5*(U(3,IM,J,K)+U(3,IM,J,KP))

                  U1JP=0.5*(U(1,I,JP,K)+U(1,IP,JP,K))
                  U1JC=0.5*(U(1,I,J ,K)+U(1,IP,J ,K))
                  U1JM=0.5*(U(1,I,JM,K)+U(1,IP,JM,K))
                  U12=0.5/H(JP)*(DY(JP)*U1JC+DY(J)*U1JP)
                  U11=0.5/H(J )*(DY(J)*U1JM+DY(JM)*U1JC)
                  U12=U12*JUP+(1-JUP)*U(1,I,N2,K)
                  U11=U11*JUM+(1-JUM)*U(1,I,0,K)

                  U3JP=0.5*(U(3,I,JP,K)+U(3,I,JP,KP))
                  U3JC=0.5*(U(3,I,J ,K)+U(3,I,J ,KP))
                  U3JM=0.5*(U(3,I,JM,K)+U(3,I,JM,KP))
                  U32=0.5/H(JP)*(DY(JP)*U3JC+DY(J)*U3JP)
                  U31=0.5/H(J )*(DY(J)*U3JM+DY(JM)*U3JC)
                  U32=U32*JUP+(1-JUP)*U(3,I,N2,K)
                  U31=U31*JUM+(1-JUM)*U(3,I,0,K)

                  DV3DX2=(U32-U31)/DY(J)
                  DV2DX3=(U2KP-U2KM)*0.5*DX3
                  DV3DX1=(U3IP-U3IM)*0.5*DX1
                  DV1DX3=(U1KP-U1KM)*0.5*DX3
                  DV2DX1=(U2IP-U2IM)*0.5*DX1
                  DV1DX2=(U12-U11)/DY(J)

                  VOR(1,I,J,K)=DV3DX2-DV2DX3
                  VOR(2,I,J,K)=DV1DX3-DV3DX1
                  VOR(3,I,J,K)=DV2DX1-DV1DX2
   10 CONTINUE

      DO 15 I=1,N1M
          IP=IPA(I)
          IM=IMA(I)
          DO 15 K=1,N3M
              UHC=0.5*(U(1,I,1,K)+U(1,IP,1,K))     ! Bottom Wall
              UHW=0.5*(U(1,I,0,K)+U(1,IP,0,K))
              WS(I,1,K)=(UHC-UHW)/(DY(1)*0.5)

              UHC=0.5*(U(1,I,N2M,K)+U(1,IP,N2M,K)) ! Top Wall
              UHW=0.5*(U(1,I,N2 ,K)+U(1,IP,N2 ,K))
              WS(I,2,K)=(UHC-UHW)/(DY(N2M)*0.5)
   15 CONTINUE

C     YZ Plane view 
      Write(20,201) N3M,N2M
      NX1=N1M/4*0+1
      NX2=N1M/4*1+1
      NX3=N1M/4*2+1
      NX4=N1M/4*3+1
      DO 20 J=1,N2M
          X2=Y(J)+DY(J)*0.5
          DO 20 K=1,N3M
              X3=(Real(K-1)+0.5)/DX3
              WRITE(20,202) X3,X2,VOR(1,NX1,J,K),VOR(1,NX2,J,K)
     >        ,VOR(1,NX3,J,K),VOR(1,NX4,J,K)
   20 CONTINUE

C     XZ Plane view 
      Write(21,211) N3M,N1M
      Write(23,211) N3M,N1M
      NY1=5
      NY2=19
      NY3=33
      NY4=49
      DO 30 I=1,N1M
          X1=(Real(I-1)+0.5)/DX1
          DO 30 K=1,N3M
              X3=(Real(K-1)+0.5)/DX3
              WRITE(21,212) X1,X3,VOR(1,I,NY1,K),VOR(1,I,NY2,K)
     >        ,VOR(1,I,NY3,K),VOR(1,I,NY4,K)
              WRITE(23,232) X1,X3,WS(I,1,K),WS(I,2,K)
     >        ,P(I,1,K),P(I,N2M,K)
   30 CONTINUE

C     XY Plane View
      Write(22,221) N1M,N2M
      NZ1=N3M/4*0+1
      NZ2=N3M/4*1+1
      NZ3=N3M/4*2+1
      NZ4=N3M/4*3+1
      DO 40 J=1,N2M
          X2=Y(J)+DY(J)*0.5
          DO 40 I=1,N1M
              X1=(Real(I-1)+0.5)/DX1
              WRITE(22,222) X1,X2,VOR(1,I,J,NZ1),VOR(1,I,J,NZ2)
     >        ,VOR(1,I,J,NZ3),VOR(1,I,J,NZ4)
   40 CONTINUE

  201 FORMAT('Zone i=',I4,'j=',I4,'f=point')
  202 FORMAT(6(E12.4,2X))
  211 FORMAT('Zone i=',I4,'j=',I4,'f=point')
  212 FORMAT(6(E12.4,2X))
  221 FORMAT('Zone i=',I4,'j=',I4,'f=point')
  222 FORMAT(6(E12.4,2X))
  231 FORMAT('Zone i=',I4,'j=',I4,'f=point')
  232 FORMAT(6(E12.4,2X))

      CLOSE(20)
      CLOSE(21)
      CLOSE(22)
      CLOSE(23)

      RETURN
      END

C  ****************************** ENERGY **********************
C     CALCULATE TOTAL ENERGY = (V1**2+V2**2+V3**2)*0.5.
C     CALCULATE MEAN AND RMS VALUES AT THE CENTER OF CELL.

      SUBROUTINE ENERGY(U,P)

      Include 'channel.h'
      Common/DIM/N1,N2,N3,N1M,N2M,N3M
      Common/PARA/RE
      Common/MESH1/DX1,DX1Q,DX3,DX3Q
      Common/MESH2/Y(0:M2)
      Common/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      Common/NSTEP/NTST
      COMMON/TSTEP/DT,CFLMAX
      Common/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      Common/INDX2/JPA(M2),JMU(M2),JMV(M2)
      COMMON/VMEAN/VM(3,M2),VRMS(6,M2)
      COMMON/WMEAN/WM(3,M2),WRMS(3,M2)
      COMMON/PMEAN/PM(M2),PRMS(M2)
      COMMON/WSMEM/WSM(3),WSMO(3),ENEJ,ENEJO
      Common/SIZE/ALX,ALY,ALZ,VOL

      Real U(3,0:M1,0:M2,0:M3),P(M1,M2,M3)
      Real VOR(3,M1,M2,M3)
      REAL NXZ

      NXZ=1.0/DBLE(N1M*N3M)

C     CALCULATE MEAN Velocity at Cell Center Points.
      DO 10 J=1,N2M
          JP=J+1
          V1M=0.0
          V2M=0.0
          V3M=0.0
          PMM=0.0
          DO 20 I=1,N1M
              IP=IPA(I)
              DO 20 K=1,N3M
                  KP=KPA(K)
                  V1M=V1M+(U(1,I,J,K)+U(1,IP,J,K))*0.5
                  V2M=V2M+(U(2,I,J,K)+U(2,I,JP,K))*0.5
                  V3M=V3M+(U(3,I,J,K)+U(3,I,J,KP))*0.5
                  PMM=PMM+P(I,J,K)
   20     CONTINUE
          VM(1,J)=V1M*NXZ
          VM(2,J)=V2M*NXZ
          VM(3,J)=V3M*NXZ
          PM(J)=PMM*NXZ
   10 CONTINUE

C     CALCULATE Mean Square VALUES at Cell Center Points.
      ENEJ=0.
      DO 30 J=1,N2M
          JP=J+1
          U11MS=0.0
          U22MS=0.0
          U33MS=0.0
          U12MS=0.0
          U13MS=0.0
          U23MS=0.0
          PMMMS=0.0
          DO 40 I=1,N1M
              IP=IPA(I)
              DO 40 K=1,N3M
                  KP=KPA(K)
                  V1=(U(1,I,J,K)+U(1,IP,J,K))*0.5
                  V2=(U(2,I,J,K)+U(2,I,JP,K))*0.5
                  V3=(U(3,I,J,K)+U(3,I,J,KP))*0.5
                  U11MS=U11MS+(V1-VM(1,J))**2
                  U22MS=U22MS+(V2-VM(2,J))**2
                  U33MS=U33MS+(V3-VM(3,J))**2
                  U12MS=U12MS+(V1-VM(1,J))*(V2-VM(2,J))
                  U13MS=U13MS+(V1-VM(1,J))*(V3-VM(3,J))
                  U23MS=U23MS+(V2-VM(2,J))*(V3-VM(3,J))
                  PMMMS=PMMMS+(P(I,J,K)-PM(J))**2
   40     CONTINUE
          VRMS(1,J)=U11MS*NXZ
          VRMS(2,J)=U22MS*NXZ
          VRMS(3,J)=U33MS*NXZ
          VRMS(4,J)=-U12MS*NXZ
          VRMS(5,J)=-U13MS*NXZ
          VRMS(6,J)=-U23MS*NXZ
          PRMS(J)=PMMMS*NXZ
          ENEJ=ENEJ+(VRMS(1,J)+VRMS(2,J)+VRMS(3,J))*DY(J)
   30 CONTINUE
      ENEJ=ENEJ/ALY

C     CALCULATE VORTICITY at Cell Center Points
      DO 100 K=1,N3M
          KP=KPA(K)
          KM=KMA(K)
          DO 100 J=1,N2M
              JP=J+1
              JM=J-1
              JUM=J-JMU(J)
              JUP=JPA(J)-J
              DO 100 I=1,N1M
                  IP=IPA(I)
                  IM=IMA(I)

                  U1KP=0.5*(U(1,I,J,KP)+U(1,IP,J,KP))
                  U1KM=0.5*(U(1,I,J,KM)+U(1,IP,J,KM))
                  U2KP=0.5*(U(2,I,J,KP)+U(2,I,JP,KP))
                  U2KM=0.5*(U(2,I,J,KM)+U(2,I,JP,KM))
                  U2IP=0.5*(U(2,IP,J,K)+U(2,IP,JP,K))
                  U2IM=0.5*(U(2,IM,J,K)+U(2,IM,JP,K))
                  U3IP=0.5*(U(3,IP,J,K)+U(3,IP,J,KP))
                  U3IM=0.5*(U(3,IM,J,K)+U(3,IM,J,KP))

                  U1JP=0.5*(U(1,I,JP,K)+U(1,IP,JP,K))
                  U1JC=0.5*(U(1,I,J ,K)+U(1,IP,J ,K))
                  U1JM=0.5*(U(1,I,JM,K)+U(1,IP,JM,K))
                  U12=0.5/H(JP)*(DY(JP)*U1JC+DY(J)*U1JP)
                  U11=0.5/H(J )*(DY(J)*U1JM+DY(JM)*U1JC)
                  U12=U12*JUP+(1-JUP)*U(1,I,N2,K)
                  U11=U11*JUM+(1-JUM)*U(1,I,0,K)

                  U3JP=0.5*(U(3,I,JP,K)+U(3,I,JP,KP))
                  U3JC=0.5*(U(3,I,J ,K)+U(3,I,J ,KP))
                  U3JM=0.5*(U(3,I,JM,K)+U(3,I,JM,KP))
                  U32=0.5/H(JP)*(DY(JP)*U3JC+DY(J)*U3JP)
                  U31=0.5/H(J )*(DY(J)*U3JM+DY(JM)*U3JC)
                  U32=U32*JUP+(1-JUP)*U(3,I,N2,K)
                  U31=U31*JUM+(1-JUM)*U(3,I,0,K)

                  DV3DX2=(U32-U31)/DY(J)
                  DV2DX3=(U2KP-U2KM)*0.5*DX3
                  DV3DX1=(U3IP-U3IM)*0.5*DX1
                  DV1DX3=(U1KP-U1KM)*0.5*DX3
                  DV2DX1=(U2IP-U2IM)*0.5*DX1
                  DV1DX2=(U12-U11)/DY(J)

                  VOR(1,I,J,K)=DV3DX2-DV2DX3
                  VOR(2,I,J,K)=DV1DX3-DV3DX1
                  VOR(3,I,J,K)=DV2DX1-DV1DX2
  100 CONTINUE

C     VORTICITY Tilting and Stretching must be considered!!!

C     CALCULATE MEAN Vorticity at Cell Center Points.
      DO 110 J=1,N2M
          JP=J+1
          VM1M=0.0
          VM2M=0.0
          VM3M=0.0
          DO 120 I=1,N1M
              DO 120 K=1,N3M
                  VM1M=VM1M+VOR(1,I,J,K)
                  VM2M=VM2M+VOR(2,I,J,K)
                  VM3M=VM3M+VOR(3,I,J,K)
  120     CONTINUE
          WM(1,J)=VM1M*NXZ
          WM(2,J)=VM2M*NXZ
          WM(3,J)=VM3M*NXZ
  110 CONTINUE

C     CALCULATE MS Voticity VALUES at Cell Center Points.
      DO 130 J=1,N2M
          VOR1MS=0.0
          VOR2MS=0.0
          VOR3MS=0.0
          DO 140 I=1,N1M
              DO 140 K=1,N3M
                  VOR1MS=VOR1MS+(VOR(1,I,J,K)-WM(1,J))**2
                  VOR2MS=VOR2MS+(VOR(2,I,J,K)-WM(2,J))**2
                  VOR3MS=VOR3MS+(VOR(3,I,J,K)-WM(3,J))**2
  140     CONTINUE
          WRMS(1,J)=VOR1MS*NXZ
          WRMS(2,J)=VOR2MS*NXZ
          WRMS(3,J)=VOR3MS*NXZ
  130 CONTINUE

      RETURN
      END


C  ****************************** TAVER **********************
C     CALCULATE TIME AVERAGED MEAN QUANTITIES :

      SUBROUTINE TAVER
      Include 'channel.h'
      Common/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/VMEAN/VM(3,M2),VRMS(6,M2)
      COMMON/WMEAN/WM(3,M2),WRMS(3,M2)
      COMMON/PMEAN/PM(M2),PRMS(M2)
      COMMON/VMEANO/VMO(3,M2),VRMSO(6,M2)
      COMMON/WMEANO/WMO(3,M2),WRMSO(3,M2)
      COMMON/PMEANO/PMO(M2),PRMSO(M2)
      COMMON/WSMEM/WSM(3),WSMO(3),ENEJ,ENEJO

      ENEJO=ENEJO+ENEJ

      DO 10 J=1,N2M
          VMO(1,J)=(VMO(1,J)+VM(1,J))
          VMO(2,J)=(VMO(2,J)+VM(2,J))
          VMO(3,J)=(VMO(3,J)+VM(3,J))
          VRMSO(1,J)=(VRMSO(1,J)+VRMS(1,J))
          VRMSO(2,J)=(VRMSO(2,J)+VRMS(2,J))
          VRMSO(3,J)=(VRMSO(3,J)+VRMS(3,J))
          VRMSO(4,J)=(VRMSO(4,J)+VRMS(4,J))
          VRMSO(5,J)=(VRMSO(5,J)+VRMS(5,J))
          VRMSO(6,J)=(VRMSO(6,J)+VRMS(6,J))
          WMO(1,J)=(WMO(1,J)+WM(1,J))
          WMO(2,J)=(WMO(2,J)+WM(2,J))
          WMO(3,J)=(WMO(3,J)+WM(3,J))
          WRMSO(1,J)=(WRMSO(1,J)+WRMS(1,J))
          WRMSO(2,J)=(WRMSO(2,J)+WRMS(2,J))
          WRMSO(3,J)=(WRMSO(3,J)+WRMS(3,J))
          PMO(J)=(PMO(J)+PM(J))
          PRMSO(J)=(PRMSO(J)+PRMS(J))
   10 CONTINUE

      DO 20 L=1,3
          WSMO(L)=WSMO(L)+WSM(L)
   20 CONTINUE
      RETURN
      END

C  ****************************** SAVER **********************
C     THIS ROUTINE IS ONLY FOR FIRST TIME STEP

      SUBROUTINE SAVER
      Include 'channel.h'
      Common/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/VMEAN/VM(3,M2),VRMS(6,M2)
      COMMON/WMEAN/WM(3,M2),WRMS(3,M2)
      COMMON/PMEAN/PM(M2),PRMS(M2)
      COMMON/VMEANO/VMO(3,M2),VRMSO(6,M2)
      COMMON/WMEANO/WMO(3,M2),WRMSO(3,M2)
      COMMON/PMEANO/PMO(M2),PRMSO(M2)
      COMMON/WSMEM/WSM(3),WSMO(3),ENEJ,ENEJO

      ENEJO=ENEJ

      DO 10 J=1,N2M
          VMO(1,J)=VM(1,J)
          VMO(2,J)=VM(2,J)
          VMO(3,J)=VM(3,J)
          VRMSO(1,J)=VRMS(1,J)
          VRMSO(2,J)=VRMS(2,J)
          VRMSO(3,J)=VRMS(3,J)
          VRMSO(4,J)=VRMS(4,J)
          VRMSO(5,J)=VRMS(5,J)
          VRMSO(6,J)=VRMS(6,J)
          WMO(1,J)=WM(1,J)
          WMO(2,J)=WM(2,J)
          WMO(3,J)=WM(3,J)
          WRMSO(1,J)=WRMS(1,J)
          WRMSO(2,J)=WRMS(2,J)
          WRMSO(3,J)=WRMS(3,J)
          PMO(J)=PM(J)
          PRMSO(J)=PRMS(J)
   10 CONTINUE
      DO 20 L=1,3
          WSMO(L)=WSM(L)
   20 CONTINUE

      RETURN
      END

C  ****************************** WALLSS **********************
C     CALCULATE WALL SHEAR STRESS and Pressure  at Top and Bottom WALL

      SUBROUTINE WALLSS(U,P,TIME)
      Include 'channel.h'

      Common/DIM/N1,N2,N3,N1M,N2M,N3M
      Common/MESH1/DX1,DX1Q,DX3,DX3Q
      Common/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      Common/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      Common/INDX2/JPA(M2),JMU(M2),JMV(M2)
      COMMON/WSMEM/WSM(3),WSMO(3),ENEJ,ENEJO

      Real U(3,0:M1,0:M2,0:M3),P(M1,M2,M3)
      Real WSB(3,M1,M3),WST(3,M1,M3),PWB(M1,M3),PWT(M1,M3)
      REAL NXZ

      NXZ=1.0/DBLE(N1M*N3M)

C     Bottom and Top WALL at I+1/2,J,K+1/2 Points
      WSM(1)=0.0
      WSM(2)=0.0
      WSM(3)=0.0

      DO 10 I=1,N1M
          IP=IPA(I)
          IM=IMA(I)
          DO 10 K=1,N3M
              KP=KPA(I)
              KM=KMA(I)

              UHC=0.5*(U(1,I,1,K)+U(1,IP,1,K))     ! Bottom Wall
              UHW=0.5*(U(1,I,0,K)+U(1,IP,0,K))
              VHC=0.5*(U(2,I,2,K)+U(2,IP,2,K))
              VHW=0.5*(U(2,I,1,K)+U(2,IP,1,K))
              WHC=0.5*(U(3,I,1,K)+U(3,IP,1,K))
              WHW=0.5*(U(3,I,0,K)+U(3,IP,0,K))
              WSB(1,I,K)=(UHC-UHW)/(DY(1)*0.5)
              WSB(2,I,K)=(VHC-VHW)/DY(1)
              WSB(3,I,K)=(WHC-WHW)/(DY(1)*0.5)
C      PWB(I,K)=P(I,1,K)

              UHC=0.5*(U(1,I,N2M,K)+U(1,IP,N2M,K)) ! Top Wall
              UHW=0.5*(U(1,I,N2 ,K)+U(1,IP,N2 ,K))
              VHC=0.5*(U(2,I,N2M,K)+U(2,IP,N2M,K))
              VHW=0.5*(U(2,I,N2 ,K)+U(2,IP,N2 ,K))
              WHC=0.5*(U(3,I,N2M,K)+U(3,IP,N2M,K))
              WHW=0.5*(U(3,I,N2 ,K)+U(3,IP,N2 ,K))
              WST(1,I,K)=(UHC-UHW)/(DY(N2M)*0.5)
              WST(2,I,K)=(VHC-VHW)/DY(N2M)
              WST(3,I,K)=(WHC-WHW)/(DY(N2M)*0.5)
C      PWT(I,K)=P(I,N2M,K)

              WSM(1)=WSM(1)+(WSB(1,I,K)+WST(1,I,K))
              WSM(2)=WSM(2)+(WSB(2,I,K)+WST(2,I,K))
              WSM(3)=WSM(3)+(WSB(3,I,K)+WST(3,I,K))
   10 CONTINUE
      WSM(1)=WSM(1)*NXZ*0.5
      WSM(2)=WSM(2)*NXZ*0.5
      WSM(3)=WSM(3)*NXZ*0.5

      WRITE(31,100) TIME,WSM(1),WSM(2),WSM(3)
  100 FORMAT(4(E12.5,2X))

      RETURN
      END

C  ****************************** OUTPUT  **********************
C     OUTPUT Flow Quantities As functions of Y,Y^+
C     Wall Quatity Scaling is imposed. 
      SUBROUTINE OUTPUT(TIME,NAVG)

      Include 'channel.h'
      Common/DIM/N1,N2,N3,N1M,N2M,N3M
      Common/PARA/RE
      Common/MESH2/Y(0:M2)
      COMMON/VMEANO/VMO(3,M2),VRMSO(6,M2)
      COMMON/WMEANO/WMO(3,M2),WRMSO(3,M2)
      COMMON/PMEANO/PMO(M2),PRMSO(M2)
      COMMON/WSMEM/WSM(3),WSMO(3),ENEJ,ENEJO
      Real VMP(3),VRMSP(6),WMP(3),WRMSP(3),PMP,PRMSP,WSMP(3)

      OPEN(10,FILE='General.avg',STATUS='UNKNOWN')
      OPEN(11,FILE='Vmean.plt',STATUS='UNKNOWN')
      OPEN(12,FILE='Vrms.plt',STATUS='UNKNOWN')
      OPEN(13,FILE='Wavg.plt',STATUS='UNKNOWN')
      OPEN(14,FILE='Pavg.plt',STATUS='UNKNOWN')

      RNAVG=DBLE(NAVG)
      ENEJP=ENEJO/RNAVG

      DO 10 L=1,3
          WSMP(L)=WSMO(L)/RNAVG
   10 CONTINUE

      UTAU=SQRT(WSMP(1)/RE)

      WRITE(10,100)
      WRITE(10,110) TIME,UTAU,ENEJP,NAVG
  100 FORMAT('TIME',2X,'UTAU',2X,'KINETIC ENERGY',2X,
     >'AVERAGE NUMBER')
  110 FORMAT(3(E12.5,2X),I5)

      DO 20 J=1,N2M
          DO 30 L=1,3
              VMP(L)=VMO(L,J)/RNAVG
   30     VRMSP(L)=SQRT(VRMSO(L,J)/RNAVG)
          DO 40 L=4,6
   40     VRMSP(L)=VRMSO(L,J)/RNAVG
          PMP=PMO(J)/RNAVG
          PRMSP=SQRT(PRMSO(J)/RNAVG)
          DO 50 L=1,3
              WMP(L)=WMO(L,J)/RNAVG
   50     WRMSP(L)=SQRT(WRMSO(L,J)/RNAVG)
          X2P=0.5*(Y(J)+Y(J+1))
          X2PLUS=RE*UTAU*X2P
          WRITE(11,111) X2P,X2PLUS,VMP(1),(VMP(L)/UTAU,L=1,3)
          WRITE(12,112) X2P,X2PLUS,(VRMSP(L)/UTAU,L=1,3)
     >    ,(VRMSP(L)/UTAU**2,L=4,6)
          WRITE(13,113) X2P,X2PLUS,(WMP(L)/UTAU,L=1,3)
     >    ,(WRMSP(L)/UTAU**2/RE,L=1,3)
          WRITE(14,114) X2P,X2PLUS,PMP/UTAU**2,PRMSP/UTAU**2
   20 CONTINUE
  111 FORMAT(6(E12.5,2X))
  112 FORMAT(8(E12.5,2X))
  113 FORMAT(8(E12.5,2X))
  114 FORMAT(4(E12.5,2X))

      CLOSE(10)
      CLOSE(11)
      CLOSE(12)
      CLOSE(13)
      CLOSE(14)

      RETURN
      END


C     DEFINE THE MAIN COEFFICIENT MATRIX 
C
C     | A11  A12  A13 |    1  | I+DT*M11  DT*M12   DT*M13  |
C     | A21  A22  A23 | = --- |  DT*M21  I+DT*M22  DT*M23  |
C     | A31  A32  A33 |   DT  |  DT*M31   DT*M32  I+DT*M33 |
C

C     M11U^N      
c      V2=0.5*(U(2,I,JP,K)+U(2,IM,JP,K))      
c      V1=0.5*(U(2,I,J ,K)+U(2,IM,J ,K))      
c      APJ=JUP*(
c     >      -0.5*HP(J)/RE
c     >      +0.5/DY(J)*V2/H(JP)*DY(J)/2.0
c     >      )
c      ACJ=   0.5*HC(J)/RE
c     >      +0.5/DY(J)*(JUP*V2/H(JP)*DY(JP)/2.0
c     >                 -JUM*V1/H(J )*DY(JM)/2.0)
c      AMJ=JUM*(
c     >      -0.5*HM(J)/RE
c     >      -0.5/DY(J)*V1/H(J)*DY(J)/2.0
c     >      )
c      U2=0.5*(U(1,IP,J,K)+U(1,I ,J,K))
c      U1=0.5*(U(1,I ,J,K)+U(1,IM,J,K)) 
c     API=  -0.5*DX1Q/RE
c     >      +DX1*U2*0.5
c      ACI=   DX1Q/RE
c     >      +DX1*(U2*0.5-U1*0.5)
c      AMI=  -0.5*DX1Q/RE
c     >      -DX1*U1*0.5
c      W2=0.5*(U(3,IM,J,KP)+U(3,I,J,KP))
c      W1=0.5*(U(3,IM,J,K )+U(3,I,J,K ))
c      APK=   -0.5*DX3Q/RE
c     >       +0.5*DX3*W2*0.5
c      ACK=    DX3Q/RE
c     >       +0.5*DX3*(W2*0.5-W1*0.5)
c      AMK=   -0.5*DX3Q/RE
c     >       -0.5*DX3*W1*0.5
c      M11U^N=APJ*U(1,I,JP,K)
c     >      +ACJ*U(1,I,J ,K)
c     >      +AMJ*U(1,I,JM,K)
c     >      +API*U(1,IP,J,K)
c     >      +ACI*U(1,I ,J,K)
c     >      +AMI*U(1,IM,J,K)
c     >      +APK*U(1,I,J,KP)
c     >      +ACK*U(1,I,J,K )
c     >      +AMK*U(1,I,J,KM)
c     M12V^N
c      U2=1.0/H(JP)*(DY(J )/2.*U(1,I,JP,K)+DY(JP)/2.*U(1,I,J ,K))
c      U1=1.0/H(J )*(DY(JM)/2.*U(1,I,J ,K)+DY(J )/2.*U(1,I,JM,K))
c      V2=0.5*(U(2,I,JP,K)+U(2,IM,JP,K))      
c      V1=0.5*(U(2,I,J ,K)+U(2,IM,J ,K))      
c      M12V^N=JUP*0.5/DY(J)*U2*V2
c            -JUM*0.5/DY(J)*U1*V1
C     M13W^N
c      U2=0.5*(U(1,I,J,KP)+U(1,I,J,K ))
c      U1=0.5*(U(1,I,J,K )+U(1,I,J,KM))
c      W2=0.5*(U(3,IM,J,KP)+U(3,I,J,KP)) 
c      W1=0.5*(U(3,IM,J,K )+U(3,I,J,K )) 
c      M13W^N=0.5*DX3*(U2*W2-U1*W1)
c
c
C     M22V^N
c
c      V2=0.5*(U(2,I,JP,K)+U(2,I,J ,K))
c      V1=0.5*(U(2,I,J ,K)+U(2,I,JM,K))
c      APJ=JUP*(
c     >      -0.5*DYP(J)/RE
c     >      +1.0/H(J)*V2*0.5
c     >      )
c      ACJ=
c     >      +0.5*DYC(J)/RE
c     >      +1.0/H(J)*(V2*0.5-V1*0.5)
c      AMJ=JUM*(
c     >      -0.5*DYM(J)/RE
c     >      -1.0/H(J)*V1*0.5
c     >      )
c      U2=1.0/H(J)*(DY(J)/2.0*U(1,IP,JM,K)+DY(JM)/2.0*U(1,IP,J,K))
c      U1=1.0/H(J)*(DY(J)/2.0*U(1,I ,JM,K)+DY(JM)/2.0*U(1,I ,J,K))
c      API=
c     >      -0.5*DX1Q/RE
c     >      +0.5*DX1*U2*0.5
c      ACI=
c     >      +DX1Q/RE
c     >      +0.5*DX1*(U2*0.5-U1*0.5)
c      AMI=
c     >      -0.5*DX1Q/RE
c     >      -0.5*DX1*U1*0.5
c      W2=1.0/H(J)*(DY(J)/2.0*U(3,I,JM,KP)+DY(JM)/2.0*U(3,I,J,KP))
c      W1=1.0/H(J)*(DY(J)/2.0*U(3,I,JM,K )+DY(JM)/2.0*U(3,I,J,K ))
c      APK=
c     >      -0.5*DX3Q/RE
c     >      +0.5*DX3*W2/2.0
c      ACK=
c     >      +DX3Q/RE
c     >      +0.5*DX3*(W2/2.0-W1/2.0)
c      AMK=
c     >      -0.5*DX3Q/RE
c     >      -0.5*DX3*W1/2.0
c      M22V^N=APJ*U(2,I,JP,K)
c     >      +ACJ*U(2,I,J ,K)
c     >      +AMJ*U(2,I,JM,K)
c     >      +API*U(2,IP,J,K)
c     >      +ACI*U(2,I ,J,K)
c     >      +AMI*U(2,IM,J,K)
c     >      +APK*U(2,I,J,KP)
c     >      +ACK*U(2,I,J,K )
c     >      +AMK*U(2,I,J,KM)
C     M21U^N
c      V2=0.5*(U(2,IP,J,K)+U(2,I ,J,K))
c      V1=0.5*(U(2,I ,J,K)+U(2,IM,J,K))
c      U2=1.0/H(J)*(DY(J)/2.*U(1,IP,JM,K)+DY(JM)/2.*U(1,IP,J,K))
c      U1=1.0/H(J)*(DY(J)/2.*U(1,I ,JM,K)+DY(JM)/2.*U(1,I ,J,K))
c      M21U^N=0.5*DX1*(V2*U2-V1*U1)
C     M23W^N
c      V2=0.5*(U(2,I,J,KP)+U(2,I,J,K ))
c      V1=0.5*(U(2,I,J,K )+U(2,I,J,KM))
c      W2=1.0/H(J)*(DY(J)/2.*U(3,I,JM,KP)+DY(JM)/2.*U(3,I,J,KP))
c      W1=1.0/H(J)*(DY(J)/2.*U(3,I,JM,K )+DY(JM)/2.*U(3,I,J,K ))
c      M23W^N=0.5*DX3*(V2*W2-V1*W1)
c
C     M33W^N
c
c      V2=0.5*(U(2,I,JP,K)+U(2,I,JP,KM))
c      V1=0.5*(U(2,I,J ,K)+U(2,I,J ,KM))
c      APJ=JUP*(
c     >      -0.5*HP(J)/RE
c     >      +0.5/DY(J)*V2/H(JP)*DY(J)/2.0 
c     >      )
c      ACJ=
c     >      +0.5*HC(J)/RE
c     >      +0.5/DY(J)*(JUP*V2/H(JP)*DY(JP)/2.0
c     >                 -JUM*V1/H(J)*DY(JM)/2.0) 
c      AMJ=JUM*(
c     >      -0.5*HM(J)/RE
c     >      -0.5/DY(J)*V1/H(J)*DY(J)/2.0 
c     >      )
c      W2=0.5*(U(3,I,J,KP)+U(3,I,J,K ))
c      W1=0.5*(U(3,I,J,K )+U(3,I,J,KM))
c      APK=
c     >      -0.5*DX3Q/RE
c     >      +DX3*W2/2.0
c      ACK=
c     >      +DX3Q/RE
c     >      +DX3*(W2/2.0-W1/2.0)
c      AMK=
c     >      -0.5*DX3Q/RE
c     >      -DX3*W1/2.0
c      U2=0.5*(U(1,IP,J,K)+U(1,IP,J,KM))
c      U1=0.5*(U(1,I ,J,K)+U(1,I ,J,KM))
c      API=
c     >      -0.5*DX1Q/RE
c     >      +0.5*DX1*U2/2.0 
c      ACI=
c     >      +DX1Q/RE
c     >      +0.5*DX1*(U2/2.0-U1/2.0)
c      AMI=
c     >      -0.5*DX1Q/RE
c     >      -0.5*DX1*U1/2.0 
c      M33W^N=APJ*U(3,I,JP,K)
c     >      +ACJ*U(3,I,J ,K)
c     >      +AMJ*U(3,I,JM,K)
c     >      +API*U(3,IP,J,K)
c     >      +ACI*U(3,I ,J,K)
c     >      +AMI*U(3,IM,J,K)
c     >      +APK*U(3,I,J,KP)
c     >      +ACK*U(3,I,J,K )
c     >      +AMK*U(3,I,J,KM)
C     M31U^N
c      W2=0.5*(U(3,IP,J,K)+U(3,I ,J,K))
c      W1=0.5*(U(3,I ,J,K)+U(3,IM,J,K))
c      U2=0.5*(U(1,IP,J,K)+U(1,IP,J,KM))
c      U1=0.5*(U(1,I ,J,K)+U(1,I ,J,KM))
c      M31U^N=0.5*DX1*(W2*U2-W1*U1)
C     M32V^N
c      W2=1.0/H(JP)*(DY(J )/2.*U(3,I,JP,K)+DY(JP)/2.*U(3,I,J ,K))
c      W1=1.0/H(J )*(DY(JM)/2.*U(3,I,J ,K)+DY(J )/2.*U(3,I,JM,K))
c      V2=0.5*(U(2,I,JP,K)+U(2,I,JP,KM))      
c      V1=0.5*(U(2,I,J ,K)+U(2,I,J ,KM))      
c      M32V^N=JUP*0.5/DY(J)*W2*V2
c            -JUM*0.5/DY(J)*W1*V1
