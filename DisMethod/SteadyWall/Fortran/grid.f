C--------------------------------------------------------------
C     GRID GENERATION IN CHANNEL FLOW
C     X2 Clustering only
C     JICHOI
C--------------------------------------------------------------

      PROGRAM MAIN
      PARAMETER (M2=129)
      DIMENSION Y(0:M2)
      PI=ACOS(-1.0)
      N2=129
      N2M=N2-1

      ALY=1.0     ! CHANNEL HALF WIDTH
      EPS=1.0E-30 ! SMALL VALUES
      Y(0)=0.0

      GAMMA0=.1  ! INITIAL GUESS

      Y(2)=1.0E-3 ! FIRST NON-ZERO Y GRID POINTS

      YTILDA=REAL(1)/REAL(N2M)*2
      OLD=ALY*(1.0-TANH(GAMMA0*(ALY-YTILDA))/TANH(GAMMA0*ALY))

      GAMMA=1. ! SECOND GUESS

      DO 5 I=1,100
          TEMP=ALY*(1.0-TANH(GAMMA*(ALY-YTILDA))/TANH(GAMMA*ALY))
          DYDG=(TEMP-OLD)/(GAMMA-GAMMA0+EPS)
          GAMMA0=GAMMA
          OLD=TEMP
          IF(DYDG.EQ.0) GOTO 200
          GAMMA=(Y(2)-TEMP)/(DYDG)+GAMMA0
          IF(ABS(TEMP-Y(2)).LE.1E-30) GOTO 200
    5 CONTINUE

  200 WRITE(*,*) GAMMA

      DO 10 J=1,N2
          YTILDA  =REAL(J-1)/REAL(N2-1)*2
          Y(J)=ALY*(1.0-TANH(GAMMA*(ALY-YTILDA))/TANH(GAMMA*ALY))
   10 CONTINUE

      OPEN(2,FILE='channel.grd',STATUS='UNKNOWN')
      WRITE(2,*) (Y(J),J=0,N2)
      CLOSE(2)

      DO J=0,N2
          WRITE(*,*)J,Y(J)
      ENDDO

      STOP
      END
