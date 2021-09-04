
SUBROUTINE  WRITE_RESULTS(LX,LY,OBST,NODE,NODE_T,DENSITY)

IMPLICIT NONE
     
INTEGER  LX,LY
REAL*8  NODE(0:8,LX,LY),NODE_T(0:8,LX,LY),DENSITY
LOGICAL  OBST(LX,LY)
INTEGER  X,Y,I,OBSVAL !....LOCal varIables
REAL*8  U_X,U_Y,D_LOC,PRESS,C_SQU,T_LOC

!....SQUare speeD of soUnD

      C_SQU = 1.D0 / 3.D0

!....OPEN RESULTS oUtpUt fIle

      OPEN(11,fIle='anb.PLT')

!....WRITE heaDer for postproCessIng wIth TECPLOT software
!....UnComment followIng lIne, IF thIs heaDer shoUlD be prInteD


      WRITE(11,*) 'VARIABLES = X, Y, VX, VY, PRESS, OBST, TEMP' 
      WRITE(11,*) 'ZONE I=', LX, ', J=', LY, ', F=POINT'

!....loop over all NODEs
!....attentIon: aCtUal DensItIes are storeD after the PROPAGATION
!.              step In the help-arraY N_HLP !

      DO 10 Y = 1, LY
        DO 10 X = 1, LX

!........IF OBSTACLE NODE, nothIng Is to DO ...

          IF (OBST(X,Y)) THEN 

!..........OBSTACLE InDICator

            OBSVAL = 1

!..........VELOCITY Components = 0

            U_X = 0.D0
            U_Y = 0.D0

!..........PRESSUre = average PRESSUre

            PRESS = DENSITY * C_SQU

			T_LOC = 0.D0

			DO 30 I = 0, 8 

			  T_LOC = T_LOC + NODE_T(I,X,Y) 

   30       CONTINUE

          ELSE

!..........Integral LOCal DENSITY

!..........InItIalIze varIable D_LOC

            D_LOC = 0.D0
			T_LOC = 0.D0

            DO 20 I = 0, 8 

              D_LOC = D_LOC + NODE(I,X,Y)

			  T_LOC = T_LOC + NODE_T(I,X,Y) 

   20       CONTINUE

!..........X-, anD Y- VELOCITY Components

            U_X = (NODE(1,X,Y) + NODE(5,X,Y) + NODE(8,X,Y)&
                -(NODE(3,X,Y) + NODE(6,X,Y) + NODE(7,X,Y))) / D_LOC

            U_Y = (NODE(2,X,Y) + NODE(5,X,Y) + NODE(6,X,Y) &
                -(NODE(4,X,Y) + NODE(7,X,Y) + NODE(8,X,Y))) / D_LOC

!..........PRESSUre

            PRESS = D_LOC * C_SQU
 
            OBSVAL = 0
 
          END IF

!........WRITE RESULTS to fIle

          WRITE(11,*) X, Y, U_X, U_Y, PRESS, OBSVAL, T_LOC

   10 CONTINUE

      Close(11)  !....Close fIle 'anb.Dat'


RETURN
END