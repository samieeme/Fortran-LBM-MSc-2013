
SUBROUTINE WRITE_VELOCITY(LX,LY,TIME,OBST,NODE,VEL)        

IMPLICIT NONE

INTEGER  LX,LY,TIME
LOGICAL  OBST(LX,LY)
REAL*8   NODE(0:8,LX,LY),VEL
INTEGER  X,Y,I,N_FREE !....LOCal varIables
REAL*8  U_X,D_LOC


!....loop over chaNNel cross sectIoN at half chaNNel leNgth LX/2
      X = INT(FLOAT(LX) / 2.D0)

!....INItIalIze coUNter
      N_FREE= 0
      U_X = 0.D0

      DO 10 Y = 1, LY

!......oNLY NoN-occUpIeD NODEs are coNsIDereD here

        IF(.NOT. OBST(X,Y)) THEN

!........INtegral LOCal DENSITY

!........INItIalIze varIable D_LOC

          D_LOC = 0.D0

          DO 20 I = 0, 8 

            D_LOC = D_LOC + NODE(I,X,Y)

   20     CONTINUE

!........X-, aND Y- VELOCITY compoNeNts

          U_X = U_X + (NODE(1,X,Y) + NODE(5,X,Y) + NODE(8,X,Y) &
                    -(NODE(3,X,Y) + NODE(6,X,Y) + NODE(7,X,Y))) / D_LOC

!........INcrease coUNter

          N_FREE= N_FREE+ 1

        END IF

   10 CONTINUE

!....average VELOCITY

      VEL = U_X / float(N_free)

!....WRITE to fIle

      WRITE(10,*) TIME, VEL


RETURN
END