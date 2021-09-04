
SUBROUTINE BOUNCEBACK(OBSCOUNT,OBSCOORD,LX,LY,OBST,NODE,N_HLP,NODE_T,N_HLP_T,THOT,TCOLD)

IMPLICIT NONE

INTEGER  LX,LY,OBSCOUNT,OBSCOORD(1000,2)
LOGICAL  OBST(LX,LY)
REAL*8   T_P,NODE(0:8,LX,LY),N_HLP(0:8,LX,LY),NODE_T(0:8,LX,LY)
REAL*8   THOT,TCOLD,T_0,T_1,T_2,T_WALL,T_I1,T_I2,T_EQ,N_HLP_T(0:8,LX,LY)
INTEGER  I,X,Y !.....local variables

open(22,file='obstacleout.txt')
write(22,*) OBSCOUNT

DO  I=1,OBSCOUNT
	
			X=OBSCOORD(I,1)
			Y=OBSCOORD(I,2)
			!write(22,*) OBSCOORD(I,1),OBSCOORD(I,2)
!...........rotate all ensities and WRITE back to NODE

!...........east

            NODE(1,X,Y) = N_HLP(3,X,Y)

!...........north

            NODE(2,X,Y) = N_HLP(4,X,Y)

!...........west

            NODE(3,X,Y) = N_HLP(1,X,Y)

!...........south

            NODE(4,X,Y) = N_HLP(2,X,Y)

!...........north-east

            NODE(5,X,Y) = N_HLP(7,X,Y)

!...........north-west

            NODE(6,X,Y) = N_HLP(8,X,Y)

!...........south-west

            NODE(7,X,Y) = N_HLP(5,X,Y)

!...........south-east

            NODE(8,X,Y) = N_HLP(6,X,Y)



END DO

T_0 = 4.D0 /9.D0
T_1 = 1.D0 /  9.D0
T_2 = 1.D0 / 36.D0


Y=1
DO  X=1,LX
	
	T_I1 = 0.D0
	T_I2 = 0.D0

	DO I = 0,8
	
	T_I1 = T_I1 + N_HLP_T(I,X,2) 
	T_I2 = T_I2 + N_HLP_T(I,X,3)

	END DO

	T_P = N_HLP_T(0,X,Y) + N_HLP_T(1,X,Y) + N_HLP_T(3,X,Y) + N_HLP_T(4,X,Y) + N_HLP_T(7,X,Y) + N_HLP_T(8,X,Y)
	T_WALL = THOT !1.D0. / 3.D0. * ( 4.D0 * T_I1 - T_I2) ADIABATIC
	T_EQ = 6.D0 * ( T_WALL - T_P )
	
	!N_HLP_T(0,X,1) =		!T_0 * T_WALL  ADIABATIC
	!N_HLP_T(1,X,1) =		!T_1 * T_WALL
	N_HLP_T(2,X,1) = T_1 * T_EQ		!T_1 * T_WALL
    !N_HLP_T(3,X,1) =		!T_1 * T_WALL
	!N_HLP_T(4,X,1) =		!T_1 * T_WALL
	N_HLP_T(5,X,1) = T_2 * T_EQ		!T_2 * T_WALL
	N_HLP_T(6,X,1) = T_2 * T_EQ		!T_2 * T_WALL
	!N_HLP_T(7,X,1) =		!T_2 * T_WALL
    !N_HLP_T(8,X,1) =		!T_2 * T_WALL

END DO

Y=LY
!T_WALL=300
DO  X=1,LX
	

	T_I1 = 0.D0
	T_I2 = 0.D0

	DO I = 0,8
	
	T_I1 = T_I1 + N_HLP_T(I,X,LY-1) 
	T_I2 = T_I2 + N_HLP_T(I,X,LY-2)

	END DO

	T_P = N_HLP_T(0,X,Y) + N_HLP_T(1,X,Y) + N_HLP_T(2,X,Y) + N_HLP_T(5,X,Y) + N_HLP_T(6,X,Y) + N_HLP_T(3,X,Y)
	T_WALL = TCOLD ! T_WALL = 1.D0 / 3.D0 * ( 4.D0 * T_I1 - T_I2)  ! ADIABATIC
	T_EQ = 6.D0 * ( T_WALL - T_P )

	!N_HLP_T(0,X,LY) = T_0 * T_WALL
	!N_HLP_T(1,X,LY) = T_1 * T_WALL
	!N_HLP_T(2,X,LY) = T_1 * T_WALL
    !N_HLP_T(3,X,LY) = T_1 * T_WALL
	N_HLP_T(4,X,LY) = T_1 * T_EQ
	!N_HLP_T(5,X,LY) = T_2 * T_WALL
	!N_HLP_T(6,X,LY) = T_2 * T_WALL
	N_HLP_T(7,X,LY) = T_2 * T_EQ
    N_HLP_T(8,X,LY) = T_2 * T_EQ

END DO



X=1
DO  Y=2,LY-1
	
	T_I1 = 0.D0
	T_I2 = 0.D0

	DO I = 0,8
	
	T_I1 = T_I1 + N_HLP_T(I,2,Y) 
	T_I2 = T_I2 + N_HLP_T(I,3,Y)

	END DO

	T_P = N_HLP_T(0,X,Y) + N_HLP_T(2,X,Y) + N_HLP_T(3,X,Y) + N_HLP_T(4,X,Y) + N_HLP_T(6,X,Y) + N_HLP_T(7,X,Y)
	T_WALL = 1.D0 / 3.D0 * ( 4.D0 * T_I1 - T_I2) !ADIABATIC
	T_EQ = 6.D0 * ( T_WALL - T_P )
	
	!N_HLP_T(0,X,1) =		!T_0 * T_WALL  ADIABATIC
	!N_HLP_T(1,X,1) =		!T_1 * T_WALL
	N_HLP_T(1,X,Y) = T_1 * T_EQ		!T_1 * T_WALL
    !N_HLP_T(3,X,1) =		!T_1 * T_WALL
	!N_HLP_T(4,X,1) =		!T_1 * T_WALL
	N_HLP_T(5,X,Y) = T_2 * T_EQ		!T_2 * T_WALL
	N_HLP_T(8,X,Y) = T_2 * T_EQ		!T_2 * T_WALL
	!N_HLP_T(7,X,1) =		!T_2 * T_WALL
    !N_HLP_T(8,X,1) =		!T_2 * T_WALL

END DO

X=LX
!T_WALL=300
DO  Y=2,LY-1
	

	T_I1 = 0.D0
	T_I2 = 0.D0

	DO I = 0,8
	
	T_I1 = T_I1 + N_HLP_T(I,LX-1,Y) 
	T_I2 = T_I2 + N_HLP_T(I,LX-2,Y)

	END DO

	T_P = N_HLP_T(0,X,Y) + N_HLP_T(1,X,Y) + N_HLP_T(2,X,Y) + N_HLP_T(4,X,Y) + N_HLP_T(5,X,Y) + N_HLP_T(8,X,Y)
	T_WALL = 1.D0 / 3.D0 * ( 4.D0 * T_I1 - T_I2)  ! ADIABATIC
	T_EQ = 6.D0 * ( T_WALL - T_P )

	!N_HLP_T(0,X,LY) = T_0 * T_WALL
	!N_HLP_T(1,X,LY) = T_1 * T_WALL
	!N_HLP_T(2,X,LY) = T_1 * T_WALL
    !N_HLP_T(3,X,LY) = T_1 * T_WALL
	N_HLP_T(3,X,Y) = T_1 * T_EQ
	!N_HLP_T(5,X,LY) = T_2 * T_WALL
	!N_HLP_T(6,X,LY) = T_2 * T_WALL
	N_HLP_T(6,X,Y) = T_2 * T_EQ
    N_HLP_T(7,X,Y) = T_2 * T_EQ

END DO


!DO  X = 1, LX                !.....loop over all NODEs
!	DO  Y = 1, LY

!.........consider onLY OBSTACLE NODEs

!		IF (OBST(X,Y)) THEN


!		END IF

!	END DO
!END DO


RETURN
END