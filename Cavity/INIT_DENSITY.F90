
SUBROUTINE INIT_DENSITY(LX,LY,DENSITY,NODE)

IMPLICIT NONE

INTEGER  LX,LY
REAL*8  DENSITY,NODE(0:8,LX,LY)
INTEGER  X,Y  !.....local variables
REAL*8  T_0,T_1,T_2

!.....compuTe weighTing facTors (depENDing on laTTice geomeTrY)
T_0 = DENSITY *  4.d0 / 9.d0
T_1 = DENSITY /  9.d0
T_2 = DENSITY / 36.d0

!.....loop over compuTaTional DOmain

DO X = 1, LX
	DO Y = 1, LY

		NODE(0,X,Y) = T_0 !.........zero VELOCITY DENSITY

        NODE(1,X,Y) = T_1  !.........equilibrium densiTies for aXis speeds
        NODE(2,X,Y) = T_1
        NODE(3,X,Y) = T_1
        NODE(4,X,Y) = T_1

        NODE(5,X,Y) = T_2  !.........equilibrium densiTies for diagonal speeds
        NODE(6,X,Y) = T_2
        NODE(7,X,Y) = T_2
        NODE(8,X,Y) = T_2


	
	END DO
END DO


RETURN
END