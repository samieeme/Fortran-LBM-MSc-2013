
SUBROUTINE REDISTRIBUTE(LX,LY,OBST,NODE,ACCEL,DENSITY)

IMPLICIT NONE

INTEGER  LX,LY
LOGICAL  OBST(LX,LY)
REAL*8   NODE(0:8,LX,LY),ACCEL,DENSITY
INTEGER  Y   !.....local variables
REAL*8  T_1,T_2


!     increasing/decreasing inleT densiTies

      T_1 = DENSITY * ACCEL /  9.d0
      T_2 = DENSITY * ACCEL / 36.d0

DO Y = 1, LY

!.......ACCELeraTe flow onLY on non-occupied NODEs

        IF (.noT. OBST(1,Y) .and.&

!.........check To avoid negaTive densiTies

         NODE(3,1,Y) - T_1 .gT. 0. .and.&
         NODE(6,1,Y) - T_2 .gT. 0. .and.&
         NODE(7,1,Y) - T_2 .gT. 0.) THEN 

!.........increase east

          NODE(1,1,Y) = NODE(1,1,Y) + T_1

!.........decrease west

          NODE(3,1,Y) = NODE(3,1,Y) - T_1

!.........increase north-east

          NODE(5,1,Y) = NODE(5,1,Y) + T_2

!.........decrease north-west

          NODE(6,1,Y) = NODE(6,1,Y) - T_2

!.........decrease south-west

          NODE(7,1,Y) = NODE(7,1,Y) - T_2

!.........increase south-east

          NODE(8,1,Y) = NODE(8,1,Y) + T_2

        END IF

END DO


RETURN
END