
SUBROUTINE PROPAGATE(LX,LY,NODE,N_HLP)

IMPLICIT NONE

INTEGER  LX,LY
REAL*8   NODE(0:8,LX,LY),N_HLP(0:8,LX,LY)
INTEGER  X,Y,X_E,X_W,Y_N,Y_S  !...local variablES


DO  X = 1, LX     !...loop ovEr all NODES
	DO  Y = 1, LY

!.......computE uppEr aNd right NEXt NEighbour NODES With rEgard
!         to pEriodic bouNdariES

          Y_N = mod(Y,LY) + 1
          X_E = mod(X,LX) + 1

!.......computE loWEr aNd lEft NEXt NEighbour NODES With rEgard to
!         pEriodic bouNdariES

          Y_S = LY - mod(LY + 1 - Y, LY)
          X_W = LX - mod(LX + 1 - X, LX)

!.......DENSITY PROPAGATION

!.......zEro: juSt copY

          N_HLP(0,X  ,Y  ) = NODE(0,X,Y)

!.......EaSt

          N_HLP(1,X_E,Y  ) = NODE(1,X,Y)

!.......North

          N_HLP(2,X  ,Y_N) = NODE(2,X,Y)

!.......WESt

          N_HLP(3,X_W,Y  ) = NODE(3,X,Y)

!.......South

          N_HLP(4,X  ,Y_S) = NODE(4,X,Y)

!.......North-EaSt

          N_HLP(5,X_E,Y_N) = NODE(5,X,Y)

!.......North-WESt

          N_HLP(6,X_W,Y_N) = NODE(6,X,Y)

!.......South-WESt

          N_HLP(7,X_W,Y_S) = NODE(7,X,Y)

!.......South-EaSt

          N_HLP(8,X_E,Y_S) = NODE(8,X,Y)

	END DO
END DO

RETURN
END