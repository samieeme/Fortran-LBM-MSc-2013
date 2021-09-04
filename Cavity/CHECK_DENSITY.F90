
SUBROUTINE CHECK_DENSITY(LX,LY,NODE,NODE_T,TIME,T_SUM)  

IMPLICIT NONE
     
INTEGER  LX,LY,TIME
REAL*8  NODE(0:8,LX,LY),NODE_T(0:8,LX,LY)
INTEGER  X,Y,N   !.....local variables
REAL*8 N_SUM,T_SUM

N_SUM = 0.d0
T_SUM = 0.d0

DO Y = 1, LY   !.....loop over computatioNal DOmaiN
	DO X = 1, LX
		DO N = 0, 8  !...........loop over all deNsities
			
			N_SUM = N_SUM + NODE(N,X,Y)  !.............sum up densities
			T_SUM = T_SUM + NODE_T(N,X,Y)

		END DO
	END DO
END DO

WRITE(6,*) '*** IteratioN Number = ', TIME
WRITE(6,*) '*** INtegral DENSITY = ', N_sum
WRITE(6,*) '*** INtegral TEMP = ', T_SUM
WRITE(6,*) '***'


RETURN
END