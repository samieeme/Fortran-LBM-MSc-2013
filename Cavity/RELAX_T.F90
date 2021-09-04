
SUBROUTINE RELAX_T(OMEGA_T,LX,LY,NODE_T,N_HLP,N_HLP_T,NODID,F_X,F_Y,OMEGA_2)

IMPLICIT NONE

INTEGER  LX,LY,NODID(LX,LY)
LOGICAL  OBST(LX,LY)
REAL*8   DENSITY,OMEGA,NODE(0:8,LX,LY),N_HLP(0:8,LX,LY)
REAL*8   OMEGA_T,NODE_T(0:8,LX,LY),N_HLP_T(0:8,LX,LY)
INTEGER  X,Y,I  !....LOCAL VARIABLES
REAL*8  C_SQU,T_0,T_1,T_2,U_X,U_Y,U_N(8),N_EQU(0:8),U_SQU,D_LOC,T_LOC
REAL*8  F_X,F_Y,OMEGA_2


!....WEIGHTING FACTORS (DEPENDING ON LATTICE GEOMETRY)

T_0 = 4.D0 /9.D0
T_1 = 1.D0 /  9.D0
T_2 = 1.D0 / 36.D0


!....SQUARE SPEED OF SOUND

      C_SQU = 1.D0

!....LOOP OVER ALL NODES
!....ATTENTION: ACTUAL DENSITIES ARE STORED AFTER THE PROPAGATION
!                STEP IN THE HELP-ARRAY N_HLP !
F_X = 0.D0

F_Y = 0.D0


DO X=1,LX

	!DO WHILE (Y .LE. LY)
	DO Y=1,LY	
		

	!	DO WHILE (NODID(X,Y)==1)

			D_LOC = 0.D0
			T_LOC = 0.D0

            DO  I = 0, 8 

              D_LOC = D_LOC + N_HLP(I,X,Y)

			  T_LOC = T_LOC + N_HLP_T(I,X,Y)
			  
			END DO

			
!..........X-, AND Y- VELOCITY COMPONENTS

            U_X = (N_HLP(1,X,Y) + N_HLP(5,X,Y) + N_HLP(8,X,Y) - (N_HLP(3,X,Y) + N_HLP(6,X,Y) + N_HLP(7,X,Y))) / D_LOC

            U_Y = (N_HLP(2,X,Y) + N_HLP(5,X,Y) + N_HLP(6,X,Y) - (N_HLP(4,X,Y) + N_HLP(7,X,Y) + N_HLP(8,X,Y))) / D_LOC

			U_X = U_X + F_X / ( D_LOC * OMEGA_2 )  
			
			U_Y = U_Y + F_Y / ( D_LOC * OMEGA_2 ) 


!..........SQUARE VELOCITY

!IF (NODID(X,Y)==1) THEN

IF (NODID(X,Y)==2) THEN
            U_X = 0.D0
			U_Y=0.D0
END IF


            U_SQU = U_X * U_X + U_Y * U_Y

!..........N- VELOCITY COMPNENTS (N = LATTICE NODE CONNECTION VECTORS)
!..........THIS IS ONLY NECESSARY FOR CLEARENCE, AND ONLY 3 SPEEDS WOULD
!..........BE NECESSARY

            U_N(1) =   U_X
            U_N(2) =         U_Y
            U_N(3) = - U_X
            U_N(4) =       - U_Y
            U_N(5) =   U_X + U_Y
            U_N(6) = - U_X + U_Y
            U_N(7) = - U_X - U_Y
            U_N(8) =   U_X - U_Y



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!..........EQUILIBRIUM TEMPEERATURE
!..........THIS CAN BE REWRITTEN TO IMPROVE COMPUTATIONAL PERFORMANCE
!..........CONSIDERABELY !

!..........ZERO VELOCITY DENSITY
!
            N_EQU(0) = T_0 * T_LOC !* (U_SQU / C_SQU)

!..........AXIS SPEEDS (FACTOR: T_1)

            N_EQU(1) = T_1 * T_LOC * ( 1.D0 + 3.D0 / C_SQU * U_N(1) ) !* (3.D0 / 2.D0 + 3.D0 * U_N(1) / ( 2.D0 * C_SQU ) + 9.D0 / 2.D0 * (U_N(1) / C_SQU) * (U_N(1) / C_SQU) - 3.D0 / 2.D0 * (U_SQU / C_SQU) * (U_SQU / C_SQU)) 

            N_EQU(2) = T_1 * T_LOC * ( 1.D0 + 3.D0 / C_SQU * U_N(2) ) !* (3.D0 / 2.D0 + 3.D0 * U_N(2) / ( 2.D0 * C_SQU ) + 9.D0 / 2.D0 * (U_N(2) / C_SQU) * (U_N(2) / C_SQU) - 3.D0 / 2.D0 * (U_SQU / C_SQU) * (U_SQU / C_SQU))

            N_EQU(3) = T_1 * T_LOC * ( 1.D0 + 3.D0 / C_SQU * U_N(3) ) !* (3.D0 / 2.D0 + 3.D0 * U_N(3) / ( 2.D0 * C_SQU ) + 9.D0 / 2.D0 * (U_N(3) / C_SQU) * (U_N(3) / C_SQU) - 3.D0 / 2.D0 * (U_SQU / C_SQU) * (U_SQU / C_SQU))

            N_EQU(4) = T_1 * T_LOC * ( 1.D0 + 3.D0 / C_SQU * U_N(4) ) !* (3.D0 / 2.D0 + 3.D0 * U_N(4) / ( 2.D0 * C_SQU ) + 9.D0 / 2.D0 * (U_N(4) / C_SQU) * (U_N(4) / C_SQU) - 3.D0 / 2.D0 * (U_SQU / C_SQU) * (U_SQU / C_SQU))

		
!..........DIAGONAL SPEEDS (FACTOR: T_2)

            N_EQU(5) = T_2 * T_LOC * ( 1.D0 + 3.D0 / C_SQU * U_N(5) ) !* (3.D0 + 6.D0 * U_N(5) / C_SQU + 9.D0 / 2.D0 * (U_N(5) / C_SQU) * (U_N(5) / C_SQU) - 3.D0 / 2.D0 * (U_SQU / C_SQU) * (U_SQU / C_SQU))

            N_EQU(6) = T_2 * T_LOC * ( 1.D0 + 3.D0 / C_SQU * U_N(6) ) !* (3.D0 + 6.D0 * U_N(6) / C_SQU + 9.D0 / 2.D0 * (U_N(6) / C_SQU) * (U_N(6) / C_SQU) - 3.D0 / 2.D0 * (U_SQU / C_SQU) * (U_SQU / C_SQU))

            N_EQU(7) = T_2 * T_LOC * ( 1.D0 + 3.D0 / C_SQU * U_N(7) ) !* (3.D0 + 6.D0 * U_N(7) / C_SQU + 9.D0 / 2.D0 * (U_N(7) / C_SQU) * (U_N(7) / C_SQU) - 3.D0 / 2.D0 * (U_SQU / C_SQU) * (U_SQU / C_SQU))

            N_EQU(8) = T_2 * T_LOC * ( 1.D0 + 3.D0 / C_SQU * U_N(8) ) !* (3.D0 + 6.D0 * U_N(8) / C_SQU + 9.D0 / 2.D0 * (U_N(8) / C_SQU) * (U_N(8) / C_SQU) - 3.D0 / 2.D0 * (U_SQU / C_SQU) * (U_SQU / C_SQU))
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$




!..........RELAXATION STEP

		DO  I = 0, 8

			NODE_T(I,X,Y) = N_HLP_T(I,X,Y) + OMEGA_T * (N_EQU(I) - N_HLP_T(I,X,Y))

		END DO
		!Y=Y+1
		!Y=Y+1
	END DO
END DO

RETURN
END