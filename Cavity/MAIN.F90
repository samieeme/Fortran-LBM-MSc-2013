!                                     IN THE NAME OF GOD
      PROGRAM SAL_SAM
      IMPLICIT NONE    
      INTEGER  LX,LY 
      PARAMETER(LX=95,LY=61)  
!*************************************************************************************************
      !.....FLOW VARIABLES 
      REAL*8  DENSITY !fluid DENSITY per link
      REAL*8  OMEGA  !RELAXATION parameter
      REAL*8  ACCEL  !ACCELleration 
      INTEGER  T_MAX,TK,OBSCOUNT,FLOWCOUNT,OBSCOORD(1000,2),OBSNUM    !maximum number of iterations
      REAL*8   R_REY,TSTART,TEND,RO,UX,F_X,F_Y  !linear dimension for Reynolds number
      INTEGER  TIME,TIMEK,NODID(LX,LY),TIMEMOD  !iteration counter
!     the densities are numbered as follows:
!
!              6   2   5
!                \ | /
!              3 - 0 - 1
!                / | \
!              7   4   8

	  REAL*8  NODE(0:8,LX,LY)  
      REAL*8  N_HLP(0:8,LX,LY) !help array for temporareLY storage of fluid densities
      REAL*8  VEL  !average VELOCITY, computed by SUBROUTINE 'WRITE_VELOCITY'
!*************************************************************************************************
	  !.....TEMP VARIABLES 
      REAL*8  T_INIT,OMEGA_T,T_SUM,NEQ_TEST(LX,LY)
	  REAL*8  NODE_T(0:8,LX,LY),G_X,G_Y,BET,T_REF,THOT,TCOLD,G,RAY,NOU,PR,DT_MAX  
      REAL*8  N_HLP_T(0:8,LX,LY)


      LOGICAL  ERROR  !ERROR flag
      LOGICAL  OBST(LX,LY)  !OBSTACLE array



	  CALL CPU_TIME(TSTART)
!=======================================================================
!     begin initialisation
!=======================================================================
!
!.....initialize ERROR flag
!
      ERROR = .false.
	  NODID=1

      CALL READ_PARAMETERS(ERROR,T_MAX,DENSITY,ACCEL,OMEGA,R_REY,RO,UX) !read parameter file EXCEPT lattice size must be fixed befor compilation
	  CALL READ_TEMP(ERROR,T_INIT,OMEGA_T)		
		
	  THOT = 1.D0

	  TCOLD = 0.D0
	  
      T_REF = 0.D0

	  G_X = 1.D0
	  
	  G_Y = 0.D0

	  G = SQRT( G_X * G_X + G_Y * G_Y ) 

      RAY = 200000.D0

      PR = 0.7D0

	  DT_MAX = THOT -TCOLD

	  NOU = (1.D0 / OMEGA - 0.5D0) /3.D0

      BET = RAY * NOU * NOU / (PR * G * DT_MAX * (REAL (LY) )**3.D0)

      IF (ERROR) goto 990

	  !OBSNUM=10
      CALL MAKING_OBSTACLE(ERROR,OBST,NODID,OBSCOUNT,OBSCOORD,LX,LY)
	! CALL READ_OBSTACLES(ERROR,OBST,NODID,OBSCOUNT,OBSCOORD,LX,LY) !read OBSTACLE file
	  
	 ! CALL SEPERATINGNODE(FLOWCOUNT,NODID,LX,LY)

      IF (ERROR) goto 990 
	!IF an I/O ERROR occurs reading while the OBSTACLE file, the "ERROR"-flag is set "true" and the program stops.

      CALL INIT_DENSITY(LX,LY,DENSITY,NODE)
      CALL INIT_TEMP(LX,LY,T_INIT,NODE_T,N_HLP_T)
	  
      OPEN(10,file='anb_qxm.out')   !mean sqaure flow is monitored and stored in the file'anb_qxm.out', one value for each iteration. The file is OPENed here.   

!=======================================================================
!     END initialisation
!=======================================================================


!***********************************************************************
!     begin iterations
!***********************************************************************
	TK=(T_MAX/10)
	TIME=1


      DO  TIMEK = 1, TK+1    !main loop
	
	
	!The integral fluid DENSITY is checked each T_MAX/10 iteration. The integral fluid DENSITY should be constant all TIME.

      CALL CHECK_DENSITY(LX,LY,NODE,NODE_T,TIME,T_SUM)
	       

	DO TIMEMOD=1,10

	TIME=(TIMEK-1)*10+TIMEMOD


	!directed flow is induced by DENSITY REDISTRIBUTE in the first lattice column.    !!! WHAT'S THIS?  
!	CALL REDISTRIBUTE(LX,LY,OBST,NODE,ACCEL,DENSITY)  
	    
	
		

	!DENSITY PROPAGATION: all fluid densities are PROPAGATED from non-occupied NODEs along the lattice connection lines to their next neighbours.
	CALL PROPAGATE(LX,LY,NODE,N_HLP)
	CALL PROPAGATE_T(LX,LY,NODE_T,N_HLP_T)


!	CALL INFLOW(N_HLP,NODE,LX,LY,UX)
	!CALL INFLOW_T(N_HLP_T,N_HLP,LX,LY,UX,NODE)
	
!	CALL OUTFLOW(N_HLP,LX,LY,RO)
	!CALL OUTFLOW_T(LX,LY,NODE,N_HLP_T)

	!bounc back from OBSTACLEs: this is the no-slip boundary-condition.
	CALL BOUNCEBACK(OBSCOUNT,OBSCOORD,LX,LY,OBST,NODE,N_HLP,NODE_T,N_HLP_T,THOT,TCOLD)
	


	!DENSITY RELAXATION: a single TIME RELAXATION with RELAXATION parameter OMEGA is applied here.
    CALL RELAX_NCONVEC(DENSITY,OMEGA,LX,LY,NODE,N_HLP,NODID,N_HLP_T,G_X,G_Y,BET,T_REF)
!	CALL RELAXATION(DENSITY,OMEGA,LX,LY,NODE,N_HLP,NODID,F_X,F_Y)
	CALL RELAX_T(OMEGA_T,LX,LY,NODE_T,N_HLP,N_HLP_T,NODID,F_X,F_Y,OMEGA) 

	!average flow VELOCITY is computed at a cross section in the  middle of the channel and written to the file "anb_qx.out" ever iteration.
    CALL WRITE_VELOCITY(LX,LY,TIME,OBST,NODE,VEL)        


	END DO
	END DO  !100   END of the main loop


	!compute fluid-VELocities u,v and pressure from VELOCITY distribution, and WRITE to file anb_rs.out.
      CALL WRITE_RESULTS (LX,LY,OBST,NODE,NODE_T,DENSITY)
	                    
	
	!compute reynolds number
      CALL COMP_REY(LX,LY,OBST,NODE,TIME,OMEGA,DENSITY,R_REY)

      goto 999
	
	
	!Here we get onLY, IF the "ERROR" flag was set "true".
990   WRITE (6,*) '!!! ERROR: program stopped during iteration =', TIME   !!!!!!!!WHAT'S 6?
      WRITE (6,*) '!!!'
      
  999 continue

      close(10) !(anb_qx.out) is closed. it tells you, IF you reached steady state and IF you got reasonable RESULTS.

      CALL CPU_TIME(TEND)
	  WRITE (6,*) '********************    END     ********************'  , OBSCOUNT,TSTART-TEND
	
	  
      stop
      END PROGRAM SAL_SAM
 
	INCLUDE 'READ_PARAMETERS.F90'
	INCLUDE	'MAKING_OBSTACLE.F90'
	INCLUDE 'READ_OBSTACLES.F90'
!	INCLUDE 'SEPERATINGNODE.F90'
	INCLUDE 'INIT_DENSITY.F90'
	INCLUDE 'CHECK_DENSITY.F90'
	INCLUDE 'REDISTRIBUTE.F90'
	INCLUDE 'BOUNCEBACK.F90'
	INCLUDE 'PROPAGATE.F90'
	INCLUDE 'INFLOW.F90' !(N_HELP,LX,LY,UX)
	INCLUDE 'OUTFLOW.F90' !(N_HELP,LX,LY,RO)
	INCLUDE 'RELAXATION.F90'
	INCLUDE 'WRITE_VELOCITY.F90'
	INCLUDE 'WRITE_RESULTS.F90'
	INCLUDE 'COMP_REY.F90'
	INCLUDE 'READ_TEMP.F90'
	INCLUDE 'INIT_TEMP.F90'
	INCLUDE 'PROPAGATE_T.F90'
	INCLUDE 'RELAX_T.F90'
	INCLUDE 'INFLOW_T.F90'
	INCLUDE 'OUTFLOW_T.F90'
	INCLUDE 'RELAX_NCONVEC.F90' 