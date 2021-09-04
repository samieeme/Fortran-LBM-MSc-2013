      SUBROUTINE COMP_REY(LX,LY,OBST,NODE,TIME,OMEGA,DENSITY,R_REY)

      implicit none

      integer  LX,LY,TIME

      REAL*8  DENSITY,NODE(0:8,LX,LY),OMEGA,R_REY

      logical  OBST(LX,LY)

!....local variables

      REAL*8  VEL,visc,rey

!....compute average VELOCITY

      CALL WRITE_VELOCITY(LX,LY,TIME,OBST,NODE,VEL)        

!....compute viscosity

      visc = 1.d0 / 6.d0 * (2.d0 / OMEGA - 1.d0)

!....compute Reynolds number

      rey = VEL * R_REY / visc


!....messages

       WRITE (6,*) '*** Calculations finished, RESULTS:'
       WRITE (6,*) '***'
       WRITE (6,*) '*** viscosity = ', visc
       WRITE (6,*) '*** average VELOCITY = ', VEL
       WRITE (6,*) '*** Reynolds number = ', rey
       WRITE (6,*) '***'
       WRITE (6,*) '*** In the file anb.dat, you can find local'
       WRITE (6,*) '*** information about the simulated flow.'
       WRITE (6,*) '***'
       WRITE (6,*) '*** In the file anb_qx.out, you can find the average &'
       WRITE (6,*) '*** flow VELOCITY plotted as a function of TIME.'
       WRITE (6,*) '***'


      return
      END