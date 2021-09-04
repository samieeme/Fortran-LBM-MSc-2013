
SUBROUTINE READ_PARAMETERS(ERROR,T_MAX,DENSITY,ACCEL,OMEGA,R_REY,RO,UX)

IMPLICIT NONE
REAL*8  DENSITY,ACCEL,OMEGA,R_REY,RO,UX
INTEGER  T_MAX
LOGICAL  ERROR


OPEN(10,file='anb.par')   !.....OPEN parameter file



READ(10,*,err=900) T_MAX  !.......line 1: number of iterations

READ(10,*,err=900) DENSITY  !.......line 2: fluid DENSITY per link

READ(10,*,err=900) ACCEL  !.......line 3: DENSITY REDISTRIBUTE

READ(10,*,err=900) OMEGA   !.......line 4: RELAXATION parameter

READ(10,*,err=900) R_REY  !.......line 5: linear dimension (for reynolds number)  $$$$$$$$$$$$$ WHAT'S THIS?

READ(10,*,err=900) UX

CLOSE(10)      !.....close parameter file

WRITE (6,*) '*** Paramters READ from file anb.par.'  !.....information message
WRITE (6,*) '***'

RO=DENSITY

goto 999

900   WRITE (6,*) '!!! ERROR READing file anb.par'   !.....ERROR message: file READ ERROR
      WRITE (6,*) '!!!'

goto 990

990   ERROR = .true.

999   CONTINUE

RETURN
END