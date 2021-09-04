
SUBROUTINE READ_TEMP(ERROR,T_INIT,OMEGA_T)

IMPLICIT NONE
REAL*8  T_INIT,OMEGA_T
LOGICAL  ERROR


OPEN(10,file='ANB_T.PAR')  

READ(10,*,err=900) T_INIT  

READ(10,*,err=900) OMEGA_T   

CLOSE(10)     

WRITE (6,*) '*** Paramters READ from file ANB_T.PAR.'  
WRITE (6,*) '***'

goto 999

900   WRITE (6,*) '!!! ERROR READing file ANB_T.PAR'   
      WRITE (6,*) '!!!'

goto 990

990   ERROR = .true.

999   CONTINUE

RETURN
END