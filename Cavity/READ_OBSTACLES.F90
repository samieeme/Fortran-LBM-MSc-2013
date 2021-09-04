
SUBROUTINE READ_OBSTACLES(ERROR,OBST,NODID,OBSCOUNT,OBSCOORD,LX,LY)

IMPLICIT NONE
     
INTEGER  LX,LY,OBSCOUNT,IO,IOSTAT,OBSCOORD(1000,2),NODID(LX,LY)  
LOGICAL  ERROR,OBST(LX,LY)
INTEGER  X,Y !.....local variables


DO Y = 1, LY                           !.....no OBSTACLEs in OBSTACLE arraY  
    DO X = 1, LX

      OBST(X,Y) = .FALSE.

	END DO
END DO


!.....OPEN OBSTACLE file
!.....the OBSTACLE file is defined bY the X- and Y- coordinates of the
!     OBSTACLEs. Each OBSTACLE has a line of its own. Also boundaries
!     have to be defined here.
!

OPEN(10,file='SALSAM.OBS')

OBSCOUNT=0
!IO=0
!DO WHILE ((IO.NE.(-1)) .AND. (IO.NE.(1)))

!	READ(10,*,IOSTAT=IO) X,Y
!	OBSCOUNT=OBSCOUNT+1;

!END DO

!ALLOCATE(OBSCOORD(OBSCOUNT,2))


20 CONTINUE

READ(10,*,END=50,err=900) X,Y    !.....READ OBSTACLE coordinates

IF (X .le. LX .and. Y .le. LY) THEN  !.......check IF OBSTACLE inside DOmain boundaries

	OBST(X,Y) = .TRUE.  !.......define OBSTACLE
	OBSCOUNT=OBSCOUNT+1
	OBSCOORD(OBSCOUNT,1)=X   
	OBSCOORD(OBSCOUNT,2)=Y
	NODID(X,Y)=2
ELSE

	WRITE(6,*) '!!! OBSTACLE out of range, skipped'
	WRITE(6,*) '!!! LX = ', X, ' , LY = ', Y
	WRITE(6,*) '!!!'

END IF

GOTO 20

50   CONTINUE

CLOSE(10)  !.....close OBSTACLE file

     WRITE (6,*) '*** GeometrY information READ from file anb.obs.'
     WRITE (6,*) '***' , OBSCOUNT

GOTO 999

900 WRITE (6,*) '!!! ERROR READing file anb.obs'   !.....ERROR message: file READ ERROR
    WRITE (6,*) ' '

GOTO 990

990 ERROR = .TRUE.

999 CONTINUE

RETURN
END