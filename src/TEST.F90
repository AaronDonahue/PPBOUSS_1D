      PROGRAM TEST
      
!       USE GLOBALS, ONLY : GRIDFILE
      
      IMPLICIT NONE
      LOGICAL            :: NAFound,TRUTH
      CHARACTER(LEN=100) :: CSTRING,MYCASE,JC1,JC2
      
      CHARACTER(LEN=100) :: GRIDFILE
      
      TRUTH = .true.
      
!.....Check for control file
      INQUIRE(FILE='fort.15',EXIST=NAFound)
      ! If it doesn't exist stop computation
      IF (.NOT.NAFound) THEN
        WRITE(*,*) 'Control File (fort.15) not found'
        WRITE(*,*) 'Terminating Program'
        CALL EXIT
      ELSE
        OPEN(UNIT=15,FILE='fort.15',ACTION='READ')
      END IF
      
!.....Process control file
      DO WHILE (TRUTH)
        READ(15,'(A100)',END=100) CSTRING
        READ(CSTRING,'(A)') MYCASE
        print*, CSTRING, MYCASE
        SELECT CASE (TRIM(ADJUSTL(MYCASE)))
!.......Gridfile
        CASE('grid')
          READ(CSTRING,'(A,A,A)') JC1, JC2, GRIDFILE
          print*, 'hello'
          print*, GRIDFILE
        END SELECT        
      END DO
100   CONTINUE      
      
      END PROGRAM TEST