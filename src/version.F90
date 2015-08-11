      SUBROUTINE version()

      USE READ_DGINP, ONLY : out_direc
      USE GLOBALS,    ONLY : FORT16
      
      IMPLICIT NONE
      
      PRINT*, "Version Information"
      PRINT*, "  Branch: master" 
      PRINT*, "  SHA: 77060b6daaac05577eea69051ff5178d80d99be4     +" 
      PRINT*, " "
      
      WRITE(FORT16,*) 'DG_WASUPP'
      WRITE(FORT16,*) ' '
      WRITE(FORT16,*) "Version Information"
      WRITE(FORT16,*) "  Branch: master" 
      WRITE(FORT16,*) "  Branch: master" 
      WRITE(FORT16,*) "  SHA: 77060b6daaac05577eea69051ff5178d80d99be4     +" 
 
      END SUBROUTINE version