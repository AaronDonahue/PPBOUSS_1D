      SUBROUTINE version()

      USE READ_DGINP, ONLY : out_direc
      USE GLOBALS,    ONLY : FORT16
      
      IMPLICIT NONE
      
      PRINT*, "Version Information"
      PRINT*, "  Branch: master" 
      PRINT*, "  SHA: be0fba366d4485a030113da97ffc39b1a13b8879     +" 
      PRINT*, " "
      
      WRITE(FORT16,*) 'DG_WASUPP'
      WRITE(FORT16,*) ' '
      WRITE(FORT16,*) "Version Information"
      WRITE(FORT16,*) "  Branch: master" 
      WRITE(FORT16,*) "  Branch: master" 
      WRITE(FORT16,*) "  SHA: be0fba366d4485a030113da97ffc39b1a13b8879     +" 
 
      END SUBROUTINE version