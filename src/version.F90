      SUBROUTINE version()

      USE READ_DGINP, ONLY : out_direc
      USE GLOBALS,    ONLY : FORT16
      
      IMPLICIT NONE
      
      PRINT*, "Version Information"
      PRINT*, "  Branch: master" 
      PRINT*, "  SHA: b48d3531bf2461079ffe4ebe5f3b04999fdabfd2     +" 
      PRINT*, " "
      
      WRITE(FORT16,*) 'DG_WASUPP'
      WRITE(FORT16,*) ' '
      WRITE(FORT16,*) "Version Information"
      WRITE(FORT16,*) "  Branch: master" 
      WRITE(FORT16,*) "  SHA: b48d3531bf2461079ffe4ebe5f3b04999fdabfd2     +" 
 
      END SUBROUTINE version