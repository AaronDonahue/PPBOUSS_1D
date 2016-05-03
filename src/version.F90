      SUBROUTINE version()

      USE READ_DGINP, ONLY : out_direc
      USE GLOBALS,    ONLY : FORT16
      
      IMPLICIT NONE
      
      PRINT*, "Version Information"
      PRINT*, "  Branch: Eddy_Visc" 
      PRINT*, "  SHA: 71b46ca14be41149af1ef9cc0a5c1267efcd5a9c     +" 
      PRINT*, " "
      
      WRITE(FORT16,*) 'DG_WASUPP'
      WRITE(FORT16,*) ' '
      WRITE(FORT16,*) "Version Information"
      WRITE(FORT16,*) "  Branch: Eddy_Visc" 
      WRITE(FORT16,*) "  SHA: 71b46ca14be41149af1ef9cc0a5c1267efcd5a9c     +" 
 
      END SUBROUTINE version