      SUBROUTINE version()

      USE READ_DGINP, ONLY : out_direc
      USE GLOBALS,    ONLY : FORT16
      
      IMPLICIT NONE
      
      PRINT*, "Version Information"
      PRINT*, "  Branch: Eddy_Visc" 
      PRINT*, "  SHA: 7d576654f4bc4a3fdd13202671de5cd21400b8cd     +" 
      PRINT*, " "
      
      WRITE(FORT16,*) 'DG_WASUPP'
      WRITE(FORT16,*) ' '
      WRITE(FORT16,*) "Version Information"
      WRITE(FORT16,*) "  Branch: Eddy_Visc" 
      WRITE(FORT16,*) "  SHA: 7d576654f4bc4a3fdd13202671de5cd21400b8cd     +" 
 
      END SUBROUTINE version