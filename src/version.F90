      SUBROUTINE version()

      USE READ_DGINP, ONLY : out_direc
      USE GLOBALS,    ONLY : FORT16
      
      IMPLICIT NONE
      
      PRINT*, "Version Information"
      PRINT*, "  Branch: Eddy_Visc" 
      PRINT*, "  SHA: 0375275a4f9182b46e649b5c8a021549b618d195     +" 
      PRINT*, " "
      
      WRITE(FORT16,*) 'DG_WASUPP'
      WRITE(FORT16,*) ' '
      WRITE(FORT16,*) "Version Information"
      WRITE(FORT16,*) "  Branch: Eddy_Visc" 
      WRITE(FORT16,*) "  SHA: 0375275a4f9182b46e649b5c8a021549b618d195     +" 
 
      END SUBROUTINE version