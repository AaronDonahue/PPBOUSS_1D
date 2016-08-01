      SUBROUTINE version()

      USE READ_DGINP, ONLY : out_direc
      USE GLOBALS,    ONLY : FORT16
      
      IMPLICIT NONE
      
      PRINT*, "Version Information"
      PRINT*, "  Branch: Eddy_Visc" 
      PRINT*, "  SHA: 3928f8a312699d8f49dd40ada8cbfafd042fa958     +" 
      PRINT*, " "
      
      WRITE(FORT16,*) 'DG_WASUPP'
      WRITE(FORT16,*) ' '
      WRITE(FORT16,*) "Version Information"
      WRITE(FORT16,*) "  Branch: Eddy_Visc" 
      WRITE(FORT16,*) "  SHA: 3928f8a312699d8f49dd40ada8cbfafd042fa958     +" 
 
      END SUBROUTINE version