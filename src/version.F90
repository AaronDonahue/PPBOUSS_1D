      SUBROUTINE version()

      USE READ_DGINP, ONLY : out_direc
      USE GLOBALS,    ONLY : FORT16
      
      IMPLICIT NONE
      
      PRINT*, "Version Information"
      PRINT*, "  Branch: master" 
      PRINT*, "  SHA: a7c9e4071782dad4717a799f0d4746db04bcc41e     +" 
      PRINT*, " "
      
      WRITE(FORT16,*) 'DG_WASUPP'
      WRITE(FORT16,*) ' '
      WRITE(FORT16,*) "Version Information"
      WRITE(FORT16,*) "  Branch: master" 
      WRITE(FORT16,*) "  SHA: a7c9e4071782dad4717a799f0d4746db04bcc41e     +" 
 
      END SUBROUTINE version