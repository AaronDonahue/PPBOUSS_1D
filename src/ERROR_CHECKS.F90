!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!
      SUBROUTINE STAB_CHK
      
      USE GLOBALS,    ONLY : NE,ZE,QE,TIME,SIMERROR,FORT16
      USE READ_DGINP, ONLY : P
      
      IMPLICIT NONE
      INTEGER :: L,I
      INTEGER :: ERROR,ERRLOC(2)
      
      SIMERROR = 0
      DO L = 1,NE
        DO I = 1,P+1
          IF (ISNAN(ZE(I,L,1)).OR.ABS(ZE(I,L,1)).GT.9999.D0) THEN
            SIMERROR  = 1
            ERRLOC(1) = L
            ERRLOC(2) = I
            GOTO 100
          END IF
          IF (ISNAN(QE(I,L,1)).OR.ABS(QE(I,L,1)).GT.9999.D0) THEN
            SIMERROR  = 2
            ERRLOC(1) = L
            ERRLOC(2) = I
            GOTO 100
          END IF
        END DO
      END DO
      
100   CONTINUE

      IF (SIMERROR.NE.0) THEN
        WRITE(*,"(A,F16.8,A)"), "*** ERROR: Solution has gone unstable at t = ", TIME, " ***"  
        WRITE(*,"(A,I8,A)"), "***        Unstable at element # ",ERRLOC(1),"             ***"
        PRINT("(A)"), "!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!"
        WRITE(FORT16,"(A,F16.8,A)") "*** ERROR: Solution has gone unstable at t = ", TIME, " ***"  
        WRITE(FORT16,"(A,I8,A)") "***        Unstable at element # ",ERRLOC(1),"             ***"
        WRITE(FORT16,"(A)") "!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!"
      END IF 
                  
      END SUBROUTINE STAB_CHK
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------! 


! List of possible things for Error Checks to look for
! 1) If H0 is <= 0 then calculate H0 based on max initial wave height
