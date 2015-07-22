!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!
      SUBROUTINE STAB_CHK
      
      USE GLOBALS,    ONLY : NE,ZE,QE,TIME
      USE READ_DGINP, ONLY : P
      
      IMPLICIT NONE
      INTEGER :: L,I
      INTEGER :: ERROR,ERRLOC(2)
      
      ERROR = 0
      DO L = 1,NE
        DO I = 1,P+1
          IF (ISNAN(ZE(I,L,1)).OR.ABS(ZE(I,L,1)).GT.9999.D0) THEN
            ERROR     = 1
            ERRLOC(1) = L
            ERRLOC(2) = I
            GOTO 100
          END IF
          IF (ISNAN(QE(I,L,1)).OR.ABS(QE(I,L,1)).GT.9999.D0) THEN
            ERROR     = 2
            ERRLOC(1) = L
            ERRLOC(2) = I
            GOTO 100
          END IF
        END DO
      END DO
      
100   CONTINUE

      IF (ERROR.NE.0) THEN
        WRITE(*,"(A,F16.8,A)"), "*** ERROR: Solution has gone unstable at t = ", TIME, " ***"  
        WRITE(*,"(A,I8,A)"), "***        Unstable at element # ",ERRLOC(1),"             ***"
        PRINT("(A)"), "!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!"
        
        CALL WRITE_63
        STOP
      END IF 
                  
      END SUBROUTINE STAB_CHK
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!      