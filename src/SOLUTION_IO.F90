!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!
      SUBROUTINE WRITE_63
      
      USE READ_DGINP, ONLY : P,NEGP
      USE GLOBALS,    ONLY : ZE,NE,MCG,CGPIV,LE,PSI,PHI,WEGP,TIME,WDFLG
      USE SIZES,      ONLY : SZ,C12
      
      IMPLICIT NONE
      INTEGER  :: I,II,J,K,L,LOCI,INFO
      REAL(SZ) :: ZE_L2(NE*P+1),ZE_IN

      PRINT '(A,F16.8)', 'Writing output for t = ', TIME
!.....Write DG output raw
      OPEN(UNIT=631,FILE='dg.63',ACCESS='APPEND')
      WRITE(631,'(F16.8)') TIME
      DO L = 1,NE
        WRITE(631,'(I,<P+1>F32.16)') WDFLG(L), (ZE(I,L,1),I=1,P+1)
      END DO
      CLOSE(631)
      
!.....Project solution into continuous framework
!.....1) Build right-hand-side vector from DG solution
!       ZE_L2(:) = 0.D0
!       DO L = 1,NE
!         J = (L-1)*P+1
!         DO K = 1,NEGP
!           ZE_IN = 0.D0
!           DO II = 1,P+1
!             ZE_IN = ZE_IN + ZE(II,L,1)*PHI(II,K)
!           END DO
!           DO I = 1,P+1
!             LOCI = J + I - 1
!             ZE_L2(LOCI) = ZE_L2(LOCI) + WEGP(K)*PSI(I,K)*ZE_IN*C12*LE(L)
!           END DO
!         END DO       
!       END DO
! !.....2) L2 projection using global Mass matrix
!       CALL DGBTRS('No transpose',NE*P+1,P,P,1,MCG,3*P+1,CGPIV,ZE_L2,NE*P+1,INFO)
! !.....3) Write solution to file
!       OPEN(UNIT=63,FILE='fort.63',ACCESS='APPEND')
!       WRITE(63,'(F16.8)') TIME
!       DO L = 1,NE*P+1,P
!         IF (ABS(ZE_L2(L)).GT.999.D0) THEN
!           ZE_L2(L) = 9999.D0
!         END IF
!         WRITE(63,'(I,F16.8)') L, ZE_L2(L)
!       END DO      
!       CLOSE(63)
      
      RETURN
      END SUBROUTINE WRITE_63
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!
      SUBROUTINE WRITE_64
      
      USE READ_DGINP, ONLY : P,NEGP
      USE GLOBALS,    ONLY : QE,NE,MCG,CGPIV,LE,PSI,PHI,WEGP,TIME
      USE SIZES,      ONLY : SZ,C12
      
      IMPLICIT NONE
      INTEGER  :: I,II,J,K,L,LOCI,INFO
      REAL(SZ) :: QE_L2(NE*P+1),QE_IN

!.....Write DG output raw
      OPEN(UNIT=641,FILE='dg.64',ACCESS='APPEND')
      WRITE(641,'(F16.8)') TIME
      DO L = 1,NE
        WRITE(641,'(I,<P+1>F32.16)') L, (QE(I,L,1),I=1,P+1)
      END DO
      CLOSE(641)
      
!.....Project solution into continuous framework
!.....1) Build right-hand-side vector from DG solution
!       QE_L2(:) = 0.D0
!       DO L = 1,NE
!         J = (L-1)*P+1
!         DO K = 1,NEGP
!           QE_IN = 0.D0
!           DO II = 1,P+1
!             QE_IN = QE_IN + QE(II,L,1)*PHI(II,K)
!           END DO
!           DO I = 1,P+1
!             LOCI = J + I - 1
!             QE_L2(LOCI) = QE_L2(LOCI) + WEGP(K)*PSI(I,K)*QE_IN*C12*LE(L)
!           END DO
!         END DO       
!       END DO
! !.....2) L2 projection using global Mass matrix
!       CALL DGBTRS('No transpose',NE*P+1,P,P,1,MCG,3*P+1,CGPIV,QE_L2,NE*P+1,INFO)
! !.....3) Write solution to file
!       OPEN(UNIT=64,FILE='fort.64',ACCESS='APPEND')
!       WRITE(64,'(F16.8)') TIME
!       DO L = 1,NE*P+1,P
!         IF (ABS(QE_L2(L)).GT.999.D0) THEN
!           QE_L2(L) = 9999.D0
!         END IF
!         WRITE(64,'(I,F16.8)') L, QE_L2(L)
!       END DO      
!       CLOSE(64)
      
      RETURN
      END SUBROUTINE WRITE_64  
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!
      SUBROUTINE INITIALIZE_OUTPUT
      
      USE READ_DGINP, ONLY : P
      USE GLOBALS,    ONLY : NE
      
      IMPLICIT NONE
      
      OPEN(UNIT=63,FILE='fort.63',STATUS='REPLACE')
      WRITE(63,'(I)') NE+1
      CLOSE(63)
      OPEN(UNIT=631,FILE='dg.63',STATUS='REPLACE')
      WRITE(631,'(I,I)') NE,P
      CLOSE(631)
      
      OPEN(UNIT=64,FILE='fort.64',STATUS='REPLACE')
      WRITE(64,'(I)') NE+1
      CLOSE(64)
      OPEN(UNIT=641,FILE='dg.64',STATUS='REPLACE')
      WRITE(641,'(I,I)') NE,P
      CLOSE(641)
      
      RETURN
      END SUBROUTINE INITIALIZE_OUTPUT
      
      
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!