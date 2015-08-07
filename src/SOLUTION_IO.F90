!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!
      SUBROUTINE WRITE_OUTPUT_GLOBAL
      
      USE GLOBALS, ONLY : TIME,FORT16
      
      IMPLICIT NONE
      
      PRINT '(A,F32.16)', 'Writing global output for t = ', TIME
      WRITE(16,'(A,F32.16)') 'Writing global output for t = ', TIME
      
      CALL WRITE_63
      CALL WRITE_64
      CALL WRITE_73
      CALL WRITE_74
      
      RETURN
      END SUBROUTINE WRITE_OUTPUT_GLOBAL
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!
      SUBROUTINE WRITE_OUTPUT_STATION
      
      USE GLOBALS, ONLY : TIME,FORT16
      
      IMPLICIT NONE
      
      PRINT '(A,F32.16)', 'Writing station output for t = ', TIME
      WRITE(16,'(A,F32.16)') 'Writing station output for t = ', TIME
      
      CALL WRITE_61
!       CALL WRITE_62
!       CALL WRITE_71
!       CALL WRITE_72
      
      RETURN
      END SUBROUTINE WRITE_OUTPUT_STATION
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!
      SUBROUTINE WRITE_61
      
      USE READ_DGINP, ONLY : P
      USE GLOBALS,    ONLY : ZE,NE,TIME,FORT611,PHISTN,NUMSTNs,STNELEM
      USE SIZES,      ONLY : SZ
      
      IMPLICIT NONE
      INTEGER  :: I,K,L
      REAL(SZ) :: ZE_IN
      
!.....Write DG output raw
      WRITE(FORT611,'(F32.16)') TIME
      DO L = 1,NUMSTNs
        ZE_IN = 0.D0
        DO I = 1,P+1        
          ZE_IN = ZE_IN + ZE(I,STNELEM(L),1)*PHISTN(L,I)
        END DO
        WRITE(FORT611,'(I,F32.16)') L, ZE_IN
      END DO

      RETURN
      END SUBROUTINE WRITE_61
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!
      SUBROUTINE WRITE_63
      
      USE READ_DGINP, ONLY : P
      USE GLOBALS,    ONLY : ZE,NE,TIME,FORT631
      
      IMPLICIT NONE
      INTEGER  :: I,L
      
!.....Write DG output raw
      WRITE(FORT631,'(F32.16)') TIME
      DO L = 1,NE
        WRITE(FORT631,'(I,<P+1>F32.16)') L, (ZE(I,L,1),I=1,P+1)
      END DO

      RETURN
      END SUBROUTINE WRITE_63
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!
      SUBROUTINE WRITE_64
      
      USE READ_DGINP, ONLY : P
      USE GLOBALS,    ONLY : QE,NE,TIME,FORT641
      
      IMPLICIT NONE
      INTEGER  :: I,L

!.....Write DG output raw
      WRITE(FORT641,'(F32.16)') TIME
      DO L = 1,NE
        WRITE(FORT641,'(I,<P+1>F32.16)') L, (QE(I,L,1),I=1,P+1)
      END DO
      
      RETURN
      END SUBROUTINE WRITE_64
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!
      SUBROUTINE WRITE_73
      
      USE READ_DGINP, ONLY : P
      USE GLOBALS,    ONLY : PD,NE,TIME,FORT731
      
      IMPLICIT NONE
      INTEGER  :: I,L

!.....Write DG output raw
      WRITE(FORT731,'(F32.16)') TIME
      DO L = 1,NE
        WRITE(FORT731,'(I,<P+1>F32.16)') L, (PD(I,L,1),I=1,P+1)
      END DO
      
      RETURN
      END SUBROUTINE WRITE_73
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!
      SUBROUTINE WRITE_74
      
      USE READ_DGINP, ONLY : P
      USE GLOBALS,    ONLY : PB,NE,TIME,FORT741
      
      IMPLICIT NONE
      INTEGER  :: I,L

!.....Write DG output raw
      WRITE(FORT741,'(F32.16)') TIME
      DO L = 1,NE
        WRITE(FORT741,'(I,<P+1>F32.16)') L, (PB(I,L,1),I=1,P+1)
      END DO
      
      RETURN
      END SUBROUTINE WRITE_74
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!
      SUBROUTINE INITIALIZE_OUTPUT
      
      USE READ_DGINP, ONLY : P,OUT_DIREC,station_file,dgbasis
      USE GLOBALS,    ONLY : NE,CPU_START,FORT16,FORT631,FORT641,FORT731,  &
     &                       FORT741,FORT611,NUMSTNS,PHISTN,STNELEM,X,LE
      USE SIZES,      ONLY : SZ
      
      IMPLICIT NONE
      INTEGER  :: I,S,L
      REAL(SZ) :: STNX,STNXI,GETMODAL,GETNODAL
      
!.....Track how long simulation takes to run
      CALL CPU_TIME(CPU_START)  
!.....Output files      
      FORT16 = 16
      OPEN(UNIT=FORT16,FILE=TRIM(out_direc)//'fort.16',STATUS='REPLACE')
      
      FORT631 = 631
      OPEN(UNIT=FORT631,FILE=TRIM(out_direc)//'dg.63',STATUS='REPLACE')
      WRITE(FORT631,'(I,I)') NE, P
      
      FORT641 = 641
      OPEN(UNIT=FORT641,FILE=TRIM(out_direc)//'dg.64',STATUS='REPLACE')
      WRITE(FORT641,'(I,I)') NE,P
      
      FORT731 = 731
      OPEN(UNIT=FORT731,FILE=TRIM(out_direc)//'dg.73',STATUS='REPLACE')
      WRITE(FORT731,'(I,I)') NE,P
      
      FORT741 = 741
      OPEN(UNIT=FORT741,FILE=TRIM(out_direc)//'dg.74',STATUS='REPLACE')
      WRITE(FORT741,'(I,I)') NE,P
      
      IF (STATION_FILE.NE."none") THEN
        OPEN(UNIT=99,FILE=TRIM(station_file),ACTION='READ')
        READ(99,*) NUMSTNS
        ALLOCATE(PHISTN(NUMSTNS,P+1),STNELEM(NUMSTNS))
        DO S = 1,NUMSTNS
          READ(99,*) STNX
          DO L = 1,NE           
            IF (X(L).LT.STNX.AND.X(L+1).GT.STNX) THEN
              STNELEM(S) = L
              STNXI = -1.D0 + 2.D0*(STNX-X(L))/LE(L)
              IF (ADJUSTL(TRIM(DGBASIS)).EQ.'NODAL'.OR.                            &
     &             ADJUSTL(TRIM(DGBASIS)).EQ.'nodal') THEN
                DO I = 1,P+1
                  PHISTN(S,I) = GETNODAL(P,I,STNXI,0)
                END DO
              ELSE
                DO I = 1,P+1
                  PHISTN(S,I) = GETMODAL(P,I,STNXI,0)
                END DO
              END IF
              GOTO 6100
            END IF
          END DO
          print*, "Station # ", S, "is not in the mesh"
          STNELEM(S) = 1
          DO I = 1,P+1
            PHISTN(S,I) = 0.D0
          END DO
6100      CONTINUE
        END DO
      END IF
      
      RETURN
      END SUBROUTINE INITIALIZE_OUTPUT
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!
      SUBROUTINE FINISH_OUTPUT
      
      USE GLOBALS, ONLY : CPU_START,CPU_FINISH,FORT16,FORT631,FORT641,     &
     &                    FORT731,FORT741
      
      IMPLICIT NONE
      
      CALL CPU_TIME(CPU_FINISH)
      PRINT '(A,F32.16,A)', 'Simulation took ',CPU_FINISH-CPU_START,' sec.'
      
      WRITE(16,*) " "
      WRITE(16,*) "Simulation finished..."
      WRITE(16,*) "Computation took ",CPU_FINISH-CPU_START," sec."
      CLOSE(FORT16)
      
      CLOSE(FORT631)
      CLOSE(FORT641)
      CLOSE(FORT731)
      CLOSE(FORT741)
      
      RETURN
      END SUBROUTINE FINISH_OUTPUT
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!

! Things to do for I/O
! 2) Setup Calculate Runup subroutine
