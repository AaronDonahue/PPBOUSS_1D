!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!
      SUBROUTINE WRITE_OUTPUT_GLOBAL
      
      USE GLOBALS,    ONLY : TIME,FORT16
      USE READ_DGINP, ONLY : INONHYDRO
      
      IMPLICIT NONE
      
      PRINT '(A,F32.16)', 'Writing global output for t = ', TIME
      WRITE(16,'(A,F32.16)') 'Writing global output for t = ', TIME
      
      CALL WRITE_63
      CALL WRITE_64
      CALL WRITE_RUNUP
      IF (INONHYDRO.GT.0) THEN
        CALL WRITE_73
        CALL WRITE_74
      END IF
      
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
      USE GLOBALS,    ONLY : ZE,NE,TIME,FORT611,PHISTN,NUMSTNS,STNELEM,    &
     &                       WDFLG
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
        WRITE(FORT611,'(I,F32.16,I)') L, ZE_IN, WDFLG(L)
      END DO

      RETURN
      END SUBROUTINE WRITE_61
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!
      SUBROUTINE WRITE_63
      
      USE READ_DGINP, ONLY : P
      USE GLOBALS,    ONLY : ZE,NE,TIME,FORT631,WDFLG
      
      IMPLICIT NONE
      INTEGER  :: I,L
      
!.....Write DG output raw
      WRITE(FORT631,'(F32.16)') TIME
      DO L = 1,NE
        WRITE(FORT631,'(I,<P+1>F32.16,I)') L, (ZE(I,L,1),I=1,P+1), WDFLG(L)
      END DO

      RETURN
      END SUBROUTINE WRITE_63
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!
      SUBROUTINE WRITE_64
      
      USE READ_DGINP, ONLY : P
      USE GLOBALS,    ONLY : QE,NE,TIME,FORT641,WDFLG
      
      IMPLICIT NONE
      INTEGER  :: I,L

!.....Write DG output raw
      WRITE(FORT641,'(F32.16)') TIME
      DO L = 1,NE
        WRITE(FORT641,'(I,<P+1>F32.16,I)') L, (QE(I,L,1),I=1,P+1),WDFLG(L)
      END DO
      
      RETURN
      END SUBROUTINE WRITE_64
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!
      SUBROUTINE WRITE_73
      
      USE READ_DGINP, ONLY : P
      USE GLOBALS,    ONLY : PD,NE,TIME,FORT731,DISPFLG
      
      IMPLICIT NONE
      INTEGER  :: I,L

!.....Write DG output raw
      WRITE(FORT731,'(F32.16)') TIME
      DO L = 1,NE
        WRITE(FORT731,'(I,<P+1>F32.16,I)') L, (PD(I,L,1),I=1,P+1), DISPFLG(L)
      END DO
      
      RETURN
      END SUBROUTINE WRITE_73
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!
      SUBROUTINE WRITE_74
      
      USE READ_DGINP, ONLY : P
      USE GLOBALS,    ONLY : PB,NE,TIME,FORT741,DISPFLG
      
      IMPLICIT NONE
      INTEGER  :: I,L

!.....Write DG output raw
      WRITE(FORT741,'(F32.16)') TIME
      DO L = 1,NE
        WRITE(FORT741,'(I,<P+1>F32.16,I)') L, (PB(I,L,1),I=1,P+1), DISPFLG(L)
      END DO
      
      RETURN
      END SUBROUTINE WRITE_74
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!
      SUBROUTINE WRITE_RUNUP
      
      USE GLOBALS,    ONLY : WDFLG,NE,ZE,PHIB,MAXRUNUP,TIME,FORTRUNUP,X
      USE READ_DGINP, ONLY : P
      USE SIZES,      ONLY : SZ
      
      IMPLICIT NONE
      INTEGER  :: I,L,XS
      REAL(SZ) :: ZE_IN,RUNLOC
      
!.....First find shoreline
      XS = 0
      DO L = 2,NE
        IF (WDFLG(L-1).EQ.1.AND.WDFLG(L).EQ.0) THEN
          XS = L
          EXIT
        END IF
      END DO
!.....If the shoreline is found then calculate runup
      RUNLOC = 0.D0
      IF (XS.NE.0) THEN
        ZE_IN = 0.D0
        DO I = 1,P+1
          ZE_IN = ZE_IN + ZE(I,XS,1)*PHIB(I,2)
        END DO
        MAXRUNUP = MAX(MAXRUNUP,ZE_IN)
        RUNLOC = MAX(RUNLOC,ZE_IN)
      ELSE
        ZE_IN = -9999.D0
      END IF
!.....Write runup results
      WRITE(FORTRUNUP,'(3F32.16)') TIME, X(XS+1), RUNLOC
      
      RETURN
      END SUBROUTINE
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!
      SUBROUTINE INITIALIZE_OUTPUT
      
      USE READ_DGINP, ONLY : P,OUT_DIREC,station_file,dgbasis,INONHYDRO
      USE GLOBALS,    ONLY : NE,CPU_START,FORT631,FORT641,FORT731,FORT611, &
     &                       FORT741,FORT611,NUMSTNS,PHISTN,STNELEM,X,LE,  &
     &                       FORT16,FORTRUNUP,MAXRUNUP
      USE SIZES,      ONLY : SZ
      
      IMPLICIT NONE
      INTEGER  :: I,S,L
      REAL(SZ) :: STNX,STNXI,GETMODAL,GETNODAL
      
!.....Track how long simulation takes to run
      CALL CPU_TIME(CPU_START)  
!.....Output files      
      FORT631 = 631
      OPEN(UNIT=FORT631,FILE=TRIM(out_direc)//'dg.63',STATUS='REPLACE')
      WRITE(FORT631,'(I,I)') NE, P
      
      FORT641 = 641
      OPEN(UNIT=FORT641,FILE=TRIM(out_direc)//'dg.64',STATUS='REPLACE')
      WRITE(FORT641,'(I,I)') NE,P
      ! Pressure Files
      IF (INONHYDRO.GT.0) THEN
        FORT731 = 731
        OPEN(UNIT=FORT731,FILE=TRIM(out_direc)//'dg.73',STATUS='REPLACE')
        WRITE(FORT731,'(I,I)') NE,P
            
        FORT741 = 741
        OPEN(UNIT=FORT741,FILE=TRIM(out_direc)//'dg.74',STATUS='REPLACE')
        WRITE(FORT741,'(I,I)') NE,P
      END IF
      ! Station File
      IF (STATION_FILE.NE."none") THEN        
        PRINT*, ' '
        PRINT*, 'Reading Station input file'
        WRITE(FORT16,*) ' '
        WRITE(FORT16,*) 'Reading Station input file'
        OPEN(UNIT=99,FILE=TRIM(station_file),ACTION='READ')
        READ(99,*) NUMSTNS
        ALLOCATE(PHISTN(NUMSTNS,P+1),STNELEM(NUMSTNS))
        DO S = 1,NUMSTNS
          READ(99,*) STNX
          DO L = 1,NE           
            IF (X(L).LE.STNX.AND.X(L+1).GT.STNX) THEN
              STNELEM(S) = L
              STNXI = -1.D0 + 2.D0*(STNX-X(L))/LE(L)
              IF (ADJUSTL(TRIM(DGBASIS)).EQ.'NODAL'.OR.                    &
     &             ADJUSTL(TRIM(DGBASIS)).EQ.'nodal') THEN
                DO I = 1,P+1
                  PHISTN(S,I) = GETNODAL(P,I,STNXI,0)
                END DO
              ELSE
                DO I = 1,P+1
                  PHISTN(S,I) = GETMODAL(P,I,STNXI,0)
                END DO
              END IF
              print*, "Station # ", S," is in element ",STNELEM(S)
              WRITE(FORT16,"(A,I,A,I)") "Station # ", S," is in element ",STNELEM(S)
              GOTO 6100
            END IF
          END DO
          print*, "Station # ", S, "is not in the mesh"
          WRITE(FORT16,*) "Station # ", S, "is not in the mesh"
          STNELEM(S) = 1
          DO I = 1,P+1
            PHISTN(S,I) = 0.D0
          END DO
6100      CONTINUE
        END DO
        FORT611 = 61
        OPEN(UNIT=FORT611,FILE=TRIM(out_direc)//'dg.61',STATUS='REPLACE')
        WRITE(FORT611,'(I)') NUMSTNS
      END IF
      ! Runup file
      MAXRUNUP = -9999.D0
      FORTRUNUP = 51
      OPEN(UNIT=FORTRUNUP,FILE=TRIM(out_direc)//'runup.51',STATUS='REPLACE')
      
      RETURN
      END SUBROUTINE INITIALIZE_OUTPUT
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!
      SUBROUTINE FINISH_OUTPUT
      
      USE GLOBALS,    ONLY : CPU_START,CPU_FINISH,FORT16,FORT631,FORT641,  &
     &                    FORT731,FORT741,FORT611,FORTRUNUP,MAXRUNUP
      USE READ_DGINP, ONLY : INONHYDRO,STATION_FILE
      
      IMPLICIT NONE
      
      CALL CPU_TIME(CPU_FINISH)
      PRINT '(A,F32.16,A)', 'Simulation took ',CPU_FINISH-CPU_START,' sec.'
      
      CLOSE(FORT631)
      CLOSE(FORT641)
      
      IF (INONHYDRO.GT.0) THEN
        CLOSE(FORT731)
        CLOSE(FORT741)
      END IF
      
      IF (STATION_FILE.NE."none") THEN
        CLOSE(FORT611)
      END IF
      
      IF (MAXRUNUP.GT.-9999.D0) THEN
        WRITE(FORT16,'(A,F32.16)') 'Maximum runup = ',MAXRUNUP
      END IF
      CLOSE(FORTRUNUP)
      
      RETURN
      END SUBROUTINE FINISH_OUTPUT
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!

! Things to do for I/O

