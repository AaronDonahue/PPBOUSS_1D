      PROGRAM MAIN
      
      USE GLOBALS
      USE SIZES
      USE READ_DGINP
      
      IMPLICIT NONE
      INTEGER :: L,M
      INTEGER :: SNAP
      REAL(SZ) :: CPU_START,CPU_FINISH
      
!.....Track how long simulation takes to run
      CALL CPU_TIME(CPU_START)
      
!.....Prep the DG run
      CALL PREP_DG
!.....Deal with output files
      CALL INITIALIZE_OUTPUT
!.....Run through timesteps
      TIME = 0.D0
      CALL WRITE_63
      DO SNAP = 1,TIMESTEPS
        ! Advance the solution in time
        CALL DG_TIMESTEP
        TIME = TIME + DT
        ! Write Solution
        IF (MOD(SNAP,TSNAP).EQ.0) THEN
          CALL WRITE_63
          CALL WRITE_64
        END IF
        ! Stability Check
        CALL STAB_CHK
      END DO
      
!.....Write final timestep in solution
      IF (TSNAP.NE.1.AND.MOD(TIMESTEPS,TSNAP).NE.0) THEN
        CALL WRITE_63
      END IF
      
!.....Deallocate all memory usage
      CALL VARI_DEALLOCATE
      
      
!.....Calculate time length of simulation
      CALL CPU_TIME(CPU_FINISH)
      PRINT '(A,F8.4,A)', 'Simulation took ',CPU_FINISH-CPU_START,' sec.'
      
      END PROGRAM MAIN