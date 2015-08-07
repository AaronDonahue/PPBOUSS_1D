      PROGRAM MAIN
      
      USE GLOBALS
      USE SIZES
      USE READ_DGINP
      
      IMPLICIT NONE
      INTEGER :: L,M
      INTEGER :: SNAP
      
!.....Prep the DG run
      CALL PREP_DG
!.....Deal with output files
      CALL INITIALIZE_OUTPUT
!.....Print version information
      CALL version 
!.....Run through timesteps
      TIME = 0.D0
      CALL WRITE_63
      DO SNAP = 1,TIMESTEPS
        ! Advance the solution in time
        CALL DG_TIMESTEP
        TIME = TIME + DT
        ! Write Solution
        IF (MOD(SNAP,TSNAP).EQ.0) THEN
          CALL WRITE_OUTPUT_GLOBAL
        END IF
        ! Stability Check
        CALL STAB_CHK
        IF (SIMERROR.NE.0) THEN
          GOTO 9999
        END IF
      END DO

9999  CONTINUE
!.....Write final timestep in solution
      IF (MOD(TIMESTEPS,TSNAP).NE.0) THEN
        CALL WRITE_OUTPUT_GLOBAL
      END IF

!.....Deallocate all memory usage
      CALL VARI_DEALLOCATE      

!.....Close all output files, finish run
      CALL FINISH_OUTPUT
      
      END PROGRAM MAIN