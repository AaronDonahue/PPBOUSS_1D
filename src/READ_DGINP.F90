!==========================================================================!
!==========================================================================!
      MODULE READ_DGINP
      
      USE SIZES, ONLY : SZ
      
      TYPE :: key_val
        CHARACTER(15)           :: key      ! keyword
        REAL(SZ), POINTER       :: rptr     ! pointer to real target
        INTEGER, POINTER        :: iptr     ! pointer to integer target
        CHARACTER(100), POINTER :: cptr     ! pointer to character target      
        
        INTEGER                 :: vartype  ! target type indicator: 1=integer, 2=real, 3=character
        
        LOGICAL                 :: required ! required/optional flag
        
        INTEGER                 :: flag     ! successful read flag
      END TYPE key_val
      
      INTEGER, PARAMETER :: maxopt = 100      ! maximum allowable fort.dg options
      TYPE(key_val), DIMENSION(maxopt) :: dginp
      
      
      INTEGER                    :: nopt      ! number of valid options in dginp structure
      INTEGER, DIMENSION(maxopt) :: dginp_ind ! indicies of valid options in dginp structure
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !          Input file variables          !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !CHARACTER(100), TARGET  :: 
      !INTEGER, TARGET         ::
      !REAL(SZ), TARGET        ::
      
      ! Characters
      CHARACTER(100), TARGET  :: grid_file
      CHARACTER(100), TARGET  :: hotstart_file
      CHARACTER(100), TARGET  :: output_file,out_direc
      CHARACTER(100), TARGET  :: dgbasis
      CHARACTER(100), TARGET  :: boundtype
      CHARACTER(100), TARGET  :: breakmodel
      CHARACTER(100), TARGET  :: station_file
      CHARACTER(100), TARGET  :: nodalattr_file
      ! Integers
      INTEGER, TARGET         :: P
      INTEGER, TARGET         :: NEGP
      INTEGER, TARGET         :: IHOT
      INTEGER, TARGET         :: NRK
      INTEGER, TARGET         :: IWET
      INTEGER, TARGET         :: ISLP
      INTEGER, TARGET         :: INONHYDRO
      INTEGER, TARGET         :: NWP
      INTEGER, TARGET         :: IBREAK
      INTEGER, TARGET         :: NONHYDRO_EXT
      INTEGER, TARGET         :: EDDY_VISCOSITY
      ! Real Numbers
      REAL(SZ), TARGET        :: MaxTime
      REAL(SZ), TARGET        :: CFL_Adj
      REAL(SZ), TARGET        :: TIMESNAP
      REAL(SZ), TARGET        :: H0
      REAL(SZ), TARGET        :: SLOPEM
      REAL(SZ), TARGET        :: STNSNAP

      CONTAINS
!..........................................................................!
!==========================================================================!
!..........................................................................!
      SUBROUTINE read_input()
      
      USE GLOBALS, ONLY : FORT16

      IMPLICIT NONE
      
      INTEGER :: i,j,opt
      INTEGER :: read_stat
      INTEGER :: opt_read
      INTEGER :: comment,blank
      INTEGER :: eqind,exind
      LOGICAL :: found      
      CHARACTER(100) :: temp,line,empty
      CHARACTER(15)  :: test_opt
      CHARACTER(100) :: test_val   
      
      LOGICAL :: file_exists
      INTEGER :: skipped
      INTEGER :: read_file
      
      
      ! initialize the dginp option structure
      CALL dginp_setup()
      
      INQUIRE(FILE='./'//'fort.wasupp', EXIST = file_exists)
      IF(file_exists == .FALSE.) THEN
        PRINT*, "fort.wasupp file does not exist"
        WRITE(FORT16,*) "fort.wasupp file does not exist"
        CALL EXIT
      ENDIF
      OPEN(UNIT=15,FILE='./'//'fort.wasupp')
      opt_read = 0
      skipped = 0 
      read_file = 0
      DO WHILE (0 < 1)
      
        READ(15,"(A100)",IOSTAT=read_stat) temp
        IF(read_stat /= 0) THEN                    ! check for end-of-file
          EXIT
        ENDIF
        
        IF ( INDEX(temp,"!") == 1 .or. LEN(TRIM(temp)) == 0) THEN
            skipped = skipped + 1
        ELSE  
          read_file = 1
          EXIT
        ENDIF
        
      ENDDO
      CLOSE(15)
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "             Input Information               "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "
      WRITE(FORT16,"(A)") "---------------------------------------------"
      WRITE(FORT16,"(A)") "             Input Information               "
      WRITE(FORT16,"(A)") "---------------------------------------------"
      WRITE(FORT16,"(A)") " "
      IF (read_file.ne.1) THEN 
        PRINT*, "ERROR: fort.wasupp does not contain any information"
        WRITE(FORT16,*) "ERROR: fort.wasupp does not contain any information"
        CALL EXIT
      ENDIF
      
      
      opt_read = 0
      comment = 0 
      blank = 0
      
      
      OPEN(25,FILE='./'//'fort.wasupp',POSITION="rewind")        
     
      
      DO WHILE (opt_read < nopt)
      
        READ(25,"(A100)",IOSTAT=read_stat) temp
        IF(read_stat /= 0) THEN                    ! check for end-of-file
          EXIT
        ENDIF
        
        line = ADJUSTL(temp)
        
        IF(INDEX(line,"!") == 1) THEN              ! lines beginning with ! are skipped
        
          comment = comment + 1
          
        ELSE IF (LEN(TRIM(line)) == 0) THEN        ! blank lines are skipped
        
          blank = blank + 1
          
        ELSE
  
          ! determine keyword and assignment value
          eqind = INDEX(line,"=")
          exind = INDEX(line,"!")
          test_opt = line(1:eqind-1)
          IF (exind > 0) THEN                          ! handle trailing comment 
            test_val = ADJUSTL(line(eqind+1:exind-1))  ! (only necessary if there is no space between value and the !)
          ELSE
            test_val = ADJUSTL(line(eqind+1:))
          ENDIF         
          
          ! Look for a match for the keyword
          found = .false.
    test: DO opt = 1,nopt
    
            i = dginp_ind(opt)    
    
            IF (test_opt == dginp(i)%key) THEN
              
              ! Set variables equal to value from fort.dg through pointer using an internal read
              SELECT CASE (dginp(i)%vartype) 
                CASE(1)
                  READ(test_val,*) dginp(i)%iptr
                  WRITE(FORT16,"(A,I3,A,A,A,I8)") "              ",i,") ",dginp(i)%key," = ",dginp(i)%iptr
                CASE(2)
                  READ(test_val,*) dginp(i)%rptr
                  WRITE(FORT16,"(A,I3,A,A,A,E21.8)") "              ",i,") ",dginp(i)%key," = ",dginp(i)%rptr                  
                CASE(3)
                  dginp(i)%cptr = TRIM(test_val)
                  WRITE(FORT16,"(A,I3,A,A,A,A)") "              ",i,") ",dginp(i)%key," = ",dginp(i)%cptr                  
              END SELECT

              found = .true.          ! flag match
              opt_read = opt_read + 1
              dginp(i)%flag = 1      ! flag option as found
              
              EXIT test
              
            ENDIF
          ENDDO test
                    
          IF (found == .false. .and. eqind > 0) THEN
            ! unmatched lines with an equal sign are either incorrect or no longer supported
            PRINT("(3A)"),"*** WARNING: ",test_opt, " is an incorrect or depreciated value ***"
            WRITE(FORT16,"(3A)"),"*** WARNING: ",test_opt, " is an incorrect or depreciated value ***"            
          ELSE IF (found == .false.) THEN
            ! unmatched lines without an equal sign are ignored
            PRINT("(A)"), "*** WARNING: non-comment line does not contain a keyword assignment***"
            WRITE(FORT16,"(A)"), "*** WARNING: non-comment line does not contain a keyword assignment***"
          ENDIF
          
        ENDIF
      ENDDO 
      
      PRINT*, ""
      WRITE(FORT16,*) ""
     
      CALL check_errors(opt_read)
      
      PRINT*, ""
      WRITE(FORT16,*) ""
      CLOSE(25)
            
      END SUBROUTINE read_input
      
!..........................................................................!
!==========================================================================!
!..........................................................................!
     SUBROUTINE check_errors(opt_read)
     
     USE GLOBALS, ONLY : FORT16
     
     IMPLICIT NONE
     
     INTEGER :: i,j,opt
     INTEGER :: opt_read
     INTEGER :: quit
     
     IF(opt_read /= nopt) THEN

       ! check for required options that are unspecifed 
       quit = 0
       DO opt = 1,nopt
         i = dginp_ind(opt)
         IF (dginp(i)%flag == 0 .and. dginp(i)%required == .true.) THEN
           quit = 1   ! flag fatal error
         ENDIF
       ENDDO
        
       IF (quit == 1) THEN
        
          PRINT("(A)"), "*** ERROR: There are missing required options in the fort.wasupp file ***"  
          PRINT("(A)"), "           The following options must be specified: "
          WRITE(FORT16,"(A)") "*** ERROR: There are missing required options in the fort.wasupp file ***"  
          WRITE(FORT16,"(A)") "           The following options must be specified: "
          j = 0        
          DO opt = 1,nopt
            i = dginp_ind(opt)
            IF (dginp(i)%flag == 0 .and. dginp(i)%required == .true.) THEN
              j = j+1
              PRINT "(A,I3,2A)", "              ",j,") ",dginp(i)%key
              WRITE(FORT16,"(A,I3,2A)") "              ",j,") ",dginp(i)%key
            ENDIF
          ENDDO          
          
          PRINT("(A)"), "!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!"
          WRITE(FORT16,"(A)"), "!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!"
          STOP
          
       ELSE
        
          PRINT("(A)"), "*** WARNING: There are missing optional options in the fort.wasupp file ***"
          PRINT("(A)"), "             The following default values will be used: "
          WRITE(FORT16,"(A)") "*** WARNING: There are missing optional options in the fort.wasupp file ***"
          WRITE(FORT16,"(A)") "             The following default values will be used: " 
          j = 0        
          DO opt = 1,nopt
            i = dginp_ind(opt)
            IF (dginp(i)%flag == 0 .and. dginp(i)%required == .false.) THEN
              
              j = j+1
              SELECT CASE (dginp(i)%vartype) 
                CASE(1)
                  PRINT("(A,I3,A,A,A,I8)"),     "              ",j,") ",dginp(i)%key," = ",dginp(i)%iptr
                  WRITE(FORT16,"(A,I3,A,A,A,I8)") "              ",j,") ",dginp(i)%key," = ",dginp(i)%iptr
                CASE(2)
                  PRINT("(A,I3,A,A,A,E21.8)"),  "              ",j,") ",dginp(i)%key," = ",dginp(i)%rptr
                  WRITE(FORT16,"(A,I3,A,A,A,E21.8)") "              ",j,") ",dginp(i)%key," = ",dginp(i)%rptr                  
                CASE(3)
                  PRINT("(A,I3,A,A,A,A)"),      "              ",j,") ",dginp(i)%key," = ",dginp(i)%cptr
                  WRITE(FORT16,"(A,I3,A,A,A,A)") "              ",j,") ",dginp(i)%key," = ",dginp(i)%cptr                  
              END SELECT
              
            ENDIF
          ENDDO 
          
          PRINT("(A)"), '!!!!!! EXECUTION WILL CONTINUE !!!!!!!!'
          WRITE(FORT16,"(A)") '!!!!!! EXECUTION WILL CONTINUE !!!!!!!!'
          
       ENDIF       
                  
     ENDIF        
     
     
     RETURN
     END SUBROUTINE check_errors

!..........................................................................!
!==========================================================================!
!..........................................................................!
      SUBROUTINE dginp_setup()
      
      USE GLOBALS, ONLY : FORT16
      
      IMPLICIT NONE
      
      INTEGER :: i
      INTEGER :: ncheck
      CHARACTER(15) :: empty
      
!..... Initialize dginp structure
      DO i = 1,maxopt
        NULLIFY(dginp(i)%iptr)
        NULLIFY(dginp(i)%rptr)
        NULLIFY(dginp(i)%cptr)
        
        dginp(i)%key = empty
        dginp(i)%flag = 0        
      ENDDO
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Configure fort.dg options here:
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !    keywords                         target variables                      requirement                 default values
      dginp(1)%key = "grid_file";       dginp(1)%cptr => grid_file ;      dginp(1)%required = .true.;      dginp(1)%cptr = ""
      dginp(2)%key = "hotstart_file";   dginp(2)%cptr => hotstart_file ;  dginp(2)%required = .false.;     dginp(2)%cptr = ""
      dginp(3)%key = "p";               dginp(3)%iptr => P;               dginp(3)%required = .true.;      dginp(3)%iptr = 1
      dginp(4)%key = "negp";            dginp(4)%iptr => NEGP;            dginp(4)%required = .false.;     dginp(4)%iptr = 8
      dginp(5)%key = "ihot";            dginp(5)%iptr => IHOT;            dginp(5)%required = .false.;     dginp(5)%iptr = 0
      dginp(6)%key = "nrk";             dginp(6)%iptr => NRK;             dginp(6)%required = .false.;     dginp(6)%iptr = 3
      dginp(7)%key = "maxtime";         dginp(7)%rptr => MAXTIME;         dginp(7)%required = .true.;      dginp(7)%rptr = 10.d0
      dginp(8)%key = "cfl";             dginp(8)%rptr => CFL_ADJ;         dginp(8)%required = .false.;     dginp(8)%rptr = 1.d0
      dginp(9)%key = "out_direc";       dginp(9)%cptr => out_direc;       dginp(9)%required = .false.;     dginp(9)%cptr = "./"
      dginp(10)%key = "dgbasis";        dginp(10)%cptr => dgbasis;        dginp(10)%required = .false.;    dginp(10)%cptr = "nodal"
      dginp(11)%key = "timesnap";       dginp(11)%rptr => TIMESNAP;       dginp(11)%required = .false.;    dginp(11)%rptr = -999.d0
      dginp(12)%key = "iwet";           dginp(12)%iptr => IWET;           dginp(12)%required = .false.;    dginp(12)%iptr = 0
      dginp(13)%key = "h0";             dginp(13)%rptr => H0;             dginp(13)%required = .false.;    dginp(13)%rptr = 0.d0
      dginp(14)%key = "boundarytype";   dginp(14)%cptr => BOUNDTYPE;      dginp(14)%required = .false.;    dginp(14)%cptr = "reflective"
      dginp(15)%key = "islp";           dginp(15)%iptr => ISLP;           dginp(15)%required = .false.;    dginp(15)%iptr = 0
      dginp(16)%key = "islpconstant";   dginp(16)%rptr => SLOPEM;         dginp(16)%required = .false.;    dginp(16)%rptr = 0.d0
      dginp(17)%key = "inonhydro";      dginp(17)%iptr => INONHYDRO;      dginp(17)%required = .false.;    dginp(17)%iptr = 0
      dginp(18)%key = "nwp";            dginp(18)%iptr => NWP;            dginp(18)%required = .false.;    dginp(18)%iptr = 0
      dginp(19)%key = "ibreak";         dginp(19)%iptr => IBREAK;         dginp(19)%required = .false.;    dginp(19)%iptr = 0 
      dginp(20)%key = "breakingmodel";  dginp(20)%cptr => BREAKMODEL;     dginp(20)%required = .false.;    dginp(20)%cptr = "none" 
      dginp(21)%key = "station_file";   dginp(21)%cptr => station_file;   dginp(21)%required = .false.;    dginp(21)%cptr = "none"
      dginp(22)%key = "station_timestep"; dginp(22)%rptr => STNSNAP;      dginp(22)%required = .false.;    dginp(22)%rptr = "-999.d0"
      dginp(23)%key = "nodalattr_file"; dginp(23)%cptr => nodalattr_file; dginp(23)%required = .false.;    dginp(23)%cptr = "fort.13"
      dginp(24)%key = "nonhydro_extrapolation"; dginp(24)%iptr => NONHYDRO_EXT; dginp(24)%required = .false.;    dginp(24)%iptr = 0
      dginp(25)%key = "eddy_viscosity"; dginp(25)%iptr => EDDY_VISCOSITY; dginp(25)%required = .false.;    dginp(25)%iptr = 0
      
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! End configuration
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      nopt = 0
      ncheck = 0
      DO i = 1,maxopt
      
        ! find and keep track of populated indicies
        IF (dginp(i)%key .ne. empty) THEN      
          nopt = nopt + 1      
          dginp_ind(nopt) = i
        ENDIF
        
        ! determine target variable type by checking association status
        dginp(i)%vartype = 0    
        
        IF (ASSOCIATED(dginp(i)%iptr)) THEN  ! integer
          ncheck = ncheck + 1   
          dginp(i)%vartype = 1
        ENDIF
        
        IF (ASSOCIATED(dginp(i)%rptr)) THEN ! real
          ncheck = ncheck + 1
          dginp(i)%vartype = 2
        ENDIF
        
        IF (ASSOCIATED(dginp(i)%cptr)) THEN ! character
          ncheck = ncheck + 1        
          dginp(i)%vartype = 3        
        ENDIF
      ENDDO
      
      ! ensure user has associated each keyword pointer
      IF (nopt /= ncheck) THEN
        PRINT("(A)"), "*** ERROR: fort.wasupp option pointer association error ***"
        PRINT("(A)"), "           check keyword configuration in dginp_setup subroutine"
        WRITE(FORT16,"(A)") "*** ERROR: fort.wasupp option pointer association error ***"
        WRITE(FORT16,"(A)") "           check keyword configuration in dginp_setup subroutine"
        STOP
      ENDIF
      
      RETURN
      END SUBROUTINE dginp_setup
!..........................................................................!
!==========================================================================!
!..........................................................................!

      END MODULE READ_DGINP
!==========================================================================!
!==========================================================================!
