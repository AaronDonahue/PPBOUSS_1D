!==========================================================================!
!==========================================================================!
      SUBROUTINE DG_TIMESTEP
      
      USE GLOBALS,    ONLY : IRK,ZE_RHS,QE_RHS,ZE,QE,ATVD,BTVD,DT,NE,      &
     &                       WDFLG,TTVD,DT,TIME,TIME_RK,DISPFLG,PD,PB,     &
     &                       PDLVL,PBLVL
      USE SIZES,      ONLY : SZ
      USE READ_DGINP, ONLY : NRK,P,IWET,ISLP,INONHYDRO,IBREAK,NONHYDRO_EXT

      IMPLICIT NONE
      INTEGER	:: L,I,K

!---------------------------------------------------------------------------!
!.....Timestep through each of the Runge-Kutta stages
      DO IRK = 1,NRK
!.......Zero out the Right Hand Side vectors for ZE and QE for this stage
        ZE_RHS(:,:,IRK) = 0.D0
        QE_RHS(:,:,IRK) = 0.D0
!.......Time at this stage
        TIME_RK = TIME+TTVD(IRK)*DT
!.......Build RHS vectors for DG solution
        CALL DG_AREA_INTEGRALS() ! Calculate area integral component

        CALL DG_INTERNAL_EDGES() ! Calculate internal edge integral (fluxes)
  
        CALL DG_BOUNDARY_EDGES() ! Calculate boundary integrals of domain (fluxes)
        
        CALL SWE_RHS_VARIABLES() ! Calculate the addition to the RHS of the SWE
     
!.......Finish RHS vectors by multiplying by the inverse of the mass matrix
        CALL COMPLETE_RHS()
        
!.......Build solution to ZE and QE at this stage using previous stages
!.......and RHS vectors
        DO L = 1,NE
          DO I = 1,P+1
            ZE(I,L,IRK+1) = 0.D0
            QE(I,L,IRK+1) = 0.D0
            DO K = 1,IRK
              ZE(I,L,IRK+1) = ZE(I,L,IRK+1) + ATVD(IRK,K)*ZE(I,L,K) +      &
     &                          DT*BTVD(IRK,K)*ZE_RHS(I,L,K)
              QE(I,L,IRK+1) = QE(I,L,IRK+1) + ATVD(IRK,K)*QE(I,L,K) +      &
     &                          DT*BTVD(IRK,K)*QE_RHS(I,L,K)
            END DO
          END DO
        END DO
        
!.......Before next stage, update status of solution...        
!.......Update WET/DRY status of solution
        IF (IWET.GT.0) THEN
          CALL WETDRY
        END IF
!.......Apply slope limiter if needed
        IF (ISLP.GT.0) THEN          
          CALL SLOPELIMIT
        END IF
!.......Update NODISP flag is needed
        IF (INONHYDRO.NE.0.AND.IBREAK.NE.0) THEN
          CALL BREAKING
          DO L = 1,NE
            DISPFLG(L) = MIN(DISPFLG(L),WDFLG(L))
          END DO
        END IF
        
      END DO
!---------------------------------------------------------------------------!


!.....Update current solution and zero out the solution at the other stages
      DO I = 1,P+1
        ZE(I,:,1) = ZE(I,:,NRK+1)
        QE(I,:,1) = QE(I,:,NRK+1)
        ZE(I,:,2:NRK+1) = 0.D0
        QE(I,:,2:NRK+1) = 0.D0
        
        PD(I,:,1) = PD(I,:,NRK)
        PB(I,:,1) = PB(I,:,NRK)
        
        DO K = 1,NONHYDRO_EXT
          PDLVL(I,:,K) = PDLVL(I,:,K+1)
          PBLVL(I,:,K) = PBLVL(I,:,K+1)
        END DO
        PDLVL(I,:,NONHYDRO_EXT+1) = PD(I,:,1)
        PBLVL(I,:,NONHYDRO_EXT+1) = PB(I,:,1)
      END DO

      RETURN
      END SUBROUTINE DG_TIMESTEP
!==========================================================================!
!==========================================================================!
!..........................................................................!
      SUBROUTINE DG_AREA_INTEGRALS()

      USE GLOBALS,    ONLY : NE,ZE,QE,PHI,DPHI,DE_IN,DX_IN,G,ZE_RHS,QE_RHS,&
     &                       IRK,WEGP,LE,WDFLG
      USE SIZES,      ONLY : SZ,C12
      USE READ_DGINP, ONLY : P,NEGP

      IMPLICIT NONE
      INTEGER	:: L,K,I
      REAL(SZ)  :: ZE_IN,QE_IN,HE_IN,ZE_FI,QE_FI,ZE_SI,QE_SI
            
!.....Calculate the area integral for each element
      DO L=1,NE
!.......If the element is dry skip area integral
        IF (WDFLG(L).NE.0) THEN
!.........Cycle through the Gauss integration points to build area integral
          DO K = 1,NEGP
            ZE_IN = 0.D0
            QE_IN = 0.D0
!...........Cycle through the basis functions to build internal solution
            DO I = 1,P+1
              ZE_IN = ZE_IN + ZE(I,L,IRK)*PHI(I,K)
              QE_IN = QE_IN + QE(I,L,IRK)*PHI(I,K)
            END DO
            HE_IN = ZE_IN + DE_IN(L,K)
!...........Using the internal solutions build the source and flux area integrals
            ZE_FI = QE_IN ! ZE flux integral
            QE_FI = (QE_IN)**2/HE_IN+C12*G*ZE_IN*(2.D0*DE_IN(L,K)+ZE_IN) ! QE flux integral
          
            ZE_SI = 0.D0               ! ZE source term
            QE_SI = G*ZE_IN*DX_IN(L,K) ! QE source term
!...........Include the flux and source integrals in the RHS terms
            DO I = 1,P+1
              ZE_RHS(I,L,IRK) = ZE_RHS(I,L,IRK) + (ZE_FI*DPHI(I,K) +         &
     &                            ZE_SI*PHI(I,K)*LE(L)*C12)*WEGP(K)
              QE_RHS(I,L,IRK) = QE_RHS(I,L,IRK) + (QE_FI*DPHI(I,K) +         &
     &                            QE_SI*PHI(I,K)*LE(L)*C12)*WEGP(K)
            END DO          
          END DO
        END IF        
      END DO            
      
      RETURN
      END SUBROUTINE DG_AREA_INTEGRALS
!..........................................................................!
!==========================================================================!
!..........................................................................!
      SUBROUTINE DG_INTERNAL_EDGES
      
      USE GLOBALS,    ONLY : NE,ZE,QE,PHIB,DE_ED,G,ZE_RHS,QE_RHS,IRK,GHAT, &
     &                    FHAT,ZE_LT,ZE_RT,QE_LT,QE_RT,HE_LT,HE_RT,FH_LT,  &
     &                    FH_RT,GH_LT,GH_RT,WDFLG,PD,PSIB
      USE SIZES,      ONLY : SZ,C12
      USE READ_DGINP, ONLY : P
      
      IMPLICIT NONE
      
      INTEGER	:: I,L
      REAL(SZ)  :: FHAT_LT,FHAT_RT,GHAT_LT,GHAT_RT
      REAL(SZ)  :: PD_LT,PD_RT
      REAL(SZ)  :: GTMP
      
!.....Cycle through all the internal edges, 2...NE
      DO L=2,NE ! Note that there are NE+1 edges in 1-D
!.......If the element edge is bounded by two dry elements skip calculation
        IF (WDFLG(L-1).NE.0.OR.WDFLG(L).NE.0) THEN
!.........Determine left and right values about the edge
!.........         -----|-----
!.........         LT-> | <- RT
          ZE_LT = 0.D0
          ZE_RT = 0.D0
          QE_LT = 0.D0
          QE_RT = 0.D0
          
          PD_LT = 0.D0
          PD_RT = 0.D0
          DO I = 1,P+1
            ZE_LT = ZE_LT + ZE(I,L-1,IRK)*PHIB(I,2)
            ZE_RT = ZE_RT + ZE(I,L,IRK)*PHIB(I,1)
          
            QE_LT = QE_LT + QE(I,L-1,IRK)*PHIB(I,2)
            QE_RT = QE_RT + QE(I,L,IRK)*PHIB(I,1)
            
!             PD_LT = PD_LT + PD(I,L-1,IRK)*PSIB(I,2)
!             PD_RT = PD_RT + PD(I,L,IRK)*PSIB(I,1)
          END DO
          HE_LT = ZE_LT + DE_ED(L)
          HE_RT = ZE_RT + DE_ED(L)
!.........Using left and right values determine left and right flux terms
          FH_LT = QE_LT ! ZE left flux term
          FH_RT = QE_RT ! ZE right flux term
        
          GH_LT = (QE_LT**2)/HE_LT + C12*G*ZE_LT*(2.D0*DE_ED(L)+ZE_LT)+PD_LT ! QE left flux term
          GH_RT = (QE_RT**2)/HE_RT + C12*G*ZE_RT*(2.D0*DE_ED(L)+ZE_RT)+PD_RT ! QE right flux term
!.........Determine numerical flux based on left and right flux terms
          CALL GET_FLUX
!.........Assign initial flux calculations
          FHAT_LT = FHAT
          FHAT_RT = FHAT
          GHAT_LT = GHAT
          GHAT_RT = GHAT
!.........With wet/dry fronts we must check to see if mass flux is coming from a dry element
          IF (WDFLG(L-1).EQ.0.OR.WDFLG(L).EQ.0) THEN
            IF (ABS(FHAT).GT.0.D0) THEN ! There is some mass flux, check
              IF (WDFLG(L-1).EQ.0) THEN ! Left element is dry, make sure mass flux isn't positive
                IF (FHAT.GT.0.D0) THEN  ! Make edge reflective
                  QE_LT = -QE_RT
                  HE_LT =  HE_RT
                  FH_LT = -FH_RT
                  GH_LT =  GH_RT
                  CALL GET_FLUX()
                  FHAT_LT = FHAT
                  FHAT_RT = FHAT
                  GHAT_LT = GHAT
                  GHAT_RT = GHAT
                END IF
!                 ELSE                    ! Allow mass flux, but zero out gravity terms
                  GTMP  = G
                  G     = 0.D0
                  GH_LT = (QE_LT**2)/HE_LT + C12*G*ZE_LT*(2.D0*DE_ED(L)+ZE_LT)
                  GH_RT = (QE_RT**2)/HE_RT + C12*G*ZE_RT*(2.D0*DE_ED(L)+ZE_RT)
                  CALL GET_FLUX()
                  G       = GTMP
                  GHAT_LT = GHAT
!                 END IF
              ELSE                      ! Right element is dry, make sure mass flux isn't negative
                IF (FHAT.LT.0.D0) THEN  ! Make edge reflective
                  QE_RT = -QE_LT
                  HE_RT =  HE_LT
                  FH_RT = -FH_LT
                  GH_RT =  GH_LT
                  CALL GET_FLUX()
                  FHAT_LT = FHAT
                  FHAT_RT = FHAT
                  GHAT_LT = GHAT
                  GHAT_RT = GHAT
                END IF
!                 ELSE                    ! Allow mass flux, but zero out gravity terms
                  GTMP  = G
                  G     = 0.D0
                  GH_LT = (QE_LT**2)/HE_LT + C12*G*ZE_LT*(2.D0*DE_ED(L)+ZE_LT)
                  GH_RT = (QE_RT**2)/HE_RT + C12*G*ZE_RT*(2.D0*DE_ED(L)+ZE_RT)
                  CALL GET_FLUX()
                  G     = GTMP
                  GHAT_RT = GHAT
!                 END IF
              END IF
            END IF
          END IF
!.........Update RHS vectors
!......... left element contribution
          DO I = 1,P+1
            ZE_RHS(I,L-1,IRK) = ZE_RHS(I,L-1,IRK) - PHIB(I,2)*FHAT_LT
            QE_RHS(I,L-1,IRK) = QE_RHS(I,L-1,IRK) - PHIB(I,2)*GHAT_LT
          END DO
!......... right element contribution        
          DO I = 1,P+1
            ZE_RHS(I,L,IRK) = ZE_RHS(I,L,IRK) + PHIB(I,1)*FHAT_RT
            QE_RHS(I,L,IRK) = QE_RHS(I,L,IRK) + PHIB(I,1)*GHAT_RT
          END DO
        END IF
      END DO
      
      RETURN
      END SUBROUTINE DG_INTERNAL_EDGES
!..........................................................................!
!==========================================================================!
!..........................................................................!     

!
!     SUBROUTINE DG_BOUNDARY_EDGES
!      
!....................................................................
      SUBROUTINE DG_BOUNDARY_EDGES
      
      USE GLOBALS,    ONLY : ZE,QE,PHIB,DE_ED,G,IRK,NE,NN,FHAT,GHAT,       &
     &                    ZE_RHS,QE_RHS,ZE_LT,ZE_RT,QE_LT,QE_RT,HE_LT,     &
     &                    HE_RT,FH_LT,FH_RT,GH_LT,GH_RT,WDFLG,TIME_RK,     &
     &                    SPNG_GEN,X,PD,PSIB
      USE SIZES,      ONLY : SZ,C12
      USE READ_DGINP, ONLY : P,BOUNDTYPE
      
      IMPLICIT NONE
      INTEGER	:: I
      REAL(SZ)  :: PD_RT,PD_LT
      REAL(SZ)  :: ZE_IMPOSED,QE_IMPOSED
      
!----------------------------------------------
!.....Left-hand boundary "edge" at node 1
!----------------------------------------------
!.....If first element is dry ignore calculation
      IF (WDFLG(1).EQ.1) THEN
!.......Build internal edge solution for node 1
        ZE_RT = 0.D0
        QE_RT = 0.D0
        PD_RT = 0.D0
      
        DO I = 1,P+1
          ZE_RT = ZE_RT + ZE(I,1,IRK)*PHIB(I,1)
          QE_RT = QE_RT + QE(I,1,IRK)*PHIB(I,1)
!           PD_RT = PD_RT + PD(I,1,IRK)*PSIB(I,1)
        END DO
!.......Left hand side boundary condition
        IF (SPNG_GEN(1,1).GT.0.D0) THEN
          ZE_LT = ZE_IMPOSED(X(1),TIME_RK)
          QE_LT = QE_IMPOSED(X(1),TIME_RK,ZE_LT)
          PD_LT = PD_RT
        ELSEIF (ADJUSTL(TRIM(BOUNDTYPE)).EQ.'radiant') THEN
          ZE_LT = ZE_RT
          QE_LT = QE_RT
          PD_LT = PD_RT
        ELSE             ! Else reflective boundary
          ZE_LT =  ZE_RT
          QE_LT = -QE_RT
          PD_LT =  PD_RT
        END IF
        
        HE_RT = ZE_RT+DE_ED(1)
        HE_LT = ZE_LT+DE_ED(1)
!.......Calculate left and right flux values
        FH_LT = QE_LT
        FH_RT = QE_RT
        
        GH_LT = (QE_LT**2)/HE_LT + C12*G*ZE_LT*(2.D0*DE_ED(1)+ZE_LT)+PD_LT
        GH_RT = (QE_RT**2)/HE_RT + C12*G*ZE_RT*(2.D0*DE_ED(1)+ZE_RT)+PD_RT
!.......Calculate flux value across boundary
        CALL GET_FLUX
!.......Update RHS vectors
        DO I = 1,P+1
          ZE_RHS(I,1,IRK) = ZE_RHS(I,1,IRK) + PHIB(I,1)*FHAT
          QE_RHS(I,1,IRK) = QE_RHS(I,1,IRK) + PHIB(I,1)*GHAT
        END DO
      END IF
!----------------------------------------------
!.....Right-hand boundary "edge" at node NE
!----------------------------------------------
!.....If last element is dry ignore calculation
      IF (WDFLG(NE).EQ.1) THEN
!.......Build internal edge solution for node NN
        ZE_LT = 0.D0
        QE_LT = 0.D0
        PD_LT = 0.D0
      
        DO I = 1,P+1
          ZE_LT = ZE_LT + ZE(I,NE,IRK)*PHIB(I,2)
          QE_LT = QE_LT + QE(I,NE,IRK)*PHIB(I,2)
!           PD_LT = PD_LT + PD(I,NE,IRK)*PHIB(I,2)
        END DO
!.......Right hand side boundary condition
        IF (SPNG_GEN(NE,1).GT.0.D0) THEN
          ZE_RT = ZE_IMPOSED(X(NE+1),TIME_RK)
          QE_RT = QE_IMPOSED(X(NE+1),TIME_RK,ZE_RT)
          PD_RT = PD_LT
        ELSEIF (ADJUSTL(TRIM(BOUNDTYPE)).EQ.'radiant') THEN
          ZE_RT = ZE_LT
          QE_RT = QE_LT
          PD_RT = PD_LT
        ELSE            ! Else reflective boundary 
          ZE_RT =  ZE_LT
          QE_RT = -QE_LT
          PD_RT =  PD_LT
        END IF
        
        
      
        HE_RT = ZE_RT+DE_ED(NN)
        HE_LT = ZE_LT+DE_ED(NN)
!.......Calculate left and right flux values
        FH_LT = QE_LT
        FH_RT = QE_RT
        
        GH_LT = (QE_LT**2)/HE_LT + C12*G*ZE_LT*(2.D0*DE_ED(1)+ZE_LT)+PD_LT
        GH_RT = (QE_RT**2)/HE_RT + C12*G*ZE_RT*(2.D0*DE_ED(1)+ZE_RT)+PD_RT
!.......Calculate flux value across boundary
        CALL GET_FLUX
!.......Update RHS vectors
        DO I = 1,P+1
          ZE_RHS(I,NE,IRK) = ZE_RHS(I,NE,IRK) - PHIB(I,2)*FHAT
          QE_RHS(I,NE,IRK) = QE_RHS(I,NE,IRK) - PHIB(I,2)*GHAT
        END DO
      END IF
      
      RETURN
      END SUBROUTINE DG_BOUNDARY_EDGES
!..........................................................................!
!==========================================================================!
!..........................................................................!
      SUBROUTINE SWE_RHS_VARIABLES
      
      USE READ_DGINP, ONLY : NEGP,P,INONHYDRO,NONHYDRO_EXT
      USE GLOBALS,    ONLY : PD,PB,WDFLG,LE,PSI,DPSI,PHI,DPHI,WEGP,ZE_RHS, &
     &                       QE_RHS,DE_IN,DX_IN,IRK,NE,G,ZE,QE,MANN,X,XEGP,&
     &                       TIME_RK,SPNG_GEN,SPNG_ABS,DISPFLG
      USE SIZES,      ONLY : SZ,C12
      
      IMPLICIT NONE
      INTEGER  :: L,K,I
      REAL(SZ) :: ZE_FI,ZE_SI,QE_FI,QE_SI
      REAL(SZ) :: QE_IN,ZE_IN,HE_IN,UE_IN
      REAL(SZ) :: PDX_IN,PB_IN,PD_IN
      REAL(SZ) :: ZE_OUT,QE_OUT,ZE_IMPOSED,QE_IMPOSED,XL,XIN
      REAL(SZ) :: PDLOC(P+1,NE),PBLOC(P+1,NE)
      
!.....Determine Nonhydrostatic pressure contribution
      ! Check to see if the pressure solution will be extrapolated.
!       IF (NONHYDRO_EXT.EQ.0.AND.IRK.EQ.1) THEN	! No extrapolation used
        CALL NONHYDRO
        PDLOC = PD(:,:,IRK)
        PBLOC = PB(:,:,IRK)
!       END IF
      
!.....Add extra variables to RHS vectors
      DO L = 1,NE
        IF (WDFLG(L).NE.0) THEN
          DO K = 1,NEGP
            ZE_IN = 0.D0
            QE_IN = 0.D0
            
            PDX_IN = 0.D0
            PB_IN = 0.D0
            PD_IN = 0.D0
            DO I = 1,P+1
              ZE_IN = ZE_IN + ZE(I,L,IRK)*PHI(I,K)
              QE_IN = QE_IN + QE(I,L,IRK)*PHI(I,K)
              
!               PDX_IN = PDX_IN + PD(I,L,IRK)*DPSI(I,K)/(C12*LE(L))
!               PB_IN  = PB_IN  + PB(I,L,IRK)*PSI(I,K)
!               PD_IN  = PD_IN  + PD(I,L,IRK)*PSI(I,K)
              
              PDX_IN = PDX_IN + PDLOC(I,L)*DPSI(I,K)/(C12*LE(L))
              PB_IN  = PB_IN  + PBLOC(I,L)*PSI(I,K)
              PD_IN  = PD_IN  + PDLOC(I,L)*PSI(I,K)
            END DO
            HE_IN = ZE_IN + DE_IN(L,K)
            UE_IN = QE_IN/HE_IN
            
            ZE_FI = 0.D0
            ZE_SI = 0.D0
            QE_FI = 0.D0
            QE_SI = 0.D0
            
            ! Pressure Terms
            IF (INONHYDRO.NE.0.AND.DISPFLG(L).NE.0) THEN
!               QE_FI = QE_FI + PD_IN
              QE_SI = QE_SI + DX_IN(L,K)*PB_IN - PDX_IN
            END IF
            ! Bottom Friction (Mannings N)
            IF (MANN(L,K).GT.0.D0) THEN
              QE_SI = QE_SI - G*MANN(L,K)**2*UE_IN*ABS(UE_IN)/(HE_IN**(1.D0/3.D0))
            END IF
            ! Sponge Layer Absorbing
            IF (SPNG_ABS(L,K).GT.0.D0) THEN
              ZE_SI = ZE_SI - SPNG_ABS(L,K)*ZE_IN
              QE_SI = QE_SI - SPNG_ABS(L,K)*QE_IN
            END IF
            ! Sponge Layer Generation
            IF (SPNG_GEN(L,K).GT.0.D0) THEN              
              XL    = X(L)
              XIN   = XL + C12*LE(L)*(XEGP(K)+1.D0)
              ZE_OUT = ZE_IMPOSED(XIN,TIME_RK)
              QE_OUT = QE_IMPOSED(XIN,TIME_RK,ZE_OUT)
              ZE_SI = ZE_SI + SPNG_GEN(L,K)*(ZE_OUT-ZE_IN)
              QE_SI = QE_SI + SPNG_GEN(L,K)*(QE_OUT-QE_IN)
            END IF
                  
            DO I = 1,P+1
              ZE_RHS(I,L,IRK) = ZE_RHS(I,L,IRK) + (ZE_FI*DPHI(I,K) +       &
     &                            ZE_SI*PHI(I,K)*LE(L)*C12)*WEGP(K)
              QE_RHS(I,L,IRK) = QE_RHS(I,L,IRK) + (QE_FI*DPHI(I,K) +       &
     &                            QE_SI*PHI(I,K)*LE(L)*C12)*WEGP(K)
            END DO
            
          END DO
        END IF
      END DO
      
      END SUBROUTINE SWE_RHS_VARIABLES
!..........................................................................!
!==========================================================================!
!..........................................................................!
      SUBROUTINE GET_FLUX
      
      USE GLOBALS, ONLY : G,QE_LT,ZE_LT,HE_LT,QE_RT,ZE_RT,HE_RT,FH_LT,     &
     &                    FH_RT,GH_LT,GH_RT,FHAT,GHAT
      USE SIZES,   ONLY : C12,SZ
      
      IMPLICIT NONE
      REAL(SZ)	:: VL,VR
      REAL(SZ)	:: ALPHAP,ALPHAN,ALPHA
      
!.....Calculate ALPHA
      VL = QE_LT/HE_LT
      VR = QE_RT/HE_RT            
      
      ALPHAP = MAX(ABS(VR+DSQRT(G*HE_RT)),ABS(VL+DSQRT(G*HE_LT)))
      ALPHAN = MAX(ABS(VR-DSQRT(G*HE_RT)),ABS(VL-DSQRT(G*HE_LT)))
      ALPHA  = MAX(ALPHAP,ALPHAN)
!.....Calculate LLF Flux based on left and right flux values and alpha      
      FHAT = C12*(FH_LT+FH_RT) - C12*ALPHA*(ZE_RT-ZE_LT)
      GHAT = C12*(GH_LT+GH_RT) - C12*ALPHA*(QE_RT-QE_LT)
      
      RETURN
      END SUBROUTINE GET_FLUX
!..........................................................................!
!==========================================================================!
!..........................................................................!
      SUBROUTINE COMPLETE_RHS
      
      USE GLOBALS,    ONLY : NE,MDG,DGPIV,ZE_RHS,QE_RHS,IRK,LE
      USE SIZES,      ONLY : C12
      USE READ_DGINP, ONLY : P
      
      IMPLICIT NONE
      INTEGER	:: I,L,INFO
      
!.....Check to see if DG method is "Finite Volume" (p=0) or not
      IF (P.GT.0) THEN
!.......Not finite volume, use LAPACK to solve inverse problem
        CALL DGETRS('No transpose',P+1,NE,MDG,P+1,DGPIV,ZE_RHS(:,:,IRK),P+1,INFO)
        CALL DGETRS('No transpose',P+1,NE,MDG,P+1,DGPIV,QE_RHS(:,:,IRK),P+1,INFO)
!.......Divide each elemental component of the solution by half of element length
        DO I = 1,P+1
          DO L = 1,NE
            ZE_RHS(I,L,IRK) = ZE_RHS(I,L,IRK)/(C12*LE(L))
            QE_RHS(I,L,IRK) = QE_RHS(I,L,IRK)/(C12*LE(L))
          END DO
        END DO
      ELSE
!......."Finite Volume" like (p=0)
!.......Much simpler to solve inverse problem, simply divide by MDG(1,1)
        DO L = 1,NE
          ZE_RHS(1,L,IRK) = ZE_RHS(1,L,IRK)/(C12*LE(L)*MDG(1,1))
          QE_RHS(1,L,IRK) = QE_RHS(1,L,IRK)/(C12*LE(L)*MDG(1,1))
        END DO
      END IF
      
      RETURN
      END SUBROUTINE COMPLETE_RHS
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!
      SUBROUTINE WETDRY
      
      USE GLOBALS,    ONLY : NE,ZE,QE,PHI,PHIB,DE_IN,DE_ED,PHIB,IRK,WDFLG, &
     &                       WEGP
      USE READ_DGINP, ONLY : NEGP,P,H0,DGBASIS
      USE SIZES,      ONLY : SZ,C12,ZERO
      
      IMPLICIT NONE
      INTEGER :: I,L,K
      INTEGER :: ELEM_CHECK,NDRYNODE
      
      INTEGER  :: HE_MAX_LOC
      REAL(SZ) :: HE_MAX,HE_AVG,H00
      REAL(SZ) :: ZE_IN(NEGP+2),QE_IN(NEGP+2)
      REAL(SZ) :: ZE_BND(2),QE_BND(2),DE_AVG,QE_SUM
      REAL(SZ) :: HE_IN(2),HHAT(2)
            
      H00 = H0
!.....Cycle through all the elements to check wet/dry status of each
      DO L = 1,NE
        ELEM_CHECK = 1
!.......Determine the value of the solution at each quadrature
!.......point and endpoint.
        NDRYNODE = 0      ! The initial number of dry nodes is zero
        HE_MAX = ZERO     ! Initialize max water level in element
        HE_MAX_LOC = 1
        ! Cycle through Gauss integration points
        DO K = 1,NEGP     
          ZE_IN(K) = 0.D0
          QE_IN(K) = 0.D0
          DO I = 1,P+1
            ZE_IN(K) = ZE_IN(K) + PHI(I,K)*ZE(I,L,IRK+1)
            QE_IN(K) = QE_IN(K) + PHI(I,K)*QE(I,L,IRK+1)
          END DO          
          IF ( (ZE_IN(K)+DE_IN(L,K)).LE.H00 ) THEN
            NDRYNODE = NDRYNODE+1
          END IF
          IF ( (ZE_IN(K)+DE_IN(L,K)).GT.HE_MAX ) THEN
            HE_MAX = (ZE_IN(K)+DE_IN(L,K))
            HE_MAX_LOC = K
          END IF
        END DO        
        ! Cycle through element edges
        DO K = 1,2
          ZE_IN(NEGP+K) = 0.D0
          QE_IN(NEGP+K) = 0.D0
          DO I = 1,P+1
            ZE_IN(NEGP+K) = ZE_IN(NEGP+K) + PHIB(I,K)*ZE(I,L,IRK+1)
            QE_IN(NEGP+K) = QE_IN(NEGP+K) + PHIB(I,K)*QE(I,L,IRK+1)
          END DO
          IF ( (ZE_IN(NEGP+K)+DE_ED(L+K-1)).LE.(H00) ) THEN
            NDRYNODE = NDRYNODE+1
          END IF
          IF ( (ZE_IN(NEGP+K)+DE_ED(L+K-1)).GT.HE_MAX ) THEN
            HE_MAX = (ZE_IN(NEGP+K)+DE_ED(L+K-1))
            HE_MAX_LOC = NEGP+K
          END IF
        END DO        
        ! Setup matrix of elemental boundary values for all elements
        DO I = 1,2
          ZE_BND(I) = ZE_IN(NEGP+I)
          QE_BND(I) = QE_IN(NEGP+I)
        END DO        
!         DE_AVG = C12*(DE_ED(L+1)+DE_ED(L))     

!.......If the element is dry at all nodes designate the element
!.......as dry.  Set the water depth parallel to the bathymetry
!.......at the minimum H0 value and zero out the velocity.
        IF (NDRYNODE.EQ.NEGP+2) THEN
          WDFLG(L) = 0
          ZE_BND(1) = H0-DE_ED(L)
          ZE_BND(2) = H0-DE_ED(L+1)
          IF (ADJUSTL(TRIM(DGBASIS)).EQ.'nodal') THEN
            DO I = 1,P+1
              ZE(I,L,IRK+1) = ZE_BND(1) + (ZE_BND(2)-ZE_BND(1))/REAL(P)*REAL(I-1)
              QE(I,L,IRK+1) = QE_BND(1) + (QE_BND(2)-QE_BND(1))/REAL(P)*REAL(I-1)
            END DO
          ELSE
            ZE(1,L,IRK+1) = C12*(ZE_BND(2)+ZE_BND(1))
            ZE(2,L,IRK+1) = C12*(ZE_BND(2)-ZE_BND(1))
            DO I = 3,P+1
              ZE(I,L,IRK+1) = 0.D0
            END DO
          END IF
          
          QE(:,L,IRK+1) = 0.D0                    
          
          ELEM_CHECK = 0
          GOTO 9100
!.......Otherwise if the element is fully wet designate the element
!.......as wet and check element.
        ELSEIF (NDRYNODE.EQ.0) THEN
          IF (WDFLG(L).EQ.1) THEN
            ELEM_CHECK = 0
          ELSE
            ELEM_CHECK = 1
          END IF
          GOTO 9200
!.......Otherwise Redistribute mass for partially wet elements
        ELSE
          HE_AVG = 0.D0 ! Keep track of the average water depth over element
          DO K = 1,NEGP
            HE_AVG = HE_AVG + WEGP(K)*(ZE_IN(K)+DE_IN(L,K))
          END DO
          HE_AVG = C12*HE_AVG
          ! If the average water depth over element is less than H0 element is dry
          IF (HE_AVG.LE.(H00)) THEN
            WDFLG(L) = 0
            HE_AVG = H0
            ZE_BND(1) = HE_AVG-DE_ED(L)
            ZE_BND(2) = HE_AVG-DE_ED(L+1)
            IF (ADJUSTL(TRIM(DGBASIS)).EQ.'nodal') THEN
              DO I = 1,P+1
                ZE(I,L,IRK+1) = ZE_BND(1) + (ZE_BND(2)-ZE_BND(1))/REAL(P)*REAL(I-1)
                QE(I,L,IRK+1) = QE_BND(1) + (QE_BND(2)-QE_BND(1))/REAL(P)*REAL(I-1)
              END DO
            ELSE
              ZE(1,L,IRK+1) = C12*(ZE_BND(2)+ZE_BND(1))
              ZE(2,L,IRK+1) = C12*(ZE_BND(2)-ZE_BND(1))
              DO I = 3,P+1
                ZE(I,L,IRK+1) = 0.D0
              END DO
            END IF
            QE(:,L,IRK+1) = 0.D0
            
            ELEM_CHECK = 0
            GOTO 9100
          ! Otherwise redistribute water depth over element
          ELSE
            HE_IN(1) = ZE_BND(1)+DE_ED(L)
            HE_IN(2) = ZE_BND(2)+DE_ED(L+1)
            ! If HE_IN(1) > HE_IN(2) then left side of element is higher
            IF (HE_IN(1).GT.HE_IN(2)) THEN
              HHAT(2) = MAX(H0,HE_IN(2))
              HHAT(1) = MAX(H0,HE_IN(1)-(HHAT(2)-HE_IN(2)))
            ! Otherwise the right side of the element is higher            
            ELSE
              HHAT(1) = MAX(H0,HE_IN(1))
              HHAT(2) = MAX(H0,HE_IN(2)-(HHAT(1)-HE_IN(1)))
            END IF
            HE_MAX = ZERO ! Redefine max water depth based on new values (for elem_check
            HE_MAX_LOC = 1
            DO I = 1,2
              IF (HHAT(I).GT.HE_MAX) THEN
                HE_MAX = HHAT(I)
                HE_MAX_LOC = I
              END IF
            END DO
            ! Determine ZE_BND based on new values
            ZE_BND(1) = HHAT(1)-DE_ED(L)
            ZE_BND(2) = HHAT(2)-DE_ED(L+1)
            
            ! Reevaluate the wet/dry status of the element based on redistibuted values
            NDRYNODE = 0
            QE_SUM   = 0.D0
            DO I = 1,2
              IF (HHAT(I).LE.H00) THEN
                NDRYNODE = NDRYNODE+1
                QE_SUM   = QE_SUM + QE_BND(I)               
              END IF
            END DO
            ! Element is fully dry
            IF (NDRYNODE.EQ.2) THEN
              QE(:,L,IRK+1) = 0.D0
            ! Element could possibly be wet, redistribute the momentum
            ELSE
              QE_SUM = QE_SUM/(2.D0-REAL(NDRYNODE))
              DO I = 1,2
                IF (HHAT(I).LE.H00) THEN ! This element side is "dry"
                  QE_BND(I) = 0.D0
                ELSE               ! This element side is not "dry"
                  IF (QE_SUM*QE_BND(I).LT.0.D0) THEN
                    QE_BND(I) = QE_BND(I)+QE_SUM
                  ELSE
                    QE_BND(I) = QE_BND(I)                  
                  END IF
                END IF
              END DO            
            END IF
            
            ! Using the redistributed values reconstruct ze and qe
            IF (ADJUSTL(TRIM(DGBASIS)).EQ.'nodal') THEN
              DO I = 1,P+1
                ZE(I,L,IRK+1) = ZE_BND(1) + (ZE_BND(2)-ZE_BND(1))/REAL(P)*REAL(I-1)
                QE(I,L,IRK+1) = QE_BND(1) + (QE_BND(2)-QE_BND(1))/REAL(P)*REAL(I-1)
              END DO
            ELSE
              ZE(1,L,IRK+1) = C12*(ZE_BND(2)+ZE_BND(1))
              ZE(2,L,IRK+1) = C12*(ZE_BND(2)-ZE_BND(1))
              QE(1,L,IRK+1) = C12*(QE_BND(2)+QE_BND(1))
              QE(2,L,IRK+1) = C12*(QE_BND(2)-QE_BND(1))
              DO I = 3,P+1
                ZE(I,L,IRK+1) = 0.D0
                QE(I,L,IRK+1) = 0.D0
              END DO
            END IF
            ZE_IN(1) = ZE_BND(1)
            ZE_IN(2) = ZE_BND(2)
            ELEM_CHECK = 1
            GOTO 9200
          END IF          
        END IF
!.......For certain elements check for change in wet dry status.
9200    IF (ELEM_CHECK.EQ.1) THEN
          IF (ZE_IN(HE_MAX_LOC).GT.(-MIN(DE_ED(L),DE_ED(L+1))+H00)) THEN
            WDFLG(L) = 1
          ELSE
            WDFLG(L) = 0
            QE(:,L,IRK+1) = 0.D0
          END IF
        END IF
9100    CONTINUE

      END DO      
      
      RETURN
      END SUBROUTINE WETDRY
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!
      SUBROUTINE SLOPELIMIT
      
      USE READ_DGINP, ONLY : P,NEGP,SLOPEM,DGBASIS,H0
      USE GLOBALS,    ONLY : ZE,QE,NE,WEGP,IRK,PHIB,PHI,LE,WDFLG,          &
     &                       DE_ED,G
      USE SIZES,      ONLY : SZ,C12
      
      IMPLICIT NONE
      
      INTEGER  :: L,I,K,ind
      REAL(SZ) :: ZESLP(NE),ZEAVG(NE)
      REAL(SZ) :: QESLP(NE),QEAVG(NE)
      REAL(SZ) :: A(3),SLP
      
      
!.....Cycle through elements to assess cell averages
      DO L = 1,NE
        ZEAVG(L) = 0.D0
        ZESLP(L) = 0.D0
        
        QEAVG(L) = 0.D0
        QESLP(L) = 0.D0
        DO I = 1,P+1
          DO K = 1,NEGP
            ZEAVG(L) = ZEAVG(L) + ZE(I,L,IRK+1)*PHI(I,K)*WEGP(K)
            QEAVG(L) = QEAVG(L) + QE(I,L,IRK+1)*PHI(I,K)*WEGP(K)
          END DO
          ZESLP(L) = ZESLP(L) + ZE(I,L,IRK+1)*PHIB(I,2)
          QESLP(L) = QESLP(L) + QE(I,L,IRK+1)*PHIB(I,2)
        END DO
        ZEAVG(L) = C12*ZEAVG(L)
        ZESLP(L) = ZESLP(L)-ZEAVG(L)
        
        QEAVG(L) = C12*QEAVG(L)
        QESLP(L) = QESLP(L)-QEAVG(L)
      END DO
!.....Apply slopelimiter if necessary (TVB)
      DO L = 2,NE-1
        IND = 999
        A(:) = 999.d0
        SLP = 999.d0
        IF (ABS(ZESLP(L)).GT.SLOPEM*4.d0) THEN  ! Slope limit
          A(1) = ZESLP(L);
          A(2) = ZEAVG(L+1)-ZEAVG(L)
          A(3) = ZEAVG(L)-ZEAVG(L-1)
          CALL MINMOD(A,IND,SLP)
          IF (IND.NE.1) THEN ! Adjust solution
            ! Solution is zeavg(l) + A(ind)*x
            IF (ADJUSTL(TRIM(DGBASIS)).EQ.'nodal') THEN
              DO I = 1,P+1
                ZE(I,L,IRK+1) = ZEAVG(L) + SLP*(-1.D0+2.D0*REAL(I-1)/REAL(P))
              END DO
            ELSE
              ZE(1,L,IRK+1) = ZEAVG(L)
              ZE(2,L,IRK+1) = SLP
              DO I = 3,P+1
                ZE(I,L,IRK+1) = 0.D0
              END DO
            END IF
          END IF
        END IF
        
        IND = 999
        A(:) = 999.d0
        SLP = 999.d0
        IF (ABS(QESLP(L)).GT.SLOPEM*4.d0*10.d0) THEN  ! Slope limit
          A(1) = QESLP(L);
          A(2) = QEAVG(L+1)-QEAVG(L)
          A(3) = QEAVG(L)-QEAVG(L-1)
          CALL MINMOD(A,IND,SLP)
          IF (IND.NE.1) THEN ! Adjust solution
            ! Solution is qeavg(l) + A(ind)*x
            IF (ADJUSTL(TRIM(DGBASIS)).EQ.'nodal') THEN
              DO I = 1,P+1
                QE(I,L,IRK+1) = QEAVG(L) + SLP*(-1.D0+2.D0*REAL(I-1)/REAL(P))
              END DO
            ELSE
              QE(1,L,IRK+1) = QEAVG(L)
              QE(2,L,IRK+1) = SLP
              DO I = 3,P+1
                QE(I,L,IRK+1) = 0.D0
              END DO
            END IF
          END IF
        END IF
      END DO
      
      RETURN
      END SUBROUTINE SLOPELIMIT
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!
      SUBROUTINE MINMOD(A,IND,SLP)
      
      USE SIZES, ONLY : SZ
      
      IMPLICIT NONE
      REAL(SZ),INTENT(IN)  :: A(3)
      INTEGER,INTENT(OUT)  :: IND
      REAL(SZ),INTENT(OUT) :: SLP
      INTEGER              :: I
      REAL(SZ)             :: SGN
      
      IF (ABS(A(1)).LT.1D-12) THEN
        IND = 1
        GOTO 1101
      ELSEIF (A(1).GT.0.D0) THEN
        SGN = 1.D0
      ELSE
        SGN = -1.D0
      END IF
      
      IND = 999
      IF ( (SGN*A(2).GT.0.D0).AND.(SGN*A(3).GT.0.D0) ) THEN
        IF ( (ABS(A(1)).LE.ABS(A(2))).AND.(ABS(A(1)).LE.ABS(A(3))) ) THEN
          IND = 1
        ELSE IF ( (ABS(A(2)).LE.ABS(A(1))).AND.(ABS(A(2)).LE.ABS(A(3))) ) THEN
          IND = 2
        ELSE IF ( (ABS(A(3)).LE.ABS(A(1))).AND.(ABS(A(3)).LE.ABS(A(2))) ) THEN
          IND = 3
        END IF
        SLP = A(IND)
      ELSE
        IND = 0
        SLP = 0.D0
      END IF
      
      
1101  RETURN
      END SUBROUTINE MINMOD
      
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!
      FUNCTION ZE_IMPOSED(XIN,TIN) RESULT(YOUT)
      
      USE GLOBALS, ONLY : NUM_FREQ,SPNG_ZAMP,SPNG_QAMP,SPNG_K,SPNG_SIG,    &
     &                    SPONGE_TYPE,SPNG_DIMP,G
      USE SIZES,   ONLY : SZ
      
      IMPLICIT NONE
      INTEGER  :: FREQ
      REAL(SZ) :: XIN,TIN,YOUT
      
      SELECT CASE (TRIM(SPONGE_TYPE))
        CASE ('SINUSOIDAL')
          YOUT = 0.D0
          DO FREQ = 1,NUM_FREQ
            YOUT = YOUT + SPNG_ZAMP(FREQ)*DCOS(SPNG_K(FREQ)*XIN            &
     &                                             -SPNG_SIG(FREQ)*TIN)
          END DO
        CASE ('CNOIDAL')
          YOUT = SPNG_ZAMP(1)          
          DO FREQ = 2,NUM_FREQ
            YOUT = YOUT + SPNG_ZAMP(FREQ)*DCOS(SPNG_K(FREQ)*(              &
     &               XIN/SPNG_DIMP-SPNG_SIG(FREQ)*TIN*DSQRT(G/SPNG_DIMP) ))
          END DO
          YOUT = SPNG_DIMP*YOUT
        CASE DEFAULT
          YOUT = 0.D0
      END SELECT
      
      RETURN
      END FUNCTION ZE_IMPOSED
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!
      FUNCTION QE_IMPOSED(XIN,TIN,ZIN) RESULT(YOUT)
      
      USE GLOBALS, ONLY : NUM_FREQ,SPNG_QAMP,SPNG_QAMP,SPNG_K,SPNG_SIG,    &
     &                    SPONGE_TYPE,SPNG_DIMP,G
      USE SIZES, ONLY : SZ
      
      IMPLICIT NONE
      INTEGER  :: FREQ
      REAL(SZ) :: XIN,TIN,ZIN,YOUT
      
      SELECT CASE (TRIM(SPONGE_TYPE))
        CASE ('SINUSOIDAL')
          YOUT = 0.D0
          DO FREQ = 1,NUM_FREQ
            YOUT = YOUT + SPNG_QAMP(FREQ)*DCOS(SPNG_K(FREQ)*XIN            &
     &                                             -SPNG_SIG(FREQ)*TIN)
          END DO
          YOUT = YOUT*(SPNG_DIMP+ZIN)
        CASE ('CNOIDAL')
          YOUT = SPNG_QAMP(1)*(1.D0+ZIN/SPNG_DIMP)
          DO FREQ = 2,NUM_FREQ
            YOUT = YOUT + SPNG_QAMP(FREQ)/DCOSH(SPNG_K(FREQ))*             &
     &             (DSINH(SPNG_K(FREQ))+DSINH(SPNG_K(FREQ)*ZIN/SPNG_DIMP)) &
     &             *DCOS(SPNG_K(FREQ)*(                                    &
     &               XIN/SPNG_DIMP-SPNG_SIG(FREQ)*TIN*DSQRT(G/SPNG_DIMP) ))
          END DO
          YOUT = YOUT*DSQRT(G*SPNG_DIMP)*SPNG_DIMP
        CASE DEFAULT
          YOUT = 0.D0
      END SELECT
      
      RETURN
      END FUNCTION QE_IMPOSED

!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!
      SUBROUTINE BREAKING
     
      USE READ_DGINP, ONLY : BREAKMODEL

      IMPLICIT NONE

      SELECT CASE(BREAKMODEL)
        CASE("duran")
          CALL BREAKING_DURAN
        CASE("tonelli")
          CALL BREAKING_TONELLI
        CASE("duran_adj")
          CALL BREAKING_DURAN_ADJ
        CASE DEFAULT
          CALL BREAKING_DURAN
      END SELECT

      RETURN
      END SUBROUTINE BREAKING
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!
      SUBROUTINE BREAKING_DURAN
      
      USE GLOBALS,    ONLY : WDFLG,NE,DISPFLG,ZE,QE,IRK,PHIB,DE_ED,LE
      USE SIZES,      ONLY : SZ,C12
      USE READ_DGINP, ONLY : P
      
      IMPLICIT NONE
      INTEGER     :: L,I
      REAL(SZ)    :: HMAX,TAUJ
      REAL(SZ)    :: QE_IN
      INTEGER     :: MINWIN,MAXWIN
      REAL(SZ)    :: ZE_LT,ZE_RT
      
      DO L = 2,NE-1
        IF (WDFLG(L).EQ.0) THEN
          dispflg(L) = 0
        ELSE
          TAUJ = 0.D0
          HMAX = 0.D0
          ! LEFT NODE
          ZE_LT = 0.D0
          ZE_RT = 0.D0
          QE_IN = 0.D0
          DO I = 1,P+1
            ZE_LT = ZE_LT + ZE(I,L-1,IRK)*PHIB(I,2)
            ZE_RT = ZE_RT + ZE(I,L,IRK)*PHIB(I,1)
            QE_IN = QE_IN + QE(I,L,IRK)*PHIB(I,1)
          END DO
          HMAX = MAX(HMAX,(DE_ED(L)+ZE_RT))
          IF (QE_IN.GE.0.D0) THEN
            TAUJ = TAUJ + (ZE_RT-ZE_LT)
          END IF
          ! RIGHT NODE
          ZE_LT = 0.D0
          ZE_RT = 0.D0
          QE_IN = 0.D0
          DO I = 1,P+1
            ZE_LT = ZE_LT + ZE(I,L,IRK)*PHIB(I,2)
            ZE_RT = ZE_RT + ZE(I,L+1,IRK)*PHIB(I,1)
            QE_IN = QE_IN + QE(I,L,IRK)*PHIB(I,2)
          END DO
          HMAX = MAX(HMAX,(DE_ED(L+1)+ZE_LT))
          IF (QE_IN.LE.0.D0) THEN
            TAUJ = TAUJ + (ZE_LT-ZE_RT)
          END IF
          
          TAUJ = TAUJ/((LE(L)**(P+1))**(C12)*ABS(HMAX))
          IF (TAUJ.GT.1.D0) THEN ! Troubled Cell
            MINWIN = MAX(L-5,1)
            MAXWIN = MIN(L+5,NE)
            DISPFLG(MINWIN:MAXWIN) = 0
          END IF
          
        END IF
      END DO
      
      DO L = 2,NE-1
        IF (DISPFLG(L).NE.0) THEN
          IF (DISPFLG(L-1).EQ.0.AND.DISPFLG(L+1).EQ.0) THEN
            DISPFLG(L) = 0
          END IF
        END IF
      END DO
      
      RETURN
      END SUBROUTINE BREAKING_DURAN
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!
      SUBROUTINE BREAKING_DURAN_ADJ
      
      USE GLOBALS,    ONLY : WDFLG,NE,DISPFLG,ZE,QE,IRK,PHIB,LE
      USE SIZES,      ONLY : SZ,C12
      USE READ_DGINP, ONLY : P
      
      IMPLICIT NONE
      INTEGER     :: L,I
      REAL(SZ)    :: ALP,ALPMAX,QIN
      INTEGER     :: MINWIN,MAXWIN
      REAL(SZ)    :: ZE_LT,ZE_RT,QE_LT,QE_RT
      
      ALPMAX = -0.577350269189626d0 ! 30 degree slope
      DO L = 1,NE
        IF (WDFLG(L).EQ.0) THEN
          dispflg(L) = 0
        ELSE
          ! LEFT NODE
          ZE_LT = 0.D0
          QE_LT = 0.D0
          DO I = 1,P+1
            ZE_LT = ZE_LT + ZE(I,L,IRK)*PHIB(I,1)
            QE_LT = QE_lt + QE(I,L,IRK)*PHIB(I,1)
          END DO

          ! RIGHT NODE
          ZE_RT = 0.D0
          QE_RT = 0.D0
          DO I = 1,P+1
            ZE_RT = ZE_RT + ZE(I,L,IRK)*PHIB(I,2)
            QE_RT = QE_RT + QE(I,L,IRK)*PHIB(I,2)
          END DO
                    
          IF (QE_LT*QE_RT.GE.0.D0) THEN
            QIN = QE_LT/ABS(QE_LT)
          ELSE
            QIN = 0.D0
          END IF 
          
          ALP = (ZE_RT-ZE_LT)/LE(L) * QIN

          IF (ALP.LT.ALPMAX) THEN ! Troubled Cell
            MINWIN = MAX(L-5,1)
            MAXWIN = MIN(L+5,NE)
            DISPFLG(MINWIN:MAXWIN) = 0
          END IF
          
        END IF
      END DO
      
      DO L = 2,NE-1
        IF (DISPFLG(L).NE.0) THEN
          IF (DISPFLG(L-1).EQ.0.AND.DISPFLG(L+1).EQ.0) THEN
            DISPFLG(L) = 0
          END IF
        END IF
      END DO
      
      RETURN
      END SUBROUTINE BREAKING_DURAN_ADJ
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!
      SUBROUTINE BREAKING_TONELLI
      
      USE GLOBALS,    ONLY : WDFLG,NE,DISPFLG,ZE,QE,IRK,PHIB,LE,DE_IN,     &
     &                       DE_ED,PHI
      USE SIZES,      ONLY : SZ,C12
      USE READ_DGINP, ONLY : P,NEGP
      
      IMPLICIT NONE
      INTEGER     :: L,I,K,BRKFLG
      REAL(SZ)    :: HBRK
      INTEGER     :: MINWIN,MAXWIN
      REAL(SZ)    :: ZE_IN

      DISPFLG(:) = 1
      HBRK = 0.8D0
      DO L = 1,NE        
        IF (WDFLG(L).EQ.0) THEN
          DISPFLG(L) = 0
        ELSE
          BRKFLG = 0
          ! LEFT NODE
          ZE_IN = 0.D0
          DO I = 1,P+1
            ZE_IN = ZE_IN + ZE(I,L,IRK)*PHIB(I,1)
          END DO
          IF (ZE_IN.GT.(HBRK*DE_ED(L))) THEN
            BRKFLG = 1
            GOTO 5111
          END IF

          ! INTERNAL NODES
          DO K = 1,NEGP
            ZE_IN = 0.D0
            DO I = 1,P+1
              ZE_IN = ZE_IN + ZE(I,L,IRK)*PHI(I,K)
            END DO
            IF (ZE_IN.GT.(HBRK*DE_IN(L,K))) THEN
              BRKFLG = 1
              GOTO 5111
            END IF
          END DO          
          
          ! RIGHT NODE
          ZE_IN = 0.D0
          DO I = 1,P+1
            ZE_IN = ZE_IN + ZE(I,L,IRK)*PHIB(I,2)
          END DO
          IF (ZE_IN.GT.(HBRK*DE_ED(L+1))) THEN
            BRKFLG = 1
            GOTO 5111
          END IF
          
5111      CONTINUE
          IF (BRKFLG.EQ.1) THEN ! Troubled Cell
            MINWIN = MAX(L-5,1)
            MAXWIN = MIN(L+5,NE)
            DISPFLG(MINWIN:MAXWIN) = 0
          END IF
          
        END IF
      END DO
      
      DO L = 2,NE-1
        IF (DISPFLG(L).NE.0) THEN
          IF (DISPFLG(L-1).EQ.0.AND.DISPFLG(L+1).EQ.0) THEN
            DISPFLG(L) = 0
          END IF
        END IF
      END DO
      
      RETURN
      END SUBROUTINE BREAKING_TONELLI
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!

! To do for DG_Timestep
