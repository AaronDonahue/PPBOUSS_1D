!==========================================================================!
!==========================================================================!
      SUBROUTINE DG_TIMESTEP
      
      USE GLOBALS,    ONLY : IRK,ZE_RHS,QE_RHS,ZE,QE,ATVD,BTVD,DT,NE
      USE SIZES,      ONLY : SZ
      USE READ_DGINP, ONLY : NRK,P,IWET

      IMPLICIT NONE
      INTEGER	:: L,I,K

!---------------------------------------------------------------------------!
!.....Timestep through each of the Runge-Kutta stages
      DO IRK = 1,NRK
!.......Zero out the Right Hand Side vectors for ZE and QE for this stage
        ZE_RHS(:,:,IRK) = 0.D0
        QE_RHS(:,:,IRK) = 0.D0

!.......Build RHS vectors for DG solution
        CALL DG_AREA_INTEGRALS() ! Calculate area integral component

        CALL DG_INTERNAL_EDGES() ! Calculate internal edge integral (fluxes)
  
        CALL DG_BOUNDARY_EDGES() ! Calculate boundary integrals of domain (fluxes)
     
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
        
!.......Update WET/DRY status of solution
        IF (IWET.GT.0) THEN
          CALL WETDRY
        END IF
        
      END DO
!---------------------------------------------------------------------------!


!.....Update current solution and zero out the solution at the other stages
      DO I = 1,P+1
        ZE(I,:,1) = ZE(I,:,NRK+1)
        QE(I,:,1) = QE(I,:,NRK+1)
        ZE(I,:,2:NRK+1) = 0.D0
        QE(I,:,2:NRK+1) = 0.D0                
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
     &                    FH_RT,GH_LT,GH_RT,WDFLG
      USE SIZES,      ONLY : SZ,C12
      USE READ_DGINP, ONLY : P
      
      IMPLICIT NONE
      
      INTEGER	:: I,L
      REAL(SZ)  :: FHAT_LT,FHAT_RT,GHAT_LT,GHAT_RT
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
          DO I = 1,P+1
            ZE_LT = ZE_LT + ZE(I,L-1,IRK)*PHIB(I,2)
            ZE_RT = ZE_RT + ZE(I,L,IRK)*PHIB(I,1)
          
            QE_LT = QE_LT + QE(I,L-1,IRK)*PHIB(I,2)
            QE_RT = QE_RT + QE(I,L,IRK)*PHIB(I,1)
          END DO
          HE_LT = ZE_LT + DE_ED(L)
          HE_RT = ZE_RT + DE_ED(L)
!.........Using left and right values determine left and right flux terms
          FH_LT = QE_LT ! ZE left flux term
          FH_RT = QE_RT ! ZE right flux term
        
          GH_LT = (QE_LT**2)/HE_LT + C12*G*ZE_LT*(2.D0*DE_ED(L)+ZE_LT) ! QE left flux term
          GH_RT = (QE_RT**2)/HE_RT + C12*G*ZE_RT*(2.D0*DE_ED(L)+ZE_RT) ! QE right flux term
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
     &                    HE_RT,FH_LT,FH_RT,GH_LT,GH_RT,WDFLG
      USE SIZES,      ONLY : SZ,C12
      USE READ_DGINP, ONLY : P,BOUNDTYPE
      
      IMPLICIT NONE
      INTEGER	:: I
      
!----------------------------------------------
!.....Left-hand boundary "edge" at node 1
!----------------------------------------------
!.....If first element is dry ignore calculation
      IF (WDFLG(1).EQ.1) THEN
!.......Build internal edge solution for node 1
        ZE_RT = 0.D0
        QE_RT = 0.D0
      
        DO I = 1,P+1
          ZE_RT = ZE_RT + ZE(I,1,IRK)*PHIB(I,1)
          QE_RT = QE_RT + QE(I,1,IRK)*PHIB(I,1)
        END DO
!.......Left hand side boundary condition (reflection)
        IF (ADJUSTL(TRIM(BOUNDTYPE)).EQ.'radiant') THEN
          ZE_LT = ZE_RT
          QE_LT = QE_RT
        ELSE        ! Else reflective boundary
          ZE_LT =  ZE_RT
          QE_LT = -QE_RT
        END IF
        
        
      
        HE_RT = ZE_RT+DE_ED(1)
        HE_LT = ZE_LT+DE_ED(1)
!.......Calculate left and right flux values
        FH_LT = QE_LT
        FH_RT = QE_RT
        
        GH_LT = (QE_LT**2)/HE_LT + C12*G*ZE_LT*(2.D0*DE_ED(1)+ZE_LT)
        GH_RT = (QE_RT**2)/HE_RT + C12*G*ZE_RT*(2.D0*DE_ED(1)+ZE_RT)
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
      
        DO I = 1,P+1
          ZE_LT = ZE_LT + ZE(I,NE,IRK)*PHIB(I,2)
          QE_LT = QE_LT + QE(I,NE,IRK)*PHIB(I,2)
        END DO
!.......Right hand side boundary condition (reflection)
        IF (ADJUSTL(TRIM(BOUNDTYPE)).EQ.'radiant') THEN
          ZE_RT = ZE_LT
          QE_RT = QE_LT
        ELSE
          ZE_RT =  ZE_LT
          QE_RT = -QE_LT
        END IF
        
        
      
        HE_RT = ZE_RT+DE_ED(NN)
        HE_LT = ZE_LT+DE_ED(NN)
!.......Calculate left and right flux values
        FH_LT = QE_LT
        FH_RT = QE_RT
        
        GH_LT = (QE_LT**2)/HE_LT + C12*G*ZE_LT*(2.D0*DE_ED(1)+ZE_LT)
        GH_RT = (QE_RT**2)/HE_RT + C12*G*ZE_RT*(2.D0*DE_ED(1)+ZE_RT)
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
!             DO I = 1,2
!               HE_IN(I) = MAX(HE_IN(I),ZERO)
!             END DO
            ! If HE_IN(1) > HE_IN(2) then left side of element is higher
            IF (HE_IN(1).GT.HE_IN(2)) THEN
              HHAT(2) = MAX(H0,HE_IN(2))
              HHAT(1) = HE_IN(1)-(HHAT(2)-HE_IN(2))
            ! Otherwise the right side of the element is higher            
            ELSE
              HHAT(1) = MAX(H0,HE_IN(1))
              HHAT(2) = HE_IN(2)-(HHAT(1)-HE_IN(1))
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