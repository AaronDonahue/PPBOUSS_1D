!==========================================================================!
!==========================================================================!
      SUBROUTINE PREP_DG
      
      USE READ_DGINP, ONLY : READ_INPUT,IHOT,NWP
      USE GLOBALS,    ONLY : IRK
      
      IMPLICIT NONE
      
!.....Read control file
      CALL READ_INPUT
!.....Setup the variables for Runge-Kutta timestepping
      CALL GETRK
      CALL GAUSSPOINTS
!.....Read in grid
      CALL READGRID
!.....Allocate Variables
      CALL VARI_ALLOCATE
!.....Setup the DG-Basis Functions
      CALL DG_BASIS_FUNCTIONS
!.....Using basis functions setup Mass Matrix
      CALL MASS_MATRIX
!.....Read in hotstart file, if applicable
      IF (IHOT.EQ.1) THEN
        CALL READHOTSTART
      END IF
!.....Timestep
      CALL SETUP_TIME
!.....Initialize the wet/dry status of the solution
      IRK = 0
      CALL WETDRY
!.....Setup Nodal Attributes, if needed
      IF (NWP.GT.0) THEN
        CALL GET_NODAL_ATTR
      END IF

      RETURN
      END SUBROUTINE PREP_DG
!==========================================================================!
!==========================================================================!
!--------------------------------------------------------------------------!
      SUBROUTINE VARI_ALLOCATE

      USE GLOBALS,    ONLY : PHI,DPHI,PHIB,DPHIB,PSI,DPSI,PSIB,DPSIB,      &
     &                    ATVD,BTVD,TTVD,MDG,DGPIV,MCG,CGPIV,NE,ZE,QE,     &
     &                    ZE_RHS,QE_RHS,WDFLG,PD,PB,MANN,SPNG_GEN,SPNG_ABS
      USE READ_DGINP, ONLY : P,NRK,NEGP
      
      IMPLICIT NONE
      
!.....Basis functions
      ALLOCATE(PHI(P+1,NEGP),DPHI(P+1,NEGP))
      ALLOCATE(PHIB(P+1,2),DPHIB(P+1,2))
      ALLOCATE(PSI(P+1,NEGP),DPSI(P+1,NEGP))
      ALLOCATE(PSIB(P+1,2),DPSIB(P+1,2))      
!.....Mass matrices
      ALLOCATE(MDG(P+1,NEGP),DGPIV(P+1),MCG(3*P+1,NE*P+1),CGPIV(NE*P+1))
!.....Solution Matrices
      ALLOCATE(ZE(P+1,NE,NRK+1),QE(P+1,NE,NRK+1))
      ALLOCATE(ZE_RHS(P+1,NE,NRK),QE_RHS(P+1,NE,NRK))
!.....Solution Status Matrices
      ALLOCATE(WDFLG(NE))
!.....Pressure Matrices
      ALLOCATE(PD(P+1,NE,NRK+1),PB(P+1,NE,NRK+1))
!.....Nodal Attributes
      ALLOCATE(MANN(NE))
      ALLOCATE(SPNG_GEN(NE),SPNG_ABS(NE))

!.....Initialize certain variables
      PHI(:,:)      = 0.D0
      DPHI(:,:)     = 0.D0
      PHIB(:,:)     = 0.0D0
      DPHIB(:,:)    = 0.D0
      PSI(:,:)      = 0.D0
      DPSI(:,:)     = 0.D0
      PSIB(:,:)     = 0.D0
      DPSIB(:,:)    = 0.D0
      MDG(:,:)      = 0.D0
      MCG(:,:)      = 0.D0
      DGPIV(:)      = 1
      CGPIV(:)      = 1
      ZE(:,:,:)     = 0.D0
      QE(:,:,:)     = 0.D0
      ZE_RHS(:,:,:) = 0.D0
      QE_RHS(:,:,:) = 0.D0      
      WDFLG(:)      = 1
      PB(:,:,:)     = 0.D0
      PD(:,:,:)     = 0.D0
      MANN(:)       = 0.D0
      SPNG_GEN(:)   = 0.D0
      SPNG_ABS(:)   = 0.D0
      
      RETURN
      END SUBROUTINE VARI_ALLOCATE
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!      
      SUBROUTINE VARI_DEALLOCATE

      USE GLOBALS,    ONLY : PHI,DPHI,PHIB,DPHIB,PSI,DPSI,PSIB,DPSIB,      &
     &                    ATVD,BTVD,TTVD,MDG,DGPIV,MCG,CGPIV,ZE,QE,        &
     &                    ZE_RHS,QE_RHS,WEGP,XEGP,ATVD,BTVD,TTVD,X,LE,     &
     &                    DE_ED,DE_IN,DX_IN,WDFLG,PD,PB,MANN,SPNG_GEN,     &
     &                    SPNG_ABS
      
      IMPLICIT NONE
      
!.....Basis functions
      DEALLOCATE(PHI,DPHI)
      DEALLOCATE(PHIB,DPHIB)
      DEALLOCATE(PSI,DPSI)
      DEALLOCATE(PSIB,DPSIB)      
!.....Mass matrices
      DEALLOCATE(MDG,DGPIV,MCG,CGPIV)
!.....Solution Matrices
      DEALLOCATE(ZE,QE)
      DEALLOCATE(ZE_RHS,QE_RHS)
!.....Gauss Points vectors
      DEALLOCATE(WEGP,XEGP)
!.....Runge-Kutta vectors
      DEALLOCATE(ATVD,BTVD,TTVD)
!.....Mesh variables
      DEALLOCATE(X,DE_ED,LE)
      DEALLOCATE(DE_IN,DX_IN)
!.....Status variables
      DEALLOCATE(WDFLG)
!.....Pressure Matrices
      DEALLOCATE(PD,PB)
!.....Nodal Attributes
      DEALLOCATE(MANN)
      DEALLOCATE(SPNG_GEN,SPNG_ABS)
      
      RETURN
      END SUBROUTINE VARI_DEALLOCATE
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!
      SUBROUTINE DG_BASIS_FUNCTIONS
      
      USE GLOBALS,    ONLY : XEGP,PHI,DPHI,PHIB,DPHIB,PSI,DPSI,PSIB,DPSIB
      USE SIZES,      ONLY : SZ
      USE READ_DGINP, ONLY : P,NEGP,DGBASIS
      
      IMPLICIT NONE
      INTEGER	:: I,K
      REAL(SZ)  :: GETNODAL,GETMODAL
      
!.....Initialize basis function values at integration points and boundaries      
      PHI(:,:)   = 0.d0
      DPHI(:,:)  = 0.d0
      PHIB(:,:)  = 0.d0
      DPHIB(:,:) = 0.d0
      
      PSI(:,:)   = 0.d0
      DPSI(:,:)  = 0.d0
      PSIB(:,:)  = 0.d0
      DPSIB(:,:) = 0.d0
!.....Construct Nodal Basis Functions
      DO I = 1,P+1
        DO K = 1,NEGP
          PSI(I,K)  = GETNODAL(P,I,XEGP(K),0)
          DPSI(I,K) = GETNODAL(P,I,XEGP(K),1)
        END DO
        PSIB(I,1)  = GETNODAL(P,I,-1.D0,0)
        PSIB(I,2)  = GETNODAL(P,I, 1.D0,0)
        DPSIB(I,1) = GETNODAL(P,I,-1.D0,0)
        DPSIB(I,2) = GETNODAL(P,I, 1.D0,0)
      END DO
!.....Construct DG Basis Functions
      IF (ADJUSTL(TRIM(DGBASIS)).EQ.'NODAL'.OR.                            &
     &             ADJUSTL(TRIM(DGBASIS)).EQ.'nodal') THEN
        DO I = 1,P+1
          DO K = 1,NEGP
            PHI(I,K)  = GETNODAL(P,I,XEGP(K),0)
            DPHI(I,K) = GETNODAL(P,I,XEGP(K),1)
          END DO
          PHIB(I,1)  = GETNODAL(P,I,-1.D0,0)
          PHIB(I,2)  = GETNODAL(P,I, 1.D0,0)
          DPHIB(I,1) = GETNODAL(P,I,-1.D0,0)
          DPHIB(I,2) = GETNODAL(P,I, 1.D0,0)
        END DO
      ELSE
        DO I = 1,P+1
          DO K = 1,NEGP
            PHI(I,K)  = GETMODAL(I,XEGP(K),0)
            DPHI(I,K) = GETMODAL(I,XEGP(K),1)
          END DO
          PHIB(I,1)  = GETMODAL(I,-1.D0,0)
          PHIB(I,2)  = GETMODAL(I, 1.D0,0)
          DPHIB(I,1) = GETMODAL(I,-1.D0,0)
          DPHIB(I,2) = GETMODAL(I, 1.D0,0)
        END DO
      END IF
      
      PRINT*, '----- DG Basis Functions -----'      
      DO I = 1,P+1
        PRINT "(F,<NEGP>F,F)", PHIB(I,1), (PHI(I,K),K=1,NEGP), PHIB(I,2)
      END DO
      
      RETURN
      END SUBROUTINE DG_BASIS_FUNCTIONS
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!
      SUBROUTINE MASS_MATRIX
      
      USE GLOBALS,    ONLY : NE,MDG,WEGP,PHI,DGPIV,MCG,CGPIV,PSI,PHI,LE
      USE SIZES,      ONLY : SZ,C12
      USE READ_DGINP, ONLY : P,NEGP
      
      IMPLICIT NONE
      INTEGER	:: L,J,I,II,K,INFO
      INTEGER   :: LOCI,LOCJ
      REAL(SZ)  :: BDG(P+1)
      REAL(SZ)  :: MLOC(P+1,P+1),BCG(NE*P+1)
      
      INTEGER   :: LDAB,N
      
!.....Build local mass matrix (MDG)      
      DO I = 1,P+1
        DO II = 1,P+1
          MDG(I,II) = 0.D0
          DO K = 1,NEGP
            MDG(I,II) = MDG(I,II) + WEGP(K)*PHI(I,K)*PHI(II,K)
          END DO
        END DO
        BDG(I) = 1.D0
      END DO
      
      PRINT "(A)", 'Local DG Mass Matrix:'
      DO I = 1,P+1
      PRINT "(<p+1>f8.4)", (MDG(I,J),J=1,P+1)
      END DO
      
      CALL DGESV(P+1,1,MDG,P+1,DGPIV,BDG,P+1,INFO)
      
!.....Build global mass matrix (MCG)
      LDAB = 3*P+1
      N    = NE*P+1
      ! Build local nodal mass matrix
      DO I = 1,P+1
        DO II = 1,P+1
          MLOC(I,II) = 0.D0
          DO K = 1,NEGP
            MLOC(I,II) = MLOC(I,II) + WEGP(K)*PSI(I,K)*PSI(II,K)
          END DO
        END DO
      END DO
      ! Using local mass matrix build global matrix
      MCG(:,:) = 0.D0
      DO L = 1,NE
        J = (L-1)*P+1
        DO I = 1,P+1
          LOCI = J + I - 1
          DO II = 1,P+1
            LOCJ = J + II - 1
            MCG(2*P+1+LOCI-LOCJ,LOCJ) = MCG(2*P+1+LOCI-LOCJ,LOCJ) +        &
     &           MLOC(I,II)*C12*LE(L)
          END DO
          BCG(LOCI) = 1.D0
        END DO
      END DO
      
      CALL DGBTRF(N,N,P,P,MCG,LDAB,CGPIV,INFO)    
      
      
      RETURN
      END SUBROUTINE MASS_MATRIX
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!
      SUBROUTINE GAUSSPOINTS()
      
      USE GLOBALS,    ONLY : WEGP,XEGP
      USE READ_DGINP, ONLY : NEGP
 
      IMPLICIT NONE

!.....Allocate Gauss Points vectors
      ALLOCATE(WEGP(NEGP),XEGP(NEGP))

!.....Assign Gauss Points
      IF (NEGP.EQ.2) THEN

         XEGP(1) = -0.57735026918963D0
         XEGP(2) =  0.57735026918963D0
         
         WEGP(1) = 1.D0
         WEGP(2) = 1.D0
         
      ELSEIF (NEGP.EQ.3) THEN
         
         XEGP(1) = -0.77459666924148D0
         XEGP(2) =  0.D0
         XEGP(3) =  0.77459666924148D0
         
         WEGP(1) =  0.55555555555556D0
         WEGP(2) =  0.88888888888888D0
         WEGP(3) =  0.55555555555556D0
         
      ELSEIF (NEGP.EQ.4) THEN

         XEGP(1) = -0.86113631159405D0
         XEGP(2) = -0.33998104358486D0
         XEGP(3) =  0.33998104358486D0
         XEGP(4) =  0.86113631159405D0
         
         WEGP(1) =  0.34785484513745D0
         WEGP(2) =  0.65214515486255D0
         WEGP(3) =  0.65214515486255D0
         WEGP(4) =  0.34785484513745D0

      ELSEIF (NEGP.EQ.5) THEN

         XEGP(1) = -0.90617984593866D0
         XEGP(2) = -0.53846931010568D0
         XEGP(3) =  0.D0
         XEGP(4) =  0.53846931010568D0
         XEGP(5) =  0.90617984593866D0
         
         WEGP(1) =  0.23692688505619D0
         WEGP(2) =  0.47862867049937D0
         WEGP(3) =  0.56888888888889D0
         WEGP(4) =  0.47862867049937D0
         WEGP(5) =  0.23692688505619D0

      ELSEIF (NEGP.EQ.6) THEN

         XEGP(1) = -0.93246951420315D0
         XEGP(2) = -0.66120938646626D0
         XEGP(3) = -0.23861918608320D0
         XEGP(4) =  0.23861918608320D0
         XEGP(5) =  0.66120938646626D0
         XEGP(6) =  0.93246951420315D0
         
         WEGP(1) =  0.17132449237917D0
         WEGP(2) =  0.36076157304814D0
         WEGP(3) =  0.46791393457269D0
         WEGP(4) =  0.46791393457269D0
         WEGP(5) =  0.36076157304814D0
         WEGP(6) =  0.17132449237917D0

      ELSEIF (NEGP.EQ.7) THEN
         
         XEGP(1) = -0.94910791234276D0
         XEGP(2) = -0.74153118559939D0
         XEGP(3) = -0.40584515137740D0
         XEGP(4) =  0.D0
         XEGP(5) =  0.40584515137740D0
         XEGP(6) =  0.74153118559939D0
         XEGP(7) =  0.94910791234276D0
         
         WEGP(1) =  0.12948496616887D0
         WEGP(2) =  0.27970539148928D0
         WEGP(3) =  0.38183005050512D0
         WEGP(4) =  0.41795918367347D0
         WEGP(5) =  0.38183005050512D0
         WEGP(6) =  0.27970539148928D0
         WEGP(7) =  0.12948496616887D0

      ELSEIF (NEGP.EQ.8) THEN
         
         XEGP(1) = -0.96028985649754D0
         XEGP(2) = -0.79666647741363D0
         XEGP(3) = -0.52553240991633D0
         XEGP(4) = -0.18343464249565D0
         XEGP(5) =  0.18343464249565D0
         XEGP(6) =  0.52553240991633D0
         XEGP(7) =  0.79666647741363D0
         XEGP(8) =  0.96028985649754D0
         
         WEGP(1) =  0.10122853629038D0
         WEGP(2) =  0.22238103445337D0
         WEGP(3) =  0.31370664587789D0
         WEGP(4) =  0.36268378337836D0
         WEGP(5) =  0.36268378337836D0
         WEGP(6) =  0.31370664587789D0
         WEGP(7) =  0.22238103445337D0
         WEGP(8) =  0.10122853629038D0

      ELSEIF (NEGP.EQ.9) THEN
         
         XEGP(1) = -0.96816023950763D0
         XEGP(2) = -0.83603110732664D0
         XEGP(3) = -0.61337143270059D0
         XEGP(4) = -0.32425342340381D0
         XEGP(5) =  0.D0
         XEGP(6) =  0.32425342340381D0
         XEGP(7) =  0.61337143270059D0
         XEGP(8) =  0.83603110732664D0
         XEGP(9) =  0.96816023950763D0
         
         WEGP(1) = 0.08127438836163D0
         WEGP(2) = 0.18064816069483D0
         WEGP(3) = 0.26061069640294D0
         WEGP(4) = 0.31234707704000D0
         WEGP(5) = 0.33023935500126D0
         WEGP(6) = 0.31234707704000D0
         WEGP(7) = 0.26061069640294D0
         WEGP(8) = 0.18064816069483D0
         WEGP(9) = 0.08127438836163D0
         
      ELSEIF (NEGP.EQ.10) THEN
         
         XEGP(1)  = -0.97390652851717D0
         XEGP(2)  = -0.86506336668898D0
         XEGP(3)  = -0.67940956829902D0
         XEGP(4)  = -0.43339539412925D0
         XEGP(5)  = -0.14887433898163D0
         XEGP(6)  =  0.14887433898163D0
         XEGP(7)  =  0.43339539412925D0
         XEGP(8)  =  0.67940956829902D0
         XEGP(9)  =  0.86506336668898D0
         XEGP(10) =  0.97390652851717D0

         WEGP(1)  =  0.06667134430869D0
         WEGP(2)  =  0.14945134915058D0
         WEGP(3)  =  0.21908636251598D0
         WEGP(4)  =  0.26926671931000D0
         WEGP(5)  =  0.29552422471475D0
         WEGP(6)  =  0.29552422471475D0
         WEGP(7)  =  0.26926671931000D0
         WEGP(8)  =  0.21908636251598D0
         WEGP(9)  =  0.14945134915058D0
         WEGP(10) =  0.06667134430869D0

      ELSEIF (NEGP.EQ.11) THEN
         
         XEGP(1)  = -0.97822865814606D0
         XEGP(2)  = -0.88706259976810D0
         XEGP(3)  = -0.73015200557405D0
         XEGP(4)  = -0.51909612920681D0
         XEGP(5)  = -0.26954315595234D0
         XEGP(6)  =  0.D0
         XEGP(7)  =  0.26954315595234D0
         XEGP(8)  =  0.51909612920681D0
         XEGP(9)  =  0.73015200557405D0
         XEGP(10) =  0.88706259976810D0
         XEGP(11) =  0.97822865814606D0
         
         WEGP(1)  =  0.05566856711627D0
         WEGP(2)  =  0.12558036946485D0
         WEGP(3)  =  0.18629021092774D0
         WEGP(4)  =  0.23319376459199D0
         WEGP(5)  =  0.26280454451025D0
         WEGP(6)  =  0.27292508677790D0
         WEGP(7)  =  0.26280454451025D0
         WEGP(8)  =  0.23319376459199D0
         WEGP(9)  =  0.18629021092774D0
         WEGP(10) =  0.12558036946485D0
         WEGP(11) =  0.05566856711627D0
         
      ELSEIF (NEGP.EQ.12) THEN
         
         XEGP(1)  = -0.98156063424672D0
         XEGP(2)  = -0.90411725637047D0
         XEGP(3)  = -0.76990267419430D0
         XEGP(4)  = -0.58731795428662D0
         XEGP(5)  = -0.36783149899818D0
         XEGP(6)  = -0.12523340851147D0
         XEGP(7)  =  0.12523340851147D0
         XEGP(8)  =  0.36783149899818D0
         XEGP(9)  =  0.58731795428662D0
         XEGP(10) =  0.76990267419430D0
         XEGP(11) =  0.90411725637047D0
         XEGP(12) =  0.98156063424672D0
         
         WEGP(1)  =  0.04717533638677D0
         WEGP(2)  =  0.10693932599520D0
         WEGP(3)  =  0.16007832854334D0
         WEGP(4)  =  0.20316742672308D0
         WEGP(5)  =  0.23349253653835D0
         WEGP(6)  =  0.24914704581340D0
         WEGP(7)  =  0.24914704581340D0
         WEGP(8)  =  0.23349253653835D0
         WEGP(9)  =  0.20316742672308D0
         WEGP(10) =  0.16007832854334D0
         WEGP(11) =  0.10693932599520D0
         WEGP(12) =  0.04717533638677D0

      ELSEIF (NEGP.EQ.13) THEN
         
         XEGP(1)  = -0.98418305471859D0
         XEGP(2)  = -0.91759839922298D0
         XEGP(3)  = -0.80157809073331D0
         XEGP(4)  = -0.64234933944034D0
         XEGP(5)  = -0.44849275103645D0
         XEGP(6)  = -0.23045831595513D0
         XEGP(7)  =  0.D0
         XEGP(8)  =  0.23045831595513D0
         XEGP(9)  =  0.44849275103645D0
         XEGP(10) =  0.64234933944034D0
         XEGP(11) =  0.80157809073331D0
         XEGP(12) =  0.91759839922298D0
         XEGP(13) =  0.98418305471859D0

         WEGP(1)  =  0.04048400476532D0
         WEGP(2)  =  0.09212149983773D0
         WEGP(3)  =  0.13887351021979D0
         WEGP(4)  =  0.17814598076195D0
         WEGP(5)  =  0.20781604753689D0
         WEGP(6)  =  0.22628318026290D0
         WEGP(7)  =  0.23255155323087D0
         WEGP(8)  =  0.22628318026290D0
         WEGP(9)  =  0.20781604753689D0
         WEGP(10) =  0.17814598076195D0
         WEGP(11) =  0.13887351021979D0
         WEGP(12) =  0.09212149983773D0
         WEGP(13) =  0.04048400476532D0

      ENDIF


      RETURN
      END
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!
      SUBROUTINE GETRK()

      USE GLOBALS,    ONLY : ATVD,BTVD,TTVD
      USE SIZES,      ONLY : SZ
      USE READ_DGINP, ONLY : NRK

      IMPLICIT NONE
      INTEGER I,J
      REAL(SZ) :: CTVD(NRK,NRK),DTVD(NRK)

      ALLOCATE(ATVD(NRK,NRK),BTVD(NRK,NRK),TTVD(NRK))
      
      TTVD(:)   = 0.D0

!.....The forward Euler method

      IF (NRK.EQ.1) THEN

         ATVD(:,:) = 0.D0
         BTVD(:,:) = 0.D0
         CTVD(:,:) = 0.D0
         DTVD(:)   = 0.D0

         ATVD(1,1) = 1.D0
         BTVD(1,1) = 1.D0

!.....SSP(s,2) schemes

      ELSEIF (NRK.EQ.2) THEN

         ATVD(:,:) = 0.D0
         BTVD(:,:) = 0.D0
         CTVD(:,:) = 0.D0
         DTVD(:)   = 0.D0

         DO I = 1,NRK
            DO J = 0,NRK-1

               IF ((J.EQ.(I-1)).AND.(I.LT.NRK)) THEN
                  ATVD(I,J+1) = 1.D0
                  BTVD(I,J+1) = 1.D0/(NRK-1)
               ELSEIF ((J.EQ.0).AND.(I.EQ.NRK)) THEN
                  ATVD(I,J+1) = 1.D0/NRK
               ELSEIF ((J.EQ.(NRK-1)).AND.(I.EQ.NRK)) THEN
                  ATVD(I,J+1) = (NRK-1.D0)/NRK
                  BTVD(I,J+1) = 1.D0/NRK
               ENDIF

            ENDDO
         ENDDO

         TTVD(2)   = 1.D0

!.....SSP(3,3) scheme

      ELSEIF (NRK.EQ.3) THEN

         ATVD(:,:) = 0.D0
         BTVD(:,:) = 0.D0
         CTVD(:,:) = 0.D0
         DTVD(:)   = 0.D0

         ATVD(1,1) = 1.D0
         ATVD(2,1) = 3.D0/4.D0
         ATVD(2,2) = 1.D0/4.D0
         ATVD(3,1) = 1.D0/3.D0
         ATVD(3,3) = 2.D0/3.D0

         BTVD(1,1) = 1.D0
         BTVD(2,2) = 1.D0/4.D0
         BTVD(3,3) = 2.D0/3.D0

         TTVD(2)   = 1.D0
         TTVD(3)   = 1.D0/2.D0

!.....SSP(4,3) scheme

      ELSEIF (NRK.EQ.4) THEN

         ATVD(:,:) = 0.D0
         BTVD(:,:) = 0.D0
         CTVD(:,:) = 0.D0
         DTVD(:)   = 0.D0

         ATVD(1,1) = 1.D0
         ATVD(2,2) = 1.D0
         ATVD(3,1) = 2.D0/3.D0
         ATVD(3,3) = 1.D0/3.D0
         ATVD(4,4) = 1.D0

         BTVD(1,1) = 1.D0/2.D0
         BTVD(2,2) = 1.D0/2.D0
         BTVD(3,3) = 1.D0/6.D0
         BTVD(4,4) = 1.D0/2.D0

!.....SSP(5,3) scheme

      ELSEIF (NRK.EQ.5) THEN

         ATVD(:,:) = 0.D0
         BTVD(:,:) = 0.D0
         CTVD(:,:) = 0.D0
         DTVD(:)   = 0.D0

         ATVD(1,1) = 1.D0
         ATVD(2,2) = 1.D0
         ATVD(3,1) = 0.355909775063327D0
         ATVD(3,3) = 0.644090224936674D0
         ATVD(4,1) = 0.367933791638137D0
         ATVD(4,4) = 0.632066208361863D0
         ATVD(5,3) = 0.237593836598569D0
         ATVD(5,5) = 0.762406163401431D0

         BTVD(1,1) = 0.377268915331368D0
         BTVD(2,2) = 0.377268915331368D0
         BTVD(3,3) = 0.242995220537396D0
         BTVD(4,4) = 0.238458932846290D0
         BTVD(5,5) = 0.287632146308408D0

!.....SSP(6,3) scheme

      ELSEIF (NRK.EQ.6) THEN

         ATVD(:,:) = 0.D0
         BTVD(:,:) = 0.D0
         CTVD(:,:) = 0.D0
         DTVD(:)   = 0.D0

         ATVD(1,1) = 1.D0
         ATVD(2,2) = 1.D0
         ATVD(3,3) = 1.D0
         ATVD(4,1) = 0.476769811285196D0
         ATVD(4,2) = 0.098511733286064D0
         ATVD(4,4) = 0.424718455428740D0
         ATVD(5,5) = 1.D0
         ATVD(6,3) = 0.155221702560091D0
         ATVD(6,6) = 0.844778297439909D0

         BTVD(1,1) = 0.284220721334261D0
         BTVD(2,2) = 0.284220721334261D0
         BTVD(3,3) = 0.284220721334261D0
         BTVD(4,4) = 0.120713785765930D0
         BTVD(5,5) = 0.284220721334261D0
         BTVD(6,6) = 0.240103497065900D0
       END IF



      RETURN
      END SUBROUTINE GETRK
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!
      FUNCTION GETNODAL(PIN,ORD,XIN,DER) RESULT(YOUT)
      
      USE SIZES,      ONLY : SZ,C12
      
      IMPLICIT NONE
      INTEGER  :: ORD,DER,PIN
      REAL(SZ) :: XIN,YOUT
      INTEGER  :: I,J
      REAL(SZ) :: NUM,DEN,XOR,XI
      
      XOR = -1.D0 + ((ORD-1)*2.D0)/PIN
      ! Function
      IF (DER.EQ.0) THEN
        NUM = 1.D0
        DEN = 1.D0        
        DO I = 1,PIN+1        
          IF (I.NE.ORD) THEN
            XI = -1.D0 + ((I-1)*2.D0)/PIN
            NUM = NUM * (XIN-XI)
            DEN = DEN * (XOR-XI)
          END IF
        END DO
        YOUT = NUM/DEN
      ! Derivative
      ELSEIF (DER.EQ.1) THEN
        YOUT = 0.D0
        DEN  = 1.D0
        DO I = 1,PIN+1
          IF (I.NE.ORD) THEN
            XI = -1.D0 + ((I-1)*2.D0)/PIN
            DEN = DEN * (XOR-XI)
          END IF
        END DO
        
        DO J = 1,PIN+1
          IF (J.NE.ORD) THEN
            NUM = 1.D0
            DO I = 1,PIN+1
              IF (I.NE.ORD.AND.I.NE.J) THEN
                XI = -1.D0 + ((I-1)*2.D0)/PIN
                NUM = NUM * (XIN-XI)
              END IF
            END DO
            YOUT = YOUT + NUM/DEN
          END IF
        END DO
      END IF
      
      END FUNCTION GETNODAL
!--------------------------------------------------------------------------!
!==========================================================================!
!--------------------------------------------------------------------------!
      FUNCTION GETMODAL(ORD,XIN,DER) RESULT(YOUT)
      
      USE SIZES, ONLY : SZ,C12
      
      IMPLICIT NONE
      INTEGER  :: ORD,DER
      REAL(SZ) :: XIN,YOUT
      
      IF (DER.EQ.0) THEN
        IF (ORD.EQ.1) THEN
          YOUT = 1.D0
        ELSEIF (ORD.EQ.2) THEN
          YOUT = XIN
        ELSEIF (ORD.EQ.3) THEN
          YOUT = C12*(3*XIN**2-1.D0)
        ELSEIF (ORD.EQ.4) THEN
          YOUT = C12*(5*XIN**3-3*XIN)
        END IF
      ELSEIF (DER.EQ.1) THEN
        IF (ORD.EQ.1) THEN
          YOUT = 0.D0
        ELSEIF (ORD.EQ.2) THEN
          YOUT = 1.D0
        ELSEIF (ORD.EQ.3) THEN
          YOUT = C12*(6*XIN)
        ELSEIF (ORD.EQ.4) THEN
          YOUT = C12*(15*XIN**2-3.D0)
        END IF
      END IF
      
      END FUNCTION GETMODAL
!..........................................................................!
!==========================================================================!
!..........................................................................!      
      SUBROUTINE READGRID
      
      USE READ_DGINP, ONLY : GRID_FILE,NEGP
      USE GLOBALS,    ONLY : RUNNAME,X,DE_ED,LE,NE,NN,XEGP,DE_IN,DX_IN
      USE SIZES,      ONLY : SZ,C12
      
      IMPLICIT NONE
      INTEGER   :: I,L,M,NODE
      REAL(SZ)  :: D1,D2
      REAL(SZ)  :: NODAL1,NODAL2,GETNODAL
      
      LOGICAL   :: file_exists
      
      
      INQUIRE(FILE='./'//ADJUSTL(TRIM(grid_file)), EXIST = file_exists)
      IF(file_exists == .FALSE.) THEN
        PRINT*, "grid file file does not exist"
        CALL EXIT
      ENDIF
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "              Read grid file                 "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "
      OPEN(UNIT=14,FILE='./'//ADJUSTL(TRIM(grid_file)),ACTION='READ')
      READ(14,"(A100)") RUNNAME
      READ(14,*)  NE,NN
      
!.....Allocate mesh variables
      ALLOCATE(X(NN),DE_ED(NN),LE(NE))
      ALLOCATE(DE_IN(NE,NEGP),DX_IN(NE,NEGP))
!.....Read nodal information from the grid file
      DO I = 1,NN
        READ(14,*) NODE, X(I), DE_ED(I)
      END DO
!.....Determine the element size for each element
      DO L = 1,NE
        LE(L) = X(L+1)-X(L)
      END DO
!.....Setup bathymetry matrices using grid information
      DO L = 1,NE
        DO M = 1,NEGP
          NODAL1 = GETNODAL(1,1,XEGP(M),0)
          NODAL2 = GETNODAL(1,2,XEGP(M),0)
          DE_IN(L,M) = DE_ED(L)*NODAL1+DE_ED(L+1)*NODAL2
          
          NODAL1 = GETNODAL(1,1,XEGP(M),1)
          NODAL2 = GETNODAL(1,2,XEGP(M),1)
          DX_IN(L,M) = (DE_ED(L)*NODAL1+DE_ED(L+1)*NODAL2)/(C12*LE(L))
        END DO
      END DO
      
      
      CLOSE(14)
      
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "        Finished reading grid file           "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "
      
      
      RETURN
      END SUBROUTINE READGRID
!..........................................................................!
!==========================================================================!
!..........................................................................!
      SUBROUTINE READHOTSTART
      
      USE READ_DGINP, ONLY : P,HOTSTART_FILE,DGBASIS
      USE GLOBALS,    ONLY : NE,ZE,QE
      USE SIZES,      ONLY : SZ
      
      IMPLICIT NONE
      
      INTEGER              :: NNHOT,READ_STAT,INFO
      INTEGER              :: I,II,L,NN,IPT
      REAL(SZ)             :: PT,DPT
      REAL(SZ),ALLOCATABLE :: XLOC(:),ZLOC(:),QLOC(:)
      INTEGER,ALLOCATABLE  :: NDLOC(:)
      
      LOGICAL              :: FILE_EXISTS
      CHARACTER(LEN=100)   :: JC1,HOTTYPE
      
      REAL(SZ)             :: VANDER(P+1,P+1)
      REAL(SZ)             :: BZ(P+1),BQ(P+1)
      INTEGER              :: VPIV(P+1)
      REAL(SZ)             :: GETMODAL
      
      
      INQUIRE(FILE='./'//ADJUSTL(TRIM(hotstart_file)), EXIST = file_exists)
      IF(file_exists == .FALSE.) THEN
        PRINT*, "hotstart file does not exist"
        CALL EXIT
      ENDIF
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "           Read hotstart file                "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "
      OPEN(UNIT=63,FILE='./'//ADJUSTL(TRIM(hotstart_file)),ACTION='READ')
      
!.....Read in number of nodes in hotstart file and local P value
      READ(63,'(A100)') JC1      
      READ(63,'(A100)') HOTTYPE
      READ(63,*) NNHOT
      ALLOCATE(XLOC(NNHOT),ZLOC(NNHOT),QLOC(NNHOT),NDLOC(NNHOT))
!.....Based on HOTTYPE read in hotstart file appropriately
      SELECT CASE (ADJUSTL(TRIM(HOTTYPE)))
        ! If the hotstart is nodal
        CASE('NODAL')
          ! Check to make sure there is enough nodal information for the P
          IF (NNHOT.NE.(P*NE+1)) THEN
            PRINT("(A)"), "*** ERROR: The hotstart file does not contain enough ***"  
            PRINT("(A)"), "           information for the simulation order (P).    "  
            PRINT("(A)"), "!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!"
            CALL EXIT
          ELSE
            ! Read in hotstart information
            DO NN = 1,NNHOT
              READ(63,*,IOSTAT=READ_STAT) XLOC(NN),ZLOC(NN),QLOC(NN)
              IF (READ_STAT/=0.AND.NN.NE.NNHOT) THEN
                PRINT("(A)"), "*** ERROR: Number of nodes in hotstart file does    "  
                PRINT("(A)"), "           does not match value on line 1 of file.  "
                PRINT("(A)"), "!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!"
                CALL EXIT
              END IF
            END DO
            
            ! Determine dg coefficients from nodal data
            ! If DG basis set is nodal then simply read in values,
            IF (ADJUSTL(TRIM(DGBASIS)).EQ.'nodal') THEN
              DO L = 1,NE
                DO I = 1,P+1
                  IPT = (L-1)*P+I
                  ZE(I,L,1) = ZLOC(IPT)
                  QE(I,L,1) = QLOC(IPT)
                END DO
              END DO
            ! Otherwise it is modal, use the Vandermonde matrix to get modal values
            ELSE
              DPT = 2.D0/P
              DO I = 1,P+1
                PT = -1.D0 + (I-1)*DPT
                DO II = 1,P+1
                  VANDER(I,II) = GETMODAL(II,PT,0)
                END DO
              END DO
              CALL DGESV( P+1,1,VANDER,P+1,VPIV,BZ,P+1,INFO ) ! Solve Vandermonde matrix once
              DO L = 1,NE
                DO I = 1,P+1
                  IPT   = (L-1)*P+I
                  BZ(I) = ZLOC(IPT)
                  BQ(I) = QLOC(IPT)
                END DO
                CALL DGETRS( 'No transpose',P+1,1,VANDER,P+1,VPIV,BZ,P+1,INFO )
                CALL DGETRS( 'No transpose',P+1,1,VANDER,P+1,VPIV,BQ,P+1,INFO )
                DO I = 1,P+1
                  ZE(I,L,1) = BZ(I)
                  QE(I,L,1) = BQ(I)
                END DO
              END DO
            END IF
          END IF
        ! If the hotstart is modal
        CASE('MODAL')
          ! Check to make sure that there is enough modal information for the P
          IF (NNHOT.NE.((P+1)*NE)) THEN
            PRINT("(A)"), "*** ERROR: The hotstart file does not contain enough ***"  
            PRINT("(A)"), "           information for the simulation order (P).    "  
            PRINT("(A)"), "!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!"
            CALL EXIT
          ELSE
            ! Read in hotstart information
            DO NN = 1,NNHOT
              READ(63,*,IOSTAT=READ_STAT) NDLOC(NN),ZLOC(NN),QLOC(NN)
              IF (READ_STAT/=0.AND.NN.NE.NNHOT) THEN
                PRINT("(A)"), "*** ERROR: Number of nodes in hotstart file does    "  
                PRINT("(A)"), "           does not match value on line 1 of file.  "
                PRINT("(A)"), "!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!"
                CALL EXIT
              END IF
            END DO
            ! If DG basis set is modal then simply read in values
            IF (ADJUSTL(TRIM(DGBASIS)).EQ.'MODAL') THEN
              DO L = 1,NE
                DO I = 1,P+1
                  IPT = (L-1)*(P+1) + I
                  ZE(I,L,1) = ZLOC(IPT)
                  QE(I,L,1) = QLOC(IPT)
                END DO
              END DO
            ! Otherwise it is nodal, project solution onto nodal points
            ELSE
              DO L = 1,NE
                DO I = 1,P+1
                  IPT = (L-1)*(P+1) + I
                  PT  = -1.D0+(I-1)*2.D0/P
                  ZE(I,L,1) = 0.D0
                  QE(I,L,1) = 0.D0
                  DO II = 1,P+1
                    ZE(I,L,1) = ZE(I,L,1) + ZLOC(IPT)*GETMODAL(II,PT,0)
                    QE(I,L,1) = QE(I,L,1) + QLOC(IPT)*GETMODAL(II,PT,0)
                  END DO
                END DO
              END DO
            END IF
          END IF
        
        CASE DEFAULT
          PRINT("(A)"), "*** ERROR: The hotstart file does not have a ***"  
          PRINT("(A)"), "           valid type, i.e. NODAL or MODAL.             "  
          PRINT("(A)"), "           Check line 2 for spelling and uppercase      "  
          PRINT("(A)"), "!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!"
          CALL EXIT
      END SELECT
        
      
      
      DEALLOCATE(XLOC,ZLOC,QLOC,NDLOC)
      CLOSE(63)
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "     Finished reading hotstart file          "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "
      
      RETURN
      END SUBROUTINE READHOTSTART
!..........................................................................!
!==========================================================================!
!..........................................................................!
      SUBROUTINE SETUP_TIME
      
      USE READ_DGINP, ONLY : MAXTIME,CFL_ADJ,P,TIMESNAP
      USE GLOBALS,    ONLY : TIMESTEPS,DT,LE,NE,G,DE_ED,TSNAP,ZE,QE,PHIB
      USE SIZES,      ONLY : SZ
      
      IMPLICIT NONE
      INTEGER  :: L,I,K
      REAL(SZ) :: MINLE,LAMMAX,UEMAX
      REAL(SZ) :: ZI,QI,HI
      
!.....Determine DT based on gridsize and CFL condition
      MINLE  = 10**16      ! Gridsize
      LAMMAX = 1.d0        ! Maximum Depth
      
      DO L = 1,NE
        DO K = 1,2
          ZI = 0.D0
          QI = 0.D0
          DO I = 1,P+1
            ZI = ZI + ZE(I,L,1)*PHIB(I,K)
            QI = QI + QE(I,L,1)*PHIB(I,K)
          END DO
          HI = ZI + DE_ED(L+K-1)
          LAMMAX = MAX(LAMMAX,ABS(QI/HI+DSQRT(G*HI)))
        END DO
        MINLE = MIN(MINLE,LE(L))        
      END DO
      
      DT        = MINLE/(LAMMAX*(2*P+1))*CFL_ADJ
      TIMESTEPS = CEILING(MAXTIME/DT)
      DT        = MAXTIME/DBLE(TIMESTEPS)
      
!.....If a timesnap is designated then determine frequency of output
      IF (TIMESNAP.LE.-999D0) THEN
        TSNAP = TIMESTEPS+1 ! No output
      ELSE IF (TIMESNAP.LE.DT) THEN
        TSNAP = 1           ! Output at every step
      ELSE
        TSNAP = CEILING(TIMESNAP/DT)
        DT    = TIMESNAP/DBLE(TSNAP)
        TIMESTEPS = CEILING(MAXTIME/DT)
      END IF
            
      PRINT "(A)", " "
      PRINT "(A14,F16.8,A14)", "         DT = ",DT,"            "      
      PRINT "(A14,I,A14)", "Total Steps = ",TIMESTEPS,"            "
      PRINT "(A)", " "
      
      END SUBROUTINE SETUP_TIME
!..........................................................................!
!==========================================================================!
!..........................................................................!
      SUBROUTINE GET_NODAL_ATTR
      
      USE READ_DGINP, ONLY : NWP
      USE GLOBALS,    ONLY : MANN,NE,SPNG_GEN,SPNG_ABS
      USE SIZES,      ONLY : SZ
      
      IMPLICIT NONE
      LOGICAL              :: file_exists
      CHARACTER(LEN=100)   :: JC1,ATTR_NAME,ATTR_UNIT
      INTEGER              :: I,NA,NE_TST,NN_TST,NWP_TST,ELEM
      INTEGER              :: ATTR_SIZE,ATTR_NONDEF
      REAL(SZ),ALLOCATABLE :: ATTR_DEFAULT(:)
      REAL(SZ)             :: VALUE
      
!..... Make sure nodal attributes file exists
      INQUIRE(FILE='fort.13', EXIST = file_exists)
      IF(file_exists == .FALSE.) THEN
        PRINT*, "nodal attributes file (fort.13) does not exist"
        CALL EXIT
      ENDIF
!..... If so open it and read in nodal attributes
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "       Read nodal attributes file            "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "
      OPEN(UNIT=13,FILE='fort.13',ACTION='READ')
      
!.....Begin reading in nodal attributes
      READ(13,'(A)') JC1
      READ(13,*) NE_TST, NN_TST
      READ(13,*) NWP_TST
      ! Check to make sure the number of elements in file matches grid
      IF (NE_TST.NE.NE) THEN
        PRINT("(A)"), "*** ERROR: Number of elements in attributes file does not ***"  
        PRINT("(A)"), "           match the gridfile. "
        PRINT("(A)"), "!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!"
        STOP
      END IF
      ! Check to make sure the number of attributes matches NWP
      IF (NWP_TST.NE.NWP) THEN
        PRINT("(A)"), "*** ERROR: Number of attributes in file does not ***"  
        PRINT("(A)"), "           match the fort.wasupp file. "
        PRINT("(A)"), "!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!"
        STOP
      END IF
      ! Read in default values
      DO NA = 1,NWP
        READ(13,'(A)') ATTR_NAME
        READ(13,'(A)') ATTR_UNIT
        READ(13,*) ATTR_SIZE
        ALLOCATE(ATTR_DEFAULT(ATTR_SIZE))
        READ(13,'(<ATTR_SIZE>F)') (ATTR_DEFAULT(I),I=1,ATTR_SIZE)
        SELECT CASE (TRIM(ATTR_NAME))
          CASE ('mannings_n_at_sea_floor')
            DO I = 1,NE
              MANN(I) = ATTR_DEFAULT(1)
            END DO
          CASE ('sponge_generation_layer')
            DO I = 1,NE
              SPNG_GEN(I) = ATTR_DEFAULT(1)
            END DO
          CASE ('sponge_absorbing_layer')
            DO I = 1,NE
              SPNG_ABS(I) = ATTR_DEFAULT(1)
            END DO
          CASE DEFAULT
            PRINT("(A,A,A)"), "*** ERROR: Attribute ",TRIM(ATTR_NAME)," does not match available attributes ***"  
            PRINT("(A)"), "!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!"
            STOP
        END SELECT
        DEALLOCATE(ATTR_DEFAULT)
      END DO
      ! Adjust nondefault values
      DO NA = 1,NWP
        READ(13,'(A)') ATTR_NAME
        READ(13,*) ATTR_NONDEF
        IF (ATTR_NONDEF.GT.0) THEN
          SELECT CASE (TRIM(ATTR_NAME))
            CASE ('mannings_n_at_sea_floor')
              DO I = 1,ATTR_NONDEF
                READ(13,*) ELEM, VALUE
                MANN(ELEM) = VALUE
              END DO              
            CASE ('sponge_generation_layer')
              DO I = 1,ATTR_NONDEF
                READ(13,*) ELEM, VALUE
                SPNG_GEN(ELEM) = VALUE
              END DO
            CASE ('sponge_absorbing_layer')
              DO I = 1,ATTR_NONDEF
                READ(13,*) ELEM, VALUE
                SPNG_ABS(ELEM) = VALUE
              END DO
            CASE DEFAULT
              PRINT("(A,A,A)"), "*** ERROR: Attribute ",TRIM(ATTR_NAME)," does not match available attributes ***"  
              PRINT("(A)"), "!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!"
              STOP
          END SELECT          
        END IF
        PRINT '(A,A)','Finished loading ',ATTR_NAME
      END DO
      
      CLOSE(13)
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "   Finished reading nodal attributes file    "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "
      
      RETURN
      END SUBROUTINE GET_NODAL_ATTR
!..........................................................................!
!==========================================================================!
!==========================================================================!      