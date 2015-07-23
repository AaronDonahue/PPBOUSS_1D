      MODULE GLOBALS
      
      USE SIZES, ONLY : SZ
      
      IMPLICIT NONE
      
!.....Solution state variables
      INTEGER                 :: NE,NN
      REAL(SZ),ALLOCATABLE    :: LE(:)
      CHARACTER(LEN=100)      :: RUNNAME
      REAL(SZ),ALLOCATABLE    :: X(:)
!       CHARACTER(LEN=100)      :: grid_file ! In READ_DGINP
      
      REAL(SZ)                :: G=10.d0!9.80665D0
!.....DG control variables
!       INTEGER                 :: P ! In READ_DGINP
      
      REAL(SZ)                :: GHAT,FHAT
      REAL(SZ)                :: FH_LT,FH_RT,GH_LT,GH_RT
      REAL(SZ)                :: ZE_LT,ZE_RT,QE_LT,QE_RT,HE_LT,HE_RT
!.....Integration and Timestepping Variables
!       INTEGER                 :: NEGP ! In READ_DGINP
      REAL(SZ),ALLOCATABLE    :: WEGP(:),XEGP(:)
      
      REAL(SZ)                :: DT,TIME,TIME_RK ! MAXTIME  ! In READ_DGINP
      INTEGER                 :: TIMESTEPS,TSNAP
      REAL(SZ),ALLOCATABLE    :: ATVD(:,:),BTVD(:,:),TTVD(:)
      INTEGER                 :: IRK !NRK  ! In READ_DGINP
!.....DG Basis Functions
!       CHARACTER(LEN=100)      :: DGBASIS ! In READ_DGINP
      REAL(SZ),ALLOCATABLE    :: PHI(:,:),DPHI(:,:),PHIB(:,:),DPHIB(:,:)
      REAL(SZ),ALLOCATABLE    :: PSI(:,:),DPSI(:,:),PSIB(:,:),DPSIB(:,:)
      
      REAL(SZ),ALLOCATABLE    :: MDG(:,:),MCG(:,:)
      INTEGER,ALLOCATABLE     :: DGPIV(:),CGPIV(:)
!.....DG Solution variables
      REAL(SZ),ALLOCATABLE    :: ZE(:,:,:),QE(:,:,:)
      REAL(SZ),ALLOCATABLE    :: ZE_RHS(:,:,:),QE_RHS(:,:,:)
      REAL(SZ),ALLOCATABLE    :: DE_IN(:,:),DX_IN(:,:),DE_ED(:)
!.....Solution status variables
      INTEGER,ALLOCATABLE     :: WDFLG(:)
!.....Pressure Variables
      REAL(SZ),ALLOCATABLE    :: PD(:,:,:),PB(:,:,:)
!.....Nodal Attributes
      REAL(SZ),ALLOCATABLE    :: MANN(:)
      REAL(SZ),ALLOCATABLE    :: SPNG_GEN(:),SPNG_ABS(:)

      END MODULE GLOBALS