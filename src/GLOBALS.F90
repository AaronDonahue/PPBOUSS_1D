      MODULE GLOBALS
      
      USE SIZES, ONLY : SZ
      
      IMPLICIT NONE
      
!.....Solution state variables
      INTEGER                 :: NE,NN
      REAL(SZ),ALLOCATABLE    :: LE(:)
      CHARACTER(LEN=100)      :: RUNNAME
      REAL(SZ),ALLOCATABLE    :: X(:)
      
      REAL(SZ)                :: G=9.80665D0
!.....Output File Variables
      INTEGER                 :: SIMERROR
      REAL(SZ)                :: CPU_START,CPU_FINISH
      INTEGER                 :: FORT16
      INTEGER	              :: FORT631,FORT641,FORT731,FORT741
      INTEGER                 :: FORT611
      INTEGER                 :: FORTRUNUP
      INTEGER                 :: NUMSTNS
      INTEGER,ALLOCATABLE     :: STNELEM(:)
!.....DG control variables
      REAL(SZ)                :: GHAT,FHAT
      REAL(SZ)                :: FH_LT,FH_RT,GH_LT,GH_RT
      REAL(SZ)                :: ZE_LT,ZE_RT,QE_LT,QE_RT,HE_LT,HE_RT
!.....Integration and Timestepping Variables
      REAL(SZ),ALLOCATABLE    :: WEGP(:),XEGP(:)
      
      REAL(SZ)                :: DT,TIME,TIME_RK
      INTEGER                 :: TIMESTEPS,TSNAP,SSNAP
      REAL(SZ),ALLOCATABLE    :: ATVD(:,:),BTVD(:,:),TTVD(:)
      INTEGER                 :: IRK
!.....DG Basis Functions
      REAL(SZ),ALLOCATABLE    :: PHI(:,:),DPHI(:,:),PHIB(:,:),DPHIB(:,:)
      REAL(SZ),ALLOCATABLE    :: PSI(:,:),DPSI(:,:),PSIB(:,:),DPSIB(:,:)
      REAL(SZ),ALLOCATABLE    :: PHISTN(:,:)
      
      REAL(SZ),ALLOCATABLE    :: MDG(:,:),MCG(:,:)
      INTEGER,ALLOCATABLE     :: DGPIV(:),CGPIV(:)
!.....DG Solution variables
      REAL(SZ),ALLOCATABLE    :: ZE(:,:,:),QE(:,:,:)
      REAL(SZ),ALLOCATABLE    :: ZE_RHS(:,:,:),QE_RHS(:,:,:)
      REAL(SZ),ALLOCATABLE    :: DE_IN(:,:),DX_IN(:,:),DE_ED(:)
!.....Solution status variables
      INTEGER,ALLOCATABLE     :: WDFLG(:)
!.....Pressure Variables
      INTEGER                 :: PP_DEL,PP_NEGP
      INTEGER,ALLOCATABLE     :: DISPFLG(:)
      REAL(SZ),ALLOCATABLE    :: PD(:,:,:),PB(:,:,:)
      REAL(SZ),ALLOCATABLE    :: PP_XEGP(:),PP_WEGP(:),PP_MU(:)
      REAL(SZ),ALLOCATABLE    :: PP_PHI(:,:),PP_DPHI(:,:),PP_DDPHI(:,:),PP_WEI(:,:)
      REAL(SZ),ALLOCATABLE    :: PDLVL(:,:,:),PBLVL(:,:,:)
      INTEGER		      :: PPCNT
!.....Nodal Attributes
      REAL(SZ),ALLOCATABLE    :: MANN(:,:)
      REAL(SZ),ALLOCATABLE    :: SPNG_GEN(:,:),SPNG_ABS(:,:)
!.....Sponge Generation variables
      REAL(SZ),ALLOCATABLE    :: SPNG_ZAMP(:),SPNG_QAMP(:),SPNG_K(:),SPNG_SIG(:)
      REAL(SZ),ALLOCATABLE    :: SPNG_PHASE(:)
      CHARACTER(LEN=100)      :: SPONGE_TYPE
      INTEGER                 :: NUM_FREQ
      REAL(SZ)                :: SPNG_DIMP
!.....Misc. Variables
      REAL(SZ)                :: MAXRUNUP
!.....Eddy Viscosity Variables
      REAL(SZ),ALLOCATABLE    :: EDDY_T(:), EDDY_V(:)
      REAL(SZ),ALLOCATABLE    :: EDDY_SRC(:,:,:)
      INTEGER,ALLOCATABLE     :: EDDY_B(:)
      INTEGER                 :: FORT651
      real(sz),allocatable    :: kbrwaveoldR(:,:)
      integer,allocatable     :: kbrwaveoldI(:,:)
      integer                 :: numbrwavold

      END MODULE GLOBALS
