!==========================================================================!
!==========================================================================!
!..........................................................................!
      SUBROUTINE NONHYDRO
      
      USE READ_DGINP, ONLY : INONHYDRO
      
      IMPLICIT NONE
      
      SELECT CASE (INONHYDRO)
        CASE(2)
          CALL NONHYDRO_MU2_FD
        CASE(4)
          CALL NONHYDRO_MU4_FD
        CASE DEFAULT
          CALL NONHYDRO_MU0                
      END SELECT
!       
      RETURN
      END SUBROUTINE NONHYDRO
!..........................................................................!
!==========================================================================!
!..........................................................................!
      SUBROUTINE NONHYDRO_MU0
      
      USE GLOBALS, ONLY : IRK,PD,PB
      
      IMPLICIT NONE
      
      PD(:,:,IRK) = 0.D0
      PB(:,:,IRK) = 0.D0
      
      RETURN
      END SUBROUTINE NONHYDRO_MU0
!..........................................................................!
!==========================================================================!
!..........................................................................!
      SUBROUTINE NONHYDRO_MU2_FD
      
      USE READ_DGINP, ONLY : P
      USE GLOBALS,    ONLY : ZE,NE,LE,G,PD,PB,IRK,WDFLG,DISPFLG
      USE SIZES,      ONLY : SZ,C12
      
      IMPLICIT NONE
      INTEGER  :: L,I,INFO
      REAL(SZ) :: DE_IN(NE),UE_IN(NE),ZE_IN(NE),QE_IN(NE)
      REAL(SZ) :: DE_X_IN(NE),UE_X_IN(NE),ZE_X_IN(NE),QE_X_IN(NE)
      REAL(SZ) :: DE_XX_IN(NE),ZE_XX_IN(NE)
      REAL(SZ) :: MP1XX,MP1X,MP1,MRHS
      REAL(SZ) :: P1(NE),AU(NE-1),AL(NE-1),AD(NE)
      REAL(SZ) :: PDTMP(NE),PBTMP(NE)
      REAL(SZ) :: P1X,P2,PTMP(2)
      
      REAL(SZ),DIMENSION(NE,2)   :: alp2,alp1
      REAL(SZ),DIMENSION(NE,2,1) :: bet3,bet2,bet1
      REAL(SZ),DIMENSION(NE)     :: alp0
      REAL(SZ),DIMENSION(NE,1)   :: bet0
      
      
!.....Project solution onto a finite difference mesh
!.....First project solution and first derivatives
      CALL PROJ_FE_FD(DE_IN,DE_X_IN,DE_XX_IN,ZE_IN,ZE_X_IN, &
                            ZE_XX_IN,UE_IN,UE_X_IN,QE_IN,QE_X_IN)
!.....Now setup tri-diagonal vectors for pressure poisson problem
      DO L = 2,NE-1
        !DISPFLG = WDFLG(L)
!         CALL PRESPOIS_ORDER_2_COEFF(ZE_IN(L),ZE_X_IN(L),ZE_XX_IN(L),DE_IN(L), &
!                            DE_X_IN(L),DE_XX_IN(L),UE_IN(L),UE_X_IN(L), &
!                               MRHS,MP1,MP1X,MP1XX)
        CALL WASUPP_Coeffs_0(1,DE_IN(L),DE_X_IN(L),DE_XX_IN(L),ZE_IN(L),   &
     &          ZE_X_IN(L),ZE_XX_IN(L),UE_IN(L),UE_X_IN(L),alp0(l),bet0(l,1))
        CALL WASUPP_Coeffs_1(1,DE_IN(L),DE_X_IN(L),DE_XX_IN(L),ZE_IN(L),   &
     &          ZE_X_IN(L),ZE_XX_IN(L),alp2(l,1),alp1(l,1),bet3(l,1,1),    &
     &          bet2(l,1,1),bet1(l,1,1))
        CALL WASUPP_Coeffs_N(2,1,DE_IN(L),DE_X_IN(L),DE_XX_IN(L),ZE_IN(L), &
     &          ZE_X_IN(L),ZE_XX_IN(L),alp2(l,2),alp1(l,2),bet3(l,2,1),    &
     &          bet2(l,2,1),bet1(l,2,1))
        
        MP1XX = bet3(l,1,1)
        MP1X  = (bet2(l,1,1)-alp2(l,1)/alp1(l,2)*bet1(l,2,1))
        MP1   = (bet1(l,1,1)-alp1(l,1)/alp1(l,2)*bet1(l,2,1))
        MRHS  = (bet0(l,1)-alp0(l)/alp1(l,2)*bet1(l,2,1))
        
        AL(L-1) = (MP1XX/LE(L)/LE(L)-MP1X/2.D0/LE(L))*(DISPFLG(L))
        AD(L)   = (-2.D0*MP1XX/LE(L)/LE(L)+MP1)      *(DISPFLG(L)) + REAL(1-DISPFLG(L))
        AU(L)   = (MP1XX/LE(L)/LE(L)+MP1X/2.D0/LE(L))*(DISPFLG(L))
        P1(L)   =  MRHS                              *(DISPFLG(L)) + G*ZE_IN(L)*(1-DISPFLG(L))
      END DO
      L = 1
      !DISPFLG = WDFLG(L)
!       CALL PRESPOIS_ORDER_2_COEFF(ZE_IN(L),ZE_X_IN(L),ZE_XX_IN(L),DE_IN(L), &
!                            DE_X_IN(L),DE_XX_IN(L),UE_IN(L),UE_X_IN(L), &
!                               MRHS,MP1,MP1X,MP1XX)
      CALL WASUPP_Coeffs_0(1,DE_IN(L),DE_X_IN(L),DE_XX_IN(L),ZE_IN(L),ZE_X_IN(L),ZE_XX_IN(L),UE_IN(L),UE_X_IN(L),alp0(l),bet0(l,1))
      CALL WASUPP_Coeffs_1(1,DE_IN(L),DE_X_IN(L),DE_XX_IN(L),ZE_IN(L),ZE_X_IN(L),ZE_XX_IN(L),alp2(l,1),alp1(l,1),bet3(l,1,1),bet2(l,1,1),bet1(l,1,1))
      CALL WASUPP_Coeffs_N(2,1,DE_IN(L),DE_X_IN(L),DE_XX_IN(L),ZE_IN(L),ZE_X_IN(L),ZE_XX_IN(L),alp2(l,2),alp1(l,2),bet3(l,2,1),bet2(l,2,1),bet1(l,2,1))
      
      MP1XX = bet3(l,1,1)
      MP1X  = (bet2(l,1,1)-alp2(l,1)/alp1(l,2)*bet1(l,2,1))
      MP1   = (bet1(l,1,1)-alp1(l,1)/alp1(l,2)*bet1(l,2,1))
      MRHS  = (bet0(l,1)-alp0(l)/alp1(l,2)*bet1(l,2,1))
      
      AD(L)   = (-MP1XX/LE(L)/LE(L)-MP1X/2.D0/LE(L)+MP1)*(DISPFLG(L)) + REAL(1-DISPFLG(L))
      AU(L)   = ( MP1XX/LE(L)/LE(L)+MP1X/2.D0/LE(L))    *(DISPFLG(L))
      P1(L)   =   MRHS                                  *(DISPFLG(L)) + G*ZE_IN(L)*(1-DISPFLG(L))
      
      L = NE
      !DISPFLG = WDFLG(L)
!       CALL PRESPOIS_ORDER_2_COEFF(ZE_IN(L),ZE_X_IN(L),ZE_XX_IN(L),DE_IN(L), &
!                            DE_X_IN(L),DE_XX_IN(L),UE_IN(L),UE_X_IN(L), &
!                               MRHS,MP1,MP1X,MP1XX)
      CALL WASUPP_Coeffs_0(1,DE_IN(L),DE_X_IN(L),DE_XX_IN(L),ZE_IN(L),ZE_X_IN(L),ZE_XX_IN(L),UE_IN(L),UE_X_IN(L),alp0(l),bet0(l,1))
      CALL WASUPP_Coeffs_1(1,DE_IN(L),DE_X_IN(L),DE_XX_IN(L),ZE_IN(L),ZE_X_IN(L),ZE_XX_IN(L),alp2(l,1),alp1(l,1),bet3(l,1,1),bet2(l,1,1),bet1(l,1,1))
      CALL WASUPP_Coeffs_N(2,1,DE_IN(L),DE_X_IN(L),DE_XX_IN(L),ZE_IN(L),ZE_X_IN(L),ZE_XX_IN(L),alp2(l,2),alp1(l,2),bet3(l,2,1),bet2(l,2,1),bet1(l,2,1))
      
      MP1XX = bet3(l,1,1)
      MP1X  = (bet2(l,1,1)-alp2(l,1)/alp1(l,2)*bet1(l,2,1))
      MP1   = (bet1(l,1,1)-alp1(l,1)/alp1(l,2)*bet1(l,2,1))
      MRHS  = (bet0(l,1)-alp0(l)/alp1(l,2)*bet1(l,2,1))
      
      AL(L-1) = ( MP1XX/LE(L)/LE(L)-MP1X/2.D0/LE(L))    *(DISPFLG(L))
      AD(L)   = (-MP1XX/LE(L)/LE(L)+MP1X/2.D0/LE(L)+MP1)*(DISPFLG(L)) + REAL(1-DISPFLG(L))
      P1(L)   =   MRHS                                  *(DISPFLG(L)) + G*ZE_IN(l)*(1-DISPFLG(L))
!.....Solve for the pressure profile
      CALL DGTSV(NE,1,AL,AD,AU,P1,NE,INFO)
      
!.....Use the solution for P1 to find P2 and subsequently Pd and Pb
      L = 1
      !DISPFLG = WDFLG(L)
      P1X      = (P1(2)-P1(1))/2.D0/LE(L)
!       CALL PRESPOIS_ORDER_2_P2(ZE_IN(L),DE_IN(L),DE_X_IN(L),DE_XX_IN(L), &
!                      UE_IN(L),P1(L),P1X,0,P2,PDtmp(L),PBtmp(L))
      P2    = 1.d0/alp1(l,2)*(alp0(l)-alp2(l,1)*P1X-alp1(l,1)*P1(L))*(DISPFLG(L))
      PDtmp(L) = C12*(ZE_IN(L)+DE_IN(L))*(P1(L)-G*ZE_IN(L))
      PBtmp(L) = P1(L)-1.D0/3.D0*P2-G*ZE_IN(L)
      DO L = 2,NE-1
        !DISPFLG = WDFLG(L)
        P1X      = (P1(L+1)-P1(L-1))/2.D0/LE(L)
!         CALL PRESPOIS_ORDER_2_P2(ZE_IN(L),DE_IN(L),DE_X_IN(L),DE_XX_IN(L), &
!                      UE_IN(L),P1(L),P1X,0,P2,PDtmp(L),PBtmp(L))
        P2    = 1.d0/alp1(l,2)*(alp0(l)-alp2(l,1)*P1X-alp1(l,1)*P1(L))*(DISPFLG(L))
        PDtmp(L) = C12*(ZE_IN(L)+DE_IN(L))*(P1(L)-G*ZE_IN(L))
        PBtmp(L) = P1(L)-1.D0/3.D0*P2-G*ZE_IN(L)
      END DO
      L = NE
      !DISPFLG = WDFLG(L)
      P1X      = (P1(NE)-P1(NE-1))/2.D0/LE(L)
!       CALL PRESPOIS_ORDER_2_P2(ZE_IN(L),DE_IN(L),DE_X_IN(L),DE_XX_IN(L), &
!                      UE_IN(L),P1(L),P1X,0,P2,PDtmp(L),PBtmp(L))
      P2    = 1.d0/alp1(l,2)*(alp0(l)-alp2(l,1)*P1X-alp1(l,1)*P1(L))*(DISPFLG(L))
      PDtmp(L) = C12*(ZE_IN(L)+DE_IN(L))*(P1(L)-G*ZE_IN(L))
      PBtmp(L) = P1(L)-1.D0/3.D0*P2-G*ZE_IN(L)
      
      DO I = 1,P+1
        L = 1
        PTMP(1) = C12*(PDTMP(L)+PDTMP(L))
        PTMP(2) = C12*(PDTMP(L+1)+PDTMP(L))          
        PD(I,L,IRK) = PTMP(1) + (PTMP(2)-PTMP(1))/2.D0*                   &
     &                                ( 2.D0*REAL(I-1)/REAL(P) )
        PTMP(1) = C12*(PBTMP(L)+PBTMP(L))
        PTMP(2) = C12*(PBTMP(L+1)+PBTMP(L))          
        PB(I,L,IRK) = PTMP(1) + (PTMP(2)-PTMP(1))/2.D0*                   &
     &                                ( 2.D0*REAL(I-1)/REAL(P) )
        DO L = 2,NE-1
          PTMP(1) = C12*(PDTMP(L-1)+PDTMP(L))
          PTMP(2) = C12*(PDTMP(L+1)+PDTMP(L))          
          PD(I,L,IRK) = PTMP(1) + (PTMP(2)-PTMP(1))/2.D0*                 &
     &                                ( 2.D0*REAL(I-1)/REAL(P) )
          PTMP(1) = C12*(PBTMP(L-1)+PBTMP(L))
          PTMP(2) = C12*(PBTMP(L+1)+PBTMP(L))          
          PB(I,L,IRK) = PTMP(1) + (PTMP(2)-PTMP(1))/2.D0*                 &
     &                                ( 2.D0*REAL(I-1)/REAL(P) )
        END DO
        L = NE
        PTMP(1) = C12*(PDTMP(L-1)+PDTMP(L))
        PTMP(2) = C12*(PDTMP(L)+PDTMP(L))          
        PD(I,L,IRK) = PTMP(1) + (PTMP(2)-PTMP(1))/2.D0*                   &
     &                                ( 2.D0*REAL(I-1)/REAL(P) )
        PTMP(1) = C12*(PBTMP(L-1)+PBTMP(L))
        PTMP(2) = C12*(PBTMP(L)+PBTMP(L))          
        PB(I,L,IRK) = PTMP(1) + (PTMP(2)-PTMP(1))/2.D0*                   &
     &                                ( 2.D0*REAL(I-1)/REAL(P) )
      END DO
      
      RETURN
      END SUBROUTINE NONHYDRO_MU2_FD
!..........................................................................!
!==========================================================================!
!..........................................................................!
      SUBROUTINE NONHYDRO_MU2_CG
      
      USE GLOBALS,    ONLY : IRK,PD,PB,ZE,QE,DE_IN,DX_IN,PHI,DPHI,LE,PSI,  &
     &                       DPSI,WEGP,LE,NE,G,CGPIV,MCG
      USE READ_DGINP, ONLY : P,NEGP
      USE SIZES,      ONLY : SZ,C12
      
      IMPLICIT NONE
      INTEGER  :: L,I,K,J
      REAL(SZ) :: ZE_IN,QE_IN,ZX_IN,QX_IN
      REAL(SZ) :: UE_IN,UX_IN,HE_IN,HX_IN
      REAL(SZ) :: P1_IN,P1X_IN,P2_IN
      REAL(SZ) :: P1COEFF,P1XCOEFF
      REAL(SZ) :: MLOC(P+1,P+1),B1LOC(P+1),B2LOC(P+1),BDLOC(P+1),BBLOC(P+1)
      REAL(SZ) :: AB(3*P+1,NE*P+1)
      REAL(SZ) :: P1LOC(NE*P+1),P2LOC(NE*P+1),PDLOC(NE*P+1),PBLOC(NE*P+1)
      REAL(SZ) :: P1(P+1,NE),P2(P+1,NE)
      INTEGER  :: IPIV(NE*P+1),INFO
      
      PD(:,:,IRK) = 0.D0
      PB(:,:,IRK) = 0.D0
      AB(:,:)     = 0.D0
      P1LOC(:)    = 0.D0
      P2LOC(:)    = 0.D0
      PDLOC(:)    = 0.D0
      PBLOC(:)    = 0.D0
      P1(:,:)     = 0.D0
      P2(:,:)     = 0.D0
      
!.....Build LHS matrix for CG-FEM PP problem
      DO L = 1,NE
        MLOC(:,:) = 0.D0
        B1LOC(:)  = 0.D0
        B2LOC(:)  = 0.D0
        BDLOC(:)  = 0.D0
        BBLOC(:)  = 0.D0
        DO K = 1,NEGP
          ! Build local solutions for velocity and free-surface
          ZE_IN = 0.D0
          QE_IN = 0.D0
          ZX_IN = 0.D0
          QX_IN = 0.D0
          DO I = 1,P+1
            ZE_IN = ZE_IN + ZE(I,L,IRK)*PHI(I,K)
            QE_IN = QE_IN + QE(I,L,IRK)*PHI(I,K)
            
            ZX_IN = ZX_IN + ZE(I,L,IRK)*DPHI(I,K)/(C12*LE(L))
            QX_IN = QX_IN + QE(I,L,IRK)*DPHI(I,K)/(C12*LE(L))
          END DO
          HE_IN = DE_IN(L,K)+ZE_IN
          HX_IN = DX_IN(L,K)+ZX_IN
          UE_IN = QE_IN/HE_IN
          UX_IN = (HE_IN*QX_IN-QE_IN*HX_IN)/HE_IN**2
          
          ! Build local matrix for P solver          
          DO J = 1,P+1
            P1COEFF  = 1.D0/12.D0*( 3.D0*PSI(J,K)*(5.D0+4.D0*ZX_IN**2-     &
     &                     2.D0*ZX_IN*DX_IN(L,K)-DX_IN(L,K)**2)+           &
     &                     2.D0*HE_IN*DPSI(J,K)/(C12*LE(L))*(2.D0*ZX_IN-   &
     &                     3.D0*DX_IN(L,K)) )
            P1XCOEFF = 1.D0/12.D0*HE_IN*( PSI(J,K)*(8.D0*ZX_IN+            &
     &                  3.D0*DX_IN(L,K))+6.D0*DPSI(J,K)/(LE(L)*C12)*HE_IN )
     
!             P1COEFF = PSI(J,K);
!             P1XCOEFF = DPSI(J,K)/(C12*LE(L))
            DO I = 1,P+1
              MLOC(J,I) = MLOC(J,I) + ( P1COEFF*PSI(I,K)+                  &
     &            P1XCOEFF*DPSI(I,K)/(C12*LE(L)) )*WEGP(K)*LE(L)*C12
            END DO
            B1LOC(J) = B1LOC(J) - ( -5.D0/12.D0*(6.D0*HE_IN*DX_IN(L,K)*    &
     &                       UE_IN*UX_IN*PSI(J,K)+4.D0*HE_IN**2*UX_IN**2*  &
     &                       PSI(J,K)+3.D0*UE_IN**2*DX_IN(L,K)*(HX_IN*     &
     &                       PSI(J,K)+HE_IN*DPSI(J,K)/(C12*LE(L))))        &
     &                            +1.D0/12.D0*G*(6.D0*ZE_IN**2*DX_IN(L,K)* &
     &                       DPSI(L,K)/(C12*LE(L))+2.D0*DE_IN(L,K)*ZX_IN*( &
     &                       PSI(J,K)*(6.D0*ZX_IN+DX_IN(L,K))+             &
     &                       DPSI(J,K)/(C12*LE(L))*2.D0*DE_IN(L,K))+ZE_IN* &
     &                       (PSI(J,K)*(-15.D0+8.D0*ZX_IN*DX_IN(L,K)+3.D0* &
     &                       DX_IN(L,K)**2)+2.D0*DE_IN(L,K)*(2.D0*ZX_IN+   &
     &                       3.D0*DX_IN(L,K))*DPSI(J,K)/(C12*LE(L)))) )*   &
     &                       WEGP(K)*LE(L)*C12
            B2LOC(J) = B2LOC(J) + 3.D0/4.D0*( HE_IN*UE_IN**2*DPSI(J,K)/    &
     &                       (C12*LE(L))-3.D0/4.D0*PSI(J,K)*UE_IN*(HX_IN*  &
     &                        UE_IN+2.D0*HE_IN*UX_IN)-                     &
     &                        G*(1.D0+DX_IN(L,K)**2)*ZE_IN*PSI(J,K) )*     &
     &                        WEGP(K)*LE(L)*C12
           BDLOC(J) = BDLOC(J) - C12*G*HE_IN*ZE_IN*PSI(J,K)*WEGP(K)*LE(L)*C12
           BBLOC(J) = BBLOC(J) - G*ZE_IN*PSI(J,K)*WEGP(K)*LE(L)*C12
!             BLOC(J) = BLOC(J) + ZE_IN*WEGP(K)*LE(L)*C12*PSI(J,K)
             
          END DO
        END DO
        ! Add local matrix to global matrix
        DO J = 1,P+1
          DO I = 1,P+1
            AB(2*P+1+J-I,(L-1)*P+I) = AB(2*P+1+J-I,(L-1)*P+I) + MLOC(J,I)
          END DO
          P1LOC((L-1)*P+J) = P1LOC((L-1)*P+J) + B1LOC(J)
          P2LOC((L-1)*P+J) = P2LOC((L-1)*P+J) + B2LOC(J)
          PDLOC((L-1)*P+J) = PDLOC((L-1)*P+J) + BDLOC(J)
          PBLOC((L-1)*P+J) = PBLOC((L-1)*P+J) + BBLOC(J)
        END DO
      END DO
!.....Solve for P1
      CALL DGBSV( NE*P+1,P,P,1,AB,3*P+1,IPIV,P1LOC,NE*P+1,INFO )
!.....Reconstruct P1 FEM solution from matrix solve
      DO L = 1,NE
        DO I = 1,P+1
          P1(I,L) = P1LOC((L-1)*P+I)
        END DO
      END DO
!.....Now determine P2, Pd and Pb
      DO L = 1,NE
        B2LOC(:) = 0.D0
        BDLOC(:) = 0.D0
        BBLOC(:) = 0.D0
        DO K = 1,NEGP
          ! Build local solutions for velocity and free-surface
          ZE_IN  = 0.D0
          ZX_IN  = 0.D0
          P1_IN  = 0.D0
          P1X_IN = 0.D0
          DO I = 1,P+1
            ZE_IN = ZE_IN + ZE(I,L,IRK)*PHI(I,K)            
            ZX_IN = ZX_IN + ZE(I,L,IRK)*DPHI(I,K)/(C12*LE(L))
            
            P1_IN  = P1_IN  + P1(I,L)*PSI(I,K)            
            P1X_IN = P1X_IN + P1(I,L)*DPSI(I,K)/(C12*LE(L))
          END DO
          HE_IN = ZE_IN + DE_IN(L,K)
          ! Add P1 contribution to P2, Pd and Pb solution vectors
          DO J = 1,P+1
            B2LOC(J) = B2LOC(J) + 3.D0/4.D0*( (1.D0+DX_IN(L,K)**2)*P1_IN   &
     &                -HE_IN*DX_IN(L,K)*P1X_IN )*PSI(J,K)*WEGP(K)*C12*LE(L)
            BDLOC(J) = BDLOC(J) + C12*HE_IN*P1_IN*PSI(J,K)*WEGP(K)*C12*LE(L)
            BBLOC(J) = BBLOC(J) + P1_IN*PSI(J,K)*WEGP(K)*C12*LE(L)
          END DO
        END DO
        ! Add local matrix to global matrix
        DO J = 1,P+1
          P2LOC((L-1)*P+J) = P2LOC((L-1)*P+J) + B2LOC(J)
          PDLOC((L-1)*P+J) = PDLOC((L-1)*P+J) + BDLOC(J)
          PBLOC((L-1)*P+J) = PBLOC((L-1)*P+J) + BBLOC(J)
        END DO
      END DO
!.....Solve for P2 and Pd
      CALL DGBTRS('No transpose',P*NE+1,P,P,1,MCG,3*P+1,CGPIV,P2LOC,P*NE+1,INFO)
      CALL DGBTRS('No transpose',P*NE+1,P,P,1,MCG,3*P+1,CGPIV,PDLOC,P*NE+1,INFO)
!.....Reconstruct P2 AND Pd FEM solution from matrix solve
      DO L = 1,NE
        DO I = 1,P+1
          P2(I,L)     = P2LOC((L-1)*P+I)
          PD(I,L,IRK) = PDLOC((L-1)*P+I)
        END DO
      END DO
!.....Finally add the contribution of P2 to get Pb
      DO L = 1,NE
        BBLOC(:) = 0.D0
        DO K = 1,NEGP
          P2_IN = 0.D0
          DO I = 1,P+1
            P2_IN = P2_IN + P2(I,L)*PSI(I,K)
          END DO
          DO J = 1,P+1
            BBLOC(J) = BBLOC(J) - 1.D0/3.D0*P2_IN*WEGP(K)*PSI(J,K)*C12*LE(L)
          END DO
        END DO
        DO J = 1,P+1
          PBLOC((L-1)*P+J) = PBLOC((L-1)*P+J) + BBLOC(J)
        END DO
      END DO
!.....Solve for Pb
      CALL DGBTRS('No transpose',P*NE+1,P,P,1,MCG,3*P+1,CGPIV,PBLOC,P*NE+1,INFO)
!.....Reconstruct Pb FEM solution from matrix solve
      DO L = 1,NE
        DO I = 1,P+1
          PB(I,L,IRK) = PBLOC((L-1)*P+I)
        END DO
      END DO      
      
!       CALL WRITE_73
!       STOP
      
      RETURN
      END SUBROUTINE NONHYDRO_MU2_CG
!..........................................................................!
!==========================================================================!
!..........................................................................!
      SUBROUTINE NONHYDRO_MU4_FD
      
      USE READ_DGINP, ONLY : P
      USE GLOBALS, ONLY    : IRK,PD,PB,NE,G,PP_NEGP,PP_XEGP,PP_WEGP,       &
     &                       PP_PHI,LE,WDFLG,DISPFLG
      USE SIZES,   ONLY    : SZ,C12
      
      IMPLICIT NONE
      INTEGER  :: L,I,J,M,N,K
      INTEGER  :: INFO
      REAL(SZ) :: DE_IN(NE),UE_IN(NE),ZE_IN(NE),QE_IN(NE)
      REAL(SZ) :: DE_X_IN(NE),UE_X_IN(NE),ZE_X_IN(NE),QE_X_IN(NE)
      REAL(SZ) :: DE_XX_IN(NE),ZE_XX_IN(NE)
      REAL(SZ) :: DX
      REAL(SZ) :: PDTMP(NE),PBTMP(NE)
      REAL(SZ) :: P1xx,P2xx,P1x,P2x,P3,P4,q
      REAL(SZ) :: PTMP(2)
      
      REAL(SZ),DIMENSION(NE,4)    :: alp2,alp1
      REAL(SZ),DIMENSION(NE,4,3)  :: bet3,bet2,bet1
      REAL(SZ),DIMENSION(NE)      :: alp0
      REAL(SZ),DIMENSION(NE,3)    :: bet0
      
      REAL(SZ),DIMENSION(10,2*NE) :: AB4
      REAL(SZ),DIMENSION(2*NE)    :: B4
      INTEGER,DIMENSION(2*NE)     :: IPIV4
      REAL(SZ),DIMENSION(NE,4,3)  :: bet3hat,bet2hat,bet1hat
      REAL(SZ),DIMENSION(NE,3)    :: bet0hat
      REAL(SZ),DIMENSION(NE,4,3)  :: bet3star,bet2star,bet1star
      REAL(SZ),DIMENSION(NE,3)    :: bet0star
      REAL(SZ),DIMENSION(NE)      :: P1_mu4,P2_mu4
      
      AB4(:,:) = 0.d0
      B4(:)    = 0.d0
      
!.....Project solution onto a finite difference mesh
!.....First project solution and first derivatives
      CALL PROJ_FE_FD(DE_IN,DE_X_IN,DE_XX_IN,ZE_IN,ZE_X_IN, &
                            ZE_XX_IN,UE_IN,UE_X_IN,QE_IN,QE_X_IN)
!.....Now setup tri-diagonal vectors for pressure poisson problem
      DO L = 2,NE-1
        !DISPFLG = WDFLG(L)
        DX = LE(L)
        DO M = 1,3
          CALL WASUPP_Coeffs_0(m,DE_IN(L),DE_X_IN(L),DE_XX_IN(L),ZE_IN(L), &
     &        ZE_X_IN(L),ZE_XX_IN(L),UE_IN(L),UE_X_IN(L),alp0(l),bet0(l,m))
          CALL WASUPP_Coeffs_1(m,DE_IN(L),DE_X_IN(L),DE_XX_IN(L),ZE_IN(L), &
     &        ZE_X_IN(L),ZE_XX_IN(L),alp2(l,1),alp1(l,1),bet3(l,1,m),      &
     &        bet2(l,1,m),bet1(l,1,m))
          DO N = 2,4
            CALL WASUPP_Coeffs_N(n,m,DE_IN(L),DE_X_IN(L),DE_XX_IN(L),      &
     &        ZE_IN(L),ZE_X_IN(L),ZE_XX_IN(L),alp2(l,n),alp1(l,n),         &
     &        bet3(l,n,m),bet2(l,n,m),bet1(l,n,m))
          END DO
        END DO
        
        DO M = 1,3
        DO N = 1,3      
          bet3hat(l,N,M) = bet3(l,n,m)
          bet2hat(l,N,M) = bet2(l,n,m)-bet1(l,4,m)/alp1(l,4)*alp2(l,n)
          bet1hat(l,N,M) = bet1(l,n,m)-bet1(l,4,m)/alp1(l,4)*alp1(l,n)
        END DO
          bet0hat(l,M)   = bet0(l,m)  -bet1(l,4,m)/alp1(l,4)*alp0(l)      
        END DO
      
        DO M = 2,3
        DO N = 1,2
          bet3star(l,n,m) = bet3hat(l,n,m)-bet1hat(l,3,m)/bet1hat(l,3,1)*  &
     &                      bet3hat(l,n,1)
          bet2star(l,n,m) = bet2hat(l,n,m)-bet1hat(l,3,m)/bet1hat(l,3,1)*  &
     &                      bet2hat(l,n,1)
          bet1star(l,n,m) = bet1hat(l,n,m)-bet1hat(l,3,m)/bet1hat(l,3,1)*  &
     &                      bet1hat(l,n,1)
        END DO
          bet0star(l,m)   = bet0hat(l,m)  -bet1hat(l,3,m)/bet1hat(l,3,1)*  &
     &                      bet0hat(l,1)
        END DO
        
        ! "A matrix"
        i = 2*l-1
        n = l-1        
        j = 2*n-1
        AB4(7+i-j,j) = (bet3star(l,1,2)/dx/dx-bet2star(l,1,2)/2.d0/dx)*(DISPFLG(L))
        n = l
        j = 2*n-1
        AB4(7+i-j,j) = (-2.d0*bet3star(l,1,2)/dx/dx+bet1star(l,1,2))  *(DISPFLG(L)) + REAL(1-DISPFLG(L))
        n = l+1
        j = 2*n-1
        AB4(7+i-j,j) = (bet3star(l,1,2)/dx/dx+bet2star(l,1,2)/2.d0/dx)*(DISPFLG(L))
        
        ! "B matrix"
        i = 2*l-1
        n = l-1        
        j = 2*n
        AB4(7+i-j,j) = (bet3star(l,2,2)/dx/dx-bet2star(l,2,2)/2.d0/dx)*(DISPFLG(L))
        n = l
        j = 2*n
        AB4(7+i-j,j) = (-2.d0*bet3star(l,2,2)/dx/dx+bet1star(l,2,2))  *(DISPFLG(L))
        n = l+1
        j = 2*n
        AB4(7+i-j,j) = (bet3star(l,2,2)/dx/dx+bet2star(l,2,2)/2.d0/dx)*(DISPFLG(L))
        
        ! "C matrix"
        i = 2*l
        n = l-1        
        j = 2*n-1
        AB4(7+i-j,j) = (bet3star(l,1,3)/dx/dx-bet2star(l,1,3)/2.d0/dx)*(DISPFLG(L))
        n = l
        j = 2*n-1
        AB4(7+i-j,j) = (-2.d0*bet3star(l,1,3)/dx/dx+bet1star(l,1,3))  *(DISPFLG(L))
        n = l+1
        j = 2*n-1
        AB4(7+i-j,j) = (bet3star(l,1,3)/dx/dx+bet2star(l,1,3)/2.d0/dx)*(DISPFLG(L))
        
        ! "D matrix"
        i = 2*l
        n = l-1        
        j = 2*n
        AB4(7+i-j,j) = (bet3star(l,2,3)/dx/dx-bet2star(l,2,3)/2.d0/dx)*(DISPFLG(L))
        n = l
        j = 2*n
        AB4(7+i-j,j) = (-2.d0*bet3star(l,2,3)/dx/dx+bet1star(l,2,3))  *(DISPFLG(L)) + REAL(1-DISPFLG(L))
        n = l+1
        j = 2*n
        AB4(7+i-j,j) = (bet3star(l,2,3)/dx/dx+bet2star(l,2,3)/2.d0/dx)*(DISPFLG(L))
        
        B4(2*l-1) = bet0star(l,2)*(DISPFLG(L)) + G*ZE_IN(L)*(1-DISPFLG(L))
        B4(2*l)   = bet0star(l,3)*(DISPFLG(L))

      END DO
      L = 1
      !DISPFLG = WDFLG(L)
      DX = LE(L)
      DO M = 1,3
        CALL WASUPP_Coeffs_0(m,DE_IN(L),DE_X_IN(L),DE_XX_IN(L),ZE_IN(L),   &
     &        ZE_X_IN(L),ZE_XX_IN(L),UE_IN(L),UE_X_IN(L),alp0(l),bet0(l,m))
        CALL WASUPP_Coeffs_1(m,DE_IN(L),DE_X_IN(L),DE_XX_IN(L),ZE_IN(L),   &
     &        ZE_X_IN(L),ZE_XX_IN(L),alp2(l,1),alp1(l,1),bet3(l,1,m),      &
     &        bet2(l,1,m),bet1(l,1,m))
        DO N = 2,4
          CALL WASUPP_Coeffs_N(n,m,DE_IN(L),DE_X_IN(L),DE_XX_IN(L),        &
     &          ZE_IN(L),ZE_X_IN(L),ZE_XX_IN(L),alp2(l,n),alp1(l,n),       &
     &          bet3(l,n,m),bet2(l,n,m),bet1(l,n,m))
        END DO
      END DO
        
      DO M = 1,3
      DO N = 1,3      
        bet3hat(l,N,M) = bet3(l,n,m)
        bet2hat(l,N,M) = bet2(l,n,m)-bet1(l,4,m)/alp1(l,4)*alp2(l,n)
        bet1hat(l,N,M) = bet1(l,n,m)-bet1(l,4,m)/alp1(l,4)*alp1(l,n)
      END DO
        bet0hat(l,M)   = bet0(l,m)  -bet1(l,4,m)/alp1(l,4)*alp0(l)      
      END DO
      
      DO M = 2,3
      DO N = 1,2
        bet3star(l,n,m) = bet3hat(l,n,m)-bet1hat(l,3,m)/bet1hat(l,3,1)*    &
     &                    bet3hat(l,n,1)
        bet2star(l,n,m) = bet2hat(l,n,m)-bet1hat(l,3,m)/bet1hat(l,3,1)*    &
     &                    bet2hat(l,n,1)
        bet1star(l,n,m) = bet1hat(l,n,m)-bet1hat(l,3,m)/bet1hat(l,3,1)*    &
     &                    bet1hat(l,n,1)
      END DO
        bet0star(l,m)   = bet0hat(l,m)  -bet1hat(l,3,m)/bet1hat(l,3,1)*    &
     &                    bet0hat(l,1)
      END DO
      ! "A matrix"
      i = 2*l-1
      n = l        
      j = 2*n-1
      AB4(7+i-j,j) = (-2.d0*bet3star(l,1,2)/dx/dx+bet1star(l,1,2))*(DISPFLG(L)) + REAL(1-DISPFLG(L))
      n = l+1
      j = 2*n-1
      AB4(7+i-j,j) = (2.d0*bet3star(l,1,2)/dx/dx)*(DISPFLG(L))
        
      ! "B matrix"
      i = 2*l-1
      n = l        
      j = 2*n
      AB4(7+i-j,j) = (-2.d0*bet3star(l,2,2)/dx/dx+bet1star(l,2,2))*(DISPFLG(L))
      n = l+1
      j = 2*n
      AB4(7+i-j,j) = (2.d0*bet3star(l,2,2)/dx/dx)*(DISPFLG(L))
        
      ! "C matrix"
      i = 2*l
      n = l        
      j = 2*n-1
      AB4(7+i-j,j) = (-2.d0*bet3star(l,1,3)/dx/dx+bet1star(l,1,3))*(DISPFLG(L))
      n = l+1
      j = 2*n-1
      AB4(7+i-j,j) = (2.d0*bet3star(l,1,3)/dx/dx)*(DISPFLG(L))
        
      ! "D matrix"
      i = 2*l
      n = l
      j = 2*n
      AB4(7+i-j,j) = (-2.d0*bet3star(l,2,3)/dx/dx+bet1star(l,2,3))*(DISPFLG(L))
      n = l+1
      j = 2*n
      AB4(7+i-j,j) = (2.d0*bet3star(l,2,3)/dx/dx)*(DISPFLG(L)) + REAL(1-DISPFLG(L))
      
      B4(2*l-1) = bet0star(l,2)*(DISPFLG(L)) + G*ZE_IN(L)*(1-DISPFLG(L))
      B4(2*l)   = bet0star(l,3)*(DISPFLG(L))
      
      L = NE
      !DISPFLG = WDFLG(L)
      DX = LE(L)
      DO M = 1,3
        CALL WASUPP_Coeffs_0(m,DE_IN(L),DE_X_IN(L),DE_XX_IN(L),ZE_IN(L),   &
     &         ZE_X_IN(L),ZE_XX_IN(L),UE_IN(L),UE_X_IN(L),alp0(l),bet0(l,m))
        CALL WASUPP_Coeffs_1(m,DE_IN(L),DE_X_IN(L),DE_XX_IN(L),ZE_IN(L),   &
     &         ZE_X_IN(L),ZE_XX_IN(L),alp2(l,1),alp1(l,1),bet3(l,1,m),     &
     &         bet2(l,1,m),bet1(l,1,m))
        DO N = 2,4
          CALL WASUPP_Coeffs_N(n,m,DE_IN(L),DE_X_IN(L),DE_XX_IN(L),        &
     &           ZE_IN(L),ZE_X_IN(L),ZE_XX_IN(L),alp2(l,n),alp1(l,n),      &
     &           bet3(l,n,m),bet2(l,n,m),bet1(l,n,m))
        END DO
      END DO
        
      DO M = 1,3
      DO N = 1,3      
        bet3hat(l,N,M) = bet3(l,n,m)
        bet2hat(l,N,M) = bet2(l,n,m)-bet1(l,4,m)/alp1(l,4)*alp2(l,n)
        bet1hat(l,N,M) = bet1(l,n,m)-bet1(l,4,m)/alp1(l,4)*alp1(l,n)
      END DO
        bet0hat(l,M)   = bet0(l,m)  -bet1(l,4,m)/alp1(l,4)*alp0(l)      
      END DO
      
      DO M = 2,3
      DO N = 1,2
        bet3star(l,n,m) = bet3hat(l,n,m)-bet1hat(l,3,m)/bet1hat(l,3,1)*    &
     &                    bet3hat(l,n,1)
        bet2star(l,n,m) = bet2hat(l,n,m)-bet1hat(l,3,m)/bet1hat(l,3,1)*    &
     &                    bet2hat(l,n,1)
        bet1star(l,n,m) = bet1hat(l,n,m)-bet1hat(l,3,m)/bet1hat(l,3,1)*    &
     &                    bet1hat(l,n,1)
      END DO
        bet0star(l,m)   = bet0hat(l,m)  -bet1hat(l,3,m)/bet1hat(l,3,1)*    &
     &                    bet0hat(l,1)
      END DO
      ! "A matrix"
      i = 2*l-1
      n = l-1        
      j = 2*n-1
      AB4(7+i-j,j) = (2.d0*bet3star(l,1,2)/dx/dx)*(DISPFLG(L))
      n = l
      j = 2*n-1
      AB4(7+i-j,j) = (-2.d0*bet3star(l,1,2)/dx/dx+bet1star(l,1,2))*(DISPFLG(L)) + REAL(1-DISPFLG(L))
        
      ! "B matrix"
      i = 2*l-1
      n = l-1        
      j = 2*n
      AB4(7+i-j,j) = (2.d0*bet3star(l,2,2)/dx/dx)*(DISPFLG(L))
      n = l
      j = 2*n
      AB4(7+i-j,j) = (-2.d0*bet3star(l,2,2)/dx/dx+bet1star(l,2,2))*(DISPFLG(L))
        
      ! "C matrix"
      i = 2*l
      n = l-1        
      j = 2*n-1
      AB4(7+i-j,j) = (2.d0*bet3star(l,1,3)/dx/dx)*(DISPFLG(L))
      n = l
      j = 2*n-1
      AB4(7+i-j,j) = (-2.d0*bet3star(l,1,3)/dx/dx+bet1star(l,1,3))*(DISPFLG(L))
        
      ! "D matrix"
      i = 2*n
      n = l-1        
      j = 2*n
      AB4(7+i-j,j) = (2.d0*bet3star(l,2,3)/dx/dx)*(DISPFLG(L))
      n = l
      j = 2*n
      AB4(7+i-j,j) = (-2.d0*bet3star(l,2,3)/dx/dx+bet1star(l,2,3))*(DISPFLG(L)) + REAL(1-DISPFLG(L))
      
      B4(2*l-1) = bet0star(l,2)*(DISPFLG(L)) + G*ZE_IN(L)*(1-DISPFLG(L))
      B4(2*l)   = bet0star(l,3)*(DISPFLG(L))
!.....Solve for the pressure profile
      CALL DGBSV(2*NE,3,3,1,AB4,10,IPIV4,B4,2*NE,INFO)
!.....Use the solution for P1 to find P2 and subsequently Pd and Pb
      do l = 1,ne
        p1_mu4(l) = B4(2*l-1)
        p2_mu4(l) = B4(2*l)
      end do
      
      L = 1
      !DISPFLG = WDFLG(L)
      DX = LE(L)
      P1X   = (P1_mu4(2)-P1_mu4(1))/2.D0/DX
      P2X   = (P2_mu4(2)-P1_mu4(1))/2.D0/DX
      P1xx  = (P1_mu4(2)-P1_mu4(1))/DX/DX
      P2xx  = (P2_mu4(2)-P2_mu4(1))/DX/DX
      P3    = (bet0hat(L,1)-(bet3hat(L,2,1)*P2xx+bet2hat(L,2,1)*P2x+       &
     &            bet1hat(L,2,1)*p2_mu4(l)+bet3hat(L,1,1)*P1xx+            &
     &            bet2hat(L,1,1)*P1x+bet1hat(L,1,1)*p1_mu4(l)))            &
     &            /bet1hat(L,3,1)
      P4    = (alp0(L)-alp1(L,3)*P3-(alp2(L,2)*P2x+alp1(L,2)*P2_mu4(L)+    &
     &                               alp2(L,1)*P1x+alp1(L,1)*P1_mu4(L))    &
     &           )/alp1(L,4)
      PDtmp(L) = 0.d0
      do k = 1,pp_negp
        q     = pp_xegp(k)
        PDtmp(L) = PDtmp(L) + (ze_in(l)+de_in(l))*(g*(q-1.d0)*(ze_in(l))+  &
     &                         P1_mu4(L)*pp_phi(1,k)+P2_mu4(L)*pp_phi(2,k)+&
     &                         P3*pp_phi(3,k)+P4*pp_phi(4,k))*pp_wegp(k)
      end do
      PBtmp(L) = -g*ze_in(l)+P1_mu4(L)*pp_phi(1,0)+P2_mu4(L)*pp_phi(2,0)+  &
     &               P3*pp_phi(3,0)+P4*pp_phi(4,0)
        
      PDtmp(L) = PDtmp(L)*(DISPFLG(L))
      PBtmp(L) = PBtmp(L)*(DISPFLG(L))
      
      DO L = 2,NE-1
        !DISPFLG = WDFLG(L)
        DX = LE(L)
        P1xx  = (P1_mu4(L+1)-2.d0*P1_mu4(L)+P1_mu4(L-1))/dx/dx
        P2xx  = (P2_mu4(L+1)-2.d0*P2_mu4(L)+P2_mu4(L-1))/dx/dx
        P1x   = (P1_mu4(L+1)-P1_mu4(L-1))/dx/2.d0
        P2x   = (P2_mu4(L+1)-P2_mu4(L-1))/dx/2.d0
        P3    = (bet0hat(L,1)-(bet3hat(L,2,1)*P2xx+bet2hat(L,2,1)*P2x+     &
     &                bet1hat(L,2,1)*p2_mu4(l)+bet3hat(L,1,1)*P1xx+        &
     &                bet2hat(L,1,1)*P1x+bet1hat(L,1,1)*p1_mu4(l))         &
     &                )/bet1hat(L,3,1)
        P4    = (alp0(L)-alp1(L,3)*P3-(alp2(L,2)*P2x+alp1(L,2)*P2_mu4(L)+  &
     &                alp2(L,1)*P1x+alp1(L,1)*P1_mu4(L)))/alp1(L,4)
        PDtmp(L) = 0.d0
        do k = 1,pp_negp
          q     = pp_xegp(k)
          PDtmp(L) = PDtmp(L) + (ze_in(l)+de_in(l))*(g*(q-1.d0)*(ze_in(l))+&
     &                     P1_mu4(L)*pp_phi(1,k)+P2_mu4(L)*pp_phi(2,k)+    &
     &                     P3*pp_phi(3,k)+P4*pp_phi(4,k))*pp_wegp(k)
        end do
        PBtmp(L) = -g*ze_in(l)+P1_mu4(L)*pp_phi(1,0)+P2_mu4(L)*pp_phi(2,0)+&
      &                P3*pp_phi(3,0)+P4*pp_phi(4,0)
        
        PDtmp(L) = PDtmp(L)*(DISPFLG(L))
        PBtmp(L) = PBtmp(L)*(DISPFLG(L))
      END DO
      L = NE
      !DISPFLG = WDFLG(L)
      DX = LE(L)
      P1X   = (P1_mu4(NE-1)-P1_mu4(NE))/2.D0/DX
      P2X   = (P2_mu4(NE-1)-P1_mu4(NE))/2.D0/DX
      P1xx  = (P1_mu4(NE-1)-P1_mu4(NE))/DX/DX
      P2xx  = (P2_mu4(NE-1)-P2_mu4(NE))/DX/DX
      P3    = (bet0hat(L,1)-(bet3hat(L,2,1)*P2xx+bet2hat(L,2,1)*P2x+       &
     &             bet1hat(L,2,1)*p2_mu4(l)+bet3hat(L,1,1)*P1xx+           &
     &             bet2hat(L,1,1)*P1x+bet1hat(L,1,1)*p1_mu4(l))            &
     &             )/bet1hat(L,3,1)
      P4    = (alp0(L)-alp1(L,3)*P3-(alp2(L,2)*P2x+alp1(L,2)*P2_mu4(L)+    &
     &             alp2(L,1)*P1x+alp1(L,1)*P1_mu4(L)))/alp1(L,4)
      PDtmp(L) = 0.d0
      do k = 1,pp_negp
        q     = pp_xegp(k)
        PDtmp(L) = PDtmp(L) + (ze_in(l)+de_in(l))*(g*(q-1.d0)*(ze_in(l))+  &
     &                   P1_mu4(L)*pp_phi(1,k)+P2_mu4(L)*pp_phi(2,k)+      &
     &                   P3*pp_phi(3,k)+P4*pp_phi(4,k))*pp_wegp(k)
      end do
      PBtmp(L) = -g*ze_in(l)+P1_mu4(L)*pp_phi(1,0)+P2_mu4(L)*pp_phi(2,0)+  &
     &                P3*pp_phi(3,0)+P4*pp_phi(4,0)        
      PDtmp(L) = PDtmp(L)*(DISPFLG(L))
      PBtmp(L) = PBtmp(L)*(DISPFLG(L))
      
      DO I = 1,P+1
        L = 1
        PTMP(1) = C12*(PDTMP(L)+PDTMP(L))
        PTMP(2) = C12*(PDTMP(L+1)+PDTMP(L))          
        PD(I,L,IRK) = PTMP(1) + (PTMP(2)-PTMP(1))/2.D0*                   &
     &                                ( 2.D0*REAL(I-1)/REAL(P) )
        PTMP(1) = C12*(PBTMP(L)+PBTMP(L))
        PTMP(2) = C12*(PBTMP(L+1)+PBTMP(L))          
        PB(I,L,IRK) = PTMP(1) + (PTMP(2)-PTMP(1))/2.D0*                   &
     &                                ( 2.D0*REAL(I-1)/REAL(P) )
        DO L = 2,NE-1
          PTMP(1) = C12*(PDTMP(L-1)+PDTMP(L))
          PTMP(2) = C12*(PDTMP(L+1)+PDTMP(L))          
          PD(I,L,IRK) = PTMP(1) + (PTMP(2)-PTMP(1))/2.D0*                 &
     &                                ( 2.D0*REAL(I-1)/REAL(P) )
          PTMP(1) = C12*(PBTMP(L-1)+PBTMP(L))
          PTMP(2) = C12*(PBTMP(L+1)+PBTMP(L))          
          PB(I,L,IRK) = PTMP(1) + (PTMP(2)-PTMP(1))/2.D0*                 &
     &                                ( 2.D0*REAL(I-1)/REAL(P) )
        END DO
        L = NE
        PTMP(1) = C12*(PDTMP(L-1)+PDTMP(L))
        PTMP(2) = C12*(PDTMP(L)+PDTMP(L))          
        PD(I,L,IRK) = PTMP(1) + (PTMP(2)-PTMP(1))/2.D0*                   &
     &                                ( 2.D0*REAL(I-1)/REAL(P) )
        PTMP(1) = C12*(PBTMP(L-1)+PBTMP(L))
        PTMP(2) = C12*(PBTMP(L)+PBTMP(L))          
        PB(I,L,IRK) = PTMP(1) + (PTMP(2)-PTMP(1))/2.D0*                   &
     &                                ( 2.D0*REAL(I-1)/REAL(P) )
      END DO
      
      RETURN
      END SUBROUTINE NONHYDRO_MU4_FD
!..........................................................................!
!==========================================================================!
!..........................................................................!
      SUBROUTINE PROJ_FE_FD(DI,DIX,DIXX,ZI,ZIX,ZIXX,UI,UIX,QI,QIX)
      
      USE READ_DGINP, ONLY : P
      USE GLOBALS,    ONLY : DE_ED,NE,ZE,QE,PHIB,IRK,LE
      USE SIZES,      ONLY : SZ,C12
      
      IMPLICIT NONE
      INTEGER  :: L,I,J,K
      REAL(SZ) :: QI(NE+1),QIX(NE+1)
      
      REAL(SZ),INTENT(OUT) :: DI(NE),UI(NE),ZI(NE)
      REAL(SZ),INTENT(OUT) :: DIX(NE),UIX(NE),ZIX(NE)
      REAL(SZ),INTENT(OUT) :: DIXX(NE),ZIXX(NE)
      
      REAL(SZ) :: DTMP(2),ZTMP(2),QTMP(2)
      REAL(SZ) :: DXTMP(2),ZXTMP(2),QXTMP(2)
      REAL(SZ) :: DXXTMP(2),ZXXTMP(2)
      

      DO L = 1,NE
        DO I = 1,2
!           DTMP(I) = 0.D0
          ZTMP(I) = 0.D0
          QTMP(I) = 0.D0
          DTMP(1) = DE_ED(L)
          DTMP(2) = DE_ED(L+1)
          DO J = 1,P+1
            ZTMP(I)  = ZTMP(I)  + ZE(J,L,IRK)*PHIB(J,I)
            QTMP(I)  = QTMP(I)  + QE(J,L,IRK)*PHIB(J,I)
          END DO
        END DO
        DI(L)  = C12*(DTMP(1)+DTMP(2))
        ZI(L)  = C12*(ZTMP(1)+ZTMP(2))
        QI(L)  = C12*(QTMP(1)+QTMP(2))
        DIX(L) = (DTMP(2)-DTMP(1))/LE(L)
        ZIX(L) = (ZTMP(2)-ZTMP(1))/LE(L)
        QIX(L) = (QTMP(2)-QTMP(1))/LE(L)
        UI(L)  = QI(L)/(ZI(L)+DI(L))
        UIX(L) = ( (ZI(L)+DI(L))*QIX(L)- &
                     QI(L)*(ZIX(L)+DIX(L)) )/(ZI(L)+DI(L))**2
      END DO

!                    
!.....Second, determine central difference approximation of second derivatives
      DO L = 2,NE-1
        DIXX(L) = (DI(L+1)-2.D0*DI(L)+DI(L-1))/LE(L)/LE(L)
        ZIXX(L) = (ZI(L+1)-2.D0*ZI(L)+ZI(L-1))/LE(L)/LE(L)
      END DO
      DIXX(1)  = (DI(2)-DI(1))/LE(L)/LE(L)
      DIXX(NE) = (DI(NE-1)-DI(NE))/LE(L)/LE(L)
      ZIXX(1)  = (ZI(2)-ZI(1))/LE(L)/LE(L)
      ZIXX(NE) = (ZI(NE-1)-ZI(NE))/LE(L)/LE(L)
            
      
      RETURN
      END SUBROUTINE PROJ_FE_FD
!..........................................................................!
!==========================================================================!
!==========================================================================!
!....................................................................
!
!     SUBROUTINE PRESPOIS_ORDER_2_COEFF
!
!     Subroutine which computes the coefficients for the O(mu^2) 
!     pressure-Poisson model
!
!     Written by Aaron S. Donahue 2015-01-28
!...................................................................
      
      SUBROUTINE PRESPOIS_ORDER_2_COEFF(ZI,ZIX,ZIXX,DI,DIX,DIXX,UI,UIX, &
                    MRHS,MP1,MP1X,MP1XX)
      
      USE GLOBALS, ONLY : G
      USE SIZES

      IMPLICIT NONE
      REAL(SZ),INTENT(IN)  :: ZI,ZIX,ZIXX,DI,DIX,DIXX,UI,UIX
      REAL(SZ),INTENT(OUT) :: MRHS,MP1,MP1X,MP1XX
      
      MP1XX =  (ZI+DI)**2/2.D0
      MP1X  =  (ZI+DI)*(3.D0*DIX+8.D0*ZIX)/12.D0
      MP1   = -(15.D0+3.D0*DIX**2-4.D0*DIX*ZIX+ &
                8.D0*ZIX**2+6.D0*DI*DIXX+ &
                6.D0*ZI*DIXX-4.D0*DI*ZIXX- &
                4.D0*ZI*ZIXX)/12.D0
      MRHS  = -5.D0/3.D0*(DI+ZI)**2*UIX**2+ &
               5.D0/4.D0*(DI+ZI)*DIXX*UI**2- &
               G*(ZI+DI)*DI*ZIXX/3.D0- &
               G*(2.D0*DIX*ZIX+ &
                    DI*(3.D0*DIX-2.D0*ZIX))*ZIX/3.D0- &
               G*ZI*(5.D0+DIX**2+2.D0*DI*DIXX+ &
                    2.D0*ZI*DIXX)/4.D0
                    
                    
      
      RETURN
      END SUBROUTINE

!....................................................................
!
!     SUBROUTINE PRESPOIS_ORDER_2_P2
!
!     Subroutine which computes the value of P2 for the O(mu^2) 
!     pressure-Poisson model
!
!     Written by Aaron S. Donahue 2015-01-28
!...................................................................
      
      SUBROUTINE PRESPOIS_ORDER_2_P2(ZI,DI,DIX,DIXX,UI,P1,P1X,ISBREAK, &
                    P2,PD,PB)
      
      USE GLOBALS, ONLY : G
      USE SIZES

      IMPLICIT NONE
      REAL(SZ),INTENT(IN)  :: ZI,DI,DIX,DIXX,UI,P1,P1X
      INTEGER,INTENT(IN)   :: ISBREAK
      REAL(SZ),INTENT(OUT) :: P2,PD,PB
      
      P2    = (3.D0/4.D0*(1.D0+DIX**2)*P1- &
                3.D0/4.D0*(DI+ZI)*DIX*P1X+ &
                3.D0/4.D0*(DI+ZI)*DIXX*UI**2- &
                3.D0/4.D0*G*ZI*(1.D0+DIX**2))
      PD = C12*(ZI+DI)*(P1-G*ZI)
      PB = (P1-1.D0/3.D0*P2-G*ZI)

      RETURN
      END SUBROUTINE
!..........................................................................!
!==========================================================================!
!..........................................................................!
      SUBROUTINE WASUPP_Coeffs_N(n,m,di,dix,dixx,zi,zix,zixx,alp2n,alp1n,  &
     &                           bet3n,bet2n,bet1n)
      
      USE GLOBALS, ONLY : pp_negp,pp_wegp,pp_xegp,pp_phi,pp_dphi,pp_ddphi, &
     &                    pp_mu,pp_wei,pp_del
      USE SIZES,   ONLY : SZ
      
      IMPLICIT NONE
      
      INTEGER,INTENT(IN)	:: n,m
      REAL(SZ),INTENT(IN)	:: di,dix,dixx,zi,zix,zixx
      REAL(SZ),INTENT(OUT)	:: alp2n,alp1n,bet3n,bet2n,bet1n
      INTEGER	 		:: mu0,mu2,k
      REAL(SZ)			:: q,del
      
      
      del = pp_del
      
      mu0 = pp_mu(n)
      mu2 = pp_mu(n+2)

!.....Bottom Boundary Condition coefficients
      alp2n = 0.d0
      alp1n = pp_dphi(n,0)      
!.....Build Pressure Poisson coefficients through Gaussian Integration      
      bet3n = 0.d0
      bet2n = 0.d0
      bet1n = 0.d0
      do k = 1,pp_negp
        bet1n = bet1n + pp_wegp(k)*pp_ddphi(n,k)*pp_wei(m,k)
      end do
      
!.....Add higher order terms (if applicable)      
      if (mu2.gt.0) then
!.......Bottom Boundary Condition coefficients
        alp2n = alp2n + di*dix*pp_phi(n,0)
        alp1n = alp1n + dix**2*pp_dphi(n,0)  
!.......Build Pressure Poisson coefficients through Gaussian Integration      
        do k = 1,pp_negp
          q = pp_xegp(k)
          bet3n = bet3n + pp_wegp(k)*( di**2*pp_phi(n,k) )*pp_wei(m,k)
          bet2n = bet2n + pp_wegp(k)*( -2.d0*di*dix*(q-1.d0)*pp_dphi(n,k) )&
     &            *pp_wei(m,k)
          bet1n = bet1n + pp_wegp(k)*( (dix**2*(q-1.d0)**2)*pp_ddphi(n,k)  &
     &            +(2.d0*dix**2-di*dixx)*(q-1.d0)*pp_dphi(n,k) )*pp_wei(m,k)
        end do      
      end if
      
      RETURN
      END SUBROUTINE WASUPP_Coeffs_N
!..........................................................................!
!==========================================================================!
!..........................................................................!
      SUBROUTINE WASUPP_Coeffs_1(m,di,dix,dixx,zi,zix,zixx,alp21,alp11,    &
     &                           bet31,bet21,bet11)
      
      USE GLOBALS, ONLY : pp_negp,pp_wegp,pp_xegp,pp_phi,pp_dphi,pp_ddphi, &
     &                    pp_del,pp_wei
      USE SIZES,   ONLY : SZ
      
      IMPLICIT NONE
      INTEGER,INTENT(IN)	:: m
      REAL(SZ),INTENT(IN)	:: di,dix,dixx,zi,zix,zixx
      REAL(SZ),INTENT(OUT)	:: alp21,alp11,bet31,bet21,bet11
      INTEGER			:: k,del
      REAL(SZ)			:: q
      
      del = pp_del

!.....Bottom Boundary Condition coefficients   
      alp21 = (del*zi+di)*dix*pp_phi(1,0)
      alp11 = (1.d0+dix**2)*pp_dphi(1,0)
!.....Build Pressure Poisson coefficients through Gaussian Integration     
      bet31 = 0.d0
      bet21 = 0.d0
      bet11 = 0.d0
      do k = 1,pp_negp
        q = pp_xegp(k)
        bet31 = bet31 + pp_wegp(k)*( (del*zi+di)**2*pp_phi(1,k) )*pp_wei(m,k)
        bet21 = bet21 + pp_wegp(k)*( -2.d0*(del*zi+di)*((del*zix+dix)*q-   &
     &                     dix)*pp_dphi(1,k) )*pp_wei(m,k)
        bet11 = bet11 + pp_wegp(k)*( ((del*zi+di)*dixx-2.d0*dix*(del*zix+  &
     &                     dix)+(2.d0*(del*zix+dix)**2-(del*zi+di)*        &
     &                     (del*zixx+dixx))*q)*pp_dphi(1,k) )*pp_wei(m,k)
      end do
      
      RETURN
      END SUBROUTINE WASUPP_Coeffs_1
!..........................................................................!
!==========================================================================!
!..........................................................................!
      SUBROUTINE WASUPP_Coeffs_0(m,di,dix,dixx,zi,zix,zixx,ui,uix,alp0,bet0)
      
      USE GLOBALS, ONLY : g,pp_negp,pp_wegp,pp_xegp,pp_phi,pp_dphi,        &
     &                    pp_ddphi,pp_del,pp_wei
      USE SIZES,   ONLY : SZ
      
      
      IMPLICIT NONE
      INTEGER,INTENT(IN)	:: m
      REAL(SZ),INTENT(IN)	:: di,dix,dixx,zi,zix,zixx,ui,uix
      REAL(SZ),INTENT(OUT)	:: alp0,bet0
      INTEGER			:: k,del
      REAL(SZ)			:: q
      
      del = pp_del
!.....Bottom Boundary Condition coefficient      
      alp0 = del*(del*zi+di)*dixx*ui**2-g*(1.d0+dix**2)*zi
!.....Build Pressure Poisson coefficient through Gaussian Integration   
      bet0 = 0.d0
      do k = 1,pp_negp
        q = pp_xegp(k)
   
        bet0 = bet0 + pp_wegp(k)*( -g*q*di*(del*zi+di)*zixx &
       &                +2.d0*g*(-q*del*zi*dix+di*((q-1)*dix+del*q*zix))*  &
       &                zix+g*(q-1)*((del*zi+di)*dixx-2.d0*dix**2)*zi      &
       &                -2.d0*(del*zi+di)**2*uix**2 ) *pp_wei(m,k)
      end do
      
      
      END SUBROUTINE WASUPP_Coeffs_0
!..........................................................................!
!==========================================================================!
!..........................................................................!
      SUBROUTINE SETUP_WASUPP
      
      USE READ_DGINP, ONLY : INONHYDRO,NONHYDRO_EXT,P
      USE GLOBALS,    ONLY : PP_NEGP,PP_WEGP,PP_XEGP,PP_PHI,PP_DPHI,       &
     &                       PP_DDPHI,PP_WEI,PP_MU,WDFLG,DISPFLG,NE,PPCNT, &
     &                       PDLVL,PBLVL
      
      IMPLICIT NONE
      
      pp_negp = max(2*INONHYDRO,1)
      ALLOCATE(pp_wegp(pp_negp),pp_xegp(pp_negp))
      ALLOCATE(pp_mu(INONHYDRO+2))
      
      ALLOCATE(pp_phi(INONHYDRO,0:pp_negp),pp_dphi(INONHYDRO,0:pp_negp))
      ALLOCATE(pp_ddphi(INONHYDRO,1:pp_negp))
      ALLOCATE(pp_wei(INONHYDRO,1:pp_negp))
      
      CALL WASUPP_Gauss
      CALL WASUPP_Basis
      
      ! Setup extrapolation variables
      ALLOCATE(PDLVL(NONHYDRO_EXT+1,NE,P+1),PBLVL(NONHYDRO_EXT+1,NE,P+1))
      PPCNT = 0

      ! Setup the dispersion flags, initialize through wet/dry flag
      DISPFLG(:) = WDFLG(:)
      
      RETURN
      END SUBROUTINE SETUP_WASUPP
!..........................................................................!
!==========================================================================!
!..........................................................................!
      SUBROUTINE WASUPP_Gauss
      
      USE READ_DGINP, ONLY : INONHYDRO
      USE GLOBALS,    ONLY : PP_NEGP,PP_WEGP,PP_XEGP
      
      IMPLICIT NONE
      
      if (pp_negp.eq.1) then
        pp_wegp(1) = 1.d0
        
        pp_xegp(1) = 0.5d0
        
      elseif (pp_negp.eq.3) then
        pp_wegp(3) = 0.277777777777777d0
        pp_wegp(2) = 0.444444444444444d0
        pp_wegp(1) = 0.277777777777777d0
        
        pp_xegp(3) = 0.887298334620742d0
        pp_xegp(2) = 0.500000000000000d0
        pp_xegp(1) = 0.112701665379258d0
        
      elseif (pp_negp.eq.4) then
        pp_wegp(4) = 0.173927422568727d0
        pp_wegp(3) = 0.326072577431273d0
        pp_wegp(2) = 0.326072577431273d0
        pp_wegp(1) = 0.173927422568727d0
   
        pp_xegp(4) = 0.930568155797026d0
        pp_xegp(3) = 0.669990521792428d0
        pp_xegp(2) = 0.330009478207572d0
        pp_xegp(1) = 0.069431844202974d0
                
      elseif (pp_negp.eq.5) then
        pp_wegp(5) = 0.118463442528095d0
        pp_wegp(4) = 0.239314335249683d0
        pp_wegp(3) = 0.284444444444444d0
        pp_wegp(2) = 0.239314335249683d0
        pp_wegp(1) = 0.118463442528095d0
        
        pp_xegp(5) = 0.953089922969332d0
        pp_xegp(4) = 0.769234655052841d0
        pp_xegp(3) = 0.500000000000000d0
        pp_xegp(2) = 0.230765344947158d0
        pp_xegp(1) = 0.046910077030668d0   
        
      elseif (pp_negp.eq.8) then  
        pp_wegp(8) = 0.050614268145188D0
        pp_wegp(7) = 0.111190517226687D0
        pp_wegp(6) = 0.156853322938944D0
        pp_wegp(5) = 0.181341891689181D0
        pp_wegp(4) = 0.181341891689181D0
        pp_wegp(3) = 0.156853322938944D0
        pp_wegp(2) = 0.111190517226687D0
        pp_wegp(1) = 0.050614268145188D0
        
        pp_xegp(8) = 0.980144928248768d0
        pp_xegp(7) = 0.898333238706813d0
        pp_xegp(6) = 0.762766204958164d0
        pp_xegp(5) = 0.591717321247825d0
        pp_xegp(4) = 0.408282678752175d0
        pp_xegp(3) = 0.237233795041836d0
        pp_xegp(2) = 0.101666761293187d0
        pp_xegp(1) = 0.019855071751232d0
        
      else
        write(*,*) 'Higher than fourth order pressure-Poisson not supported'
        write(*,*) 'Fix IDYNAMIC parameter and rerun simulation'
        STOP
      end if     
      
      END SUBROUTINE WASUPP_Gauss
!..........................................................................!
!==========================================================================!
!..........................................................................!
      SUBROUTINE WASUPP_Basis
      
      USE READ_DGINP, ONLY : INONHYDRO
      USE GLOBALS,    ONLY : PP_PHI,PP_DPHI,PP_DDPHI,PP_WEI,PP_DEL,PP_MU,  &
     &                       PP_NEGP,PP_XEGP
      USE SIZES,      ONLY : SZ
      
      IMPLICIT NONE
      INTEGER  :: i,j,k
      REAL(SZ) :: pp_coeff(INONHYDRO,INONHYDRO),we_coeff(INONHYDRO,INONHYDRO)
      
      pp_del = 1.d0
      
      pp_mu(1) = 1
      do i = 2,INONHYDRO+2
        if ((i+mod(i,2)).le.INONHYDRO) then
          pp_mu(i) = 1
        else
          pp_mu(i) = 0
        end if
      end do

      pp_coeff(:,:) = 0.d0
      we_coeff(:,:) = 0.d0
      if (INONHYDRO.eq.2) then        
        pp_coeff(1,1) = 1.d0
        
        pp_coeff(2,1) = -4.d0/3.d0
        pp_coeff(2,2) = 1.d0
        
        we_coeff(1,1) = -4.d0/3.d0
        we_coeff(1,2) = 1.d0        
        
      elseif (INONHYDRO.eq.4) then
        pp_coeff(1,1) =  1.d0
        
        pp_coeff(2,1) = -1.790467309d0
        pp_coeff(2,2) =  1.d0
        
        pp_coeff(3,1) = -0.1206295845d0
        pp_coeff(3,2) = -1.023252036375844d0
        pp_coeff(3,3) =  1.d0
        
        pp_coeff(4,1) =  1.d0
        pp_coeff(4,2) =  15.84705801627071d0
        pp_coeff(4,3) = -15.902695460613225d0
        pp_coeff(4,4) =  1.d0
        
        we_coeff(1,1) = -1.d0
        we_coeff(1,2) =  1.d0
        
        we_coeff(2,1) =  0.790467309000000d0
        we_coeff(2,2) = -1.790467309000000d0
        we_coeff(2,3) =  1.d0
        
        we_coeff(3,1) =  0.143881620875844d0
        we_coeff(3,2) = -0.120629584500000d0
        we_coeff(3,3) = -1.023252036375844d0
        we_coeff(3,4) =  1.d0
      end if
      
      
      pp_phi(:,:)   = 0.d0
      pp_dphi(:,:)  = 0.d0
      pp_ddphi(:,:) = 0.d0
      pp_wei(:,:)   = 0.d0
      
      do i = 1,INONHYDRO
        pp_phi(i,0)   = 0.d0
        pp_dphi(i,0)  = -pp_coeff(i,1)
        do j = 1,i
          pp_phi(i,0) = pp_phi(i,0) + pp_coeff(i,j)*(1.d0-0.d0**j)
        end do
        do k = 1,pp_negp
          pp_phi(i,k)   = 0.d0
          pp_dphi(i,k)  = 0.d0
          pp_ddphi(i,k) = 0.d0
          if (i.lt.INONHYDRO) then
            pp_wei(i,k) = we_coeff(i,1)
          end if
          do j = 1,i
            pp_phi(i,k)   = pp_phi(i,k)   + pp_coeff(i,j)*(1.d0-pp_xegp(k)**j)
            pp_dphi(i,k)  = pp_dphi(i,k)  + pp_coeff(i,j)*(-j*pp_xegp(k)**(j-1))
            if (i.lt.INONHYDRO) then
              pp_wei(i,k) = pp_wei(i,k)   + we_coeff(i,j+1)*pp_xegp(k)**j
            end if
          end do
          do j = 2,i
            pp_ddphi(i,k) = pp_ddphi(i,k) + pp_coeff(i,j)*(-j*(j-1)*pp_xegp(k)**(j-2))
          end do
        end do
      end do
      
      END SUBROUTINE WASUPP_Basis


!..........................................................................!
!==========================================================================!
!..........................................................................!
