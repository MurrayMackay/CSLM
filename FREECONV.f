!-------------------------------------- LICENCE BEGIN ------------------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END --------------------------------------
         SUBROUTINE FREECONV(LKICEH,T0,TLAK,RHOIW,NLAK,NLAKMAX,
     >                  ILG,IL1,IL2,IYEAR,IHOUR,IDAY,IMIN )
C=======================================================================
      IMPLICIT NONE
C
C ----* LAKE MODEL VARIABLES *----------------------------------------
C
      INTEGER NLAKMAX
      INTEGER,DIMENSION(ILG) :: NLAK
      REAL,DIMENSION(ILG) :: T0, LKICEH
      REAL,DIMENSION(ILG,NLAKMAX) :: TLAK
C
C ----* INPUT *-------------------------------------------
C
      INTEGER ILG,IL1,IL2,IYEAR,IDAY,IHOUR,IMIN
      REAL RHOIW
C
C ----* CLASS COMMON BLOCKS *------------------------------------------
C
      REAL DELT,TFREZ
      REAL TKECN,TKECF,TKECE,TKECS,TKECL,HDPTHMIN,
     1     TKEMIN,DELMAX,DELMIN,EMSW,DELZLK,DELSKIN,DHMAX,DUMAX
      COMMON /CLASS1/ DELT,TFREZ                                       
      COMMON /LAKECON/ TKECN,TKECF,TKECE,TKECS,HDPTHMIN,        
     2                 TKEMIN,DELMAX,DELMIN,EMSW,DELZLK,DELSKIN,
     3                 DHMAX,TKECL,DUMAX
C
C ----* LOCAL VARIABLES *---------------------------------------------
C
      INTEGER I,J,K,NMIX
      REAL ZTOP,ZBOT,TTEST,RHO1,RHO2,TC1,TC2,TBAR,XXX,ICEBOT
C=======================================================================
C
      DO 100 I=IL1,IL2
         IF (LKICEH(I) .LE. 0.0) THEN
           TC1=T0(I)-TFREZ
           TC2=TLAK(I,1)-TFREZ
           CALL EQNST(XXX,RHO1,TC1,0.05)
           CALL EQNST(XXX,RHO2,TC2,0.5)
           IF (RHO1 .GT. RHO2) THEN
             TBAR=((DELSKIN*RHO1*T0(I))+(DELZLK*RHO2*TLAK(I,1)))/
     >             ((DELSKIN*RHO1)+(DELZLK*RHO2))
             T0(I)=TBAR
             TLAK(I,1)=TBAR
           ENDIF
         ENDIF
         ICEBOT=RHOIW*LKICEH(I)
        
        NMIX=1
        DO 420, J=1,NLAK(I)-1
          ZTOP=DELSKIN + (J-1)*DELZLK
          ZBOT=ZTOP+DELZLK
          IF (ICEBOT .LE. ZTOP) THEN
           TC1=TLAK(I,J)-TFREZ
           TC2=TLAK(I,J+1)-TFREZ
Cmdm       TTEST=(TC1-3.9816)*(TC2-3.9816)
           TTEST=(TC1-3.98275)*(TC2-3.98275)
           CALL EQNST(XXX,RHO1,TC1,ZBOT)
           CALL EQNST(XXX,RHO2,TC2,ZBOT+DELZLK)
C--------- MIX LAYERS IF RHO1>RHO2 OR TEMPERATURES SPAN 
C--------- T_MAXDENSITY=3.9816 C.
          IF ((RHO1 .GT. RHO2) .OR. (TTEST .LT. 0.0)) THEN
            TBAR=((NMIX*RHO1*TLAK(I,J))+(RHO2*TLAK(I,J+1)))/
     >                      (NMIX*RHO1+RHO2)
            DO 430, K=J-NMIX+1,J+1
              TLAK(I,K)=TBAR
430         CONTINUE
            NMIX=NMIX+1
C           WRITE(6,6666) "static instability removed under ice:",
C    >                IYEAR,IDAY,IHOUR,IMIN,j*DELZLK
          ELSE
            NMIX=1
          ENDIF
         ENDIF
420     CONTINUE
100   CONTINUE
6666  FORMAT(A37,4I5,F5.1)
      RETURN
      END
