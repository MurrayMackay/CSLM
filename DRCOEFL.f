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
      SUBROUTINE DRCOEFL (CDM,CDH,CDMN,VA,T0,TA,QA,PRES,ZREFM,ZREFH,
     1                    FICE,FLS,ILG,IL1,IL2)
C=======================================================================
C     * DEC 15/16 - M.LAZARE/   - NSTEP REMOVED.
C     *             D.VERSEGHY. - BUGFIX IN ACCOUNTING FOR EFFECT
C     *                           OF FICE ON DRAG COEFFICIENTS.
C     * JAN 25/16 - M.MACKAY.   FRACTIONAL ICE COVER INCLUDED
C     *                         BUG FIX FOR NEAR SFC HUMIDITY
C     *                         THERMODYNAMIC REF HEIGHT ADDED
C     * MAY 22/15 - D.VERSEGHY. WEIGHT DRAG COEFFICIENTS FOR PRESENCE
C     *                         OF ICE.
C     * NOV 27/07 - M.MACKAY.  	TURBULENT TRANSFER COEFFICIENTS
C     *
C=======================================================================
C
      IMPLICIT NONE
C
C ----* INPUT FIELDS *------------------------------------------------
C
      INTEGER ILG,IL1,IL2
      REAL,DIMENSION(ILG) :: VA,T0,TA,QA,PRES,ZREFM,ZREFH,FICE,FLS
C
C ----* OUTPUT FIELDS *------------------------------------------------
C
      REAL,DIMENSION(ILG) ::
     1  CDH, CDM, CDMN
C
C ----* COMMON BLOCK PARAMETERS *--------------------------------------
C
      REAL RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN,DELT,TFREZ
C
C ----* CLASS COMMON BLOCKS *------------------------------------------
C
      COMMON /CLASS1/ DELT,TFREZ                                       
      COMMON /CLASS2/ RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN
C
C ----* LOCAL VARIABLES *---------------------------------------------
C
      INTEGER I,J,K,ITER,ITMAX
      REAL G,TOL,CHARN,RESID,VSQ,Cold,EASAT,QASAT,CA,CB,SHF,LHF,
     1     TVIRT,MOL,Z_L,CDHN,X,PSIM,PSIH,PI,MOLold,TOL2,
     2     CDHRAT,CDMNW,CDHNW,CDMNI,CDHNI,DENOM,XICE
C
C ----* LOCAL PARAMETERS *--------------------------------------------
C
      G=GRAV
      TOL=1.0E-8
      TOL2=1.0E-2
      ITMAX=100
      CHARN=0.0175
      PI=3.14159
C
C======================================================================
      DO 100 I=IL1,IL2
C-----------------------------------------------------------------------
C NEUTRAL DRAG COEFFICIENTS  
C  Iterative routine to solve for CDN based on Charnock relation for 
C  roughness
C
        CDHNW=1.35E-3
        CDMNW=1.0E-3        !initial trial value
        VSQ=VA(I)*VA(I)
        RESID=999.0
        DO 
          IF (ABS(RESID) .LE. TOL ) EXIT
          Cold=CDMNW
          CDMNW=VKC*VKC/
     >     (LOG(ZREFM(I)*G/(CHARN*Cold*VSQ))*LOG(ZREFM(I)*G/
     1     (CHARN*Cold*VSQ)))
          RESID=CDMNW-Cold
        END DO   
C
        CDMNI=(VKC/(LOG(ZREFM(I)/0.002)))**2
        CDHNI=(VKC/(LOG(ZREFH(I)/0.00067)))**2
        IF(FICE(I).GT.(FLS(I)+0.001))                THEN
           XICE=MAX(FICE(I)-FLS(I),0.0)
           CDMN(I)=(XICE*CDMNI+(1.0-FICE(I))*CDMNW)/(1.0-FLS(I))
           CDHN=(XICE*CDHNI+(1.0-FICE(I))*CDHNW)/(1.0-FLS(I))
        ELSE
           CDMN(I)=CDMNW
           CDHN=CDHNW
        ENDIF
C-----------------------------------------------------------------------
C INITIAL TRIAL VALUES FOR TRANSFER COEFFICIENTS: SET TO NEUTRAL VALUES
C
       CDH(I)=CDHN
       CDM(I)=CDMN(I)
C-----------------------------------------------------------------------
C ITERATIVELY COMPUTE DRAG COEFFICIENTS UNTIL M.O. LENGTH CONVERGES
C
        RESID=999.0
        MOL=9999.0
        ITER=0
        DO
          IF (ABS(RESID) .LE. TOL2 .OR. ITER .GE. ITMAX) EXIT
C-----------------------------------------------------------------------
C HEAT FLUXES
C----------------------------------------------------
C     * CALCULATION OF EASAT CONSISTENT WITH CLASSI
C     * BUT CONSTANTS DIFFER FROM ROGERS&YAU
C     * Rogers and Yau values
C         CA=17.67
C         CB=29.65
C----------------------------------------------------
          IF(T0(I).GE.TFREZ) THEN                                             
              CA=17.269                                                       
              CB=35.86                                                        
          ELSE                                                                
              CA=21.874                                                       
              CB=7.66                                                         
          ENDIF                                                               
        SHF=CDH(I)*VA(I)*(T0(I)-TA(I))
        EASAT=611.0*EXP(CA*(T0(I)-TFREZ)/(T0(I)-CB))
        QASAT=0.622*EASAT/(PRES(I)-0.378*EASAT)
        LHF=CDH(I)*VA(I)*(QASAT-QA(I))

C-----------------------------------------------------------------------
C VIRTUAL TEMPERATURE AND M.-O. LENGTH
C
        TVIRT=TA(I)*(1.0+0.61*QA(I))
        MOLold=MOL
        MOL = -VA(I)*VA(I)*VA(I)*CDM(I)*SQRT(CDM(I))*TVIRT/
     >       ( VKC*G*(SHF + 0.61*LHF*TA(I)) )
        Z_L = ZREFM(I)/MOL
        RESID=MOL-MOLold

C-----------------------------------------------------------------------
C STABILITY CORRECTIONS
C
C
C- UNSTABLE CASE
C---------------
        IF (Z_L .LT. 0.0) THEN
           X = (1.0 - (16.0*Z_L))**0.25
           PSIM = 2.0*LOG((1.0+X)/2.0) + LOG((1.0+X*X)/2.0)
     >           - 2.0*ATAN(X) + PI/2.0
           PSIH = 2.0*LOG((1.0+X*X)/2.0)
C
C- STABLE CASE
C-------------
        ELSE IF (Z_L .GE. 0 .AND. Z_L .LT. 0.5) THEN
           PSIM = -5.0*Z_L
           PSIH=PSIM
        ELSE IF (Z_L .GE. 0.5 .AND. Z_L .LT. 10.0) THEN
           PSIM = (0.5/(Z_L*Z_L)) - (4.25/Z_L) - 7.0*LOG(Z_L) - 0.852
           PSIH=PSIM
        ELSE 
           PSIM = LOG(Z_L) - 0.76*Z_L - 12.093
           PSIH=PSIM
        END IF

C-----------------------------------------------------------------------
C RECOMPUTE DRAG COEFFICIENTS WITH STABILTY CORRECTIONS
C
        DENOM = (1.0 + (CDMN(I)/(VKC*VKC))*(PSIM*PSIM 
     >          - (2.0*VKC*PSIM/SQRT(CDMN(I)))) )
        IF(DENOM.LT.1.0E-6) DENOM=1.0E-6
C        IF(ABS(DENOM).LT.1.0E-6) DENOM=SIGN(1.0E-6,DENOM)
        CDM(I) = CDMN(I)/DENOM
        DENOM = (1.0 + (CDHN/(VKC*VKC))*(PSIM*PSIH 
     >         -(VKC*PSIH/SQRT(CDMN(I)))-(VKC*PSIM*SQRT(CDMN(I))/CDHN)))
        IF(DENOM.LT.1.0E-6) DENOM=1.0E-6
C        IF(ABS(DENOM).LT.1.0E-6) DENOM=SIGN(1.0E-6,DENOM)
        CDH(I) = CDHN/DENOM

        ITER=ITER+1
        END DO

      IF(ITER.GE.ITMAX)      THEN
        CDM(I)=CDMN(I)
        CDH(I)=CDHN
      ENDIF

C     IF (ITER .GE. ITMAX) print*, "** max iters reached"
      
      IF (CDH(I) .LT. 0.0) THEN 
       CDH(I)=CDHN
      ENDIF

      IF (CDM(I) .LT. 0.0) THEN 
       CDM(I)=CDMN(I)
      ENDIF

      CDHRAT=CDH(I)/CDHN
      IF (CDHRAT .GE. 8.0) THEN
        CDH(I)=CDHN
      ENDIF

100   CONTINUE

      RETURN                                                                      
      END        
