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
      SUBROUTINE SLDIAG(SUT,SVT,STT,SQT,CDM,CDH,UA,VA,TA,QA,T0,Q0,      
     1                  Z0M,Z0E,F,ZA,ZU,ZT,ILG,IL1,IL2,JL)              
C     * JUN 23/14 - M.LAZARE.   New version for gcm18+:                 
C     *                         - Bugfix to calculation of              
C     *                           screen temperature and                
C     *                           screen specific humidity.             
C     *                         - Accumulation removed (now             
C     *                           done in classt/oiflux11) so           
C     *                           that a screen relative humidity       
C     *                           can be calculated. Therefore,        
C     *                           "instantaneous" fields are            
C     *                           calculated and passed out             
C     *                           instead.                              
C     * OCT 17/11 - D.VERSEGHY. ADD CODE TO CIRCUMVENT SPECIAL          
C     *                         CASE WHERE TA~T0 OR QA~QO, THUS         
C     *                         AVOIDING A DIVIDE BY ZERO.              
C     * NOV 04/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.            
C     * JUL 19/96 - Y. DELAGE.                                          
C     * CALCULATES NEAR SURFACE OUTPUT VARIABLES                        
c     * OUTPUT FIELDS ARE:                                              
C     *   SUT  : U COMPONENT OF THE WIND AT ZU                          
C     *   SVT  : V COMPONENT OF THE WIND AT ZU                          
C     *   STT  : TEMPERATURE AT ZT                                      
C     *   SQT  : SPECIFIC HUMIDITY AT ZT                                
C     * INPUT FIELDS ARE:                                               
C     *   CDM : DRAG COEFFICIENT                                        
C     *   CDH : TRASFER COEFFICIENT FOR HEAT AND MOISTURE               
C     *   UA  : U COMPONENT OF THE WIND AT ZA                           
C     *   VA  : V COMPONENT OF THE WIND AT ZA                           
C     *   TA  : POTENTIAL TEMPERATURE AT ZA                             
c     *   T0  : TEMPERATURE AT BOTTOM OF SURFACE LAYER                  
C     *   Q0  : SPECIFIC HUMIDITY AT BOTTOM OF SURFACE LAYER            
C     *   Z0M : ROUGHNESS LENGTH FOR MOMENTUM                           
C     *   Z0E : ROUGHNESS LENGTH FOR HEAT AND MOISTURE                  
C     *   F   : FRACTION OF GRID POINT AFFECTED BY PROCESS              
C     *   ZA  : TOP OF SURFACE LAYER                                    
C     *   ZU  : HEIGHT OF OUTPUT WIND                                   
C     *   ZT  : HEIGHT OF OUTPUT TEMPERATURE AND HUMIDITY               
C     *   ILG : NUMBER OF POINTS TO BE TREATED                          
C                                                                       
      IMPLICIT NONE                                                     
c                                                                       
      INTEGER ILG,IL1,IL2,JL,I                                          
c                                                                       
      REAL  SUT(ILG), SVT(ILG), STT(ILG), SQT(ILG), CDM(ILG), CDH(ILG), 
     1       UA(ILG),  VA(ILG),  TA(ILG),  QA(ILG), Z0M(ILG), Z0E(ILG), 
     2        F(ILG),  T0(ILG),  Q0(ILG),  ZA(ILG),  ZU(ILG),  ZT(ILG)  
c                                                                       
      REAL PR,WSPD,CM,US,TS,QS,L,UVA,RATIO,UVU,TTA,CE                   
c                                                                       
      REAL RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN                              
c                                                                       
      REAL PSM,PSE,Y,PIM,PIE,X                                          
c                                                                       
      COMMON /CLASS2/ RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN                   
c     * STABILITY FUNCTIONS FOR THE STABLE CASE                         
      PSM(X)= -X -.667*(X-5/.35)*EXP(-.35*X)                            
      PSE(X)= -(1+.667*X)**1.5 -.667*(X-5/.35)*EXP(-.35*X)              
c     * STABILITY FUNCTIONS FOR THE UNSTABLE CASE                       
      Y(X)=(1-16*X)**.25                                                
      PIM(X)= LOG((1+X)**2*(1+X**2)) -2*ATAN(X)                         
      PIE(X)= 2*LOG(1+X**2)                                             
      PR=1.0                                                            
      DO 100 I=IL1,IL2                                                  
      IF(F(I).GT.0.)                                                THEN
C     * CALCULATION OF SURFACE FLUXES AND MONIN-OBUKHOV LENGTH          
        WSPD=MAX(VMIN,SQRT(UA(I)**2+VA(I)**2))                          
        CM=SQRT(CDM(I))                                                 
        US=CM*WSPD                                                      
        IF(ABS(TA(I)-T0(I)).LT.0.01) THEN                               
            TS=-0.01*CDH(I)/CM                                          
        ELSE                                                            
            TS=CDH(I)*(TA(I)-T0(I))/CM                                  
        ENDIF                                                           
        IF(ABS(QA(I)-Q0(I)).LT.1.0E-7) THEN                             
            QS=-1.0E-7*CDH(I)/CM                                        
        ELSE                                                            
            QS=CDH(I)*(QA(I)-Q0(I))/CM                                  
        ENDIF                                                           
        L=TA(I)*US**2/(VKC*GRAV*(TS*(1+.61*QA(I))+.61*TA(I)*QS))        
C     * CALCULATE CORRECTION FACTORS TO TAKE INTO ACCOUNT THE APPROXIMATIONS
C     * IN DRCOEF                                                          
        IF(L.GT.0.)                                                 THEN
C     * STABLE CASE                                                     
         UVA=US/VKC*(LOG(ZA(I)/Z0M(I))-PSM(ZA(I)/L)+PSM(Z0M(I)/L))      
         RATIO=WSPD/UVA                                                 
         UVU=US/VKC*(LOG((ZU(I)+Z0M(I))/Z0M(I))-PSM((ZU(I)+Z0M(I))/L)   
     1       +PSM(Z0M(I)/L))*RATIO                                      
         TTA=T0(I)+TS/VKC*PR*(LOG(ZA(I)/Z0E(I))-PSE(ZA(I)/L)+           
     1          PSE(Z0E(I)/L))                                          
         RATIO=(TA(I)-T0(I))/SIGN(MAX(ABS(TTA-T0(I)),1.E-4),TTA-T0(I))  
         CE=(LOG((ZT(I)+Z0M(I))/Z0E(I))-PSE((ZT(I)+Z0M(I))/L)           
     1      +PSE(Z0E(I)/L))*RATIO*PR/VKC                                
        ELSE                                                            
C     * UNSTABLE CASE                                                   
         UVA=US/VKC*(LOG(ZA(I)/Z0M(I))-PIM(Y(ZA(I)/L))+PIM(Y(Z0M(I)/L)))
         RATIO=WSPD/UVA                                                 
         UVU=US/VKC*(LOG((ZU(I)+Z0M(I))/Z0M(I))-PIM(Y((ZU(I)+Z0M(I))/L))
     1        +PIM(Y(Z0M(I)/L)))*RATIO                                  
         TTA=T0(I)+TS/VKC*PR*(LOG(ZA(I)/Z0E(I))-PIE(Y(ZA(I)/L))+        
     1          PIE(Y(Z0E(I)/L)))                                       
         RATIO=(TA(I)-T0(I))/SIGN(MAX(ABS(TTA-T0(I)),1.E-4),TTA-T0(I))  
         CE=(LOG((ZT(I)+Z0M(I))/Z0E(I))-PIE(Y((ZT(I)+Z0M(I))/L))        
     1      +PIE(Y(Z0E(I)/L)))*RATIO*PR/VKC                             
        ENDIF                                                           
C                                                                       
        SUT(I)=UVU*UA(I)/WSPD                                           
        SVT(I)=UVU*VA(I)/WSPD                                           
        STT(I)=T0(I)+TS*CE                                              
        SQT(I)=Q0(I)+QS*CE                                              
      ENDIF                                                             
  100 CONTINUE                                                          
      RETURN                                                            
      END                                                               
