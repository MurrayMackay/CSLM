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
      SUBROUTINE SNOALBA(ALVSSN,ALIRSN,ALVSSC,ALIRSC,ALBSNO,            
     1                   TRSNOWC, ALSNO, TRSNOWG, FSDB, FSFB, RHOSNO,   
     2                   REFSN,BCSN,SNO,CSZ,ZSNOW,FSNOW,ASVDAT,ASIDAT,  
     3                   ALVSG, ALIRG,                                  
     4                   ILG,IG,IL1,IL2,JL,IALS,NBS,ISNOALB)            
C                                                                       
C     * MAR 08/16 - M.Mackay    4-BAND SOLAR RADIATION DISABLED
C     * NOV 16/13 - J.COLE.     Final version for gcm17:                
C     *                         - Fixes to get the proper BC mixing ratio in 
C     *                           snow, which required passing in and using  
C     *                           the snow density RHON.                
C     * JUN 22/13 - J.COLE.     ADD CODE FOR "ISNOALB" OPTION,          
C     *                         WHICH IS BASED ON 4-BAND SOLAR.         
C     * FEB 05/07 - D.VERSEGHY. STREAMLINE CALCULATIONS OF              
C     *                         ALVSSN AND ALIRSN.                      
C     * APR 13/06 - D.VERSEGHY. SEPARATE ALBEDOS FOR OPEN AND           
C     *                         CANOPY-COVERED SNOW.                    
C     * NOV 03/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.            
C     * MAR 18/02 - D.VERSEGHY. UPDATES TO ALLOW ASSIGNMENT OF          
C     *                         USER-SPECIFIED VALUES TO SNOW           
C     *                         ALBEDO.                                 
C     * JUN 05/97 - D.VERSEGHY. CLASS - VERSION 2.7.                    
C     *                         SPECIFY LOCATION OF ICE SHEETS          
C     *                         BY SOIL TEXTURE ARRAY RATHER            
C     *                         THAN BY SOIL COLOUR INDEX.              
C     * NOV 29/94 - M.LAZARE.   CLASS - VERSION 2.3.                    
C     *                         CALL ABORT CHANGED TO CALL XIT TO       
C     *                         ENABLE RUNNING ON PC'S.                 
C     * MAR 13/92 - M.LAZARE.   CODE FOR MODEL VERSION GCM7 -           
C     *                         DIVIDE PREVIOUS SUBROUTINE              
C     *                         "SNOALB" INTO "SNOALBA" AND             
C     *                         "SNOALBW" AND VECTORIZE.                
C     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -          
C     *                         CLASS VERSION 2.0 (WITH CANOPY).        
C     * APR 11/89 - D.VERSEGHY. DISAGGREGATE SNOW ALBEDO INTO           
C     *                         VISIBLE AND NEAR-IR PORTIONS;           
C     *                         CALCULATE TRANSMISSIVITY TO             
C     *                         SHORTWAVE RADIATION.                    
C                                                                       
      IMPLICIT NONE                                                     
C                                                                       
C     * INTEGER CONSTANTS.                                              
C                                                                       
      INTEGER ILG,IG,IL1,IL2,JL,IALS,IPTBAD,I,IB,NBS,ISNOALB            
C                                                                       
C     * OUTPUT ARRAYS.                                                  
C                                                                       
      REAL   ALSNO(ILG,NBS), TRSNOWG(ILG,NBS)                           
      REAL   ALVSSN(ILG),  ALIRSN(ILG),  ALVSSC(ILG),  ALIRSC(ILG),     
     1       ALVSG (ILG),  ALIRG (ILG),  TRSNOWC(ILG)                   
C                                                                       
C     * INPUT ARRAYS.                                                   
C                                                                       
      REAL   FSDB(ILG,NBS), FSFB(ILG,NBS)                               
      REAL   ALBSNO(ILG),  ZSNOW (ILG),  FSNOW (ILG),                   
     1       ASVDAT(ILG),  ASIDAT(ILG),  REFSN (ILG),  BCSN  (ILG),     
     2       CSZ   (ILG),  SNO   (ILG),  RHOSNO(ILG)                    
C                                                                       
C     * LOCAL ARRAYS                                                    
C                                                                       
      REAL SALBG(ILG,NBS), ALDIR(ILG,NBS), ALDIF(ILG,NBS),              
     +                     TRDIR(ILG,NBS), TRDIF(ILG,NBS)               
      REAL REFSNO(ILG), BCSNO(ILG)                                      
      INTEGER C_FLAG(ILG)                                               
C                                                                       
C     * CONSTANTS.                                                      
C                                                                       
      REAL WDIRCT, WDIFF                                                
      INTEGER SUM_C_FLAG                                                
C------------------------------------------------------------------     
      IF (ISNOALB .NE. 0) THEN   !4-BAND SW DISABLED 
         CALL XIT('SNOALBA',-9)                                         
      ENDIF                                                             

      IPTBAD=0                                                          
      DO 100 I=IL1,IL2                                                  
         IF(ALBSNO(I).LT.0.50.AND.ALBSNO(I).GT.0.499) ALBSNO(I)=0.50    
         IF(FSNOW(I).GT.0.0 .AND. IALS.EQ.0)              THEN          
             IF(ALBSNO(I).GT.0.70)                    THEN              
                 ALVSSN(I)=0.79*(ALBSNO(I)-0.70)+0.84                   
                 ALIRSN(I)=1.21*(ALBSNO(I)-0.70)+0.56                   
             ELSE                                                       
                 ALVSSN(I)=0.97*(ALBSNO(I)-0.50)+0.62                   
                 ALIRSN(I)=1.03*(ALBSNO(I)-0.50)+0.38                   
             ENDIF                                                      
             IF(ALVSSN(I).GT.0.999.OR.ALVSSN(I).LT.0.001) IPTBAD=I      
             IF(ALIRSN(I).GT.0.999.OR.ALIRSN(I).LT.0.001) IPTBAD=I      
         ELSE IF(FSNOW(I).GT.0.0 .AND. IALS.EQ.1)         THEN          
             ALVSSN(I)=ASVDAT(I)                                        
             ALIRSN(I)=ASIDAT(I)                                        
         ENDIF                                                          
         ALVSSC(I)=ALVSSN(I)                                            
         ALIRSC(I)=ALIRSN(I)                                            
         TRSNOWC(I)=EXP(-25.0*ZSNOW(I))                                 
  100 CONTINUE                                                          
C                                                                       
      IF(IPTBAD.NE.0) THEN                                              
         WRITE(6,6100) IPTBAD,JL,ALVSSN(IPTBAD),ALIRSN(IPTBAD)          
 6100    FORMAT('0AT (I,J)= (',I3,',',I3,'), ALVSSN,ALIRSN = ',2F10.5)  
         CALL XIT('SNOALBA',-1)                                         
      ENDIF                                                             
C                                                                       
      IF (ISNOALB .EQ. 0) THEN                                          
         DO I = IL1, IL2                                                
            ALSNO(I,1) = ALVSSN(I)                                      
            ALSNO(I,2) = ALIRSN(I)                                      
            ALSNO(I,3) = ALIRSN(I)                                      
            ALSNO(I,4) = ALIRSN(I)                                      
                                                                        
            TRSNOWG(I,1:NBS) = TRSNOWC(I)                               
         END DO ! I                                                     
C     ELSE IF (ISNOALB .EQ. 1) THEN                                     
C        DO IB = 1, NBS                                                 
C           DO I = IL1, IL2                                             
C              IF (IB .EQ. 1) THEN                                      
C                 SALBG(I,IB) = ALVSG(I)                                
C                 ALSNO(I,IB) = ALVSSN(I)                               
C              ELSE                                                     
C                 SALBG(I,IB) = ALIRG(I)                                
C                 ALSNO(I,IB) = ALIRSN(I)                               
C              END IF                                                   
C           END DO ! I                                                  
C        END DO ! IB                                                    
C        SUM_C_FLAG = 0                                                 
C        DO I = IL1, IL2                                                
C           IF (ZSNOW(I) .GT. 0.0) THEN                                 
C              C_FLAG(I) = 1                                            
C           ELSE                                                        
C              C_FLAG(I) = 0                                            
C           END IF                                                      
C           SUM_C_FLAG = SUM_C_FLAG + C_FLAG(I)                         
C        END DO ! I                                                     
                                                                        
C        IF (IALS .EQ. 0) THEN                                          
C           IF (SUM_C_FLAG .GT. 0) THEN                                 
! Convert the units of the snow grain size and BC mixing ratio          
! Snow grain size from meters to microns and BC from kg BC/m^3 to ng BC/kg SNOW 
C              DO I = IL1,IL2                                           
C                IF (C_FLAG(I) .EQ. 1) THEN                             
C                 REFSNO(I) = REFSN(I)*1.0E6                            
C                 BCSNO(I)  = (BCSN(I)/RHOSNO(I))*1.0E12                
C                END IF                                                 
C              END DO ! I                                               
                                                                        
C              CALL SNOW_ALBVAL(ALDIF, ! OUTPUT                         
C    +                          ALDIR,                                  
C    +                          CSZ,   ! INPUT                          
C    +                          SALBG,                                  
C    +                          BCSNO,                                  
C    +                          REFSNO,                                 
C    +                          SNO,                                    
C    +                          C_FLAG,                                 
C    +                          IL1,                                    
C    +                          IL2,                                    
C    +                          ILG,                                    
C    +                          NBS)                                    
C                                                                       
C              CALL SNOW_TRANVAL(TRDIF, ! OUTPUT                        
C    +                           TRDIR,                                 
C    +                           CSZ,   ! INPUT                         
C    +                           SALBG,                                 
C    +                           BCSNO,                                 
C    +                           REFSNO,                                
C    +                           SNO,                                   
C    +                           C_FLAG,                                
C    +                           IL1,                                   
C    +                           IL2,                                   
C    +                           ILG,                                   
C    +                           NBS)                                   
                                                                        
C              DO IB = 1, NBS                                           
C                 DO I = IL1, IL2                                       
C                    IF (C_FLAG(I) .EQ. 1) THEN                         
C                       WDIRCT = FSDB(I,IB)                             
C    +                         /(FSDB(I,IB)+FSFB(I,IB)+1.E-10)          
C                       WDIFF  = 1.0-WDIRCT                             
C                       ALSNO(I,IB) = ALDIF(I,IB)*WDIFF                 
C    +                              + ALDIR(I,IB)*WDIRCT                
C                       TRSNOWG(I,IB) = TRDIF(I,IB)*WDIFF               
C    +                                + TRDIR(I,IB)*WDIRCT              
C                    END IF ! C_FLAG                                    
C                 END DO ! I                                            
C              END DO ! IB                                              
C           ELSE ! SUM_C_FLAG .EQ. 0                                    
C              DO I = IL1, IL2                                          
C                 ALSNO(I,1)     = ALVSSN(I)                            
C                 ALSNO(I,2:NBS) = ALIRSN(I)                            
C                 TRSNOWG(I,1:NBS) = TRSNOWC(I)                         
C              END DO ! I                                               
C           ENDIF ! SUM_C_FLAG                                          
C        ELSE IF (IALS .EQ. 1) THEN                                     
C           DO I = IL1, IL2                                             
C              ALSNO(I,1) = ASVDAT(I)                                   
C              ALSNO(I,2) = ASIDAT(I)                                   
C              ALSNO(I,3) = ASIDAT(I)                                   
C              ALSNO(I,4) = ASIDAT(I)                                   
                                                                        
C              TRSNOWG(I,1:NBS) = TRSNOWC(I)                            
C           END DO ! I                                                  
C        END IF ! IALS                                                  
      END IF ! ISNOALB                                                  
                                                                        
      RETURN                                                            
      END                                                               
