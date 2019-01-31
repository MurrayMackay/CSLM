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
      SUBROUTINE EQNST(EXPW,RHO,TCEL,H)
C======================================================================
C     * AUG  2/16 - M.MACKAY.  	Equation of state moved to subroutine
C
      IMPLICIT NONE
C
C ----* INPUT FIELDS *------------------------------------------------
C
      REAL EXPW,RHO,TCEL
      REAL GRAV,H,C0,ALPHA,P,BETA,RHO0,T0,T,Z,XI,KAPPA,S
C
C SALINITY AND PRESSURE EFFECTS NOT INCLUDED 
c
      Z=0.0
      S=0.0
Cmdm  S=0.3
Cmdm  S=0.014
Cmdm  S=0.03
Cmdm  Z=MIN(10.0,H)
C-------------------------------------------------------------------------
      GRAV=9.80616
      C0=4.9388E-5
      ALPHA=3.3039E-7
      BETA=8.2545E-6
      XI=8.0477E-1
      KAPPA=0.2151
Cmdm  T0=3.9816
Cmdm  T0=3.9839
      T0=3.98275
      RHO0=999.975 + XI*S
      P=RHO0*GRAV*Z*1.0E-5
      T=TCEL-T0-kappa*S
C======================================================================
C Farmer and Carmack, 1981
C
       RHO= RHO0*(1. + P*(C0-ALPHA*T) - BETA*T*T )
       EXPW = (2.*RHO0*BETA*T + RHO0*P*ALPHA )/RHO

      RETURN
      END
