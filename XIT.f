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
      SUBROUTINE XIT(NAME,N)
C 
C     * OCT 01/92 - E.CHAN. (CHANGE STOP 1 TO STOP)
C     * JUN 10/91 - E.CHAN. (TRANSLATE HOLLERITH LITERALS AND 
C     *                      DIMENSION STRINGS) 
C 
C     * OCT 10/78 - J.D.HENDERSON.
C     * TERMINATES A PROGRAM BY PRINTING THE PROGRAM NAME AND 
C     * A LINE ACROSS THE PAGE FOLLOWED BY A NUMBER N. 
C 
C     * N.GE.0 IS FOR A NORMAL END. THE LINE IS DASHED. 
C     * NORMAL ENDS TERMINATE WITH   STOP.
C 
C     * N.LT.0 IS FOR AN ABNORMAL END. THE LINE IS DOTTED.
C     * IF N IS LESS THAN -100 THE PROGRAM SIMPLY TERMINATES. 
C     * OTHERWISE IF N IS LESS THAN ZERO THE PROGRAM ABORTS.
C 
      CHARACTER*(*) NAME
      CHARACTER*8   NAME8, DASH, STAR
C 
      DATA DASH /'--------'/, STAR /'********'/ 
C---------------------------------------------------------------------
C 
      NAME8 = NAME
      IF(N.GE.0) WRITE(6,6010) DASH,NAME8,(DASH,I=1,9),N 
C 
      IF(N.LT.0) WRITE(6,6010) STAR,NAME8,(STAR,I=1,9),N 
C 
      IF ( N.GE.0 .OR. N.LT.-100 ) THEN
	STOP 
      ELSE
	CALL ABORT
      ENDIF
C 
C---------------------------------------------------------------------
 6010 FORMAT('0',A8,'  END  ',A8,9A8,I8)
      END   