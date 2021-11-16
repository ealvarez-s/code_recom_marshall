C $Header: /csys/software/cvsroot/REcoM/MITgcm/recom/recom_init_fixed.F,v 1.4 2008/04/17 14:25:52 mlosch Exp $
C $Name:  $

#include "GCHEM_OPTIONS.h"
#include "RECOM_OPTIONS.h"

CBOP
C !ROUTINE: RECOM_INIT_FIXED

C !INTERFACE: ============================================================
      SUBROUTINE RECOM_INIT_FIXED( myThid )

C !DESCRIPTION:
C     Initialize fixed (not subject to adjoint) RECOM data structures
C     (to be called from S/R GCHEM_INIT_FIXED)

C !USES: ================================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#ifdef ALLOW_MNC
#include "MNC_PARAMS.h"
#endif
#ifdef ALLOW_EXF
#include "RECOM.h"
#include "RECOM_EXF.h"
#include "cal.h"
#include "EXF_PARAM.h"
#include "EXF_CONSTANTS.h"
#endif 
#ifdef RECOM_ATMOSPCO2_HISTORY
#include "RECOM_PARAMS.h"
#endif

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      INTEGER myThid

C !OUTPUT PARAMETERS: ==================================================
C  none

#ifdef ALLOW_RECOM

C !LOCAL VARIABLES: ====================================================
      LOGICAL useVariableK
#ifdef ALLOW_EXF
      INTEGER date_array(4), difftime(4), yearStartDate(4)
#endif
#ifdef RECOM_ATMOSPCO2_HISTORY
      INTEGER k, iUnit
#endif
CEOP

C
C     initialise additional output variables
C
#ifdef ALLOW_MNC     
      IF (useMNC) THEN
       CALL MNC_CW_ADD_VNAME('WCtave', 'Cen_xy_Hn__C__t', 4,5,myThid)
       CALL MNC_CW_ADD_VATTR_TEXT('WCtave','units','mmolC/m^2/s', 
     &      myThid)
       CALL MNC_CW_ADD_VNAME('WNtave', 'Cen_xy_Hn__C__t', 4,5, myThid)
       CALL MNC_CW_ADD_VATTR_TEXT('WNtave','units','mmolN/m^2/s',  
     &      myThid)
CCV#ifdef ALLOW_RECOM_SILICATE
       CALL MNC_CW_ADD_VNAME('WSitave', 'Cen_xy_Hn__C__t', 4,5, myThid)
       CALL MNC_CW_ADD_VATTR_TEXT('WSitave','units','mmolSi/m^2/s',
     &      myThid)
CCV#endif /* ALLOW_RECOM_SILICATE */
       useVariableK = useKPP .OR. usePP81 .OR. useMY82 .OR. useGGL90
     &      .OR. useGMredi .OR. ivdc_kappa.NE.0.
       IF ( useVariableK ) THEN
        CALL MNC_CW_ADD_VNAME('Cdiftave', 'Cen_xy_Hn__C__t',4,5,myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('Cdiftave','units','mmolC/m^2/s',
     &       myThid)
        CALL MNC_CW_ADD_VNAME('Ndiftave', 'Cen_xy_Hn__C__t',4,5,myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('Ndiftave','units','mmolN/m^2/s',
     &       myThid)
CCV#ifdef ALLOW_RECOM_SILICATE
        CALL MNC_CW_ADD_VNAME('Sidiftave', 'Cen_xy_Hn__C__t',4,5,myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('Sidiftave','units','mmolSi/m^2/s',
     &       myThid)
CCV#endif /* ALLOW_RECOM_SILICATE */

        CALL MNC_CW_ADD_VNAME('BenC','Cen_xy_Hn__-__t',3,4,myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('BenC','units','mmol/m^3', myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('BenC','description',
     &       'benthic carbon', myThid)
        CALL MNC_CW_ADD_VNAME('BenN','Cen_xy_Hn__-__t',3,4,myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('BenN','units','mmol/m^3', myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('BenN','description',
     &       'benthic nitrogen', myThid)
CCV#ifdef ALLOW_RECOM_SILICATE
        CALL MNC_CW_ADD_VNAME('BenSi','Cen_xy_Hn__-__t',3,4,myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('BenSi','units','mmol/m^3', myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('BenSi','description',
     &       'benthic silicate', myThid)
CCV#endif /* ALLOW_RECOM_SILICATE */

       ENDIF
      ENDIF
#endif /* ALLOW_MNC */

#ifdef ALLOW_CAL
      IF ( recom_pCO2file .NE. ' ' .and. pCO2period .gt. 0. ) THEN
       call cal_FullDate  ( pCO2startdate1,   pCO2startdate2,
     &      date_array, mythid )
       call cal_TimePassed(modelstartdate,date_array,difftime,mythid)
       call cal_ToSeconds ( difftime, pCO2startdate,   mythid )
       pCO2startdate    = modelstart  + pCO2startdate
      ENDIF
      IF ( recom_ironfile .NE. ' ' .and. fedustperiod .gt. 0. ) THEN
       call cal_FullDate  ( feduststartdate1, feduststartdate2,
     &      date_array, mythid )
       call cal_TimePassed(modelstartdate,date_array,difftime,mythid)
       call cal_ToSeconds ( difftime, feduststartdate, mythid )
       feduststartdate    = modelstart  + feduststartdate
      ENDIF
#endif /* ALLOW_CAL */

#ifdef ALLOW_EXF
      IF ( useEXF ) THEN
      _BEGIN_MASTER( myThid )
      IF ( recom_pCO2file .NE. ' ' .AND. pCO2period .GT. 0. ) THEN
       CALL cal_FullDate  ( pCO2startdate1, pCO2startdate2,
     &      date_array, myThid )
       IF ( useExfYearlyFields ) THEN
        yearStartDate(1) = int(date_array(1)/10000.) * 10000 + 101
        yearStartDate(2) = 0
        yearStartDate(3) = date_array(3)
        yearStartDate(4) = date_array(4)
        CALL cal_TimePassed( yearStartDate,date_array,difftime,myThid)
        CALL cal_ToSeconds ( difftime,   pCO2startdate      ,myThid)
       ELSE
        CALL cal_TimePassed(modelstartdate,date_array,difftime,myThid)
        CALL cal_ToSeconds ( difftime,   pCO2startdate      ,myThid)
        pCO2startdate  =  modelstart + pCO2startdate
       ENDIF
      ENDIF
      IF ( recom_ironfile .NE. ' ' .AND. fedustperiod .GT. 0. ) THEN
       CALL cal_FullDate  ( feduststartdate1, feduststartdate2,
     &      date_array, myThid )
       IF ( useExfYearlyFields ) THEN
        yearStartDate(1) = int(date_array(1)/10000.) * 10000 + 101
        yearStartDate(2) = 0
        yearStartDate(3) = date_array(3)
        yearStartDate(4) = date_array(4)
        CALL cal_TimePassed( yearStartDate,date_array,difftime,myThid)
        CALL cal_ToSeconds ( difftime,   feduststartdate      ,myThid)
       ELSE
        CALL cal_TimePassed(modelstartdate,date_array,difftime,myThid)
        CALL cal_ToSeconds ( difftime,   feduststartdate      ,myThid)
        feduststartdate  =  modelstart + feduststartdate
       ENDIF
      ENDIF
      _END_MASTER( myThid )
      _BARRIER
      ENDIF
#endif /* ALLOW_EXF */

#ifdef ALLOW_DIAGNOSTICS
C     Define diagnostics Names :
      IF ( useDiagnostics ) THEN
        CALL RECOM_DIAGNOSTICS_INIT( myThid )
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

#ifdef RECOM_ATMOSPCO2_HISTORY
      _BEGIN_MASTER( myThid )
C     Read in a history of atmopheric pCO2
      IF ( recom_pco2_int1.EQ.2 ) THEN
        CALL MDSFINDUNIT( iUnit, mythid )
        OPEN(UNIT=iUnit,FILE='co2atmos.dat',STATUS='old')
        DO k=1,recom_pco2_int2
          READ(iUnit,*) co2atmos(k,1),co2atmos(k,2)
          WRITE(standardMessageUnit,*) 'co2atmos',
     &         co2atmos(k,1),co2atmos(k,2)
        ENDDO
        CLOSE(iUnit)
      ENDIF
      _END_MASTER( myThid )
      _BARRIER
#endif /* RECOM_ATMOSPCO2_HISTORY */

#endif /* ALLOW_RECOM */

      RETURN
      END
