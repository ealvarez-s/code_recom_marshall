C $Header: /csys/software/cvsroot/REcoM/MITgcm/recom/recom_init_vari.F,v 1.11 2008/04/17 14:25:52 mlosch Exp $
C $Name:  $

#include "CPP_OPTIONS.h"
C#include "PTRACERS_OPTIONS.h"
#include "GCHEM_OPTIONS.h"
#include "RECOM_OPTIONS.h"
#ifdef ALLOW_EXF
# include "EXF_OPTIONS.h"
#endif

CBOP
C !ROUTINE: RECOM_INIT_VARI

C !INTERFACE: ============================================================
      SUBROUTINE RECOM_INIT_VARI( myThid )

C !DESCRIPTION:
C     Initialize RECOM data structures that could potentially be subject
C     to adjoint computations (to be called from S/R GCHEM_INIT_VARI)

C !USES: ================================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#ifdef ALLOW_MNC
#include "MNC_PARAMS.h"
#endif
#ifdef ALLOW_PTRACERS
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "RECOM.h"
#endif /* ALLOW_PTRACERS */
#ifdef ALLOW_EXF
# include "EXF_PARAM.h"
# include "RECOM_EXF.h"
#endif /* ALLOW_EXF */
CTW increased sinking speed
#ifdef ALLOW_SINK_INCREASE
#include "RECOM_PARAMS.h"
#endif /* ALLOW_SINK_INCREASE */

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      INTEGER myThid

C !OUTPUT PARAMETERS: ==================================================
C  none

#ifdef ALLOW_RECOM

C !LOCAL VARIABLES: ====================================================
C     i,j,k,bi,bj - Loop counters

      INTEGER i, j, k, bi, bj
      LOGICAL useVariableK
#ifdef ALLOW_EXF
      INTEGER interp_method
      PARAMETER (interp_method=1)
      _RL loc0, loc1
#endif /* ALLOW_EXF */
CEOP

C
C     initialise specific forcing fields
C      
      K = 1
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO J=1-Oly,sNy+Oly
         DO I=1-Olx,sNx+Olx
          dicFlux  (i,j,bi,bj)      = 0. _d 0
          windSpeed(i,j,bi,bj)      = 0. _d 0
          atmospCO2(i,j,bi,bj)      = 0. _d 0
          atmospCO2(i,j,bi,bj)      = 380. _d 0 * maskC(I,J,K,bi,bj)
          feDustDeposit(i,j,bi,bj)  = 0. _d 0
C     initialise hPlus with some sensible values here
          hPlus(i,j,bi,bj)          = 1. _d -08 * maskC(I,J,K,bi,bj)
C     
          benthicLayerC(i,j,bi,bj)  = 0. _d 0
          benthicLayerN(i,j,bi,bj)  = 0. _d 0
CCV#ifdef ALLOW_RECOM_SILICATE
          benthicLayerSi(i,j,bi,bj) = 0. _d 0
CCV#endif /* ALLOW_RECOM_SILICATE */
          benthicLayerCaCO3(i,j,bi,bj) = 0. _d 0
         ENDDO
        ENDDO
C
C     initialise total matter fields (diagnostic)
C
        DO K=1, Nr
         DO J=1,sNy
          DO I=1,sNx
           recomC(i,j,k,bi,bj)  = 0. _d 0
           recomN(i,j,k,bi,bj)  = 0. _d 0
           recomSi(i,j,k,bi,bj) = 0. _d 0
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
C
CTW  calculate increased sinking speed
C
#ifdef ALLOW_SINK_INCREASE
      DO k=1,Nr
       Vdetfast(k)= Vdet_a*abs(rF(k))+Vdet_b
       print *, 'layer & vdetfast',k,Vdetfast(k)
      ENDDO
#endif /* ALLOW_SINK_INCREASE */

#ifdef ALLOW_EXF
      loc0 = 0. _d 0
      loc1 = 1. _d 0
      IF ( recom_windFile .NE. ' '  ) THEN
       call exf_init_gen (
     &     recom_windfile, wspeedperiod, loc1, wspeedmask,
     &     wspeedconst, windSpeed, windSpeed0, windSpeed1, 
#ifdef USE_EXF_INTERPOLATION
     &     wspeed_lon0, wspeed_lon_inc, 
     &     wspeed_lat0, wspeed_lat_inc,
     &     wspeed_nlon, wspeed_nlat, xC, yC, interp_method,
#endif
     &     mythid )
      ENDIF
      IF ( recom_pCO2File .NE. ' ' ) THEN
       call exf_init_gen (
     &     recom_pCO2file, pCO2period, loc1, sfluxmask,
     &     loc0, atmospCO2, atmosp0, atmosp1, 
#ifdef USE_EXF_INTERPOLATION
     &     pCO2_lon0, pCO2_lon_inc, 
     &     pCO2_lat0, pCO2_lat_inc,
     &     pCO2_nlon, pCO2_nlat, xC, yC, interp_method,
#endif
     &     mythid )
      ENDIF
CCV#ifdef ALLOW_RECOM_FEINPUT
      IF ( recom_ironFile .NE. ' ' ) THEN
       call exf_init_gen (
     &     recom_ironFile, fedustperiod, loc1, sfluxmask,
     &     loc0, fedustDeposit, fedustDeposit0, fedustDeposit1, 
#ifdef USE_EXF_INTERPOLATION
     &     fedust_lon0, fedust_lon_inc, 
     &     fedust_lat0, fedust_lat_inc,
     &     fedust_nlon, fedust_nlat, xC, yC, interp_method,
#endif
     &     mythid )
      ENDIF
CCV#endif /* ALLOW_RECOM_FEINPUT */
#else
      _BEGIN_MASTER(myThid)
      IF ( recom_WindFile .NE. ' '  ) THEN
       CALL READ_FLD_XY_RS( recom_WindFile, ' ',windSpeed, 
     &      0, myThid )
      ENDIF
      IF ( recom_pCO2File .NE. ' '  ) THEN
       CALL READ_FLD_XY_RS( recom_pCO2File, ' ',atmospCO2,
     &      0, myThid )
      ENDIF
      IF ( recom_ironFile .NE. ' '  ) THEN
       CALL READ_FLD_XY_RS( recom_ironFile, ' ',feDustDeposit,
     &      0, myThid )
      ENDIF
      _END_MASTER(myThid)
#endif /* ALLOW_EXF */
C
      _EXCH_XY_RS(windSpeed    , myThid )
      _EXCH_XY_RS(atmospCO2    , myThid )
      _EXCH_XY_RS(feDustDeposit, myThid )
      _EXCH_XY_RS(hPlus        , myThid )

#ifdef RECOM_IRON_HYDROTHERMAL
C     if a file for hydrothermal iron input is defined, read it in; else
C     initialize it to zeros
      _BEGIN_MASTER(myThid)
      IF ( recom_hydroFile .NE. ' '  ) THEN
         CALL READ_FLD_XYZ_RL(recom_hydroFile,' ',
     &        feHydrothermal(1-OLx,1-OLy,1,1,1),0,myThid)
      ENDIF
      _END_MASTER(myThid)

C shifted the exchange routine outside of master thread!
      _EXCH_XYZ_RL(feHydrothermal , myThid)

C     Apply mask
      DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
            DO k=1,Nr
               DO j=1-OLy,sNy+OLy
                  DO i=1-OLx,sNx+OLx
                     IF (maskC(i,j,k,bi,bj).EQ.0.)
     &           feHydrothermal(i,j,k,bi,bj)=0. _d 0
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
#endif /* RECOM_IRON_HYDROTHERMAL */					 

C Read benthic layer variables from a pickup file if needed
C (analogous to ptracers_read_checkpoint)
      IF (nIter0.GT.PTRACERS_Iter0) THEN
       CALL RECOM_READ_PICKUP( nIter0,myThid )
      ENDIF
#endif /* ALLOW_RECOM */

      RETURN
      END
