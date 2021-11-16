C $Header: /csys/software/cvsroot/REcoM/MITgcm/recom/RECOM_EXF.h,v 1.4 2008/04/17 14:25:52 mlosch Exp $
C $Name:  $

#ifdef ALLOW_EXF

      NAMELIST /RECOM_PARM02/
     &     pco2startdate1, pco2startdate2, pco2period,
     &     pco2repeatperiod,
     &     feduststartdate1, feduststartdate2, fedustperiod,
     &     fedustrepeatperiod
#ifdef USE_EXF_INTERPOLATION
     &     ,pco2_lon0, pco2_lon_inc, pco2_lat0, pco2_lat_inc,
     &     fedust_lon0, fedust_lon_inc, fedust_lat0, fedust_lat_inc,
     &     pco2_nlon, pco2_nlat, fedust_nlon, fedust_nlat
#endif /* USE_EXF_INTERPOLATION */

      INTEGER pco2startdate1
      INTEGER pco2startdate2
      _RL     pco2startdate
      _RL     pco2period
      _RL     pco2repeatperiod

      INTEGER feduststartdate1
      INTEGER feduststartdate2
      _RL     feduststartdate
      _RL     fedustperiod
      _RL     fedustrepeatperiod

      COMMON /RECOM_EXF_PARMS_I/
     &     pco2startdate1, pco2startdate2,
     &     feduststartdate1, feduststartdate2
      COMMON /RECOM_EXF_PARMS_R/
     &     pco2startdate, pco2period, pco2repeatperiod,
     &     feduststartdate, fedustperiod, fedustrepeatperiod

#ifdef USE_EXF_INTERPOLATION
      _RL pco2_lon0, pco2_lon_inc
      _RL pco2_lat0, pco2_lat_inc(MAX_LAT_INC)
      INTEGER pco2_nlon, pco2_nlat
      _RL fedust_lon0, fedust_lon_inc
      _RL fedust_lat0, fedust_lat_inc(MAX_LAT_INC)
      INTEGER fedust_nlon, fedust_nlat
      COMMON /RECOM_INTERPOLATION/
     & pco2_lon0, pco2_lon_inc,
     & pco2_lat0, pco2_lat_inc,
     & pco2_nlon, pco2_nlat,
     & fedust_lon0, fedust_lon_inc,
     & fedust_lat0, fedust_lat_inc,
     & fedust_nlon, fedust_nlat
#endif /* USE_EXF_INTERPOLATION */


#endif /* ALLOW_EXF */
