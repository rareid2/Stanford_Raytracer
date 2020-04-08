c  ##########################################################################
c  #                                                                        #
c  #                GEOPACK-2005 to GEOPACK-2008 Adapter                    #
c  #                                                                        #
c  ##########################################################################
c
c  Adapter specifically for this project - made to avoid changes to the original code in order to get
c  updated igrf coefficients (up to 2020) and still use original functions
c
c  author: Riley R
c  date: 3/27/19
c----------------------------------------------------------------------------------
c
      SUBROUTINE IGRF_GSM (XGSM,YGSM,ZGSM,HXGSM,HYGSM,HZGSM)

      call IGRF_GSW_08 (XGSM,YGSM,ZGSM,HXGSM,HYGSM,HZGSM)

      RETURN
      END

      SUBROUTINE TSY_RECALC (IYEAR,IDAY,IHOUR,MIN,ISEC)

      real VGSEX, VGSEY, VGSEZ
      VGSEX = -400.0
      VGSEY = 0.0
      VGSEZ = 0.0
      call RECALC_08 (IYEAR,IDAY,IHOUR,MIN,ISEC,VGSEX,VGSEY,VGSEZ)

      RETURN
      END
