"""

function to generate config file for ngo model

inputs: kp, config file name
could make more....

"""

import os
import shutil
import datetime
import random
import string
import time
import numpy as np

# Apologies for the terrible naming conventions.
# This is a re-spin of rayMaker.c -- generates the input
# cards to the fortran raytracer (newray1.8.f), etc

# see ngo_infile_old.in for original

# Directory to write to:

# Open the card file for writing
def gen_ngo_defaults(Kp, configfile):

  print('Kp is set at = ', Kp)

  with open(configfile, "w") as card:

    # // -------------------   CARD 1:  ---------------------------
    #   // Interactive: if n>0, enter NUMRES and NSUPPR interactively and
    #   // override the existing values.  If n<=0 use existing values.
    INTERA    = 0;

    # // Maximum number of raytracing continuations at the resonant
    # // cone. If zero, then no continuation.
    NUMRES    = 10;

    # // values 1/0. If =1, then satellite trajectory in the range
    # // lambda_min-3 <= lambda <= lambda_max+3 is printed.
    NSUPPR    = 0;

    # // Special latitude: if SPELAT ~= 0, pogram tries to compute
    # // and print the output at the latitude lambda = 10 + eps
    SPELAT    = 0;

    card.write("%g %g %g %g\n"%(INTERA, NUMRES, NSUPPR, SPELAT))
    # printf(filePtr,"%g %g %g %g\n", INTERA, NUMRES, NSUPPR, SPELAT);


    # // ---------------- CARD 2:  Satellite info ------------------
    # // Geoccentric distance of Satellite
    DISTRE = 22400;

    # // Latitude of satellite (Don't confuse with LATITU of ray)
    LATITUD = 0;

    card.write("%g %g\n"%(DISTRE, LATITUD))

    # // ------------------  CARD 3: Separator------------------------
    # // This just seems to be a kind of separator card
    card.write("%g %g \n"%(-1.0, 0))


    # // ------------------- CARD 4: Model input -------------------
    # // Number of components in the plasma - number from 2-4, we
    # // have H, He, O, and el.'s
    NUM = 4;

    # // Is 0/1. If 0, integration size automatically adjusted to
    # // maintain accuracy specified by ABSB and RELB.  If 1,
    # // stepsize held constant
    KSKIP = 0;

    # // Is -1/1. If -1 proton whistler mode, if 1, el whistler mode
    MODE  = 1;

    # // Spacing of points on line printer and (if KTAPE~=0) tape.
    # // If 0 - all points outputted on printer, otherwise, only
    # // every nth point is outputted. Critical points always printed
    # // i.e <0 critical points only, =0 all points, >0 critical
    # // and spaced ponits.
    KOUNT = 2;

    # // Indicates number of ducts in density model.  We usually have
    # // 1 or 2, because plasmapause is a duct and we had another one
    # // to make the steepness of electron dropoff accurate.
    # change 2 to 4
    KDUCTS = 4;

    # // File number of results recorded on tape.  If =0, only
    # // lineprinter, if =1, the first file is the program itself.
    KTAPE = 1;

    # // No description: maybe reference altitude of some sort?
    REFALT =  200;  # (looks to be unused -- aps 12.2016)

    # // Density reference point - range in L shell
    # change 2.5 to 2
    DSRRNG  = 2;

    # // Density reference point - latitude in degrees
    DSRLAT  = 0;

    # // Density reference point - el/cc @ reference point
    # changed from 2500 to 2000
    DSRDENS = 2000 # 2000 #2500;

    card.write("%g %g %g %g %g %g %g %g %g %g\n"%
      (NUM, KSKIP, MODE, KOUNT, KDUCTS, KTAPE, REFALT,
      DSRRNG, DSRLAT, DSRDENS ))

    # // ---------------------  CARD 5: Model Input ----------------
    # // electron gyrofrequency in kHz, at surface of the earth, at
    # // equator, for dipole magnetic field model
    EGFEQ = 880;  # 30 microtesla * (Q_E/M_E) ~ 840 kHz

    # // Temperature used in the diffusive equilibrium density model
    # changed to 11600 from 1600 (not sure why)
    THERM = 11600 #1600;  1 ev = 1.602e-19 / (boltzman const, 1.38e-23)

    # // Initial integration stepsize H = 50*HM/sqrt(f), except below
    # // RDIV, where H is 1/16 of this value
    HM  = 10;

    # // Absolute error in integration.  When any of the 5 variables
    # // exceeds 14.2*ABSB and RELB, stepsize is halved.
    ABSB  = 0.0001;

    # // Relative error in integration - absolute error divided by the
    # // value. When relative error exceeds 14.2*RELB and absolute
    # // error, stepsize is halved.
    RELB  = 1e-6;

    card.write("%g %g %g %g %g\n"%(EGFEQ, THERM, HM, ABSB, RELB))

    # // ----------- CARD 6: Diffusive Equilirium input  -----------
    # // Geocentric base (in km) of DE density model
    RBASE = 6580; # changed from 7370

    # // Electron density (cm^-3) at RBASE.  Final density here could
    # // actually be different due to other things such as ducts, etc.
    ANE0  = 1

    # // relative concentration of H+ at RBASE, i.e. [H+]/[H+ + He+ + O+]
    ALPHA0_2= 0.08;

    # // Relative concnetration of He+ at RBASE
    ALPHA0_3= 0.02;

    # // Relative concentration of O+ at RBASE
    ALPHA0_4= 0.9;

    card.write("%g %g %g %g %g\n"%(RBASE, ANE0, ALPHA0_2, ALPHA0_3, ALPHA0_4))

    # // ---------------- CARD 7: Model Input ------------------------
    # // Geocentric distance of the bottom of the ionosphere, where
    # // density =0 (km)
    RZERO = 6460;

    # // Scale height at the bottom side of the lower ionosphere (km)
    SCBOT = 140;

    # // Geocentric distance below which raytracing stops (km)
    RSTOP = 6470;   # 100km

    # // Geocentric distance below which stepsize reset to 1/16 of
    # // normal starting value  (km)
    RDIV  = 6873;

    # // Minimum allowed value of integration stepsize
    HMIN  = 1e-8; #0.001;

    card.write("%g %g %g %g %g\n"%(RZERO, SCBOT, RSTOP, RDIV, HMIN))

    # // ---------------  CARD 8: Plasmapause input ---------------
    # // L-value of inner edge of the plasmapause, center of the
    # // knee is at LK + DDK
    # // (LK currently defined in consts.h -- aps 11.2015)

    # LK  = 5.55;
    # //LK = 4.61;  // KP = 2
    # //LK = 3.69;  // KP = 4
    # // Moldwin 2002 model of the plasmapause location
    LK =  5.39 - 0.382*Kp;
    # LK = 5.39;


    # // Exponential component of density decrease beyond the knee,
    # // R^{-EXPK}.  This is in addition to decrease due to DE model.
    EXPK  = 0.13 #0.13;

    # // Halfwidth in L of the knee
    DDK =0.2 #0.1 #0.2 # 0.07 #0.07;

    # // Geocentric distance (km) where density outside the knee is
    # // equal to the density inside (???).  Modeil is not good
    # // below RCONSS outside of the knee.
    RCONSN = 1e-8;

    # // Scale height of the radial density decrease above RCONSN
    # // above the knee.
    SCR = 50;

    card.write("%g %g %g %g %g\n"%(LK, EXPK, DDK, RCONSN, SCR))




    # // ---------- CARD 9: Duct input (optional, repeatable) ----------
    # // repeat for as many ducts as we want to include

    # // L-value of the center of the duct
    L0    = 1.1;

    # // Enhancement factor of the duct
    DEF   = 30.0;

    # // Duct half-width in L.  If -ve, then sinusoidal perturbation,
    # // and DD is the period in L of the sinusoidal perturbation
    # changed from 1.0 to 4.0
    DD    = 4.0;

    # // Geocentric radius (km) of the power end of the duct.  Below
    # // this height, the duct merges into the background plasma. If
    # // extends to earth, put RDUCLN = 6370 - In the NORTH
    RDUCLN  =  10370;

    # // Radial scale height of duct, with which the LOWER end merges
    # // into the background plasma
    HDUCLN  = 20000 #20;

    # // Geocentric radius of the upper end of the duct.  If it extends
    # // to equator, set = 6370*L0 - In the NORTH
    RDUCUN  = 70081;

    # // Radial scale height with which upper end of duct merges into
    # // the background plasma
    HDUCUN  = 20000;

    # // Lower radius of duct in the SOUTH
    RDUCLS  = 10370;

    # // Lower scale height of duct in the South
    HDUCLS  = 20000;

    # // Upper radius of duct in the South
    RDUCUS  = 70081;

    # // Upper scale height of duct in the South
    HDUCUS  = 20000;

    # // Allows one-sided ducts.  If >0, all lower L values below L0
    # // stay at the enhancement factor DEF, like a secondary
    # // plasmapause. If =0, then is a normal 2 sides duct, if <0,
    # // all higher L values above L0 stay at the enhancement factor
    # // DEF, like a recovery region.
    SIDEDU  = 1;

    # // These 3 ducts seem to be important 1-sided shelves
    # to make the density curve line up correctly.

    card.write("%g %g %g %g %g %g %g %g %g %g %g %g \n"%
      (L0,  DEF,  DD,  RDUCLN,  HDUCLN,
      RDUCUN,  HDUCUN,  RDUCLS,  HDUCLS,
      RDUCUS,  HDUCUS,  SIDEDU ))
    card.write("%g %g %g %g %g %g %g %g %g %g %g %g \n"%
      (L0,  DEF,  DD,  RDUCLN,  HDUCLN,
      RDUCUN,  HDUCUN,  RDUCLS,  HDUCLS,
      RDUCUS,  HDUCUS,  SIDEDU ))
    card.write("%g %g %g %g %g %g %g %g %g %g %g %g \n"%
      (LK,  1.8,  1.8,  RDUCLN,  HDUCLN,
      RDUCUN,  HDUCUN,  RDUCLS,  HDUCLS,
      RDUCUS,  HDUCUS,  SIDEDU ))

    # # Hey Austin, here's an example of a duct that works.
    # If adding more, don't forget to go back up and set KDUCTS!

    # card.write("%g %g %g %g %g %g %g %g %g %g %g %g \n"%
    #   (4,  100,  0.05,  RDUCLN,  HDUCLN,
    #   RDUCUN,  HDUCUN,  RDUCLS,  HDUCLS,
    #   RDUCUS,  HDUCUS,  0 ))

    # // ------------------  CARD 10: Profile Input --------------------
    # // Altitude step in Earth radii between points on radial profile,
    # // =0 for no profile
    # //  (This is what sets the step-size in Lprofile.dat)
    PSTALT  = 0.05;

    # // Geocentric distance (Earth radii) of lower limit of radial
    # // profile
    PALT1 = 1;

    # // Geocentric distance of upper limit of radial profile
    PALT2 = 7;

    # // Constant latitude (in degrees) of radial profile
    PLATIT  = 0;

    # // Latitude step (in degrees) between points on teh latitude
    # // profile. =0 for no profile
    PSTLAT  = 1;

    # // Lower latitude limit of profile (in degrees)
    PLAT1 = -60;

    # // Upper latitude limit of profile (in degrees)
    PLAT2 = 60;

    # // Constant altitude at which lat profile is to be computed
    PALTIT  = 1000;

    card.write("%g %g %g %g %g %g %g %g\n"%
      (PSTALT, PALT1, PALT2, PLATIT, PSTLAT, PLAT1, PLAT2, PALTIT))

  print('finished writing Ngo config file')

  return