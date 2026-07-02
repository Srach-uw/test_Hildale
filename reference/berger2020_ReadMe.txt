J/AJ/159/280 Gaia-Kepler stellar properties catalog.I. KIC stars (Berger+, 2020)
================================================================================
The Gaia-Kepler stellar properties catalog.
I. Homogeneous fundamental properties for 186301 Kepler stars.
    Berger T.A., Huber D., Van Saders J.L., Gaidos E., Tayar J., Kraus A.L.
   <Astron. J., 159, 280-280 (2020)>
   =2020AJ....159..280B    (SIMBAD/NED BibCode)
================================================================================
ADC_Keywords: Stars, giant; Infrared; Optical; Parallaxes, trigonometric;
              Abundances, [Fe/H]; Stars, masses; Stars, diameters; Stars, ages;
              Effective temperatures
Keywords: Catalogs; Fundamental parameters of stars; Exoplanet systems

Abstract:
    An accurate and precise Kepler Stellar Properties Catalog is essential
    for the interpretation of the Kepler exoplanet survey results.
    Previous Kepler Stellar Properties Catalogs have focused on reporting
    the best-available parameters for each star, but this has required
    combining data from a variety of heterogeneous sources. We present the
    Gaia-Kepler Stellar Properties Catalog, a set of stellar properties of
    186301 Kepler stars, homogeneously derived from isochrones and
    broadband photometry, Gaia Data Release 2 parallaxes, and
    spectroscopic metallicities, where available. Our photometric
    effective temperatures, derived from g to Ks colors, are calibrated on
    stars with interferometric angular diameters. Median catalog
    uncertainties are 112K for Teff, 0.05dex for logg, 4% for R_*_, 7% for
    M_*_, 13% for {rho}_*_, 10% for L_*_, and 56% for stellar age. These
    precise constraints on stellar properties for this sample of stars
    will allow unprecedented investigations into trends in stellar and
    exoplanet properties as a function of stellar mass and age. In
    addition, our homogeneous parameter determinations will permit more
    accurate calculations of planet occurrence and trends with stellar
    properties.

Description:
    In this paper, we utilize Gaia DR2 parallaxes, homogeneous stellar g
    and Ks photometry, and spectroscopic metallicities, where available,
    to improve on previous analyses and present the most accurate,
    homogeneous, and precise analysis of stars in the Kepler field. We
    re-derive stellar Teff, logg, radii, masses, densities, luminosities,
    and ages for 186301 Kepler targets, and investigate the stellar
    properties of a number of noteworthy Kepler exoplanet-hosting stars.

File Summary:
--------------------------------------------------------------------------------
 FileName    Lrecl    Records    Explanations
--------------------------------------------------------------------------------
ReadMe          80          .    This file
table1.dat      95     186548    Gaia-Kepler Stellar input parameters
table2.dat     222     186301    Gaia-Kepler Stellar output parameters
--------------------------------------------------------------------------------

See also:
 II/246  : 2MASS All-Sky Catalog of Point Sources (Cutri+ 2003)
 VII/233 : The Two Micron All Sky Survey: Extended sources (Skrutskie+, 2006)
 V/133   : Kepler Input Catalog (Kepler Mission Team, 2009)
 I/345   : Gaia DR2 (Gaia Collaboration, 2018)
 J/AJ/142/112  : KIC photometric calibration (Brown+, 2011)
 J/AJ/144/24   : The Kepler-INT survey (Greiss+, 2012)
 J/ApJS/211/2  : Revised stellar properties of Kepler targets (Huber+, 2014)
 J/ApJ/809/25  : Stellar and planet properties for K2 candidates (Montet+, 2015)
 J/A+A/588/A87 : Seismic global parameters of 6111 KIC (Vrard+, 2016)
 J/ApJS/224/2  : K2 EPIC stellar properties for 138600 targets (Huber+, 2016)
 J/AJ/154/107  : California-Kepler Survey (CKS). I. 1305 stars (Petigura+, 2017)
 J/ApJ/844/102 : KIC star plxs from asteroseismology vs Gaia (Huber+, 2017)
 J/ApJS/229/30 : Revised stellar properties of Q1-17 Kepler (Mathur+, 2017)
 J/ApJ/866/99  : Radii of Kepler stars & planets using Gaia DR2 (Huber+, 2018)

Byte-by-byte Description of file: table1.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label   Explanations
--------------------------------------------------------------------------------
   1-  8 I8     ---     KIC     Kepler Input Catalog identifier
  10- 16 F7.4   mag     gmag    [6.24/22.1] g band magnitude (1)
  18- 23 F6.4   mag   e_gmag    [0.02/0.18] Uncertainty in gmag
  25- 31 F7.4   mag     Ksmag   [4.33/21.9] Ks band magnitude (2)
  33- 37 F5.3   mag   e_Ksmag   [0.01/0.11] Uncertainty in Ksmag
  39- 45 F7.4   mas     plx     [0.001/83.6] Gaia DR2 Parallax
  47- 52 F6.4   mas   e_plx     [0.01/1.3] Uncertainty in Par
  54- 59 F6.3   [-]     [Fe/H]  [-2.56/0.75]? Spectroscopic metallicity (3)
  61- 64 F4.2   [-]   e_[Fe/H]  [0.15]? Uncertainty in [Fe/H]
  66- 72 F7.4   ---     RUWE    [0.63/73.5] Re-normalized unit-weight error (4)
  74- 74 I1     ---     Ncomp   [0/6]? Number, Gaia DR2 companions within 4"
  76- 89 A14    ---     KsCorr  Potential corrections compared to 2MASS Ks (5)
  91- 95 A5     ---     State   [RGB clump] Evolutionary State (6)
--------------------------------------------------------------------------------
Note (1): g-band photometry from either the Kepler Input Catalog
    (Brown+, 2011 J/AJ/142/112) or the Kepler-INT Survey
    (Greiss+, 2012, J/AJ/144/24) corrected as in Section 2.2.
Note (2): 2MASS Ks band photometry with corrections as defined in Section 2.2.
Note (3): Spectroscopic metallicity from the Kepler Stellar Properties
    Catalog (Mathur+, 2017, J/ApJS/229/30), the California-Kepler Survey
    (Petigura+, 2017, J/AJ/154/107), APOGEE DR14
    (Abolfathi+, 2018ApJS..235..42A), and/or LAMOST DR5
    (Ren+, 2018MNRAS.477.4641R).
Note (4): RUWE is the magnitude and color independent re-normalization
    of the astrometric {chi}2 of Gaia DR2 (unit-weight error or UWE).
Note (5): Corrections as follows:
    BinDet_NoCorr = Kepler star with one Gaia resolved companion but no
                     correction (10173 occurrences)
    BinaryCorr = Kepler star with one Gaia resolved companion with a correction
                  (19704 occurrences)
    TerDet_NoCorr = Multiple resolved companions; no correction
                     (376 occurrences)
    TerDet_BinCorr = Multiple resolved companions; correction based on 1
                      companion (1639 occurrences)
    TertiaryCorr = Multiple resolved companions; correction based on multiple
                    companions (1582 occurrences)
Note (6): Red giant evolutionary state flags (6956 RGB or 7840 clump (for Red
    Clump)), Vrard+, 2016, J/A+A/588/A87; Hon+, 2018MNRAS.476.3233H.
--------------------------------------------------------------------------------

Byte-by-byte Description of file: table2.dat
--------------------------------------------------------------------------------
   Bytes Format Units Label  Explanations
--------------------------------------------------------------------------------
   1-  8 I8     ---       KIC    Kepler Input Catalog identifier
  10- 14 F5.3   Msun      Mass   [0.1/4.99] Isochrone derived stellar mass
  16- 20 F5.3   Msun    E_Mass   [0.001/2.2] Upper error on Mass
  22- 27 F6.3   Msun    e_Mass   [-2.23/-0.001] Lower error on Mass
  29- 35 F7.1   K         Teff   [2475/19044] Isochrone derived effective
                                  temperature
  37- 42 F6.1   K       E_Teff   [1.6/4960] Upper error on Teff
  44- 50 F7.1   K       e_Teff   [-4051/-3.2] Lower error on Teff
  52- 57 F6.3   [cm/s2]   logg   [-0.85/5.24] Isochrone derived surface gravity
  59- 63 F5.3   [cm/s2] E_logg   [0.002/4.07] Upper error on logg
  65- 70 F6.3   [cm/s2] e_logg   [-3.27/-0.001] Lower error on logg
  72- 77 F6.3   [-]       [Fe/H] [-2.26/0.6] Isochrone derived surface
                                  metallicity
  79- 83 F5.3   [-]     E_[Fe/H] [0.002/0.68] Upper error on [Fe/H]
  85- 90 F6.3   [-]     e_[Fe/H] [-0.6/-0.002] Lower error on [Fe/H]
  92- 98 F7.3   Rsun      Rad    [0.13/588] Isochrone derived stellar radius
 100-106 F7.3   Rsun    E_Rad    [0.001/177] Upper error on Rad
 108-115 F8.3   Rsun    e_Rad    [-112/-0.001] Lower error on Rad
 117-122 F6.3   [-]       rho    [-8.08/1.65] Isochrone derived density
 124-129 F6.3   [-]     E_rho    [-8.37/1.92] Upper error on rho
 131-136 F6.3   [-]     e_rho    [-8.5/1.08] Lower error on rho
 138-143 F6.3   [Lsun]    Lum    [-2.87/4.1] Isochrone derived luminosity
 145-150 F6.3   [Lsun]  E_Lum    [-4/3.87] Upper error on Lum
 152-157 F6.3   [Lsun]  e_Lum    [-3.91/3.6] Lower error on Lum
 159-163 F5.2   Gyr       Age    [0.1/19.5] Isochrone derived age
 165-165 A1     ---     f_Age    [*] Age flag (1)
 167-171 F5.2   Gyr     E_Age    [0/18] Upper error on Age
 173-178 F6.2   Gyr     e_Age    [-17.5/0] Lower error on Age
 180-186 F7.1   pc        Dist   [11.9/14475] Isochrone derived distance
 188-193 F6.1   pc      E_Dist   [0.3/5385] Upper error on Dist
 195-201 F7.1   pc      e_Dist   [-4525/-0.3] Lower error on Dist
 203-207 F5.3   mag       Avmag  [0/3.26] Isochrone derived V-band extinction
 209-214 F6.4   ---       GOF    [0/1] Combined likelihood goodness-of-fit
 216-222 F7.2   Gyr       TAMS   [0.1/2842] Terminal age of the main sequence
--------------------------------------------------------------------------------
Note (1): Ages with uninformative posteriors (TAMS>20 Gyr) or unreliable ages
    (GOF<0.99) are flagged with an asterisk (26581 occurrences).
--------------------------------------------------------------------------------

History:
    From electronic version of the journal

References:
    Berger et al.   Paper II:   2020AJ....160..108B

================================================================================
(End)                          Prepared by [AAS], Coralie Fix [CDS], 29-Sep-2020
