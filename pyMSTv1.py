"""
221006 15:50
Simplified Multiple Streamtube model for VAWT performance estimation.
References:
    (1) Strickland J. H., 1975, "The Darrieus Turbine: A Performance Prediction
        Model Using Multiple Streamtubes," SAND75-0431, Sandia Laboratories energy report
    (2) Sanyer W. E., 2011, "The Development of a Wind Turbine for Residential Use,"
        M.Sc. thesis, Dept. of Mech. Eng., North Carolina State University
    (3) Ning, A., “Using Blade Element Momentum Methods with Gradient-Based
        Design Optimization,” Structural and Multidisciplinary Optimization,
        May 2021.

"""

import numpy as np
import pandas as pd
import modMST as mM
                    # input Strickland's data
dfcNcT = pd.read_csv("STRIDATA_ORIG.csv", sep=',', usecols=[0, 1, 2])
valcNcT = dfcNcT.values
alpArrUD = np.array(valcNcT[:, 0])  # (deg): angle of attack (AoA) vector
cNArr = np.array(valcNcT[:, 1])     # (non): normal force coefficients vector
cTArr = np.array(valcNcT[:, 2])     # (non): tangential force coefficient vec.
                    # remaining input data:
alpStalUD = 15.0    # (deg): AoA(deg) below which XCD0 is subtracted
switShear = 0       # (non): zero if the wind shear not accounted for
numST = 19          # (non): number of StreamTubes
numDelS_1 = 10      # (non): steps over rotor height if SHEAR = 1
numDelS_0 =  5      # (non): steps over rotor height if SHEAR NOT = 1
aErr = 0.001        # (non): induction factor iteration error
cD0 = 0.0085        # (non): cD at AoA = 0
maxIter = 100       # (non): max.permitted  # of iterations
numTSR = 14         # (non): the number of sampled Tip Speed Ratios (TSR)
iniTSR = 1.0        # (non): Initial TSR value
delTSR = 0.5        # (non): TSR increment
h2r = 1.0           # (non): height ratio: rotor height / max. rotor radius
sig = 0.27          # (non): rotor solidity
switTropo = 1       # (non): if switTropo == 1 then troposkein shape evaluated
gamUD = 0.0         # (deg): blade pitch
deg2rad = np.pi/180.
gamUR = gamUD * deg2rad # (rad): blade pitch converted
                    #
                    # LOOP OVER SAMPLED TSR VALUES
                    #
for ordTSR in range(1, numTSR + 1):
    curTSR = iniTSR + (ordTSR - 1.) * delTSR  # (non): current TSR to be
                    # considered
    if switShear:
        numDelZ = numDelS_1    # (non): the number of steps along VAWT's Z-axis
    else:
        numDelZ = numDelS_0
    cPsum = 0.0     # (non): initialize cP summation variable
                    # over the rotor height
    rLoc2rMaxSum = 0.0     # (non): initialize RLOCAL/RMAX summation
                    # variable over the rotor height
                    #
                    # LOOP OVER BLADE VERTICAL COORDINATE
                    #
    for ordZ in range(1, numDelZ):
        if switShear:
            zDh = (2. * ordZ - 1.) / (2. * numDelZ)  # (non): the ratio of
                    # the height above the base of the rotor to the overall
                    # rotor height
            TSRatZ = curTSR / ((2. * zDh) ** (1./7.))  # (non): define the TSR
                    # accounting for wind shear, so that when SHEAR = 1,
                    # TSR=UTMAX/UINF, with UINF defined as:
                    # UINF = UINF * (2. * ZH)**(1. / 7.)
        else:
            zDh = (2. * ordZ - 1.) / (4. * numDelZ)
            TSRatZ = curTSR    # (non). define the auxiliary variable TSRatZ
                    # as TSR = UTMAX/UINF
        if switTropo:
            rLoc2rMax = np.sin(np.pi * zDh) # (non): RLOCAL/RMAX. NOTE:
                    # Troposkein is approximated with a sinusoid
            betUR = np.arctan(h2r / (np.pi * np.cos(np.pi * zDh)))  # (non):
                    # slope of the blade segment under consideration
        else:
            rLoc2rMax = 1.0       # (non): RLOCAL/RMAX
            betUR = np.pi/2.   # (rad): slope
        rLoc2rMaxSum = rLoc2rMaxSum + rLoc2rMax
        SbetUR = np.sin(betUR)
                    #
                    # LOOP OVER STREAMTUBES
                    #
        for ordST in range(1, numST + 1):
            theUD = 90. * (2. * ordST - 1.) / numST # (deg): local blade azimuth
            theUR = theUD * deg2rad  # (rad): local blade azimuth
                    #
                    # ITERATE TO FIND a Induction Factor FOR THIS STREAMTUBE
                    #
            cPloc = mM.cPST(aErr,sig,curTSR,maxIter,alpStalUD,theUR,betUR,
                            rLoc2rMax,TSRatZ,gamUR,alpArrUD,cNArr,cTArr,cD0)

            cPsum = cPsum + cPloc
        # END LOOP OVER STREAMTUBES
    # END LOOP OVER BLADE LENGTH/VERTICAL COORDINATES
    cPatTSR = cPsum / (numST * rLoc2rMaxSum)
    print(curTSR, " ", cPatTSR)
# END LOOP OVER TSRs
