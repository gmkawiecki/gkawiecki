"""
221006 15:51
Contains methods used in pyMSTv1
"""
import numpy as np
from scipy import interpolate

print("methods needed to run the simplified VAWT MST model imported from modMST")


def cNcT(alpStalUD,alpUD,alpArrUD,cNArr,cTArr,cD0):
    """
    221003 05:26
    Interpolate Strickland's data (Table 1/16, for Re = 3.e5) to get cN, CT
        for local alpUD
    :param alpStalUD:   (deg):  stall angle of attack (AoA). cD at AoA=0º (cD0)
                                value is subtracted from cT for AoA<alpStalUD,
                                to adhere to Strickland's approach.
                                See p. 25/Ref. 1.
    :param alpUD:       (deg): local AoA
    :param alpArrUD:    (deg): vector of AoAs Strickland's Tab. 1 data are for
    :param cNArr:       (non): normal force coefficient vector
    :param cTArr:       (non): tangent force coefficient vector
    :param cD0:         (non): cD at AoA = 0
    :return: cN, cT     (non, non): normal and tangential force coefficients
    """
    cN, cT = 999., 99.  # (non): dummy variables, to avoid a warning message

    numRows = alpArrUD.shape[0]  # (non): the number of rows in the array with
                    # "AoA", "cT", "cT+cD0" vs AoA data
    for I_ in range(numRows):  # LOOP OVER AoAs
        for ind in range(len(alpArrUD) - 1):
            if alpArrUD[ind] < alpUD < alpArrUD[ind + 1]:
                x = np.array([alpArrUD[ind], alpArrUD[ind + 1]])
                yN = np.array([cNArr[ind], cNArr[ind + 1]])
                yT = np.array([cTArr[ind], cTArr[ind + 1]])
                fN = interpolate.interp1d(x, yN)
                fT = interpolate.interp1d(x, yT)
                cN = fN(alpUD)
                cT = fT(alpUD)
                if alpUD <= alpStalUD:
                    cT = cT - cD0
    return cN, cT


def cPST(aErr,sig,curTSR,maxIter,alpStalUD,theUR,betUR,
         rLoc2rMax,TSRatZ,gamUR,alpArrUD,cNArr,cTArr,cD0):
    """
    221002 17:22
    This method computes local cP for given azimuth, vertical coordinates
    :param aErr:        (non): induction factor iteration error
    :param sig:         (non): rotor solidity
    :param curTSR:      (non): Tip Speed Ratio
    :param maxIter:     (non): max. permitted  # of iterations
    :param alpStalUD:   (deg): AoA below which XCD0 value is subtracted
    :param theUR:       (rad): streamtube azimuth coordinate
    :param betUR:       (deg): blade inclination w.r.t. plane perp. 2 VAWT axis
    :param rLoc2rMax:   (non): (local radius)/(max radius)
    :param TSRatZ:      (non): TSR corresponding to the sampled coordinat
    :param gamUR:       (rad): blade pitch
    :param alpArrUD:    (deg): vector of AoAs Strickland's Tab. 1 data are for
    :param cNArr:       (non): normal force coefficient vector
    :param cTArr:       (non): tangent force coefficient vector
    :param cD0:         (non): cD at AoA = 0
    :return: cP         (non): power coefficient
    """
    deg2rad = np.radians(1.)
    aIndFac = 0.0   # (non):    initialize the INTERFERENCE FACTOR:
                    #           a ratio of the reduction in:
                    #           "wind velocity seen at the rotor"
                    # to
                    #           "wind velocity at infinity"
    ordIter = 0       # (non): initialize the number of iterations
    StheUR = np.sin(theUR)
    SbetUR = np.sin(betUR)
    CtheUR = np.cos(theUR)
    while True:
        u2uInf = 1. - aIndFac  # (non): velocity through the rotor (u) / U_infinity
        uTan2u = rLoc2rMax * TSRatZ / u2uInf   # (non): tangential velocity (uTan)
                    # divided by the velocity through the rotor (u)
        try:
            alpUR = np.arctan(StheUR * SbetUR/(CtheUR + uTan2u))  # (rad): AoA
        except ZeroDivisionError as error:
            print(error)
            print("for thetaUD = {} (deg)".format(theUR/deg2rad))
        alpUR = alpUR + gamUR # (rad): add blade's local pitch setting
        if alpUR < 0.0:
            alpUR += np.pi
        SalpUR = np.sin(alpUR)
        alpUD = alpUR/deg2rad  # (deg): conversion of AoA to degrees
        cN, cT = cNcT(alpStalUD,alpUD,alpArrUD,cNArr,cTArr,cD0)   # (non):
                        #  compute cN, cT for a given alpUD
        uRes2u = StheUR * SbetUR/SalpUR    # (non): ratio of the resultant
                        # velocity (uRes) to the velocity through the rotor
        uRes2uInf2 = (uRes2u * u2uInf)**2    # (non): (UR_/UINF_)**2, auxiliary
                    # variable needed to obtain fS, see below
        fS = sig * uRes2uInf2 * (cN - cT * CtheUR/(StheUR * SbetUR)
                          )/(4. * np.pi * rLoc2rMax)    # (non): nondimensional
                    # streamwise force
        aIter = aIndFac * aIndFac + fS # (non): update the INTERFERENCE FACTOR
        ordIter = ordIter + 1   # (non): iteration number
        if ordIter > maxIter:
            print("ordIter > maxIter")
            u2uInf = -u2uInf
            uRes2u = StheUR * SbetUR / SalpUR  # (non): UR_/U_
            break
        if aIter > 1.0: # Check whether physically possible
            u2uInf = 0.0   # (non): imposes aIndFac = 1, cPloc = 0
            uRes2u = StheUR * SbetUR / SalpUR  # (non): UR_/U_
            break
        if abs(aIter - aIndFac) < aErr: # Check whether the iteration goal
                    # has been achieved
            break
        else:
            aIndFac = aIter # (non): update the INTERFERENCE FACTOR
    cPloc = (cT * sig * curTSR * rLoc2rMax *
              (curTSR * u2uInf * uRes2u / TSRatZ) ** 2) / (2. * SbetUR)  # (non):
                    # local power coefficient based on the local
                    # torque and the area 2*R*del_h
    return cPloc



