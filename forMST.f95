PROGRAM forMST
      ! 221003 05:00 
      ! BASED ON Strickland's DART.
      ! cT, cT+cD0 data obtained from Strickland Tab. 1/16 and stored in
      ! STRIDATA_ORIG.txt. Cross-validated using data digitized from 
      ! Fig. 2.5/p15/Sanyer.
      ! REFERENCES:
      ! 1. Strickland J. H., 1975, "The Darrieus Turbine: A Performance 
      ! Prediction Model Using Multiple Streamtubes," SAND75-0431, 
      ! Sandia Laboratories energy report.
      ! 2. Reef, J. W., Maydew, R. C. and Blackwell, B. F., 1974, "Wind Energy
      ! Potential in New Mexico, "Sandia Laboratory Report SAND74-0071, July.
      ! 3. Ning, A., “Using Blade Element Momentum Methods with Gradient-Based 
      ! Design Optimization,” Structural and Multidisciplinary Optimization, 
      ! May 2021. 
      ! 4. Sanyer W. E., 2011, "The Development of a Wind Turbine for 
      ! Residential Use," M.Sc. thesis, Dept. of Mech. Eng., North Carolina 
      ! State University
      COMMON/TABLS/TA(50),TCN(50),TCT(50),NTBL1,XCD0
      open (unit=99, file='STRIDATA_ORIG.txt', status='old', action='read')
      ALIM = 15.0       ! (deg): AoA(deg) below which XCD0 is subtracted
      SHEAR = 0.
      NT = 19           ! (non): number of streamtubes
      ERR = 0.001       ! (non): induction factor iteration error
      XCD0 = 0.0085     ! (non): cD at AoA = 0
      MAXTRY = 100      ! (non): max. permitted # of iterations
      NTBL = 30         ! (non): the number of rows in the array with
                        ! "alpUD", "cT", "cT+cD0" vs AoA data      
      NUMCOLS = 3       ! (non): the number of columns in this array
      NTSR = 14         ! (non): the number of sampled Tip Speed Ratios (TSR)
      TSRI = 1.0        ! (non): Initial TSR value
      DTSR = 0.5        ! (non): TSR increment
      HR = 1.0          ! (non): height ratio: rotor height/max. rotor radius
      S = 0.27          ! (non): rotor solidity
      WRITE(*,200) XCD0 ! where XCD0=0.0085, see Strickland Tab. 1/16
 200  FORMAT(F7.5)
 201  FORMAT(I3)
      WRITE(*,201) NTBL
      WRITE(*,201) NUMCOLS
 203  FORMAT(I3,4F10.2)
      WRITE(*,203) NTSR, TSRI, DTSR, HR, S 

      PY = 4.*ATAN(1.)  ! (non): trig constants
      DTR = PY/180.
      NTBL1=NTBL-1      ! number of sampled AoAs, less one

      DO I=1,NTBL       ! read the array with airfoil charecterics vs AoA
        READ(99,*) TA(I),TCN(I),TCT(I)
        !WRITE(*,22) TA(I),TCN(I),TCT(I)
      END DO
  !22  FORMAT(2X,F6.0,2F8.5)
      DO 60 J = 1, NTSR         ! BEGIN LOOP OVER TSR
        X = TSRI + (J - 1.) * DTSR  ! define the TSR to be evaluated
        IF (SHEAR.EQ.1.) THEN   ! wind shear considered, if SHEAR = 1, see
                                ! Stri/S3.5/p23/e26
            GO TO 6             ! Reduce step over rotor height if SHEAR = 1
        END IF
        NZH = 5                 ! (non): # of blade elements along rotor 
                                !        vertical axis Z, if no shear 
        GO TO 7
    6   NZH = 10                ! (non): # of blade elements increased, to
                                !        accomodate shear consideration
    7   CPSUM = 0.0             ! (non): initialize cP summation variable
                                !        over the rotor height 
        RRSUM = 0.0             ! (non): initialize RLOCAL/RMAX summation
                                !        variable over the rotor height
        DO 90 I = 1, NZH    ! BEGIN LOOP OVER BLADE LENGTH
            IF(SHEAR.EQ.1.) THEN    ! wind shear considered, if SHEAR = 1, see
                                    ! Stri/S3.5/p23
                GO TO 8
            END IF
            ZH = (2. * I - 1.)/(4. * NZH)   ! (non): define the nondimensional
                                    ! coordinate Z/H along rotor vertical axis
                                    ! i.e., the ratio of the height above
                                    ! the base of the rotor to the overall
                                    ! rotor height for no shear (SHEAR=0). 
                                    ! See Eq. 22/p25/e27.
                                    ! 
            U1 = X                  ! (non): define the auxiliary var. U1 as
                                    ! TSR according to standard definition 
                                    ! TSR=UTMAX/UINF for no shear (SHEAR=0).
            GO TO 9                 ! Skip Z/H and TSR definitions for SHEAR=1 
    8       ZH = (2. * I - 1.)/(2. * NZH)
            U1 = X/((2. * ZH)**(1. / 7.))   ! (non): define the TSR accounting
                                    ! for wind shear, so that when SHEAR = 1,
                                    ! TSR=UTMAX/UINF, with UINF defined as
                                    ! UINF = UINF * (2. * ZH)**(1. / 7.). See 
                                    ! Ref. 1, above.
    9       CONTINUE
            RR = SIN(PY*ZH)         ! (non): RLOCAL/RMAX. 
                                    ! Troposkein is approximated
                                    ! with a sinusoid. So that RR=1 for ZH=1/2
            RRSUM = RRSUM + RR
            BETA = ATAN(HR/(PY*COS(PY*ZH))) ! (rad): blade segment slope
            SBETA = SIN(BETA)
            DO 89 K = 1, NT             ! BEGIN LOOP OVER STREAMTUBES
                T = 90.*(2.*K-1.)/NT    ! (deg): local blade azimuth
                THETA = T*DTR           ! (rad): local blade azimuth
                STH = SIN(THETA)
                CTH = COS(THETA)
                AA = 0.                 ! (non): initialize the INTERFERENCE 
                                        ! FACTOR, a ratio of the reduction in
                                        ! "wind velocity seen at the rotor", 
                                        ! to "wind velocity at infinity".
                NTRY = 0                ! (non): number of iterations
                    !
                    ! BEGIN LOOP OVER THE INTERFERENCE FACTOR AA
                    ! 
  100               U2 = 1. - AA    ! (non): U2 = U/UINF = 1 - (Uinf - U)/Uinf
                                    ! See Stri/E12/p10/e12
                    U3 = RR * U1/U2 ! (non): auxiliary variable needed in 
                                    ! the expression for ALPHA, below. 
                                    ! U3 = UT_LOC/U, see (*)/DMS/58
                    ALPHA = ATAN(STH*SBETA/(CTH + U3))  ! (rad): the AoA, see
                    ! Stri/E10/p10/e12
                    IF (ALPHA.LT.0.) THEN
                        ALPHA = PY + ALPHA  ! (rad): adding PI (rad) results in
                                    ! changing the sign of sin, cos(ALPHA)
                    END IF
                    SAL = SIN(ALPHA)
                    ALD = ALPHA/DTR ! (deg): conversion of AoA to degrees
                    CALL CNCT(ALIM, ALD, CN, CT)  ! (non): compute cN, cT for
                                    !  a given AoA           
                    U4 = STH*SBETA/SAL  ! U4 = UR/U, auxiliary variable needed 
                                    ! to obtain U5, derived from 
                                    ! Stri/E11/p10
                    U5 = (U4 * U2)**2   ! U5 = (UR/UINF)**2, auxiliary variable
                                    ! needed to obtain FX
                    FX = S*U5*(CN-CT*CTH/(STH*SBETA))/(4.* PY*RR)  ! (non): E9/
                                    ! p10: nondimensional streamwise force
                    ANEW = AA * AA + FX ! (non): updated INTERFERENCE FACTOR
                    NTRY = NTRY + 1 ! (non): iteration number
                    IF(NTRY.LE.MAXTRY) THEN
                        GO TO 81
                    END IF
                    U2 = -U2        ! (m/s): suggests opposite directions
                                    ! of U (through the rotor) and UINF
                                    ! (wind) velocities directions
                                    ! See also Algorithms 1, 2/p9/e10/Ref. 2
                    GO TO 130
   81               IF(ANEW.GT.1.) THEN ! Check whether physically possible
                        GO TO 70    ! Skip to the line that imposes AA = 1
                    END IF
                    IF(ABS(ANEW-AA).LT.ERR) THEN    ! Check whether iteration
                                    ! goal achieved
                        GO TO 130   ! Skip to local cP computation
                    END IF
                    AA = ANEW       ! (non): update the INTERFERENCE FACTOR
                    GO TO 100       ! END LOOP OVER AA
   70           U2 = 0.         ! (non): imposes AA = 1, cPLOC = 0
  130           CPLOC = (CT*S*X*RR*(X*U2*U4/U1)**2)/(2.*SBETA)  ! (non): cP
                                ! local power coefficient based on the local
                                ! torque and the area 2*R*del_h
                                ! See Eq. 18/p12/e16
                CPSUM = CPSUM + CPLOC   ! (non): Eq. 18/p12/e16
                FN = -CN * U5/(U1**2)   ! (non): Normal and tangent forces
                FT = -CT * FN/CN        ! (non): Eq. 7/p8/e10
   89       CONTINUE    ! END LOOP OVER STREAMTUBES (CN, CT VS AOA)
   !80       FORMAT(1X,8E12.4)
   90   CONTINUE    ! END LOOP OVER BLADE LENGTH
        CP = CPSUM/(NT * RRSUM) ! (non): Eq. 19/p12/e16
        WRITE(*,30) X, CP
   30   FORMAT(1X, 2E14.6) 
   60 CONTINUE      ! END LOOP OVER SAMPLED TSR 
      END PROGRAM forMST

      SUBROUTINE CNCT(ALIM, A, CN, CT)
        COMMON/TABLS/TA(50),TCN(50),TCT(50),NTBL1,XCD0
        DO 301 I = 1, NTBL1 ! LOOP OVER AoAs
            J = I
            IF(A.GE.TA(I) .AND. A.LE.TA(I+1)) THEN  ! (deg): Skip to interpolation
                GO TO 302
            END IF
  301   END DO
  302   X = (A - TA(J))/(TA(J + 1) - TA(J)) ! (deg): linear interpolation
        CN = TCN(J) + X * (TCN(J + 1) - TCN(J))
        CT = TCT(J) + X * (TCT(J + 1) - TCT(J))
        IF(A.LE.ALIM) THEN
            CT = CT - XCD0
        END IF
        RETURN
      END SUBROUTINE CNCT