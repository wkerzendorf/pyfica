        FUNCTION WVAC(WAIR)

!       Wavelength(A) in vacuum  computed from air wvl.
!       - see pp 83,124 in Allen AQ(1973)
!       Formula has singularities at 827.6 AND 1561.7A. Therefore set WVAC=
!       WAIR for wav below 2000A


        WVAC=WAIR                       ! 1st approx 

        IF(WAIR .LT. 2000) GO TO 9999

        WV0=WVAC/1.D4                   ! A to microns
        RW2=1.D0/WV0**2
        TM1=29498.1D-6/(146-RW2)
        TM2=255.4D-6/(41-RW2)
        XN=1+64.328D-6+TM1+TM2

        WVAC=WAIR*XN                    ! 2nd approx

9999    RETURN
        END

