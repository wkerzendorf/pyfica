        FUNCTION WAIR(WVAC)

!       Wavelength(A) in air computed from vac. wvl.
!       - see pp 83,124 in Allen AQ(1973)
!       NB This formula has singularities at 827.6 AND 1561.7A. Therefore
!       set WAIR=WVAC for wav below 2000A.

        WAIR=WVAC
        IF(WVAC .LT. 2000.) wvac=2000. !!!GO TO 9999

        WV0=WVAC/1.D4                   ! A to microns
        RW2=1.D0/WV0**2
        TM1=29498.1D-6/(146-RW2)
        TM2=255.4D-6/(41-RW2)
        XN=1+64.328D-6+TM1+TM2

        WAIR=WVAC/XN

9999    RETURN
        END

