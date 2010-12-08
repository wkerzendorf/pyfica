        FUNCTION BNU1(XNU,TMP,AP)
        SAVE
!
!       Planck fn  Bnu = 2hnu**3/c**2/(exp(hnu/kT)-1)   - single precision
!
        DIMENSION AP(50)
!
!       SHChange: ALOG10 -> LOG10
        ZNU=LOG10(XNU)
        CF=2*10**(AP(14)+3*ZNU-2*AP(16))
        UU=10**(AP(14)+ZNU-AP(17))/TMP
!
        BNU1=CF*EXP(-UU)/(1-EXP(-UU))
!
        RETURN
        END
