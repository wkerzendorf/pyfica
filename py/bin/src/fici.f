        SUBROUTINE INITIATE
        SAVE

!       Initial model for radiation field.
!       Wavelength grid for spectra.

        INCLUDE 'ficc.f'



!       Initial radiation field - use optically-thin limit

        TEFF=10**ZTEFF
        DO 30 KN=1,LN
        DO 20 J=1,JS
20      XJI(KN,J)=WD(J)*BNU1(XN2(KN),TEFF,AP)
30      CONTINUE



!       Wavelength grid for spectra - used in FICM and FICZ
!       SHChange: alog -> log
        ZWVM(1)=LOG10(WVM1)
        ZWVM(KG+1)=LOG10(WVM2)
        DZWV=(ZWVM(KG+1)-ZWVM(1))/KG

        DO 40 K=1,KG+1
        ZWVM(K)=ZWVM(1)+(K-1)*DZWV

        WVM=10**ZWVM(K)                                      ! microns - vac
        wvman = 1.D4*WVM                                     ! Ang - vac 

!        transf. to vac. wvl. - only if wv(min)>2000A to avoid discont.         
        IF (WVM1 .GE. 2000.) THEN               
          WVMA = 1.D-4*WAIR(WVMan)                           ! microns - air
        ELSE
          WVMA = 1.D-4*wvman
        END IF  

!       SHChange: alog -> log
        ZWVMA(K)=LOG10(WVMA)
        write(29,2999)k,zwvm(k),dzwv,wvm,wvman,wvma,zwvma(k)
 2999   format(i5,6f12.6)

40      CONTINUE

        DO 50 K=1,KG
        DWVM(K)=10**ZWVM(K+1)-10**ZWVM(K)                       ! bin width
        DWVMA(K)=10**ZWVMA(K+1)-10**ZWVMA(K) 
50      CONTINUE



9999    RETURN
        END

