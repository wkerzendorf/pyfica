        SUBROUTINE SPECTRUM
        SAVE

!       Calculate spectrum from formal integral. Quantity calculated is
!       luminosity density . Output wavs in air.

        INCLUDE 'ficc.f'

!       SHChange: was: DIMENSION FN(5000),ZLD(5000)
        DIMENSION FN(KZ),ZLD(KZ)

!T1=SECNDS(0.0)

        CALL SRCQ

        write(*,*)' Start Spectrum '

        WRITE(29,2009)ITER

        RE=1/XI(I1)                                ! Outer radius
        DP=0.25D0
        LL=RE/DP+1

        XLUM=0
        ZLUM=0

        DO 80 K=1,KG

        ZWVMB=(ZWVM(K)+ZWVM(K+1))/2
        XNU=10**(AP(16)+4-ZWVMB)                     ! Rest frame - vacuum

        FN(K)=0

        ZLD(K)=-99                                   ! Default

        DO 60 L=1,LL
        PP=(0.5D0+L-1)*DP
        IF(PP .GT. RE) GO TO 60

        CALL INTENS(PP,XNU,XIN)

        FN(K)=FN(K)+2*DP*PP*XIN

60      CONTINUE


        ZCF=AP(16)-2*(ZWVMB-4)

!       SHChange: ALOG -> LOG
        IF(FN(K).NE.0) ZLD(K)=2*ZRDCM+LOG10(4*AP(9)**2*FN(K))+ZCF


!       Emergent spectrum in air

        ZWVMAB=(ZWVMA(K)+ZWVMA(K+1))/2                       ! microns
        ZWVMIC(K) = ZWVMAB

        FC=DWVM(K)/DWVMA(K)                                  ! Jacobian
        ZLDFM(K) = FC*ZLD(K)                                 ! Formal
        ZLDMC(K) = FC*ZLDM(K)                                ! MC spec

        WRITE(29,1001)K,ZWVMAB,FC,ZLDFM(K),ZLDMC(K)          ! SPCT.OUD
1001    FORMAT(I9,F16.6,f16.8,2F16.5)

        XLUM=XLUM+10**(ZWVMB-4+ZLD(K)-40)*DZWV/AP(11)
!       ZLUM=ZLUM+10**(ZWVMB-4+ZLD(K))*DZWV/AP(11)

80      CONTINUE

!       SHChange: ALOG -> LOG
        ZLUM=LOG10(XLUM)+40

        CALL SPECTR


        WRITE(*,*) ' Finish formal '
        WRITE(*,1003)ZLUM,WVM1,WVM2
1003    FORMAT(/5X,'Log L =',F8.4,5X,'WVM1 =',F7.4,4X,'WVM2 =',F7.4)


!DELTA=SECNDS(T1)/60
!WRITE(*,2005)DELTA
!WRITE(39,2007)DELTA

2005    FORMAT(/35X,'DELTA =',F8.2/)
2007    FORMAT(/5X,'FICZ',50X,'Time =',F7.2)
2009    FORMAT(/5X,'Iteration',I4)

9999    RETURN
        END







        SUBROUTINE SRCQ
        SAVE

!       Calc source term S*(1-exp(-tau)) for each line in each shell from
!       XJI , the mean intensity in violet wing.


        INCLUDE 'ficc.f'


        CF=TM*10**(ZEP0-ZVOL)/4/AP(9)
        DD=4*AP(9)/TM

        DO 90 J=1,JS                                 ! Loop over shells


        DO 20 KK=1,NZ
20      EALV(KK)=0


!       Energy absorption rates to all upper levels

        DO 30 KN=1,LN

        WV=10**AP(16)/XN2(KN)

        XJI(KN,J)=CF*WV/VT(J)*XJI(KN,J)                   ! cgs

        M=IN(LN+1-KN)                                     ! Basic list
!       SHChange: JIABS -> ABS
        KKU=ABS(NUP(M))

        EALN=DD/WV*XJI(KN,J)*(1-EX(KN,J))                 ! Line

!  write(10,*)' srcq',kku
        EALV(KKU)=EALV(KKU)+EALN                          ! Level

30      CONTINUE


!       Energy emission rates for all lines

        DO 50 KN=1,LN
        WV=10**AP(16)/XN2(KN)

        M=IN(LN+1-KN)                                     ! Basic list
!       SHChange: JIABS -> ABS        
        KKU=ABS(NUP(M))

        EELN=QQ(M,J)*EALV(KKU)

        STM(KN,J)=WV/DD*EELN                              ! S*(1-exp(-tau))

!       if (m.eq.4757) then
!        write(10,*) wv,kn,j,DD,EELN,QQ(M,J),EALV(KKU)
!       end if

        XJE(KN,J)=XJI(KN,J)*EX(KN,J)+STM(KN,J)            ! J red wing

!       if (m.eq.4757) then
!        write(10,*)'z',j,XJE(KN,J),EX(KN,J),XJI(KN,J),STM(KN,J)
!       end if

50      CONTINUE


90      CONTINUE                                      ! End loop over shells


9999    RETURN
        END







        SUBROUTINE INTENS(PP,XNU,XIN)
        SAVE

!       Calc emergent intensitiy for line-of-sight with impact parameter PP.

        INCLUDE 'ficc.f'


        RE=1/XI(I1)                                ! Outer radius
!       XE is replaced by XE1 in input file
!       XE=0.200d0                      !!! MUST be same as in MODEL (fics.f)

        ZF=SQRT(RE**2-PP**2)
        ZS=-ZF
        IF(PP.LE.1) ZS=SQRT(1-PP**2)

        UPH=VPH*1.D5/10**AP(16)

        XNUS=XNU*(1-ZS*UPH)
        XNUF=XNU*(1-ZF*UPH)

        KNS=KNCAL(XNUS)
        KNF=KNCAL(XNUF)-1

!        lower boundary (Ic)
        IF(KC(3).EQ.1) THEN
          XIN0=0
          IF(PP.LE.1) XIN0=BNU1(XNU,TB,AP)        ! Static black body
        ELSE IF(KC(3).EQ.0) THEN
          XIN0=0
          IF(XNU.LT.XNUA .AND. XNU.GT.XNUB .AND. PP.LE.1) XIN0=XINB  ! Top-hat
        END IF

        XIN=XIN0

        ZP=ZS
        DO 50 KN=KNS,KNF

        ZKN=(1-XN2(KN)/XNU)/UPH
        XKN=1/SQRT(ZKN**2+PP**2)
        if (xkn.le.(XE1)) then
         write(*,*)' xkn=',xkn,' XE1=',xe1
         xkn = XE1+0.0001d0
        end if
        J=1+0.5D0*(1.00d0-XKN)/DX                        ! shell
        if (j.gt.JS) write(*,*)' J = ',j,' xkn=',xkn,' XE1=',xe1
        if (j.gt.JS) j = JS
        if (j.lt.1) write(*,*)' J = ',j,' xkn=',xkn,' XE1=',xe1
        if (j.lt.1) j = 1

!       Intensity change due to el scatt along path between consecutive
!       active lines - linearized, not centred.

        DTAUE=SG(J)*10**(ZRH(J)+ZRDCM)*(ZKN-ZP)
        XJM=0.5D0*(XJE(KN-1,J)+XJI(KN,J))
        XIN=XIN+(XJI(KN,J)-XIN)*DTAUE

!       if (kn.eq.32465) then
!        write(10,*) 'z1',kn,j,XIN,EX(KN,J),XJI(KN,J),STM(KN,J)
!       end if

!       Intensity change due to transition KN

        XIN=XIN*EX(KN,J)+STM(KN,J)

!       if (kn.eq.32465) then
!        write(10,*) 'z1',kn,j,XIN
!       end if

        ZP=ZKN
50      CONTINUE


9999    RETURN
        END




        SUBROUTINE SPECTR
        SAVE
!
!       Calculates spectrum & mags for output

        INCLUDE 'ficc.f'

!       SHChange: CHECK: This was
!       DIMENSION FEFM(1500),FEMC(1500),FEFMSL(1500),FEMCSL(1500)
!       DIMENSION DLMA(1500),WVR(200),RS(200,4)    --> replace 1500 by KG (1500 was default KG)
        DIMENSION FEFM(KG),FEMC(KG),FEFMSL(KG),FEMCSL(KG)
        DIMENSION DLMA(KG),WVR(200),RS(200,4)
!       SHChange: CHECK: This was DIMENSION BB(1500),XMGR(4),XMGF(4),XMG(4)
        DIMENSION BB(KG),XMGR(4),XMGF(4),XMG(4)
!
3       FORMAT(2A1)
4       FORMAT(F12.3,F10.2)
6       FORMAT(6F10.3,F12.3)
8       FORMAT(I5,F7.3,F7.2,F8.1,F7.4,46A1)
10      FORMAT(F10.2,2F9.3,F8.2/)
!
        write(*,*)' start Spectr '

        PI    = AP(9)
        HPL   = 10.0D0**(AP(14))
        CLIGHT= 10.0D0**(AP(16))
        BOLZK = 10.0D0**(AP(17))
        SIGMA = 10.0D0**(AP(24))

        WRITE(46,9)TDY,ZRD,ZLM
9       FORMAT(/'   DT =',F6.2,'   LOG R/RS =',F7.3,'   LOG L/LS =',F7.3)
        WRITE(46,5)ZTAUG,NP
5       FORMAT('      TAUB =',F5.1,'       NPR =',I8)

!          read in UBVR filters
        DO L=1,109
          READ (41,7) WVR(L),(RS(L,J),J=1,4)
        END DO
7       FORMAT(F9.1,4F8.3)

        ZL = ZLM + AP(26)

        NP  = 0
!XXNT = NT
!         Distance from distance modulus (m-M) (pc)
        DST = 10.D0**(1.D0 + 0.2D0*XMU)
        WRITE (46,2310) XMU,DST,EBMVGAL,EBMVCMF
 2310   FORMAT('  DISTANCE ',F8.3,F14.1,'  E(B-V)= ',F8.2,F8.2)

!         log(dist) (cm)
        ZDSN = log10(DST * 3.0857E18)
!         Flux at earth (erg/s/cm**2)
        FT    = 10.E0**(ZL-2.D0*ZDSN)/4.E0/AP(9)
        FEM   = 1.E0
        FETOT = 0.E0
!       SHChange: INIT FEFMTOT, FEMCTOT
        FEFMTOT=0.0D0
        FEMCTOT=0.0D0 
!       SHChange: CHECK: This was DO 40, K=1,1500
        DO 40, K=1,KG
!           lambda (Ang)
          WAV     = 10.E0**(ZWVMIC(K)+4.E0)
!           Delta_lambda (Angst)
          DLMA(K) = DWVMA(K) * 1.E+4
!           F_lambda (erg/s/cm**2/Ang)
          FEFM(K) = 10.E0**(ZLDFM(K)-8. -2.D0*ZDSN)/4.E0/AP(9)
          FEMC(K) = 10.E0**(ZLDMC(K)-8. -2.D0*ZDSN)/4.E0/AP(9)
!           Integrate observed flux
          FEFMTOT   = FEFMTOT + (FEFM(K) * DLMA(K))
          FEMCTOT   = FEMCTOT + (FEMC(K) * DLMA(K))

!           calculate black body spectrum
          WAVCM = WAV * 1.0E-8
          DLMCM = DLMA(K) * 1.0E-8
          BBLM1 = 2.0E0 * HPL * (CLIGHT**2) / WAVCM**5
          BBLME = HPL * CLIGHT / (WAVCM * BOLZK * TB)
          IF (BBLME .LE. 88.0E0) THEN
            BBLM  = BBLM1 / (exp(BBLME) - 1.0E0)
          ELSE
            BBLM  = 0.0E0
          END IF
          BBT   = TB * SIGMA * (TB**3) / PI
          RATIO = BBLM * DLMCM / BBT
          BB(K) = RATIO * FT / DLMA(K)

   40   CONTINUE

!       SHChange: CHECK: This was DO 50 K = 4,1497
        DO 50 K=4,KG-3

          WAV  = 10.E0**(ZWVMIC(K)+4.E0)

!         Smoothing
          FEFMS  = -2.D0*(FEFM(K-3)+FEFM(K+3))
          FEFMS  = FEFMS + 3.D0*(FEFM(K-2)+FEFM(K+2))
          FEFMS  = FEFMS + 7.D0*FEFM(K) + 6.D0*(FEFM(K-1)+FEFM(K+1))
          FEFMS0 = Amax1(FEFMS,0.E0)/21.E0

          FEMCS  = -2.D0*(FEMC(K-3)+FEMC(K+3))
          FEMCS  = FEMCS + 3.D0*(FEMC(K-2)+FEMC(K+2))
          FEMCS  = FEMCS + 7.D0*FEMC(K) + 6.D0*(FEMC(K-1)+FEMC(K+1))
          FEMCS0 = Amax1(FEMCS,0.E0)/21.E0

!           Host galaxy reddening

          AVCMF = EBMVCMF * (3.10E0 + 2.002E0*(1.E4/WAV - 1.E0/0.55E0))

          FEFMS   = FEFMS0 * 10.E0**(-0.4E0*AVCMF)
          FEMCS   = FEMCS0 * 10.E0**(-0.4E0*AVCMF)
          BBRED   = BB(K)  * 10.E0**(-0.4E0*AVCMF)

!           Apply Redshift (FLUX REDUCTION IN PLOTTING)

          WAVZ  = WAV * (1.E0 + REDSHIFT)

!           Galactic Reddening

          AVGAL = EBMVGAL * (3.10E0 + 2.002E0*(1.E4/WAVZ - 1.E0/0.55E0))

          FEFMS   = FEFMS * 10.E0**(-0.4E0*AVGAL)
          FEMCS   = FEMCS * 10.E0**(-0.4E0*AVGAL)
          BBRED   = BBRED * 10.E0**(-0.4E0*AVGAL)

!           Test for wavelength region to model

          WAVM  = WAV/1.E4
          IF (WAVM .LT. WVM2 .AND. WAVM .GT. WVM1) THEN
!             F_lambda (erg/s/cm**2/Ang)
            FEFMSL(K) = FEFMS                          ! OUTPUT FOR MONGO
            FEMCSL(K) = FEMCS
!             F_nu (J/s/m**2/Hz)
            FEFMSN = FEFMS * (WAV**2 /(CLIGHT*1.D+8)) / 1.D+3   ! O/P FOR MONGO
            FEMCSN = FEMCS * (WAV**2 /(CLIGHT*1.D+8)) / 1.D+3

            WRITE(32,4550)WAV,FEFMS0,FEFMSL(K),FEMCS0,FEMCSL(K),BBRED
 4550       FORMAT(F9.2,4E15.5,E15.5)
          END IF

   50   CONTINUE

        DO 62 J=1,4
62      XMG(J)=0.E0
!       SHChange: CHECK: This was DO 70 K=1,1500
        DO 70 K=1,KG
          WAV  = 10.E0**(ZWVMIC(K)+4.E0)
          IF(WAV.GE.8350.E0 .OR. WAV.LE.3050.E0) GO TO 70

!        interpolate filter responses
66      L=(WAV-3000.E0)/50.E0+1.E0
        DO 68 J=1,4
        DR=(RS(L+1,J)-RS(L,J))/50.E0
        AG=(RS(L+2,J)-RS(L+1,J)-RS(L,J)+RS(L-1,J))/4.E0/(50.E0)**2
        DL=WAV-WVR(L)
        RSI=RS(L,J)+DL*DR+AG*DL*(-50.E0+DL)
        XMGR(J) = XMGR(J) + DLMA(K) * FEMCSL(K) * RSI      ! redd. flux
        XMG(J) = XMG(J) + DLMA(K) * FEMC(K) * RSI          ! unredd. fl.
68      XMGF(J) = XMGF(J) + DLMA(K) * RSI                  ! filter resp.

70      CONTINUE

!        apply filter, find mags
        DO 75 J=1,4
        XMGR(J) = XMGR(J) / XMGF(J)
        XMGR(J)=-2.5E0*LOG10(XMGR(J))
        XMG(J) = XMG(J) / XMGF(J)
75      XMG(J)=-2.5E0*LOG10(XMG(J))

!       SHChange: replace filter curves, replace zeropoints
!       scheme: X0 = X_raw(Vega) - X_norm(Vega) + X_raw(unit spectrum)
!               X_raw(unit spectrum) from commenting line after DO 75 J=1,4 and comparison
!               X_norm(Vega) from recomendations STSDAS 3.5
!               X_raw(Vega) from midpoint-rule integrator

!        FU0 = 4.22E-9
!        FB0 = 6.27E-9
!        FV0 = 3.60E-9
!        FR0 = 2.16E-9
!        U0 = -2.5E0*log10(FU0)
!        B0 = -2.5E0*log10(FB0)
!        V0 = -2.5E0*log10(FV0)
!        R0 = -2.5E0*log10(FR0)
              
        U0 = 13.37846D0-0.056D0+7.562D0
        B0 = 12.23641D0-0.036D0+8.301D0
        V0 = 12.66775D0-0.026D0+8.437D0
        R0 = 12.55870D0-0.038D0+9.082D0
!       SHChange: strict parentheses to avoid compiler warnings 
        FU0 = 10**(U0/(-2.5D0))
        FB0 = 10**(B0/(-2.5D0))
        FV0 = 10**(V0/(-2.5D0))
        FR0 = 10**(R0/(-2.5D0))

!        unreddened colours
        U=XMG(1) - U0
        B=XMG(2) - B0
        V=XMG(3) - V0
        R=XMG(4) - R0
        UB=U-B
        BV=B-V
        VR=V-R
        XMB=AP(29)-2.5E0*ZLM+5.E0*LOG10(DST)-5.E0
        BC=XMB-V

!        observed colours
        UOB=XMGR(1) - U0
        BOB=XMGR(2) - B0
        VOB=XMGR(3) - V0
        ROB=XMGR(4) - R0
        UBOB=UOB-BOB
        BVOB=BOB-VOB
        VROB=VOB-ROB

        WRITE (3,3400) PW(2),FT,FEMCTOT
 3400   FORMAT(5X,'POWER: ',E12.4,4X,'INPUT FLUX: ',E12.4,4X,
     &         'INT. FLUX:',E12.4)

        WRITE(46,327)U,B,V,R
        WRITE(3,327)U,B,V,R
        WRITE(*,327)U,B,V,R
  327   FORMAT(/,'UNREDD. MAGS: U=   ',F7.3,' B=   ',F7.3,' V=   ',F7.3,' R=   ',F7.3)
        WRITE(46,328)UB,BV,VR
        WRITE(3,328)UB,BV,VR
        WRITE(*,328)UB,BV,VR
  328   FORMAT('UNREDD. COL.: U-B =',F7.3,' B-V =',F7.3,' V-R =',F7.3/)

        WRITE(46,329)UOB,BOB,VOB,ROB
        WRITE(3,329)UOB,BOB,VOB,ROB
        WRITE(*,329)UOB,BOB,VOB,ROB
  329   FORMAT('REDD. MAGS:   U=   ',F7.3,' B=   ',F7.3,' V=   ',F7.3,' R=   ',F7.3)
        WRITE(46,330)UBOB,BVOB,VROB
        WRITE(3,330)UBOB,BVOB,VROB
        WRITE(*,330)UBOB,BVOB,VROB
  330   FORMAT('REDD. COL.:   U-B =',F7.3,' B-V =',F7.3,' V-R =',F7.3/)

        UABS = U - XMU
        BABS = B - XMU
        VABS = V - XMU
        WRITE(46,332)UABS,BABS,VABS
        WRITE(3,332)UABS,BABS,VABS
        WRITE(*,332)UABS,BABS,VABS
  332   FORMAT('ABS. MAG:     U =  ',F7.3,' B =  ',F7.3,' V =  ',F7.3)

        XMBABS = XMB - XMU
        WRITE(46,331)XMB,BC,XMBABS
        WRITE(3,331)XMB,BC,XMBABS
        WRITE(*,331)XMB,BC,XMBABS
  331   FORMAT('MBol = ',F7.3,'  Bol. Corr. = ',F7.3,'  (MBol)o= ',F7.3)
!       WRITE(46,333)XMBABS
!       WRITE(3,333)XMBABS
! 333   FORMAT(/' (MBol)o =',F7.3)

        RETURN
        END
