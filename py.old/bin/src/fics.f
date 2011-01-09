!
!        NB density changed for 97ef case!!!!!! REMOVE !!!
!
        SUBROUTINE MODEL
        SAVE
        
!       Set initial model

        INCLUDE 'ficc.f'

        DIMENSION ZRD1(300),ZRH1(300)
        DIMENSION RD1(300),RH1(300),REFM(300),SHM(600)
        DIMENSION ZROUT(300),VK1(300),REFMNEW(300)


        I1=2*JS+1

!       XE is replaced by XE1 in input file
!       XE=0.200D0              !!! REMEMBER to copy this in INTENS (fiyz.f) 
        DX=(1-XE1)/(I1-1)
        DO 20 I=1,I1
20      XI(I)=1-(I-1)*DX


!       Density profile                 W7 
!        read density structure
!         LMAX = 19
!         DT1 = 31.535D0                                ! Reference day
!         ZDT1 = log10(DT1)
!         DT1SEC = DT1 * 24.D0 * 3600.D0
!         DO 22 L=1,LMAX
!           READ (19,222) LL,ZRD1(L),ZRH1(L)
!  222      FORMAT(I3,F17.13,F10.5)
!           ZRD1(L) = ZRD1(L) + AP(27)
!           zrh1(l) = zrh1(l) + log10(1.0d0)       !normal w7 density
!c if (l .eq. 17) then 
!c  zrh1(l) = zrh1(l) + log10(3.00d0)   
!c else if (l .eq. 18) then
!c  zrh1(l) = zrh1(l) + log10( 6.00d0)        !half density for bg
!c else if (l .eq. 19) then
!c  zrh1(l) = zrh1(l) + log10( 3.00d0)   
!c end if
!   22   CONTINUE

!       Density profile                 Wdd1 
!        read density structure
!       LMAX = 25
!       DT1SEC = 20.                            ! Reference time
!c       RHOC = 4.4511E-12                         
!       DT1 = DT1SEC / 24.D0 / 3600.D0
!       ZDT1 = log10(DT1)
!       READ (19,*)
!       DO 22 L=1,LMAX
!         READ (19,*) LL,rr,vcms,tt,rhorel
!         print *,LL,rhorel,vcms
!        RD1(L) = vcms * DT1SEC
!c 222     FORMAT()
!        ZRD1(L) = log10(RD1(L))
!        rh1(L) =  rhorel               !!!* RHOC
!         zrh1(l) = log10(rh1(l)) 
!  22   CONTINUE

!       Density profile                 Koichi's CO6 (tref=75sec, 176 lines)
!                                                CO138/100  100,   258 ) 
!        read density structure
!     LMAX = 19                          ! 25 for co110
!     DT1SEC = 100.D0            ! Reference day
!       DT1 =  1.d0                       !Maeda's rescaled model
!       DT1SEC = DT1*24.D0*3600.D0

!     header of density structure file
      READ(19,*)LMAX,DT1SEC ! number of points, epoch in s
      DT1 =  DT1SEC / 24.D0 / 3600.D0
      ZDT1 = log10(DT1)
      DO 22 L=1,LMAX
        READ (19,*) LL,VKS,RD1(L),RH1(L)
        ZRD1(L) = log10(RD1(L)) + log10(1.00d0) 
        zrh1(l) = log10(RH1(L)) + log10(1.00d0)   !!!change density for 97ef  
 22   CONTINUE

!         set up shell outer radii
         DO 21 L=1,LMAX
          IF (L .EQ. LMAX) THEN
            ZROUT(L) = ZRD1(L)
          ELSE
            ZROUT(L) = (ZRD1(L) + ZRD1(L+1))/2.D0
          END IF
   21   CONTINUE

!         write radii and densities
        WRITE (3,2200)
 2200   FORMAT(/5X,'  RADIUS  ','  VELOCITY ',' DENSITY  ',
     &         '  MASS  ','  TOT. MASS','   NEWMASS ','NEWTOT. MASS',
     &          ' Kin. En.')
        DO 220, L=1,LMAX
          RD1(L) = 10.D0**(ZRD1(L))
          RH1(L) = 10.D0**(ZRH1(L))
          VK1(L) = RD1(L) / DT1SEC / 1.E5
  220   CONTINUE

        SUNMASS = 10.D0**(AP(25))
        PI = AP(9)
!         masses
        ISH = 1
        REFM(1) = 4.D0/3.D0 * PI * (RD1(1)**2) * RH1(1) * RD1(1)/SUNMASS
        REFMNEW(1) = 4.D0/3.D0*PI *10.D0**(3.D0*ZROUT(1))*RH1(1)/SUNMASS
        TOTM = REFM(1)
        TOTMNEW = REFMNEW(1)
!        kinetic energy
        ENKIN  = 0.5d0*REFMNEW(1)*(VK1(1)*1.d5)**2  
        WRITE (3,2210)ISH,RD1(1),VK1(1),RH1(1),REFM(1),TOTM,
     &                 REFMNEW(1),TOTMNEW,ENKIN
        DO 221, ISH = 2,LMAX
          RSH = (RD1(ISH) + RD1(ISH-1))/2.D0
!        rsh=rd1(ish)
          DELR = RD1(ISH) - RD1(ISH-1)
          REFM(ISH) = 4.D0 * PI * (RSH**2) * RH1(ISH) / SUNMASS *DELR
          TOTM = TOTM + REFM(ISH)
          RINC = 10.D0**(ZROUT(ISH-1))
          ROUC = 10.D0**(ZROUT(ISH))
          rincsm = rinc*rinc/sunmass*rinc
          roucsm = rouc*rouc/sunmass*rouc
!          DELRCU = ROUCUB - RINCUB
           delrcu = roucsm - rincsm             !incl /sunmass
          REFMNEW(ISH) = 4.D0/3.D0 * PI * DELRCU * RH1(ISH)  !/sunmass
          TOTMNEW = TOTMNEW + REFMNEW(ISH)
!          kinetic energy
          ENKIN  = ENKIN + 0.5d0*REFMNEW(ISH)*(VK1(ISH)*1.d5)**2   
          WRITE (3,2210)ISH,RD1(ISH),VK1(ISH),RH1(ISH),REFM(ISH),TOTM,
     &                  REFMNEW(ISH),TOTMNEW,ENKIN
 2210     FORMAT(I4,D10.2,F10.2,D10.2,5D11.3)
  221   CONTINUE

!          kinetic energy
        EKLOG = log10(ENKIN) + log10(sunmass)
        WRITE (3,*)'/, Kin. En. = ', EKLOG

!         Loop over major shells (odd numbered)

!       SHChange: Comment: experiment with diff. initial temp: TST=0.7*(10**ZTEFF)
	TST=10**ZTEFF

        DO 25 J=1,JS
          XS=XI(2*J)
          ZR1 = ZRDCM - log10(XS) + log10(DT1/TDY)
          IF (ZR1 .GE. ZRD1(LMAX)) THEN
            L = LMAX - 1
          ELSE IF (ZR1 .LT. ZRD1(1)) THEN
            L = 1
          ELSE
            DO 23, LL = 1,LMAX-1
              IF (ZR1 .GE. ZRD1(LL) .AND. ZR1 .LT. ZRD1(LL+1)) THEN
                L = LL
              END IF
   23       CONTINUE
          END IF

          DRV = (ZRH1(L+1)-ZRH1(L))/(ZRD1(L+1)-ZRD1(L))
          ZRHT1 = ZRH1(L) + DRV*(ZR1-ZRD1(L))
          ZRH(J) = ZRHT1+3.D0*LOG10(DT1/TDY)
          VS(J) = VPH/XS
          TR(J) = TST
          TS(J) = 0.9D0*TR(J)
          WD(J) = 0.5D0*(1.D0-SQRT(1.D0-XS**2))

25      CONTINUE


!        compute mass above photosphere
        WRITE(3,109)
109     FORMAT(/6X,'XS',6X,'VS',8X,'LOG RH',6X,'TE',7X,'TR',7X,'W',4X,
     &         'SH M    TOTM'/)
        TOTMAB = 0.D0
        rrefcu = 10.D0**(ZRDCM)
        rrefcu = rrefcu*rrefcu/sunmass*rrefcu
        DO  J=1,JS
          XS=XI(2*J)
          XLOCU = (1.d0/XI(2*J -1))**3
          XUPCU = (1.d0/XI(2*J +1))**3
          DELRCU = (XUPCU - XLOCU) * RREFCU
          SHM(J) = 4.D0/3.d0*PI * DELRCU *10.D0**(ZRH(J))    !!!/SUNMASS
          TOTMAB = TOTMAB + SHM(J)
          WRITE(3,9)XS,VS(J),ZRH(J),TS(J),TR(J),WD(J),SHM(J),TOTMAB
        END DO
9       FORMAT(5F10.2,3F7.3)

        WRITE(3,2009)
2009    FORMAT(/8X,'XS',8X,'VS',10X,'LOG RH',7X,'TE',10X,'TR',9X,'W'/)
        DO 35 J=1,JS
          XS=XI(2*J)
35      WRITE(3,1009)XS,VS(J),ZRH(J),TS(J),TR(J),WD(J)           ! STST.OUD
1009    FORMAT(5F12.2,F9.3)


!       Volumes of shells - unit 4/3*pi*R**3
        DO 50 J=1,JS
          I=2*J
          VT(J)=1/XI(I+1)**3-1/XI(I-1)**3
50      CONTINUE


        WRITE(*,2001)
2001    FORMAT(/5X,'Structure ','calculated')


9999    RETURN
        END

