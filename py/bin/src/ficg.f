        SUBROUTINE MONTEC1(LC)
        SAVE

!       Photosphere emits pkts with frequencies uniformly distributed in 
!       interval (XNUA,XNUB).


        INCLUDE 'ficc.f'

!       SHChange: US, XL must be adjustable (see fiym.f)
        DIMENSION US(JZ*2+1),XL(JZ)

!T1=SECNDS(0.0)

        WRITE(*,2001)


        CALL ZEROS1

        U1=VPH*1.D5/10**AP(16)

        EP0=1                                     ! Rescaled below  

        NP1 = NP
        IF (LC .EQ. 0) THEN                     ! TEMP IT
          NP1 = 20000
        ELSE IF (LC .EQ. 1) THEN                     ! Trial soln for GMX
          NP1 = NP/2
        END IF
        PW(1)=NP1

        DO 25 J=1,JS
25      XL(J)=10**(-ZRH(J)-ZRDCM)

        DO 30 I=1,I1
30      US(I)=U1/XI(I)



        NN=0
        DO 90 II=1,NP1
        NN=NN+1
        I=0
        ID=1
        LE=0
        KNE=0

        XB=1
        UB=U1


!       Emission at photosphere from interval (XNUA,XNUB)


        DCN=SQRT(1-RAN(KR))

        XNR=XNUA+RAN(KR)*(XNUB-XNUA)
        XPB=XNR*(1-DCN*U1)


        EP=EP0
        

        KN=KNCAL1(XPB)                                  ! Location in line list


        GO TO 65

45      TAUR=-LOG(1-RAN(KR))
        TAUS=TAUR
        TAUQ=-TAUR
        SP=0

        XB=XI(NB)
        ARG=1-SN2*(XB/XE)**2

        IF(ARG.LE.0) THEN 
          WRITE(*,163)J,I,NB,ID,XB,XE,SN2,DC,ARG
163       FORMAT(4I4,4F10.5,D15.6)
          ARG = -ARG
        END IF


        DCN=ID*SQRT(ARG)
!DCN=ID*SQRT(1-SN2*(XB/XE)**2)
        XPB=XNR*(1-DCN*US(NB))

        SB=-DC/XE+DCN/XB

        IF(SB.LE.0) WRITE(*,2009)

        CFB=-1/U1/XNR

48      ST=SB
        if (lt.ne.1) then 
         ltold=lt
        end if
        LT=1

        SE=SP+TAUS/RLE                   
        IF(SE.LT.ST)THEN
        ST=SE
        LT=2
        END IF

        SP=CFB*(XN2(KN)-XPE)

        IF(SP.LT.0) WRITE(*,2009)


        IF(SP.LT.ST)THEN

        XJI(KN,J)=XJI(KN,J)+EP
        
        TAUQ=TAUQ+TAUL(KN,J)
        TAUS=-TAUQ-SP*RLE

        KN=KN+1                                  ! next line
        IF(TAUS .GT. 0) GO TO 48
        ST=SP
        LT=4
        KN=KN-1                                  ! Line absorbing pkt
        END IF

        IF(LT.NE.1)THEN

        RE=1/XE
        SQ=ST+DC*RE
        XE=1/SQRT(RE*(RE+DC*ST)+ST*SQ)

        DI=SQ*XE
        UE=U1/XE

        DC=-1+2*RAN(KR)
!  DC=(DC+UE)/(1+DC*UE)                         !neglect terms O(v/c)
        EF=(1-DI*UE)/(1-DC*UE)
        XNR=XNR*EF
        XPE=XNR*(1-DC*UE)
        EPI=EP
        EP=EPI*EF
        PW(5)=PW(5)+(EP-EPI)                        ! mech goes to rad energy

        IF(LT .EQ. 4) THEN

        EA(KN,J)=EA(KN,J)+EPI                              ! Absorbed energy 

        CALL NEWFRQ(KN,XPE,EP,J)                           ! FICN  

        KNE=KN                                             ! Emitted line
        KN=KN+1                                            ! next line
        XNR=XPE/(1-DC*UE)
        LE=1                                               ! Line abs-em
        END IF

        GO TO 80
        END IF

65      L=I+(ID+1)/2
        I=I+ID
        IF(I.EQ.0 .OR. I.EQ.I1) GO TO 85                   ! Leave env

        XE=XB
        DC=DCN
        XPE=XPB
        J=(I+1)/2

        IF(L .NE. 2*J) THEN
          RLE=SG(J)/XL(J)

        ELSE
          CF=(XPE/XNR)**2
          ADC=ABS(DC)
          if (ADC.eq.0.0) ADC=1.d-6
          DCC=(DC-US(L))/(1-DC*US(L))

          RQ(J,1)=RQ(J,1)+EP*CF/ADC
          RQ(J,2)=RQ(J,2)+EP*CF*DCC/ADC
          RQ(J,3)=RQ(J,3)+EP*CF*DCC**2/ADC
          RQ(J,5)=RQ(J,5)+EP*DC/ADC
          RQ(J,8)=RQ(J,8)+EP*CF/ADC*XPE
          IF(XPE .GT. XPL) RQ(J,4)=RQ(J,4)+EP*CF/ADC
          RQ(J,9) = RQ(J,9) + 1                       ! count packets
        END IF



80      SN2=1-DC**2
        DSC=-SN2+(XE/XI(I))**2
        ID=1
        IF(DC.LT.0 .AND. DSC.GT.0) ID=-1
        NB=I+(ID+1)/2
        GO TO 45

85      LB=1+(ID+1)/2                         ! IC(1)   Re-enter photosphere
        IC(LB)=IC(LB)+1                       ! IC(2)   Escape to inf  
        PW(3)=PW(3)+EP

        IF(I.EQ.0) then                       ! re-enter  
          nn = nn -1
          GO TO 90
        END IF  

        IF(LE .EQ. 1) IC(3)=IC(3)+1           ! IC(3)   LINE EVENT
        IF(LE .EQ. 0) IC(4)=IC(4)+1           ! IC(4)   NO LINE EVENT  

        PW(2)=PW(2)+EP                        ! escapees to infinity

!       Diagnostics of spectrum formation               

        IF(KNE.NE.0) THEN
          PWL(KNE)=PWL(KNE)+EP                  ! Emitted power in lines 
        END IF


!       Record the spectrum

!       SHChange: ALOG -> LOG
        ZWVMC=AP(16)+4-LOG10(XNR)
        K=(ZWVMC-ZWVM(1))/DZWV+1
        IF(K.GT.KG .OR. K.LT.1) GO TO 90
        NO(K)=NO(K)+1
        XLO(K)=XLO(K)+EP


!       Sum energy emitted in (wv1,wv2) 

        WVV=10**(AP(16)+8)/XNR                            ! Angstrom
        IF(WVV.GT.WVV1 .AND. WVV.LT.WVV2) PW(4)=PW(4)+EP 


90      CONTINUE                              ! End loop over pkts

        WRITE(*,2003)


        SM=((PW(3)-PW(1))-PW(5))/PW(1)
        WRITE(*,2007)SM
        WRITE(3,2007)SM

        WRITE(3,1006)
 1006   FORMAT(6X,'NPK',8X,'ABS',6X,'ESC',6X,'event',6X,'no event',
     &          12X,'out of rg') 
        WRITE(3,1007)NP1,(IC(N),N=1,6)
1007    FORMAT(/7I10/)

        WRITE(*,1008)NP1,IC(2),IC(6)
        WRITE(3,1008)NP1,IC(2),IC(6)
1008    FORMAT(/25X,'No. pkts   =',I8/25X,'Esc to inf =',I8/25X,'No line ev ='I8)

!        line list for spectrum, nlte
        IF (LC .EQ. 2) THEN                     !only on real MC
          do M=1,LN
            KN = IR(M)
            WLAM = 10.D0**(AP(16) + 8.D0) / XN(M)
            IF (WLAM. GT. 2000.) THEN
              WLAM = WAIR(WLAM)
            END IF
          end do            
!         if (KC(10) .EQ. 1) THEN               ! only if branching     
!           CALL IDEM                           !wvl jumps
!         end if 
          IF (ITER .eq. ITT) THEN 
            CALL STRLIN                           !ficst.f
          END IF
          
        END IF      

        write(3,*)' PW(2)= ',PW(2),' RQ(JS,5)= ',RQ(JS,5)
        CLP=(TB/TST)**4/NP1
!       CL=1/PW(2)
        CL=1/RQ(JS,5)                                   !!!PW(2)==RQ(JS,5)
        AP(41)=RQ(JS,5)/NP1
        write(3,*)' TB=',TB,' CLP=',CLP,' CL=',CL
!       SHChange: JZ+1 instead of 51
        RQ(JZ+1,5) = CLP * RQ(JS,5)

        WRITE(3,2005)
2005    FORMAT(/5X,'JP',7X,'LP',7X,'KP',8X,'KJ',9X,'L',8X,'NUJ',8X,'BT',
     &          8X,'NPK'/)
        DO 100 J=1,JS
        I=2*J
        CM=0.25*CL*XI(I)**2

!        normalization
        RQ(J,4)=RQ(J,4)/RQ(J,1)
        RQ(J,8)=RQ(J,8)/RQ(J,1)

        RQ(J,1)=CM*RQ(J,1)
        RQ(J,2)=CL*RQ(J,2)
        RQ(J,3)=CM*RQ(J,3)
        RQ(J,5)=CL*RQ(J,5)

        WRITE(3,1005)(RQ(J,K),K=1,5),RQ(J,8),US(I),J,RQ(J,9)
1005    FORMAT(3F9.5,E11.3,F10.5,E11.3,F9.4,I4,I8)

100     CONTINUE

!       SHChange: JZ+1 instead of 51
        WRITE(3,1991)RQ(JZ+1,5)
1991    FORMAT(F48.5)   

        WRITE(*,2003)


!       if (LC.EQ.2) THEN
        CALL SCALE1(NP1)


        CALL EMSPEC1
!       END If




2001    FORMAT(/5X,'Start ','Monte ','Carlo')
2003    FORMAT(/5X,'Finish ','Monte ','Carlo')
2007    FORMAT(/9X,'En consv',' check',D12.3//)

2009    FORMAT(/5X,'ERROR ','MONTEC')
2011    FORMAT(/5X,'ERROR',I5)
2015    FORMAT(/35X,'DELTA =',F8.2/)



9999    RETURN
        END






        SUBROUTINE SCALE1(NP1)
        SAVE

!       Choose EP0 to reproduce required luminosity in erg/sec.
!       Compute corresponding L+.
!       Scale XJI to physical units.

        INCLUDE 'ficc.f'

        
!       SHChange: ALOG -> LOG
        ZEP0=ZLM+AP(26)-LOG10(PW(2))                    ! Scale to L(bol)

        IF(KC(1).EQ.0) ZEP0=ZLMO+AP(26)-LOG10(PW(4))    !   "   "  L(wv1,wv2)



!       XLUM=PW(2)*10**ZEP0                              ! L(bol) 
!       SHChange: ALOG -> LOG
        ZLUM=LOG10(PW(2)) + ZEP0                              ! L(bol)

!       Corresponding L+ at photosphere

!       XLUMB=NP*10**ZEP0
        xnp = NP1
        ZLUMB=log10(xNP) + ZEP0

        WRITE(*,2017)ZLUM,ZLUMB,ZEP0
        WRITE(3,2017)ZLUM,ZLUMB,ZEP0

        XINB=XLUMB/AP(9)/(XNUA-XNUB)                     ! Intensity I+(R)


!C      Mean intensities at line frequencies
!       SHChange: Comment: ATTENTION when activating this -- WV is already used
!       CF=TM*10**(ZEP0-ZVOL)/4/AP(9)
!       DO 95 KN=1,LN
!       WV=10**AP(16)/XN2(KN)
!       DO 95 J=1,JS
!       XJI(KN,J)=CF*WV/VT(J)*XJI(KN,J)                       ! Jnu (cgs)
!95     CONTINUE


2017    FORMAT(/5X,'L =',D13.4,5X,'L+ =',D13.4,8X,'Log EP0 =',F7.3)

9999    RETURN
        END



        SUBROUTINE EMSPEC1
        SAVE

!       Output of spectrum and diagnostic file of strongly emitting lines.      
!       ZLDM also printed in FICZ into file SPCT.OUD.


        INCLUDE 'ficc.f'

        WRITE(30,2013)ITER


!       Monte Carlo spectrum - vac wavlengths

        NT=0
        DO 110 K=1,KG
        ZLDM(K)=-99
        NT=NT+NO(K)
        IF(XLO(K).NE.0) THEN 
          RAT = XLO(K)/(DWVM(K)*1.D-4)
!         SHChange: ALOG -> LOG          
          ZLDM(K)=LOG10(RAT)+ZEP0
        END IF  
        ZWVMB=(ZWVM(K)+ZWVM(K+1))/2
        WVMB=10**ZWVMB
110     WRITE(30,1011)WVMB,NO(K),XLO(K),NT,ZLDM(K)             ! SICA.OUD


!       Line emission

        DO 150 KN=1,LN
        !FRC=PWL(KN)/PW(2)
        !IF(FRC.LT.1.D-6) GO TO 150

        IF(PWL(KN) .EQ. 0) GO TO 150

        M=IN(LN+1-KN)

        NA=ELM(M)+1.D-3
        NJ=100*(ELM(M)-NA)+1+1.D-3

        WVMC=10**(AP(16)+4)/XN2(KN)     

!       SHChange: JIABS -> ABS
        KKL=ABS(NLW(M))
        KKU=ABS(NUP(M))
        write(31,*)m,kkl,kku


        WRITE(31,1013)WVMC,CHE(NA),CHI(NJ),EL(KKU),XJ(KKU),EL(KKL),XJ(KKL),PWL(KN) ! EICA
150     CONTINUE

        write(*,*)' Finish Emspec '

1011    FORMAT(F16.6,I15,F18.3,I12,F15.4)
1013    FORMAT(F10.5,4X,A3,A5,2(F12.3,F5.1),F10.2)

2013    FORMAT(//5X,'Spectrum',8X,'Iteration =',I4/)

9999    RETURN
        END





        FUNCTION KNCAL1(XPE)
        SAVE

!       Finds index KNCAL of first line in line list with freq less than XPE.

        INCLUDE 'ficc.f'


        KNL=1
        KNU=LN

        IF(XPE .GT. XN2(1)) THEN
        KNCAL1=1
        GO TO 9999
        END IF

        IF(XPE .LT. XN2(LN)) THEN
        KNCAL1=LN                                     ! No redward line
        GO TO 9999                                   ! XN(LN+1)=0 - ok 
        END IF


15      KN=(KNL+KNU)/2


        IF(XN2(KN) .GE. XPE) THEN
        KNL=KN
        ELSE
        KNU=KN
        END IF


        IF(KNU-KNL .EQ. 1) THEN
        KNCAL1=KNU                                       ! Next lwr nu
        GO TO 9999
        ELSE
        GO TO 15 
        END IF


9999    RETURN
        END







        SUBROUTINE ZEROS1
        SAVE

        INCLUDE 'ficc.f'
!       SHChange: KZ instead of 5000 (see bounds in ficc.f)
        DO 20 M=1,KZ
        XLO(M)=0
20      NO(M)=0

        DO 30 L=1,10
        IC(L)=0
!       SHChange: JZ+1 instead of 51
        DO 30 J=1,JZ+1
30      RQ(J,L)=0

        DO 40 KN=1,LN
        PWL(KN)=0
        DO 40 J=1,JS
40      XJI(KN,J)=0 

        DO 50 L=1,5
50      PW(L)=0

9999    RETURN
        END
