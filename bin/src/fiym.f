        SUBROUTINE MONTEC(LC)
!       SHChange: Declare Mersenne Twister,SAVE
        use mtprng
        SAVE

!       Black body emission at photosphere.

!REAL*4 T1,DELTA


        INCLUDE 'ficc.f'

!       SHChange: DIMENSION US(101),XL(50) [US appears also in ficst (and ficg) -- no common!]
        DIMENSION US(2*JZ+1),XL(JZ)
        DIMENSION KNIN(NP*10),KNOUT(NP*10),IQ(NP*10)
!       SHChange: sh packet counting
!        INCLUDE 'shpc-varinit.f'

!T1=SECNDS(0.0)

        WRITE(*,2001)

        CALL ZEROS
        
!       SHChange: Init LT (because readout to LTOLD fails otherwise), NBR
        LT=0
        NBR=0

        U1=VPH*1.D5/10**AP(16)

        EP0=1                                     ! Rescaled below  

        NP1 = NP
        IF (LC .EQ. 0) THEN                     ! TEMP IT
          NP1 = NP 
        ELSE IF (LC .EQ. 1) THEN                     ! Trial soln for GMX
          NP1 = NP/2
        ELSE IF (LC .EQ. 2 .or. LC .EQ. 4) THEN                     ! Trial soln for GMX
          NP1 = NP
        ELSE IF (LC .EQ. 3) THEN                     ! Trial soln for GMX
          NP1 = 1000000
        END IF
        PW(1)=NP1

        DO 25 J=1,JS
25      XL(J)=10**(-ZRH(J)-ZRDCM)

        DO 30 I=1,I1
30      US(I)=U1/XI(I)

        XM=0

        CR=10**(AP(17)-AP(14))
        TST = 10.d0**(ZTEFF)
        CB=AAA*TB*CR
        XPL=5*CR*TB
        KT=NP1/(KB-1)                            ! NP min = KB-1
!-----------------------

!-----------------------

        NN=0
        DO 90 KK=2,KB
        
        ! SHChange: Print out iter %
        IF (LC.eq.4 .and. MOD(INT(ANINT((KK-1.0)/(KB-1.0)*300))/3.0,10.0)==0) 
     $    WRITE(*,'(I5,2X,A)') INT(ANINT((KK-1.0)/(KB-1.0)*100)),'%'

        DO 90 II=1,KT
        NN=NN+1
        I=0
        ID=1
        LE=0
        KNE=0
!       SHChange: sh packet counting
!        INCLUDE 'shpc-varinit-iter.f'
        IF (LC.eq.4 .and. WHIST.eq.1) WRITE(50,'(A,I6,A)',advance='NO') 'N',NN,'; '

        XB=1
        UB=U1


!       Emission at photosphere from static black body

        KP=KK-1
	IF (USEMT==1) THEN
          DCN=SQRT(1-MTRAN(MTSTAT))
          ZP=ZK(KP)+MTRAN(MTSTAT)*(ZK(KP+1)-ZK(KP))
          IF(KP.EQ.1) ZP=ZK(2)*SQRT(1-MTRAN(MTSTAT))
          IF(KP.EQ.KB-1) ZP=1-(1-ZK(KB-1))*SQRT(1-MTRAN(MTSTAT))
	ELSE
	  DCN=SQRT(1-RAN(KR))
          ZP=ZK(KP)+RAN(KR)*(ZK(KP+1)-ZK(KP))
          IF(KP.EQ.1) ZP=ZK(2)*SQRT(1-RAN(KR))
          IF(KP.EQ.KB-1) ZP=1-(1-ZK(KB-1))*SQRT(1-RAN(KR))
	END IF

        
	
        XNR=CB*SQRT(ZP/(1-ZP))

        XPB=XNR*(1-DCN*U1)

        XM=XM+XNR*10**(AP(14)-AP(17))/TB/NP1            ! Mean freq


        EP=EP0
        
        KN=KNCAL(XPB)                                  ! Location in line list

        GO TO 65

45      IF (USEMT==1) THEN
          TAUR=-LOG(1-MTRAN(MTSTAT))
        ELSE
          TAUR=-LOG(1-RAN(KR))
        END IF
	
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

!       SHChange: more verbose error reporting, distinguish hints
        IF(SB.LT.0) WRITE(*,2009)
        IF(SB.LT.0) WRITE(*,'(A,F11.8,A,F11.8,A,F11.8,A,F11.8)') '-DC:',-DC,'  XE:',XE,' DCN:',DCN,'  XB:',XB
        IF(SB.EQ.0) WRITE(*,2019)


        CFB=-1/U1/XNR

48      ST=SB
!       ST: Strecke bis zum nächsten Event, SB: Strecke bis Border
        if (lt.ne.1) then 
         ltold=lt
        end if
        LT=1
!       DEFAULT LT = 1

!       SE: Strecke bis El. scattering
!	RLE: Freie Weglänge
        SE=SP+TAUS/RLE          ! falls ESP 48: SP_old<>0. SE=sp_old+(-tauq-sp_old*rle_old)/rle=-tauq/rle (rle nicht aktualisiert)
        IF(SE.LT.ST)THEN
          ST=SE
          LT=2
        END IF

!	SP: Strecke bis line encounter
        SP=CFB*(XN2(KN)-XPE)

        IF(SP.LT.0) WRITE(*,2099)



        IOLD = I !SHChange: additional diagnosis variable

        IF(SP.LT.ST)THEN

          XJI(KN,J)=XJI(KN,J)+EP
        
          TAUQ=TAUQ+TAUL(KN,J) !TAUQ ist tatsächliche Summe der line-Taus
          TAUS=-TAUQ-SP*RLE  !TAUS ist nur da, um line-Overflow zu testen
!         SHChange: sh packet counting
!          INCLUDE 'shpc-ctencounters.f'

! SHChange: If no line interaction (taus.gt.0), record
!          INCLUDE 'shpc-notepassthrough.f'

          KN=KN+1                                  ! next line
          IF(TAUS .GT. 0) GO TO 48                 ! proceed with next line , IF no overflow
!         PROCEED IF TAUS < 0 (line KNold caused overflow)
          ST=SP
          LT=4
          KN=KN-1                                  ! Line absorbing pkt
        END IF

        IF(LT.NE.1)THEN
!         LT.NE.default -> scattering / branching
          IF (LC.eq.4) THEN 
            IF(10**(8+AP(16))/XN2(KN).lt.100000) THEN
              IF (WHIST.eq.1) WRITE (50,'(A,I2,A,F7.5,A,I2,A,F8.5,A,F8.2,A)',advance='NO') 
     $        'A ',J,' ',XE,' ',ID,' ',DC,' ',10**(8+AP(16))/XN2(KN),';  '
            ELSE
              IF (WHIST.eq.1) WRITE (50,'(A,I2,A,F7.5,A,I2,A,F8.5,A,F9.2,A)',advance='NO') 
     $        'A ',J,' ',XE,' ',ID,' ',DC,' ',10**(8+AP(16))/XN2(KN),';  '
            END IF             
          END IF
          RE=1/XE
          SQ=ST+DC*RE
          XE=1/SQRT(RE*(RE+DC*ST)+ST*SQ)

          DI=SQ*XE
          UE=U1/XE

          IF (USEMT==1) THEN
            DC=-1+2*MTRAN(MTSTAT)
          ELSE
            DC=-1+2*RAN(KR)
          END IF

!  DC=(DC+UE)/(1+DC*UE)                         !neglect terms O(v/c)
          EF=(1-DI*UE)/(1-DC*UE)
!  SHChange: store XNR in XNROLD
          XNROLD=XNR
          XNR=XNR*EF
          XPE=XNR*(1-DC*UE)
          EPI=EP
          EP=EPI*EF
          PW(5)=PW(5)+(EP-EPI)                        ! mech goes to rad energy

!         SHChange: Record shell no. for scattering / branching event (last persists if emitted)
          JN1 = J

          IF (LT.ne.4 .and. LC.eq.4 .and. WHIST.eq.1) WRITE (50,'(A,F7.5,A,I2,A,F8.5,A)',advance='NO') 'S ',XE,' ',ID,' ',DC,';  '
	  
          IF(LT .EQ. 4) THEN
             
            KN1 = KN
            CALL NEWFRQ(KN,XPE,EP,J)                           ! FIyN, branching
	    
!           SHChange: Record shell no. for true branching event (last persists if emitted)
            JNL1=J 
!           SHChange: Write MC history
            IF (LC.eq.4) THEN
              IF (10**(8+AP(16))/XN2(KN).lt.100000) THEN
                IF (WHIST.eq.1) WRITE (50,'(A,F7.5,A,I2,A,F8.5,A,F8.2,A)',advance='NO') 
     &            'B ',XE,' ',ID,' ',DC,' ',10**(8+AP(16))/XN2(KN),';  '
              ELSE              
                IF (WHIST.eq.1) WRITE (50,'(A,F7.5,A,I2,A,F8.5,A,F9.2,A)',advance='NO') 
     &            'B ',XE,' ',ID,' ',DC,' ',10**(8+AP(16))/XN2(KN),';  '
              END IF  
            END IF  

!           store data for matrix, only for branching -- NBR (max. NP) events
!           SHChange: Increase maximum store count
            IF (KC(2).EQ.1) THEN
             IF (ITER.eq.ITT .and. NBR.lt.NP*10) THEN     
              NBR = NBR + 1 
              IF (NBR .LE. NP*10) THEN
                KNIN(NBR) = KN1
                KNOUT(NBR) = KN
              END IF
             END IF
            END IF 

            KNE=KN                                             ! Emitted line
            KN=KN+1                                            ! next line
            XNR=XPE/(1-DC*UE)

!           SHChange: sh packet counting
!            INCLUDE 'shpc-ctlineia.f'

            LE=1                                               ! Line abs-em
          END IF          !END IF LT==4

          GO TO 80
        END IF            !END IF LT<>1

!       65: Leave half-shell, new half-shell
65      L=I+(ID+1)/2
        I=I+ID
        IF(I.EQ.0 .OR. I.EQ.I1) GO TO 85                   ! Leave env

        XE=XB
        DC=DCN
        XPE=XPB
        J=(I+1)/2
!       Abrundung: for I=1,2,3,...: J=1,1,2,2,3,3,...: IOLD=1-ID,2-ID,3-ID,...: ID=1 -> L=1,2,3, ID=-1 -> L=2,3,4 

        IF(L .NE. 2*J) THEN
          RLE=SG(J)/XL(J)
!         ELSE means: L.EQ.2*J means: forward -> I(new) even, backward -> I(new) even
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
!         SHChange: sh packet counting
          IF (LC.eq.4 .and. WHIST.eq.1) WRITE (50,'(A,I6,/)') 'R',NN
!          INCLUDE 'shpc-ctreabsorbed.f'
          nn = nn -1
          GO TO 90
        END IF  

        IF(LE .EQ. 1) IC(3)=IC(3)+1           ! IC(3)   LINE EVENT
        IF(LE .EQ. 0) IC(4)=IC(4)+1           ! IC(4)   NO LINE EVENT  

        PW(2)=PW(2)+EP                        ! escapes to infinity
!       SHChange: sh packet counting
!        INCLUDE 'shpc-ctemiss.f'
        IF (LC.eq.4) THEN 
          IF (10**(8+AP(16))/XNR.lt.100000) THEN
            IF (WHIST.eq.1) WRITE (50,'(A,I6,A,F8.2,/)') 'E',NN,' ',10**(8+AP(16))/XNR
          ELSE
            IF (WHIST.eq.1) WRITE (50,'(A,I6,A,F9.2,/)') 'E',NN,' ',10**(8+AP(16))/XNR
          END IF
        END IF
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
!       SHChange: sh packet counting
!        INCLUDE 'shpc-writeout.f'

        XME=3.8322295D0
        FRE=(XM-XME)/XME
        WRITE(*,2013)XM,FRE                   ! Accuracy of Blck bdy sampling 


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
        IF (ITER .eq. ITT) THEN                 !only on real MC
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
          ! SHChange: only 1 IF ITER.eq.ITT suffices
          CALL STRLIN                           !ficst.f
            
!        GO TO 99  ! NOT build branching matrix
!        build branching matrix
!        SHChange: Routine replaced, introduced WBRAMA param. (old: !!)
          IF (KC(2).EQ.1 .and. WBRAMA==1) THEN                ! do for branching only

!!            CALL SORT1(KNOUT,NBR,IQ)

            DO K=1,NBR
              
!!              N=IQ(K)
!!              M=IN(LN+1-KNOUT(N))
!!              NA=ELM(M)+1.E-3
!!              NJ=100*(ELM(M)-NA)+1+1.E-3
!!              if (na.eq.1) then                 ! do HI only
!!              XNOUT=XN2(KNOUT(N))
!!              XNIN=XN2(KNIN(N))
!!              WVAOUT=10**(AP(16)+8)/XNOUT
!!              WVAIN=10**(AP(16)+8)/XNIN
!!              IF (WVAIN. GT. 2000.) THEN
!!               WVAIN = WAIR(WVAIN)
!!              END IF
!!              IF (WVAOUT. GT. 2000.) THEN
!!              WVAOUT = WAIR(WVAOUT)
!!              END IF
!!              WRITE(10,1111)WVAOUT,NA,NJ,WVAIN
!! 1111         FORMAT(F10.2,I10,I4,F10.2) 

              M=IN(LN+1-KNOUT(K))
              NA=ELM(M)+1.E-3
              NJ=100*(ELM(M)-NA)+1+1.E-3
              WVAOUT = 10**(AP(16)+8)/XN2(KNOUT(K))
              WVAIN = 10**(AP(16)+8)/XN2(KNIN(K))
              IF (WVAOUT.gt.2000) WVAOUT=WAIR(WVAOUT)
              IF (WVAIN.gt.2000) WVAIN=WAIR(WVAIN)
              
              IF (NA.ge.14 .and. WVAOUT.lt.100000 .and. WVAIN.lt.100000) THEN
                WRITE(10,'(1X,F9.3,1X,F9.3,1X,I2,1X,I1)') WVAOUT,WVAIN,NA,NJ
              ELSE IF (WVAOUT.lt.100000 .and. WVAIN.lt.100000) THEN
                WRITE(10,'(1X,F9.3,1X,F9.3)') WVAOUT,WVAIN
              END IF

!!            end if ! end if (na.eq.1)

            END DO             
          END IF
            
        ! SHChange: only 1 IF ITER.eq.ITT suffices
        END IF

 99     write(3,*)' PW(2)= ',PW(2),' RQ(JS,5)= ',RQ(JS,5)
        CLP=(TB/TST)**4/NP1
!       CL=1/PW(2)
        CL=1/RQ(JS,5)                                   !!!PW(2)==RQ(JS,5)
        AP(41)=RQ(JS,5)/NP1
        write(3,*)' TB=',TB,' CLP=',CLP,' CL=',CL
!       SHChange: 51 -> JZ+1
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
1005    FORMAT(3F9.5,E11.3,F10.5,E11.3,F9.4,I4,ES10.3)

100     CONTINUE

!       SHChange: 51 -> JZ+1
        WRITE(3,1991)RQ(JZ+1,5)
1991    FORMAT(F48.5)   

!       SHChange: Deallocate package counting vars
!        INCLUDE 'shpc-deallocate.f'

        WRITE(*,2003)


!       if (LC.EQ.2) THEN
        CALL SCALE(NP1)


        CALL EMSPEC
!       END If




2001    FORMAT(/5X,'Start ','Monte ','Carlo')
2003    FORMAT(/5X,'Finish ','Monte ','Carlo')
2007    FORMAT(/9X,'En consv.',' check',D12.3//)
!       SHChange: verbose reporting (was all ERROR MONTEC)
2009    FORMAT(/5X,'ERROR ','MONTEC - Strecke<0')
2019    FORMAT(/5X,'HINT ','MONTEC - Strecke=0')
2099    FORMAT(/5X,'ERROR ','MONTEC - SP<0')

2011    FORMAT(/5X,'ERROR',I5)
2013    FORMAT(/5X,'XM =',F9.6,4X,'FRE ='D10.2)
2015    FORMAT(/35X,'DELTA =',F8.2/)



9999    RETURN
        END






        SUBROUTINE SCALE(NP1)
        SAVE
        
!       Choose EP0 to reproduce required luminosity in erg/sec.
!       Compute corresponding TB and L+.
!       Scale XJI to physical units.

        INCLUDE 'ficc.f'

        
!       SHChange: ALOG -> LOG
        ZEP0=ZLM+AP(26)-LOG10(PW(2))                    ! Scale to L(bol)

        IF(KC(1).EQ.0) ZEP0=ZLMO+AP(26)-LOG10(PW(4))    !   "   "  L(wv1,wv2)



!       XLUM=PW(2)*10**ZEP0                              ! L(bol) 
!       SHChange: ALOG -> LOG
        ZLUM=log10(PW(2)) + ZEP0                              ! L(bol)

!       Corresponding Tb and L+ at photosphere

!       XLUMB=NP*10**ZEP0
        xnp = NP1
        ZLUMB=log10(xNP) + ZEP0
        ZTB=AP(28)+0.25D0*(ZLUMB-AP(26))-0.5D0*(ZRDCM-AP(27))
        TB=10**ZTB      

        WRITE(*,2017)ZLUM,ZLUMB,TB,ZEP0
        WRITE(3,2017)ZLUM,ZLUMB,TB,ZEP0
!       SHChange: Increase EP0 precision (was 7.3)
2017    FORMAT(/5X,'Scale:  L =',F7.3,5X,'L+ =',F7.3,//15X,'TB1 =',F8.1,10X,'Log EP0 =',F10.6)


9999    RETURN
        END












        SUBROUTINE EMSPEC
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


        WRITE(31,1013)WVMC,CHE(NA),CHI(NJ),EL(KKU),XJ(KKU),EL(KKL),XJ(KKL),PWL(KN)! EICA
150     CONTINUE

        write(*,*)' Finish Emspec '

1011    FORMAT(F16.6,I15,F18.3,I12,F15.4)
1013    FORMAT(F10.5,4X,A3,A5,2(F12.3,F5.1),F10.2)

2013    FORMAT(//5X,'Spectrum',8X,'Iteration =',I4/)

9999    RETURN
        END





        FUNCTION KNCAL(XPE)
        SAVE
        
!       Finds index KNCAL of first line in line list with freq less than XPE.

        INCLUDE 'ficc.f'


        KNL=1
        KNU=LN

        IF(XPE .GT. XN2(1)) THEN
        KNCAL=1
        GO TO 9999
        END IF

        IF(XPE .LT. XN2(LN)) THEN
        KNCAL=LN                                     ! No redward line
        GO TO 9999                                   ! XN(LN+1)=0 - ok 
        END IF


15      KN=(KNL+KNU)/2


        IF(XN2(KN) .GE. XPE) THEN
        KNL=KN
        ELSE
        KNU=KN
        END IF


        IF(KNU-KNL .EQ. 1) THEN
        KNCAL=KNU                                       ! Next lwr nu
        GO TO 9999
        ELSE
        GO TO 15 
        END IF


9999    RETURN
        END







        SUBROUTINE ZEROS
        SAVE

        INCLUDE 'ficc.f'

!       SHChange: 5000 -> KZ (see bounds in ficc.f)
        DO 20 M=1,KZ
        XLO(M)=0
20      NO(M)=0

        DO 30 L=1,10
        IC(L)=0
!       SHChange: 51 -> JZ+1
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


!        INCLUDE 'shpc-chkarr.f'
