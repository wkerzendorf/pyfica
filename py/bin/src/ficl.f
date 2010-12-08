        SUBROUTINE LINES
        SAVE

!       Basic (input) line list (LLST.IND) groups lines ion by ion and for each
!       ion according to upper level.
!       Ordered and censored working line lists calculated here + 
!       mappings between lists. BLKL1 and BLKL3 are ordered.

        INCLUDE 'ficc.f'


!       DIMENSION ELW1(50000),XN1(50000),NUP1(50000),NLW1(50000)
!       DIMENSION ZGF1(50000),GL1(50000),GU1(50000),ELW(50000)
        DIMENSION ELW1(LZ),XN1(LZ),NUP1(LZ),NLW1(LZ)
        DIMENSION ZGF1(LZ),GL1(LZ),GU1(LZ),ELW(LZ)
        DIMENSION ZGFU(1200),KU(1200),ZGFL(1200),KL(1200)
        DIMENSION TAUST(LZ)


        WRITE(*,2001)
2001    FORMAT(/5X,'Read & censor',' line list'/)

        IS(1)=1                                         ! Start for ion 1
        M=0                                             ! Label for basc list 
        IJ=0                                            !   "    "  ions
  

        DO 60 II=1,1000                                 ! Loop over ions
        KA=0

        READ(13,1005)NA,NJ,K1                           ! LLST - Kurucz
        !SHChange: Comment out
!        WRITE(10,1005)NA,NJ,K1
1005    FORMAT(//I8,I4,I10/)

          
        DO 20 LM=1,100000
          L1=4*(LM-1)+1
          L2=L1+3
          IF(L1 .GT. K1) GO TO 25
          READ(13,1007)(NUP1(L),NLW1(L),ZGF1(L),L=L1,L2)
20      CONTINUE

25      IF(IA(NA,NJ).EQ.0) GO TO 55       ! Ion excluded


        IF(NA.EQ.14 .AND. NJ.EQ.1) THEN            ! Temporary fix because
        K1=6795                                    ! of levels without
        NLV(NA,NJ)=803                             ! transitions
        END IF

         if (na.eq.2 .and. nj.eq.1) then           !include all He I lines
           NLV(NA,NJ)=737    
         end if
         if (na.eq.2 .and. nj.eq.2) then           !include all He II lines
           NLV(NA,NJ)=13                           
         end if




!       Locate in EL,XJ list

        KD=KS(NA,NJ)-1
        
        IF(EL(KD+1).NE.0) CALL EXIT                    ! Check on FICK

        DO 120 K=1,K1

!       SHChange: JIABS -> ABS
        KKL=ABS(NLW1(K))+KD                    ! Indeces in EL,XJ list
        KKU=ABS(NUP1(K))+KD
        GL1(K)=2*XJ(KKL)+1
        GU1(K)=2*XJ(KKU)+1

        ELW1(K)=EL(KKL)

        DE=ABS(EL(KKU))-ABS(EL(KKL))
        XN1(K)=10**AP(16)*DE
120     CONTINUE


!       Create working line lists  -  ordered + corresponding basic



!       Identify strongest up- and down-ward transition from each level

        DO 30 LV=1,NLV(NA,NJ)
        KL(LV)=0
        KU(LV)=0
        ZGFL(LV)=-99
30      ZGFU(LV)=-99

        DO 40 K=1,K1
!       SHChange: JIABS -> ABS        
        LL=ABS(NLW1(K))
        LU=ABS(NUP1(K))

        IF(LU .GT. NLV(NA,NJ)) GO TO 40                   ! Exc limit

        IF(ZGF1(K) .GT. ZGFU(LU)) THEN
        ZGFU(LU)=ZGF1(K)
        KU(LU)=K
        END IF

        IF(ZGF1(K) .GT. ZGFL(LL)) THEN
        ZGFL(LL)=ZGF1(K)
        KL(LL)=K
        END IF

40      CONTINUE

!       Eliminate weak subordinate lines but keep strongest

        TAUG=10**ZTAUG                                         ! Cut-off

        DO 50 K=1,K1
!       SHChange: JIABS -> ABS
        LL=ABS(NLW1(K))
        LU=ABS(NUP1(K))

        IF(LU .GT. NLV(NA,NJ)) GO TO 50                   ! Excitation limit   

        CALL TAUSM(NA,NJ,XN1(K),ZGF1(K),ELW1(K),TAUM)

        IF(TAUM.GT.TAUG .OR. K.EQ.KU(LU) .OR. K.EQ.KL(LL) 
     &      .or. (na.eq.2 ) ) THEN              ! fiddle to include He
        M=M+1                                      ! No. of lines   
        KA=KA+1                                    ! "    "   "   for ion II

        XN(M)=XN1(K)
        ZGF(M)=ZGF1(K)

!       SHChange: JISIGN -> SIGN
        NLW(M)=SIGN(LL+KD,NLW1(K))
        NUP(M)=SIGN(LU+KD,NUP1(K))
!       write(10,*)k,NLW1(K),NUP1(K)
!       write(10,*)m,kd,NLW(M),NUP(M)

        ELW(M)=ELW1(K)
        GL(M)=GL1(K)
        GU(M)=GU1(K)
        ELM(M)=NA+0.01D0*(NJ-1)
        TAUST(M) = TAUM                                 !store tau for LHEA.dat

        CALL ECFFS(M)

        END IF

50      CONTINUE


        WRITE(*,2003)CHE(NA),CHI(NJ),K1,KA
2003    FORMAT(8X,A3,A5,5X,'No. of lines =',I6,3X,'Included =',I6)

        IF(KA.EQ.0) GO TO 55

        IJ=IJ+1


        ION(NA,NJ)=IJ                                   ! Mapping for ions

        NAI(IJ)=NA                                      ! Inverse mapping
        NJI(IJ)=NJ                                    

        IS(IJ+1)=M+1                                    ! Start for next ion

55      IF(NA.EQ.56 .AND. NJ.EQ.2) GO TO 65             ! Last ion

60      CONTINUE                                        ! End loop over ions



65      LN=M                                            ! No. of lines
        I2=IJ                                           !  "  "  ions  
        IS(I2+1)=LN+1

        WRITE(*,2005)I2,LN
        WRITE(3,2005)I2,LN
        WRITE(46,2005)I2,LN
        IF (LN .GT. LZ) stop 'Too many lines - check TAU limit!'


!       Locate lines with same upper level

        DO 70 M=1,LN
!       SHChange: JIABS -> ABS
        LU=ABS(NUP(M))                                ! Index for EL,XJ list
        IF(MS(LU) .EQ. 0) MS(LU)=M
        IF(MF(LU) .LT. M) MF(LU)=M  
70      CONTINUE



!       Create ordered line list for Monte Carlo calculation. Because nu
!       decreases, if MC asks for line after XN2(LN) it finds zero freq - OK

!       SHChange: explicitly give borders
        CALL SORT1(XN(1:LN),LN,IN(1:LN))

        DO 80 LL=1,LN
        M=IN(LL)
        KN=LN-LL+1                                         ! nu decreasing
        IR(M)=KN                                           ! invert sort
        XN2(KN)=XN(M)
        ELW2(KN)=ELW(M)
80      CONTINUE


        do  m=1,ln
          kn=ir(m)
          WLAM = 10.D0**(AP(16) + 8.D0) / XN(M)
          IF (WLAM. GT. 2000.) THEN
            WLAM = WAIR(WLAM)
          END IF
          if (wlam.le.0.d0) then
           write(3,*)' wair ',elm(m),xn(m),wlam
          end if
          WRITE(21,1030)M,kn,XN(M),WLAM,ELM(M),ZGF(M),NLW(M),NUP(M),
     &                  TAUST(M)                                ! LHEA.dat
1030      FORMAT(2I7,E13.6,F11.3,F8.2,F9.3,2I7,E10.3)
        end do


!       Identify effective resonance lines

        LR=0
        DO 90 KN=1,LN
        IDL(KN)=1                                  ! Default - non-res ln
        M=IN(LN+1-KN)
!       SHChange: JIABS -> ABS
        LU=ABS(NUP(M))
        IF(MS(LU) .NE. MF(LU)) GO TO 90

        IF(NLW(M) .LT. 0) THEN
        IDL(KN)=0                                  ! Eff res line
        LR=LR+1
        END IF

90      CONTINUE

        WRITE(*,2011)LR

        WRITE(*,2009)
        write(*,*)' lines',nlw(1),nup(1)



1003    FORMAT(I8,E15.6,F12.2,F14.3,I12,2I7)
1007    FORMAT(4(2I6,F8.3))

2005    FORMAT(/5X,'No. of ions =',I4,12X,'Total No. ','Included =',I6/)
2009    FORMAT(/5X,'Line list',' read')

2011    FORMAT(/25X,'No. eff res',' lines ='I6/)

9999    RETURN
        END






        SUBROUTINE ECFFS(M)
        SAVE

!       Einstein coeffs for transition M in basic list

        INCLUDE 'ficc.f'

        WVMC=10**(AP(16)+4)/XN(M)
!       SHChange: ALOG -> LOG 
        CF=2*10**(AP(14)+3*LOG10(XN(M))-2*AP(16))

        AUL(M)=0.6670D8*10**ZGF(M)/GU(M)/WVMC**2
        BUL(M)=AUL(M)/CF
        BLU(M)=BUL(M)*GU(M)/GL(M)


9999    RETURN
        END






        SUBROUTINE TAUSM(NA,NJ,XN1,ZGF1,ELW1,TAUM)
        SAVE

!       Max Sobolev opt depths - stim em neglected

        INCLUDE 'ficc.f'


        BB=10**(AP(14)+AP(16)-AP(17))
        BT=0.02654D0*10**(ZGF1+AP(16))/XN1

        MT=2
        IF(ELW1 .GT. 0) MT=3

        TAUM=0
        DO 30 J=1,JS
        EF=EXP(-BB*ABS(ELW1)/TR(J))                             ! NB abs


        IF(MT .EQ. 3) EF=WD(J)*EF                              ! Normal lvl
        
        TAU=BT*EF/UU(J,NA,NJ)*FR(NJ,NA,J)*XNA(NA,J)*TM

        TAUM=DMAX1(TAU,TAUM)

30      CONTINUE


9999    RETURN
        END

