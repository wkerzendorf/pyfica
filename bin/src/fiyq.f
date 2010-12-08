        SUBROUTINE BRATS
        SAVE

!       Branching ratios for NEWFRQ. NB These branching ratios are for energy
!       not photons - whence factor nu . This routine must be within the
!       iteration loop if level populations change during iterations.

        INCLUDE 'ficc.f'

        DIMENSION SMM(1000)

        
        WRITE(*,2001)


        DO 15 KK=1,KKT
        DO 10 J=1,JS
10      GMD(KK,J)=0
15      CONTINUE


        DO 90 I=1,I2                         ! Loop over ions

        NA=NAI(I)
        NJ=NJI(I)

        KD=KS(NA,NJ)-1
        N1=NLV(NA,NJ)

        M1=IS(I)
        M2=IS(I+1)-1

        DO 80 J=1,JS                         ! Loop over shells

        CALL OCCP(NA,NJ,J)

        DO 35 M=M1,M2                               ! bsc list
!       SHChange: JIABS -> ABS
        KKU=ABS(NUP(M))                           ! Indeces in list EL,XJ   

        KN=IR(M)                                    ! ord list

        BG(KN)=BETA(TAUL(KN,J))

        GMD(KKU,J)=GMD(KKU,J)+AUL(M)*BG(KN)*XN2(KN)

35      CONTINUE



!       Branching probabilities

        DO 55 N=1,N1
55      SMM(N)=0

        DO 60 M=M1,M2                               ! bsc list

!       SHChange: JIABS -> ABS
        KKU=ABS(NUP(M))                           ! Index in list EL,XJ   
        LU=KKU-KD
        KN=IR(M)

        IF(GMD(KKU,J) .EQ. 0) GO TO 60

        QQ(M,J)=AUL(M)*BG(KN)*XN2(KN)/GMD(KKU,J)

!       if (kn.eq.32465) then
!        write(10,*) 'q',m,j,QQ(M,J)
!       end if

        SMM(LU)=SMM(LU)+QQ(M,J)

60      CONTINUE


!       Check sum rule for branching probabilities

        DO 70 LU=2,N1
        KKU=LU+KD
        IF(GMD(KKU,J) .EQ. 0) GO TO 70

        IF(ABS(SMM(LU)-1) .GT. 1.E-5) CALL EXIT          ! QQ's single prec.

70      CONTINUE

80      CONTINUE                             ! End loop over shells

90      CONTINUE                             ! End loop over ions


        WRITE(*,2003)   


1001    FORMAT(I9,I4,I10,I15,D13.3)

2001    FORMAT(/5X,'Start','  branch',' ratios',50X,'FIYQ')
2003    FORMAT(/5X,'Finish ','branch',' ratios')

9999    RETURN
        END



