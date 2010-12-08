        SUBROUTINE ATLVLS
        SAVE

!       Read atomic models - exclude absent atoms and high ions.

        INCLUDE 'ficc.f'


!       Criteria for rejection

        DO 30 NA=1,NEL
        DO 20 NJ=1,8
        IA(NA,NJ)=1
20      IF(XA(NA).EQ.0) IA(NA,NJ)=0                   ! Not in mixture
        DO 25 NJ=2,9
        IF(CH(NJ-1,NA) .GT. CHL) IA(NA,NJ)=0          ! Test on lower ion 
25      CONTINUE
30      CONTINUE


!       Read atomic models

        KK=1
        I=0
        DO 70 II=1,200                                ! Loop over ions
        READ(15,1001)NA,NJ,N1

        IF(IA(NA,NJ) .EQ. 1) THEN
        KS(NA,NJ)=KK
        I=I+1
        NAI(I)=NA
        NJI(I)=NJ
        ION(NA,NJ)=I
        END IF


        NN=0
        DO 65 N=1,N1
        READ(15,1003)ND,ED,XJD                          ! ATLV - Kurucz

        IF(IA(NA,NJ).EQ.1) THEN

        EL(KK)=ED
        XJ(KK)=XJD

        IF(ABS(ED).LT.ELL .or. (na.eq.2.and.ED.lt.425000.)) NN=NN+1  !more levels don't work ??? 

        KK=KK+1
        END IF

65      CONTINUE


        NLV(NA,NJ)=NN

        IF(NN.EQ.1) THEN
        IA(NA,NJ)=0                             ! No exc states below limit
        I=I-1
        END IF

        IF(NA.EQ.56 .AND. NJ.EQ.2) GO TO 75     ! Last ion

70      CONTINUE                                ! End loop over ions 




75      KKT=KK-1                                ! Length of EL(KK),XJ(KK)
        I2=I                                    ! No. of ions     

        WRITE(*,2003)KKT
        WRITE(*,2005)I2


!       Print atomic models

        DO 90 I=1,I2
        NA=NAI(I)
        NJ=NJI(I)
        N1=NLV(NA,NJ)
        WRITE(37,1005)CHE(NA),CHI(NJ),N1                  ! ATMD.OUD 

        KD=KS(NA,NJ)-1
        DO 85 N=1,N1
        KK=KD+N
        WRITE(37,1007)N,EL(KK),XJ(KK)
85      CONTINUE
90      CONTINUE



        WRITE(*,2001)


1001    FORMAT(///I8,I4,I10/)
1003    FORMAT(I9,F18.3,F12.1)
1005    FORMAT(//5X,A3,A5,5X,'NLV ='I4/)
1007    FORMAT(I9,F16.3,F8.1)

2001    FORMAT(/5X,'Energy lvls',' read')
2003    FORMAT(/5X,'Total no.',' enrg lvls =',I7)
2005    FORMAT(/5X,'Total no.',' ions =',I4)

9999    RETURN
        END

