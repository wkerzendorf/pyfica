        SUBROUTINE RADX
        SAVE

!       Radiative excitation rate coeffs - first calculated with model rad
!       field from FICI. Corresponding decay coeff GMDC calculated in 
!       BRATS (FICQ) outside iteration loop.

        INCLUDE 'ficc.f'

        
        WRITE(*,2001)

        write(*,*)kkt,js

        DO 15 KK=1,KKT
        DO 10 J=1,JS
10      GMX(KK,J)=0
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
        KKL=ABS(NLW(M))

        LU=KKU-KD                                   ! Indeces for ion
        LL=KKL-KD

        KN=IR(M)                                    ! ord list

        CSE=1-GG(LL)/GG(LU)*XNL(LU)/XNL(LL)         ! corr. for stim. em.

        GMX(KKL,J)=GMX(KKL,J)+BLU(M)*CSE*XJI(KN,J)*BETA(TAUL(KN,J))

35      CONTINUE



!       Pre-calculate exc ratio for NEWFRQ (FICN). NB FRX set = 1 if GMD =0 
!       even if GMX = 0 since latter may be due to absence of pkts
!        - low dyn range of MC.


        DO 40 KK=1,KKT
        FRX(KK,J)=1
        IF(GMD(KK,J) .EQ. 0) GO TO 40
        FRX(KK,J)=GMX(KK,J)/(GMX(KK,J)+GMD(KK,J))
40      CONTINUE

80      CONTINUE

90      CONTINUE


        WRITE(*,2003)   


1001    FORMAT(I9,I4,I10,I15,D13.3)

2001    FORMAT(/5X,'Start  ','rad ex rates',51X,'FICR')
2003    FORMAT(/5X,'Finish ','rad ex rates')

9999    RETURN
        END

