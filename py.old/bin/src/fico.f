        SUBROUTINE OCCP(NA,NJ,J)
        SAVE

!       Approximate level populations for ion NA,NJ in shell J.

        INCLUDE 'ficc.f'


        BB=10**(AP(14)+AP(16)-AP(17))                     ! hc/k


        K1=KS(NA,NJ)
        N1=NLV(NA,NJ)
        KD=K1-1

        IF(EL(K1) .NE. 0) CALL EXIT                       ! grd state ?

        SM=0
        DO 60 N=1,N1
        KK=KD+N
        GG(N)=2*XJ(KK)+1
        YY=BB*ABS(EL(KK))/TR(J)

        EXF=EXP(-YY)
        IF(EL(KK) .GT. 0) EXF=WD(J)*EXF                   ! Normal lvl

        XNL(N)=GG(N)*EXF
        SM=SM+XNL(N)
60      CONTINUE

        DO 65 N=1,N1
65      XNL(N)=XNL(N)/SM                                  ! Normalization


!       Print-out

        !IF(J.NE.10) GO TO 9999

        !WRITE(10,1001)NA,NJ,N1

        !DO 80 N=1,N1
        !KK=KD+N
        !WRITE(10,1003)N,EL(KK),GG(N),XNL(N)
        !80     CONTINUE


1001    FORMAT(//I9,I4,I14/)
1003    FORMAT(I7,F18.3,F10.1,D19.4)

9999    RETURN
        END

