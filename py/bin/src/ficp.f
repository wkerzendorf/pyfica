        SUBROUTINE PARTF(L)
        SAVE

!       Approx partition fns calculated using simple model for pop of excited
!       levels.

        INCLUDE 'ficc.f'

        DIMENSION IGT(56,10)


!       Read total st wts of grd terms

        IF (L .EQ. 1) THEN
          DO 20 NA=1,56
            READ(17,1001)ND,CHR,(IGT(NA,NJ),NJ=1,10)           ! NAMA.IND
20        CONTINUE
        END IF

        DO 30 NA=1,56
        NI=NA+1
        IF(NI .GT. 5) NI=5
        DO 30 NJ=1,NI
        DO 30 J=1,JS
        UU(J,NA,NJ)=IGT(NA,NJ)                            ! Default values
30      CONTINUE

        BB=10**(AP(14)+AP(16)-AP(17))                     ! hc/k
        DO 50 J=1,JS
        DO 50 NA=1,56
        NE=NA                                             ! Non-stripped ions 
        IF(NE .GT. 5) NE=5                                ! Ions I-V
        DO 50 NJ=1,NE
        IF(IA(NA,NJ).EQ.0) GO TO 50                       ! Excluded in FICK
        N1=NLV(NA,NJ)
        K1=KS(NA,NJ)
        IF(N1 .EQ. 0) GO TO 50
        UU(J,NA,NJ)=2*XJ(K1)+1                            ! Grd state
        DO 40 N=K1+1,K1+N1-1
          EF=EXP(-BB*ABS(EL(N))/TR(J))
          IF(EL(N) .GT. 0) EF=WD(J)*EF                      ! Normal lvl
          UU(J,NA,NJ)=UU(J,NA,NJ)+(2*XJ(N)+1)*EF
40      CONTINUE
50      CONTINUE


!       Print-out
!       SHChange: Print partfcs only at last iteration, but for every shell
        IF (L.EQ.(NRTEMP+1)) THEN
        DO 80 J=1,JS
        WRITE(42,'(I3,I3)')L,J
        DO 70 I=1,I2
          NA=NAI(I)
          NJ=NJI(I)
          WRITE(42,1003)I,NA,NJ,UU(J,NA,NJ),IGT(NA,NJ)      ! PTFN.OUD
70      CONTINUE
80      CONTINUE
        END IF

        WRITE(*,2001)


1001    FORMAT(I3,A3,10I4)
1003    FORMAT(I8,I14,I4,F20.3,I10)

2001    FORMAT(/5X,'Partition fns',' calculated')

9999    RETURN
        END

