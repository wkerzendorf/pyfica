        SUBROUTINE IONEQ(L)
        SAVE

!       Iterate for Ne and ionization fractions. Calc el scattering coeff.

        INCLUDE 'ficc.f'
        
        
        IF (L.EQ.1) THEN                        !densities, 1st guess
          DO 30 J=1,JS
            XNT(J)=0
            DO 20 NA=1,NEL
              XNA(NA,J)=XA(NA)*10**ZRH(J)/(AW(NA)*10**AP(13)) ! No.dens.atom NA
              XNT(J)=XNT(J)+XNA(NA,J)                         ! Total no. dens.
20          CONTINUE
30        CONTINUE
        END IF


!       Calculate Ne

        DO 80 J=1,JS

        XNEL=XNT(J)                                     ! 1st guess

        DO 40 ITR=1,10                                  ! Iterate for Ne

        CALL IONIZ(J,XNEL)

        XNELP=XNEL
        XNEL=0
        DO 35 NA=1,NEL
          DO 35 NJ=2,5
            XNEL=XNEL+(NJ-1)*FR(NJ,NA,J)*XNA(NA,J)
35      CONTINUE

        IF(ABS(XNEL-XNELP)/XNELP .LT. 1.D-4) GO TO 45
40      CONTINUE

45      XNE(J)=XNEL
!       SHChange: ALOG10 -> LOG10
        SG(J)=10**(AP(40)+LOG10(XNE(J))-ZRH(J))             ! El scatt (cm2/g) 


!       SHChange: NRTEMP now variable -> eliminate hardcoded 7
        IF (L.eq.NRTEMP+1) THEN
          WRITE(20,1001)J,XNT(J),XNE(J)                        ! YICA.OUD
1001      FORMAT(I9,5X,'Ntot =',D11.3,5X,'Ne =',D11.3)
        END IF

!       DO 60 NA=1,30
!         IF(XA(NA) .EQ. 0) GO TO 60
!         WRITE(20,1000)NA,CHE(NA),(FR(NJ,NA,J),NJ=1,5)
!1000     FORMAT(I6,2X,A3,5F12.5)
!60     CONTINUE

80      CONTINUE

!       SHChange: NRTEMP now variable -> eliminate hardcoded 7
        IF (L.eq.NRTEMP+1) THEN
         DO NA=1,56
          IF (XA(NA) .GT. 0) THEN
            WRITE(20,1002)NA,CHE(NA),(CHI(I),I=1,4)
!           SHChange: Yhea Headings: larger spaces
1002        FORMAT(/I3,2X,A3,33X,4(4X,A6,2X))
            DO 1031 J=1,JS  !,3
              I=2*J
              VELC=VPH/XI(I)
            WRITE(20,1003)XI(I),VELC,ZRH(J),XNE(J),(FR(NJ,NA,J),NJ=1,4)         
!             SHChange: Yhea: NEED better precision!
1003          FORMAT(F8.2,F9.1,F10.4,D11.3,3X,4D12.4)
1031        CONTINUE
          END IF
         END DO
        END IF

        WRITE(*,2001)
2001    FORMAT(/5X,'Ionization ','calculated')


9999    RETURN
        END








        SUBROUTINE IONIZ(J,XNEP)
        SAVE

!       Calculate ionization fractions at level J given Ne

        INCLUDE 'ficc.f'

        DIMENSION FX(8),AA(12)
        DIMENSION CF(9,60),ZEP(8,60)

!       SHChange: ALOG10 -> LOG10
        ZQ=0.5D0*LOG10(TS(J)/TR(J))
        TH=5040/TR(J)
        THE=5040/TS(J)
        BB=-1.5D0*LOG10(TH)+20.9366D0+ZQ

!       TH1=THE

        DO 90 NA=1,NEL                                   ! Start wi H

!       if (NA .EQ. 8) then  !! fudge ionization of oxygen
!         THE=TH1*.8D0
!       else
!         THE=TH1
!       endif

        IF (NA .LE. 3) THEN
          NJUP = NA
        ELSE
          NJUP = 4
        END IF

        DO 65 NJ=1,NJUP

        IF(CH(NJ,NA) .GE. CH(2,20)) THEN                 !for Ia
          CF(NJ,NA) = 10.D0**(-CH(NJ,NA)*(THE-TH))*TS(J)/TR(J)/WD(J)
        ELSE IF(CH(NJ,NA) .LT. CH(2,20)) THEN
          CF1 = 1.D0-10.D0**(-(CH(2,20)-CH(NJ,NA))*TH)
          ACF = 10.D0**(-CH(2,20)*THE+CH(NJ,NA)*TH)*TS(J)/TR(J)/WD(J)
          CF(NJ,NA)=CF1+ACF
!       ELSE
!         CF(NJ,NA)=1.D0
        END IF
!
        ZET  = 10.D0**ZZ(NJ,NA)
        ETA1 = CF(NJ,NA)*ZET
        ETA2 = WD(J)*(1.D0-ZET)
        ZEP(NJ,NA) = ETA1 + ETA2

!       SHChange: ALOG10 -> LOG10
        ZRHS = -CH(NJ,NA)*TH+BB+LOG10(2*UU(J,NA,NJ+1)/UU(J,NA,NJ))
        FX(NJ)=10**ZRHS * ZEP(NJ,NA)

!      ZRHS=ZZ(NJ,NA)-CH(NJ,NA)*TH+BB+ALOG10(2*UU(J,NA,NJ+1)/UU(J,NA,NJ))
!      FX(NJ)=10**ZRHS

65      CONTINUE

        XD=WD(J)/XNEP

        AA(1)=1.D-10
        YY=AA(1)
        DO 80 NJ=1,NJUP
          AA(NJ+1)=AA(NJ)*FX(NJ)*XD
80      YY=YY+AA(NJ+1)

        DO 85 NJ=1,njup+1
85      FR(NJ,NA,J)=AA(NJ)/YY

!        change NaI ab.
!       DO I=1,JS
!         FR2OLD = FR(2,11,I)
!         FR(1,11,I) =   7.0d-7 * (I/20.d0)**2
!         FR2ADD = fr2old - FR(1,11,I)
!       END DO

90      CONTINUE

! fudge OI->OII
!         FR1OLD = FR(1,8,J)
!         FR2OLD = FR(2,8,J)
!          FR(1,8,I) =   7.0d-7 * (I/20.d0)**2
!         FR(1,8,J) =  FR(1,8,J)*.1
!         FR(2,8,J) =  FR1OLD+FR2OLD - FR(1,8,J)


9999    RETURN
        END

