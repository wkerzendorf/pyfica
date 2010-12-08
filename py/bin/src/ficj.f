        SUBROUTINE INPUT
        SAVE

!       Input data

        INCLUDE 'ficc.f'

        DIMENSION XAIN(60),XNR(60),NAT(60)
        character*2 AD


!       Read atomic symbols

        READ(18,1009)(CHE(NA),NA=1,60)
1009    FORMAT(15A3/15A3/15A3/15A3)
        READ(18,1011)(CHI(II),II=1,9)
1011    FORMAT(9A5)


!       Read composition - mass fractions

        NEL=0
!       SHChange: INIT SM
        SM=0.0D0
        DO 50 NA=1,60
!         READ(12,1013,END=51)ND,AD,XA(NA)                  ! COMP.IND
          READ(12,*,END=51)ND,AD,XA(NA)                  ! COMP.IND
          SM=SM+XA(NA)
          NEL=NEL+1
50      CONTINUE
51      WRITE(*,'(5x,"Read abundances for",I3," elements.")')NEL

!        normalize mass to 1.0
        DO NA=1,NEL
          XAIN(NA) = XA(NA)
          XA(NA) = XAIN(NA) / SM
        END DO


!       Read atomic wts, get number Ab.

!       SHChange: INIT SN
        SN=0.0D0
        DO NA=1,NEL
          READ(16,1025)NAT(NA),AD,AW(NA)                   ! NBIA.IND
          XNR(NA) = XA(NA) / AW(NA)
          SN = SN + XNR(NA)
        END DO

!        norm. number Ab.
        DO NA=1,NEL
          XNR(NA) = XNR(NA) / SN
          IF (XA(NA) .ne. 0.0) THEN
            WRITE(46,1022)NAT(NA),CHE(NA),XAIN(NA),XA(NA),XNR(NA)
1022        FORMAT(I9,2X,A3,3F13.6)
          END IF
        END DO
        WRITE(3,1023)SM                                       ! should = 1

!       Read ionization potentials etc

        DO 60 N=1,56
60      READ(11,1015)ND,CHR,(CH(J,N),J=1,9)
        DO 70 N=1,56
70      READ(11,1019)ND,CHR,(ZZ(J,N),J=1,8)


1013    FORMAT(I3,A3,F14.6)
1015    FORMAT(I3,A3,9F8.3)
1019    FORMAT(I3,A3,8F7.2)

1021    FORMAT(I9,2X,A3,F13.4)
1023    FORMAT(/F27.4//)
1025    FORMAT(I3,A3,F9.4)


9999    RETURN
        END

