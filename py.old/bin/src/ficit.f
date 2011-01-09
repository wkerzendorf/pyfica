        SUBROUTINE TEMPIT(L)
        SAVE
!

        INCLUDE 'ficc.f'

!       SHChange: was: DIMENSION TRO(50)
        DIMENSION TRO(JZ)
!
        TST = 10.d0**ZTEFF
        UC  = 0.7D0
!
        DO 12 J=1,JS
          TRO(J) = TR(J)
          TR(J) = RQ(J,8)/3.8322295D0 * 10.D0**(AP(14)-AP(17))
          TR(J) = (1.D0-UC)*TRO(J) + UC*TR(J)
!         tr(j)=tro(j)
          WD(J) = RQ(J,1)*(TST/TR(J))**4
12      CONTINUE
!
        CALL SMOOTH(TR,JS)
        CALL SMOOTH(WD,JS)
!
        SM=0.D0
        DO 30 J=1,JS
30      SM=SM+ABS((TR(J)-TRO(J))/TRO(J))
        SM=100.D0*SM/JS
!
        DO 45 J=1,JS
45      TS(J)=0.9D0*TR(J)
!
        TB=(1.D0-UC)*TB+UC*TST/AP(41)**0.25
!
        write(3,*)' TEMP IT no:',L
        WRITE(3,782)SM,AP(40),AP(41),TST,TB
782     FORMAT(/F8.1,2F12.4,'   T* = ',F7.1,'   TB = ',F7.1)
!
        WRITE(3,2009)
2009    FORMAT(/8X,'XS',8X,'VS',10X,'LOG RH',7X,'TE',10X,'TR',9X,'W'/)
        DO 35 J=1,JS
          XS = XI(2*J)
          VS(J) = VPH/XS
35      WRITE(3,1009)XS,VS(J),ZRH(J),TS(J),TR(J),WD(J)
1009    FORMAT(5F12.2,F9.3)
!
!       SHChange: convergence criterion
        SM=0.0D0
        DO J=1,JS
          SM=MAX(SM,ABS((TR(J)-TRO(J))/TRO(J)))
        END DO  
        IF (SM.lt.1.25D-2 .and. L.gt.8) L=MAX(NRTEMP-1,L)

99      RETURN
        END



        SUBROUTINE SMOOTH(DV,JS1)
        SAVE
!
!       SHChange: was DIMENSION DV(50),WV(50) -> include ficc instead, JZ accessible
!                 and rename JS (ambiguity) -> JS1
        INCLUDE 'ficc.f'
        DIMENSION DV(JZ),WV(JZ)
!
        J2=JS1-1
        DO 20 J=2,J2
20      WV(J)=(DV(J-1)+2.D0*DV(J)+DV(J+1))/4.D0
        WV(1)=2.D0*WV(2)-WV(3)
        WV(JS1)=2.D0*WV(JS1-1)-WV(JS1-2)
!
        DO 30 J=1,JS1
30      DV(J)=WV(J)
!
        RETURN
        END
