        SUBROUTINE STRLIN
        SAVE
!
!        writes list of relevant lines and obs. wvl.


        INCLUDE 'ficc.f'
!       SHChange: Was US(101) [appears also in MC routines - no common]
        DIMENSION US(2*JZ+1)
!
!
        TAUMIN = 10.E0**(ZTAUG)
!
        U1=VPH*1.E5/10**AP(16)
        DO 20 I=1,I1
          US(I)=U1/XI(I)
20      CONTINUE
!
        DO 90 KN=1,LN

!          find wavelength
          WVA=10.E0**(AP(16)+8.E0)/XN2(KN)                      ! vacuum
!        back to air wvl's
          IF (WVA. GT. 2000.) THEN 
            XLV = WAIR(WVA)
          ELSE 
            XLV = WVA 
          END IF
!               
!          find atomic species
          M=IN(LN+1-KN)                                  ! Index in basic list
          NA=ELM(M)+1.E-3
          NJ=100*(ELM(M)-NA)+1+1.E-3
!
!          find optical depth
          IF (XLV.GT.WVM1*1.E4 .AND. XLV.LT.WVM2*1.E4) THEN !range of lines included
            SA=0.E0
            SB=0.E0
            DO 93 J=1,JS
              IM=2*J
              TAU=aMIN1(20.E0,TAUL(KN,J))
!             SHChange: replace goto with cycle and continue with end do to avoid compiler warnings 
              IF (TAU .LT. 0.D0) CYCLE
              DPT=1.D0-EXP(-TAU)
              DLM=XLV*(US(IM+1)-US(IM-1))
              WVL=XLV*(1.D0-US(IM))
              SA=SA+WVL*DPT*DLM
              SB=SB+DPT*DLM
!             if (na.eq.2.and.(4685.d0.lt.xlv.and.xlv.lt.4686.d0))then
!               write(46,3333)j,wvl,TAUL(KN,J),tau,sb
!             end if
   93       END DO
   
            IF (SB .GT. 0.0) THEN 
              XLOB = SA/SB
            END IF

!           go to 999  SHChange: deactivate go to 999, write details for all IME and Ti       
            if ((na.eq.14 .or. na.eq.16 .or. na.eq.20) .and. nj.eq.2) then   !write details of tau
!           SHChange: Comment: nj is the ionization stage (2=1+), na the element number
!                     Change odep criterion to have 4130.87 and ~5800 Ti lines in any case
              IF ((SB.GT.1.0E-2) .and. (XLV.ge.WV1 .and. XLV.le.WV2)) THEN    
                SA=0.E0
                SB=0.E0
                DO J=1,JS
                  IM=2*J
                  TAU=aMIN1(20.E0,TAUL(KN,J))
!                 SHChange: replace goto with cycle to avoid compiler warnings
                  IF (TAU .LT. 0.D0) CYCLE
                  DPT=1.D0-EXP(-TAU)
                  DLM=XLV*(US(IM+1)-US(IM-1))
                  WVL=XLV*(1.D0-US(IM))
                  SA=SA+WVL*DPT*DLM
                  SB=SB+DPT*DLM
!!                 SHChange: was: write(46,3333)j,wvl,TAUL(KN,J),tau,sb -> format(...4f...)
!                  write(46,3333)j,wvl,tau,DPT,(DPT*DLM),sb
3333              format(i4,5f13.7)
  99            END DO
              END IF
            end if
 999        continue        
            

!       print linelist to SBIB.DAT 

!        selection criteria for printing 
!           IF (SB .GT. 3.E0 .and. (XLV.ge.WV1 .and. XLV.le.WV2)) THEN
!!           SHChange: always show 6012.8 Ti II line and all Ti II lines
!            IF ((SB.GT.1.0E0 .or. abs(xlv-6012.8)<0.4 .or. (na==22 .and. nj==2)) .and. (XLV.ge.WV1 .and. XLV.le.WV2)) THEN
            IF (SB.GT.1.0E0 .and. (XLV.ge.WV1 .and. XLV.le.WV2)) THEN
!           IF (SB .GT. 1.E-1 .and. (XLV.ge.WV1 .and. XLV.le.WV2)) THEN
        
!            energy of lower level
              ENLOW = ELW2(KN) * 10.D0**(AP(14)+AP(16)-AP(21))
                       
              WRITE (46,1000) SB,XLOB,XLV,CHE(NA),CHI(NJ),
     &                  ENLOW,ELW2(KN),ZGF(M)
 1000         FORMAT(F9.3,F8.1,F12.2,3X,A3,A5,F10.2,F12.2,F9.3)
!             WRITE (46,1000) SB,XLOB,XLV,NA,NJ,
!     &                  ENLOW,MT(M),ZGF(M)
! 1000         FORMAT(F9.3,F8.1,F12.2,3X,I3,I5,F10.2,I5,F9.2)
            END IF
          END IF
   90   CONTINUE
!
        RETURN
        END    
