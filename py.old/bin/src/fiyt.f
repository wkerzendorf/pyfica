        SUBROUTINE TAUS
        SAVE

!       Sobolev optical depths.
!       This calc outside iteration loop if level pops fixed.


        INCLUDE 'ficc.f'

        write(*,*)' taus',nlw(1),nup(1)

        DO 80 I=1,I2

        NA=NAI(I)
        NJ=NJI(I)

!       write(*,*)I,I2,NA,NJ
        N1=NLV(NA,NJ)
        KD=KS(NA,NJ)-1

        M1=IS(I)
        M2=IS(I+1)-1

        DO 70 J=1,JS
!write(*,*)J,JS

        CALL OCCP(NA,NJ,J)

        CFN=10**(AP(14)+AP(16))/4/AP(9)*FR(NJ,NA,J)*XNA(NA,J)*TM

        DO 50 M=M1,M2                                   ! indx for bsc list
!write(*,*)M1,M,M2

        KN=IR(M)                                        ! indx for ord list

!       SHChange JIABS -> ABS
        LU=ABS(NUP(M))-KD                             ! Indeces for ion
        LL=ABS(NLW(M))-KD
!write(*,*)NUP(M),NLW(M),KD
!write(*,*)LL,LU

        BF=GG(LL)/GG(LU)*XNL(LU)/XNL(LL)

        TAUL(KN,J)=CFN*BLU(M)*XNL(LL)*(1-BF)
!        SHChange: Comment (unreliable?): M: line index ( ->KN ), I: Species index (-> NA, NJ, [M1, M2] = M borders)


!CCCCCCCCCC   Fiddle with tau's   ccccccccccccccccccc

!        if ((na.eq.14) .and. nj.eq.1) then
!                 TAUL(KN,J) = TAUL(KN,J) * 6.e3
!        end if

!        1994I-0404
!        if ((na.eq.2) .and. nj.eq.1) then
!                 TAUL(KN,J) = TAUL(KN,J) * 5.e0
!        end if

!     1994I - 0415
      if(.false.)then                                 
         if ((na.eq.2) .and. nj.eq.1) then
           if(J.GT.js-1)then
                  TAUL(KN,J) = TAUL(KN,J) * 1.d12
           elseif(J.GT.JS-2)then
                  TAUL(KN,J) = TAUL(KN,J) * 3.d11
           elseif(J.GT.JS-3)then
                  TAUL(KN,J) = TAUL(KN,J) * 1.d11
           elseif(J.GT.JS-4)then
                  TAUL(KN,J) = TAUL(KN,J) * 2.d10
           elseif(J.GT.JS-5)then
                  TAUL(KN,J) = TAUL(KN,J) * 1.d09
!          elseif(J.GT.JS-6)then
!                 TAUL(KN,J) = TAUL(KN,J) * 5.d9
!          elseif(J.GT.JS-7)then
!                 TAUL(KN,J) = TAUL(KN,J) * 1.d9
!          elseif(J.GT.JS-8)then
!                 TAUL(KN,J) = TAUL(KN,J) * 1.d8
!          elseif(J.GT.JS-9)then
!                 TAUL(KN,J) = TAUL(KN,J) * 1.d0
           else
                  TAUL(KN,J) = TAUL(KN,J) * 1.d0
           endif
         endif


! fudge CI for 1994I t=19d
         if ((na.eq.6) .and. nj.eq.1) then
               TAUL(KN,J) = TAUL(KN,J) * 0.5e2
         endif

        
         if ((na.eq.6) .and. nj.eq.1) then
           if(J.GT.js-1)then
                  TAUL(KN,J) = TAUL(KN,J) * 2.d03
           elseif(J.GT.JS-2)then
                  TAUL(KN,J) = TAUL(KN,J) * 2.d03
           elseif(J.GT.JS-3)then
                  TAUL(KN,J) = TAUL(KN,J) * 2.d03
           elseif(J.GT.JS-4)then
                  TAUL(KN,J) = TAUL(KN,J) * 2.d03
           elseif(J.GT.JS-5)then
                  TAUL(KN,J) = TAUL(KN,J) * 1.d03
           elseif(J.GT.JS-6)then
                  TAUL(KN,J) = TAUL(KN,J) * 1.d03
           elseif(J.GT.JS-7)then
                  TAUL(KN,J) = TAUL(KN,J) * 1.d03
!          elseif(J.GT.JS-8)then
!                 TAUL(KN,J) = TAUL(KN,J) * 1.d8
!          elseif(J.GT.JS-9)then
!                 TAUL(KN,J) = TAUL(KN,J) * 1.d0
           else
                  TAUL(KN,J) = TAUL(KN,J) * 3.d02
           endif
         endif
      endif         
!
!          if(J.GT.JS-3)THEN
!            TAUL(KN,J) = TAUL(KN,J) * 1.e4
!          elseif(J.GT.JS-5)then
!            TAUL(KN,J) = TAUL(KN,J) * 1.e4
!          elseif(J.GT.JS-10)then
!            TAUL(KN,J) = TAUL(KN,J) * 1.e3
!          endif
!        end if

!       fudge SiI 1994I-0415
!       if ((na.eq.14) .and. nj.eq.1) then
!         TAUL(KN,J) = TAUL(KN,J) * 1.e5
!       endif



!        if (na.eq.26 .and. nj.eq.3) then
!                 TAUL(KN,J) = TAUL(KN,J) * 1.e5
!        end if
!
!          WVA=10.E0**(AP(16)+8.E0)/XN2(KN)                     ! vacuum
!        back to air wvl's
!          IF (WVA. GT. 2000.) THEN
!            XLV = WAIR(WVA)
!          END IF
!        if ((xlv.gt.5875.e0 .and. xlv.lt.5876.e0)       ! fiddle HeI5876 flat
!    &           .and. (TAUL(KN,J).lt.0.1e0)) then
!             TAUL(KN,J) = 0.1e0
!          end if
!        end if

!       if (na.eq.26 .and. nj.eq.3) then
!         if (j.le.5) then                             !fudge FeIII
!           TAUL(KN,J)  = TAUL(KN,J) *  1000.
!         else ! if (j.gt.10) then
!           TAUL(KN,J)  = TAUL(KN,J) *  .010
!         end if
!       end if
!       if (na.eq.26 .and. nj.eq.2) then
!         if (j.eq.15) then                             !fudge FeII
!           TAUL(KN,J)  = TAUL(KN,J) *  10.
!         else if (j.gt.15) then
!           TAUL(KN,J)  = TAUL(KN,J) *  100.
!         end if
!       end if

        go to 2222
        if (na.eq.14 .and. nj.eq.2) then        !99ee 9 0ct lick
          TAUL(KN,J)  = TAUL(KN,J) * 5.
          if (j.eq.9) then                      !fudge Si II (99ee)
            TAUL(KN,J)  = TAUL(KN,J) * 1.6
          else if (j.eq.10) then
            TAUL(KN,J)  = TAUL(KN,J) * 6.
          else if (j.eq.11) then
            TAUL(KN,J)  = TAUL(KN,J) * 25.
          else if (j.eq.12) then
            TAUL(KN,J)  = TAUL(KN,J) * 35.
          else if (j.eq.13) then
            TAUL(KN,J)  = TAUL(KN,J) * 50.
          else if (j.ge.14) then
            TAUL(KN,J)  = TAUL(KN,J) * 50.
          end if
        end if
        if (na.eq.20 .and. nj.eq.2) then
          if (j.eq.11) then
            TAUL(KN,J)  = TAUL(KN,J) *  5.
          else if (j.gt.11 .and. j.le.13) then          !fudge Ca II (99ee)
            TAUL(KN,J)  = TAUL(KN,J) * 25.
          else if (j.gt.13) then
            TAUL(KN,J)  = TAUL(KN,J) * 400.
          end if
        end if
        if (na.eq.6 .and. nj.eq.2) then
          if (j.gt.12) then
            TAUL(KN,J)  = TAUL(KN,J) * 400.
          end if
        end if
        if (na.eq.26 .and. nj.eq.2) then
          if (j.eq.15) then                             !fudge FeII
            TAUL(KN,J)  = TAUL(KN,J) *  10.
          else if (j.gt.15) then
            TAUL(KN,J)  = TAUL(KN,J) *  100000.
          end if
        end if

        goto 2222
        if (na.eq.6 .and. nj.eq.2) then         !fudge C (90N)
          if (j.eq.1) then
            TAUL(KN,J)  = TAUL(KN,J) * 1.00
          else if (j.ge.2 .and. j.le.6) then
            TAUL(KN,J)  = TAUL(KN,J) * 1.00
          else if (j.eq.7) then
            TAUL(KN,J)  = TAUL(KN,J) * 1.0      ! 1.35
          else if (j.eq.8) then
            TAUL(KN,J)  = TAUL(KN,J) * 1.0      ! 1.75
          else if (j.eq.9) then
            TAUL(KN,J)  = TAUL(KN,J) * 1.       ! 2.5
          else if (j.eq.10) then
            TAUL(KN,J)  = TAUL(KN,J) * 1.0      ! 5
          else if (j.eq.11) then
            TAUL(KN,J)  = TAUL(KN,J) * 10.0     ! 10
          else if (j.eq.12) then
            TAUL(KN,J)  = TAUL(KN,J) * 1.0      ! 9
          else if (j.eq.13) then
            TAUL(KN,J)  = TAUL(KN,J) * 1.0      ! 9
          else if (j.eq.14) then
            TAUL(KN,J)  = TAUL(KN,J) * 1.0      ! 9
          else if (j.eq.15) then
            TAUL(KN,J)  = TAUL(KN,J) * 1.0      ! 9
          else if (j.eq.16) then
            TAUL(KN,J)  = TAUL(KN,J) * 1.0      ! 9
          else if (j.eq.17) then
            TAUL(KN,J)  = TAUL(KN,J) * 1.0      ! 9
          else if (j.ge.18) then
            TAUL(KN,J)  = TAUL(KN,J) * 1.0      ! 9
          end if
        end if
2222    continue
!       if (na.eq.56 .and. nj.eq.2) then                !fudge Ba
!         if (j.le.8) then
!           tau = 3.* tau
!         end if
!       end if


!CCCCCCCCCCCCCC  end fiddle  cccccccccccccccccccccccc

        EX(KN,J)=EXP(-TAUL(KN,J))                       ! S P

!if (na.eq.14 .and. nj.eq.2) then
!  if (xlv.gt.6347.e0 .and. xlv.lt.6371.e0) then
!    write(10,*) kn,j,XN2(KN),TAUL(KN,J),EX(KN,J)
!  end if
!end if

50      CONTINUE

70      CONTINUE

80      CONTINUE


9999    RETURN
        END





        FUNCTION BETA(TAU)

!       Escape probability - single precision.
!       Expansion for small tau necessary because beta tends to 0/0.

        IF(ABS(TAU) .LT. 0.1) THEN
          BETA=1-TAU/2+TAU**2/6-TAU**3/24+TAU**4/120
        ELSE
          BETA=(1-EXP(-TAU))/TAU
        END IF

9999    RETURN
        END


