!       Monte Carlo code with photon cascades -  conditional emissivities 
!       - imposed level populations - ie NLTE iterations excluded.

!       SHChange: Declare Mersenne Twister
        use mtprng
        SAVE

        INCLUDE 'ficc.f'

        write(*,*)' included'  


        CALL OPENING                                ! FICOP

        CALL PARAMS

        CALL OPTIONS

!       SHChange: Init Mersenne Twister
        CALL mtprng_init(KR,MTSTAT)
        IF (USEMT==1) WRITE(*,*) 'Using MERSENNE TWISTER PRNG'

        CALL INPUT                                  ! FICJ

        CALL ATLVLS                                 ! FICK 

        CALL MODEL                                  ! FICS
                
!       SHChange: NRTEMP in ficc.f -> variable #temp iter 
!                 TEMP ITER regulated by max deviation criterion
        L=1
        DO WHILE (L<=NRTEMP)
          write(*,*)'  TEMP ITER  ',L

          CALL PARTF(L)                               ! FICP

          CALL IONEQ(L)                               ! FICF

          IF (L.EQ.1) THEN
            CALL LINES                                ! FICL
          END IF

          CALL FRQSMP

          CALL INITIATE                               ! FICI
          
          CALL TAUS                                   ! FICT

          CALL BRATS                                  ! FICQ

          CALL RADX                                   ! FICR

!         IF(KC(3).EQ.1) CALL MONTEC(0)               ! FICM
!         IF(KC(3).EQ.0) CALL MONTEC1(0)              ! FICG

          CALL MONTEC(0)                          ! FIYM

          CALL TEMPIT(L)                              ! FICIT
          L=L+1 
        END DO


        CALL PARTF(L)                                 ! FICP

        CALL IONEQ(L)                                 ! FICF

!c      CALL FRQSMP

!c      CALL INITIATE                                 ! FICI

        CALL TAUS                                     ! FICT

        CALL BRATS                                    ! FICQ


        IF(KC(3).EQ.1) THEN 

          DO 50 ITER=1,ITT

            write(*,*)'  Rad. field ITER  ',iter

            CALL RADX                                    ! FICR
            if (iter.lt.itt) then 
              nn = 2 
            else if (iter.eq.itt) then  !!! use 3 only if need good MC spectrum 
              nn = 4
            end if       
            CALL MONTEC(nn)                              ! FIYM

50        CONTINUE

        ELSE IF(KC(3).EQ.0) THEN 
          
          CALL MONTEC1(2)                               ! FICG
            
        END IF  

        CALL SPECTRUM                               ! FICZ 

        CALL CLOSING                                ! FICCL

        write (*,*) ' Finish Code'

        STOP
        END










        SUBROUTINE PARAMS
        SAVE

        INCLUDE 'ficc.f'


        CALL CONSTS(AP)

        write(*,*) 'start params'
        READ(2,*)
        READ(2,*) NP4,ITT
        READ(2,*)
        READ(2,*) JS,MB,KB,XE1
        IF(XE1.EQ.0.)XE1=0.200D0

!1006    FORMAT(/I7,I4//I5,I8,I6,F7.3)
        WRITE(3,2006) NP4,ITT,JS,KB,XE1
2006    FORMAT(/5X,'NP4 =',I6,5X,'ITT =',I3,5X,'JS =',I3,5X,' KB =',I4,' XE1=',F6.3)  

        NP=10000*NP4

        READ(2,*)
        READ(2,*) ZLM,VPH,TB
!1005    FORMAT(/F8.3,2F10.1)

        READ(2,*)
        READ(2,*) TDY,XMU,REDSHIFT
!1007    FORMAT(/2F8.2,F10.4)

        READ(2,*)
        READ(2,*) EBMVGAL,EBMVCMF
!1016    FORMAT(/2F8.2)

        READ(2,*)
        READ(2,*) KR,ZTAUG
!1018    FORMAT(/I12,F8.2)

        WRITE(3,2007) TDY,ZTAUG
2007    FORMAT(/5X,'Day =',F6.2,5X,'Log Taug =',F7.3)

        READ(2,*)
        READ(2,*) CHL,NC
        WRITE(3,2017) CHL,NC

        READ(2,*)
        READ(2,*) WVM1,WVM2,KG
        WRITE(3,2019) WVM1,WVM2,KG

        READ(2,*)
        READ(2,*) WV1,WV2,ZLMO                  ! not used if KC(1)=1

        READ(2,*)
        READ(2,*) WVA,WVB                       ! not used if KC(3)=1 
        
        XNUA=10**(AP(16)+8)/WVA
        XNUB=10**(AP(16)+8)/WVB                     ! used in FICG 


        READ(2,*)
        READ(2,*) (KC(L),L=1,5)
        WRITE(*,1023) (KC(L),L=1,5)

        IF (KC(1).EQ.0) WRITE(3,2021) WV1,WV2,ZLMO
        WVV1=WVAC(WV1)                              ! in vac for MONTEC
        WVV2=WVAC(WV2)


        TM=TDY*24*3600                                    ! sec
        ZRDCM=DLOG10(VPH*1.D5*TM)                         ! Log R(cm)
        ZRD=ZRDCM-AP(27)        

        ZTEFF=AP(28)+0.25D0*ZLM-0.5D0*ZRD

        ELL=CHL*8065.46D0                                 ! Exc limit (cm-1)

        ZVOL=DLOG10(4.D0/3*AP(9))+3*ZRDCM                 ! Log Vol(cm3)

        WRITE(3,2008)
        WRITE(46,2008)
2008    FORMAT(/18X,'LOG R/RS',5X,'LOG L/LS',9X,'VEL LB',10X,'TB'/)

        WRITE(3,1008)ZRD,ZLM,VPH,TB
        WRITE(46,1008)ZRD,ZLM,VPH,TB
1008    FORMAT(F24.3,F12.3,F18.1,F14.1/)


1017    FORMAT(/F6.1,I8)
1019    FORMAT(/2F9.4,I7)
1021    FORMAT(/2F10.2,F8.3)
1023    FORMAT(/5I4)
1025    FORMAT(/2F10.2)

2017    FORMAT(/5X,'CHL(eV) =',F5.1,5X,'NC =',I3)
2019    FORMAT(/5X,'WVM1 =',F7.4,4X,'WVM2 =',F7.4,7X,'KG =',I6)
2021    FORMAT(/5X,'Obs sp ','   wv1,wv2',F11.2,F9.2,5X,'Log L(wv1,wv2) =' ,F6.3) 

9999    RETURN
        END








        SUBROUTINE OPTIONS
        SAVE

!       Print-out for the options specified in DICA.IND.

        INCLUDE 'ficc.f'

        IF(KC(1).EQ.1) WRITE(*,2001) 
        IF(KC(1).EQ.0) WRITE(*,2003)

        IF(KC(2).EQ.1) WRITE(*,2005)
        IF(KC(2).EQ.0) WRITE(*,2007)

        IF(KC(3).EQ.1) WRITE(*,2009)
        IF(KC(3).EQ.0) WRITE(*,2011)WVA,WVB



2001    FORMAT(/5X,'Scale MC sp',' to L(bol)')
2003    FORMAT(/5X,'Scale MC sp',' to L(wv1,wv2)')

2005    FORMAT(/5X,'Cascades')
2007    FORMAT(/5X,'Resonance ','scattering')

2009    FORMAT(/5X,'Black body',' em at ','lwr bdry')
2011    FORMAT(/5X,'Top hat em',' at lwr bdry',F12.2,F10.2)

9999    RETURN
        END







        SUBROUTINE FRQSMP
        SAVE

!       Bins for sampling blackbody emission.

        INCLUDE 'ficc.f'

        DIMENSION FN(5000),ZL(5000),GN(5000)

        AAA=3.5D0
        L1=4000

        L2=L1-1
        K2=KB-1
        DZ=1.D0/(L1-1)
        DO 15 L=1,L1
15      ZL(L)=(L-1)*DZ
        DO 20 L=2,L2
        XX=AAA*SQRT(ZL(L)/(1-ZL(L)))
        GN(L)=0
        IF(XX.GT.30) GO TO 20
        GN(L)=XX**2*(AAA**2+XX**2)**2/(EXP(XX)-1)
20      CONTINUE
        GN(1)=0
        GN(L1)=0

        FN(1)=0
        SM=0
        DO 22 L=2,L1
        SM=SM+0.5D0*DZ*(GN(L)+GN(L-1))
22      FN(L)=FN(L-1)+0.5D0*DZ*(GN(L)+GN(L-1))

        ZK(1)=0
        ZK(KB)=1
        K=2
        DO 50 L=2,L1
25      Z=(K-1.D0)/(KB-1.D0)
        IF(FN(L)-Z*FN(L1))50,30,30
30      ZK(K)=ZL(L-1)+DZ*(Z*FN(L1)-FN(L-1))/(FN(L)-FN(L-1))
        K=K+1
        IF(K-KB)25,9999,9999
50      CONTINUE

9999    RETURN
        END

