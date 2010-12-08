        SUBROUTINE NEWFRQ(KN,XPE,EP,J)
!       SHChange: USE MTPRNG,SAVE
        use mtprng
        SAVE
!       KN,XPE of absorbed pkt replaced by values for emitted pkt

        INCLUDE 'ficc.f'

        DIMENSION XNC(1000),LC(1000)


        IF(KC(2).EQ.0 .OR. IDL(KN).EQ.0) THEN                
        XPE=XN2(KN)                                    ! Resonance scattering
        GO TO 9999
        END IF



        M=IN(LN+1-KN)                                  ! Index in basic list
        MI=M 

        NA=ELM(M)+1.D-3
        NJ=100*(ELM(M)-NA)+1+1.D-3


        M=MI
!       SHChange: JIABS -> ABS
        LU=ABS(NUP(M))                                ! Initial uppr lvl


!       Select particular downward transition from lvl LU

        SM=0
!       SHChange: optionally use Mersenne Twister 
	
        IF (USEMT==1) THEN
	  PP=MTRAN(MTSTAT)
	ELSE
          PP=RAN(KR)
        END IF 
	 
        DO 40 M=MS(LU),MF(LU)
!write(*,*)j,m
          SM=SM+QQ(M,J)
          IF(SM .GT. PP) GO TO 50
40      CONTINUE
        M=MF(LU)

50      CONTINUE
        
        KN=IR(M)                                  ! New KN

        XPE=XN2(KN)                               ! New XPE  

        NA1=ELM(M)+1.D-3
        NJ1=100*(ELM(M)-NA1)+1+1.D-3
        IF(NA1.NE.NA .OR. NJ1.NE.NJ) THEN     ! Emitting = absorbing ion ?  
          write (*,*)' NA1=',na1,'  NJ1=',nj1,' NA=',na,'  NJ=',nj
          write (*,*)' MI=',MI,'  M=',M, 'KN = ',KN, 'XPE= ',XPE, 'J= ', J
          WRITE(*,2003)
          CALL EXIT
        END IF


2003    FORMAT(/4X,'ERROR 3 ','NEWFRQ')

9999    RETURN
        END


