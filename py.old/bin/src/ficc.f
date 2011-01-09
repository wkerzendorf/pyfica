
!       NB When this file changes all subroutines must be recompiled.

!       Parameters are upper limits to no. of elements (NEZ) , of ions (NIZ) ,
!       of ions per atom (IZ) , of energy levels (NZ) , of energy levels per 
!       ion (MZ) ,of lines (LZ) , of shells (JZ) , of no. of grid points for
!       emergent spectrum (KZ).
!
!       SHChange: Additional parameter NRTEMP: No. of temp iterations 
!       SHChange: LZ= 80000 -> 200000, JZ= 28 -> 100 (JZ, KZ need code adaption)
!       SHChange: Optionally use Mersenne Twister

        PARAMETER (USEMT=0)
        PARAMETER (WHIST=0)
        PARAMETER (WBRAMA=0)
        PARAMETER (NEZ=60,NIZ=200,IZ=5,NZ=25000,MZ=850,LZ=400000)
        PARAMETER (JZ=85,KZ=5000)
        PARAMETER (NRTEMP=30)

        COMMON /BLKNEL/NEL
        
        COMMON /BLKA1/AP(50)    

        COMMON /BLKA2/ZK(501),AAA       

        COMMON /BLKA3/KC(5)       

        COMMON /BLKA4/ZRD,ZLM,ZTEFF,TB,TDY,XMU,REDSHIFT,
     &                EBMVGAL,EBMVCMF,TM,VPH,ZRDCM,ZVOL
        
        COMMON /BLKA5/ZTAUG,CHL,ELL,NC

        COMMON /BLKA6/NP,ITT,ITER,JS,MB,KB,KR,I1,XE1

        COMMON /BLKA7/WVM1,WVM2,KG

        COMMON /BLKA8/WV1,WV2,WVV1,WVV2,ZLMO

        COMMON /BLKA9/WVA,WVB,XNUA,XNUB


        COMMON /BLKF1/FR(IZ,NEZ,JZ),XNA(NEZ,JZ),XNT(JZ),XNE(JZ)    


        COMMON /BLKG1/XINB


        COMMON /BLKI1/DZWV,ZWVM(KZ+1),DWVM(KZ),ZWVMA(KZ+1),DWVMA(KZ)

        COMMON /BLKI2/GMX(NZ,JZ),GMD(NZ,JZ)



        COMMON /BLKJ1/CHE(60),CHI(9)                               ! A3,A5

        COMMON /BLKJ2/CH(9,60),ZZ(8,60)    

        COMMON /BLKJ3/XA(60),AW(60)     



        COMMON /BLKK1/IA(NEZ,9),NLV(NEZ,9),KS(NEZ,IZ)

        COMMON /BLKK2/EL(NZ),XJ(NZ),KKT      

        COMMON /BLKK3/I2,ION(NEZ,IZ),NAI(NIZ),NJI(NIZ)       



        COMMON /BLKL1/XN2(LZ),IDL(LZ)                              ! ordered

        COMMON /BLKL2/IN(LZ),IR(LZ)    

        COMMON /BLKL3/XN(LZ),ZGF(LZ),GL(LZ),GU(LZ),LN   

        COMMON /BLKL4/ELM(LZ),NUP(LZ),NLW(LZ),ELW2(LZ)

        COMMON /BLKL5/IS(NIZ),MS(NZ),MF(NZ)

        COMMON /BLKL6/AUL(LZ),BUL(LZ),BLU(LZ)


        COMMON /BLKM1/XJI(LZ,JZ),EA(LZ,JZ),EE(LZ,JZ)

        COMMON /BLKM2/PW(5),PWL(LZ),ZEP0
!       SHChange: replace 51 by JZ+1
        COMMON /BLKM3/RQ(JZ+1,10),XLO(KZ),NO(KZ),ZLDM(KZ),IC(10)



        COMMON /BLKO1/XNL(MZ),GG(MZ)


        COMMON /BLKP1/UU(JZ,NEZ,IZ)       


        COMMON /BLKQ1/QQ(LZ,JZ)


        COMMON /BLKR1/FRX(NZ,JZ)


        COMMON /BLKS1/XI(2*JZ+1),DX,VT(JZ),VS(JZ),ZRH(JZ),SG(JZ)       

        COMMON /BLKS2/TS(JZ),WD(JZ),TR(JZ)   


        COMMON /BLKSP/ZWVMIC(KZ),ZLDFM(KZ),ZLDMC(KZ)
        

        COMMON /BLKT1/TAUL(LZ,JZ),EX(LZ,JZ),BG(LZ)                   ! ordered



        COMMON /BLKZ1/EALV(NZ),XJE(LZ,JZ),STM(LZ,JZ)
