        SUBROUTINE CONSTS
        SAVE
!
!       Single precision version of CONSTS
!
        INCLUDE 'ficc.f'        
!
!
!       MATHEMATICAL CONSTS
!
!       PI
        AP(9)=3.141592654
!       LOG PI
        AP(10)=LOG10(AP(9))
!       LOG E
        AP(11)=LOG10(EXP(1.0))
!       LN 10
        AP(12)=LOG(10.0)
!
!
!       PHYSICAL CONSTS
!
!       LOG AMU         ATOMIC MASS UNIT
        AP(13)=LOG10(1.660531E-24)
!       LOG H           PLANCK CONST
        AP(14)=LOG10(6.62620E-27)
!       LOG G           GRAV CONST
        AP(15)=LOG10(6.670E-8)
!       LOG C           VEL LIGHT
        AP(16)=LOG10(2.9979250E10)
!       LOG K           BOLTZMANN CONST
        AP(17)=LOG10(1.38062E-16)
!       LOG A           RAD CONST
        AP(18)=LOG10(7.56464E-15)
!       LOG MH          MASS H ATOM
        AP(19)=LOG10(1.67352E-24)
!       LOG MHE         MASS HE ATOM
        AP(20)=AP(19)+LOG10(4.0026/1.0080)
!       LOG 1EV         1EV IN ERGS
        AP(21)=LOG10(1.602192E-12)
!       LOG ME          MASS ELECTRON
        AP(22)=LOG10(9.10956E-28)
!       LOG E           EL CHARGE ESU
        AP(23)=LOG10(4.80325E-10)
!       LOG SG          STEFAN-BOLTZMANN
        AP(24)=LOG10(5.6696E-5) 
!
!       LOG SGE         THOMSON SCATTERING COEFF
        AP(40)=LOG10(0.66524E-24)
!
!
!       ASTRONOMICAL CONSTS
!
!       LOG MS          MASS SUN
        AP(25)=LOG10(1.989E33)
!       LOG LS          LUM SUN
        AP(26)=LOG10(3.826E33)
!       LOG RS          RAD SUN
        AP(27)=LOG10(6.9599E10)
!       LOG TS          TEFF SUN
        AP(28)=LOG10(5770.0)
!       MBOL SUN
        AP(29)=4.75
!       LOG SEC IN YR
        AP(30)=LOG10(3.15581495E7)
!       LOG CM IN PC
        AP(31)=LOG10(3.085678E18)
!
        RETURN
        END
