        SUBROUTINE OPENING
        SAVE

        ! SHChange: INCLUDE ficc
        INCLUDE 'ficc.f'

!
!
!      input files
       CHARACTER*256 home,datapath,workdir,wdlast


!       call getenv('HOME',home)
!       datapath = home(1:len_trim(home))//'/Codes/MonteCarlo/InputData/'
        call getenv('MCDATA',datapath)
        call getenv('PWD',workdir)

!        if (datapath .eq. '') stop 'Set MCDATA environment variable!'
!        datapath='/root/modelling/code/InputData'
        datapath = datapath(1:len_trim(datapath))


        write(*,*)' opening'
        OPEN (UNIT=2,FILE='dica.dat',STATUS='OLD')              ! i/p data
        OPEN (UNIT=11,FILE=datapath(1:len_trim(datapath))//'/ihea.ind',STATUS='OLD') ! ion. pot.etc
        OPEN (UNIT=12,FILE='comp.ind',STATUS='OLD')             ! Ab/mass frac.
!       OPEN (UNIT=13,FILE='/home/mazzali/lines/leon2',
!    &                                          STATUS='OLD')   ! Line list
!       OPEN (UNIT=15,FILE='/home/mazzali/lines/leon1.new',
!    &                                          STATUS='OLD')   ! Atomic lvls
        OPEN (UNIT=13,FILE=datapath(1:len_trim(datapath))//'/lines/leon2',
     &                                          STATUS='OLD')   ! Line list
        OPEN (UNIT=15,FILE=datapath(1:len_trim(datapath))//'/lines/leon1.new',
     &                                          STATUS='OLD')   ! Atomic lvls
        OPEN (UNIT=16,FILE=datapath(1:len_trim(datapath))//'/nbia.ind')	        ! Atomic wgts
        OPEN (UNIT=17,FILE=datapath(1:len_trim(datapath))//'/nama.ind')	        ! St wt gd confs
        OPEN (UNIT=18,FILE=datapath(1:len_trim(datapath))//'/nbba.ind')	        ! Atomic symbols

        OPEN (UNIT=19,FILE='hydro.dat',STATUS='OLD',IOSTAT=IO)
        IF(IO.ne.0)THEN
          WRITE(*,*)'File hydro.dat not found in current working directory.'
          WRITE(*,*)'Will use default W7 model.'
          OPEN(UNIT=19,FILE=datapath(1:len_trim(datapath))//'/rhw7.dat')
        ENDIF
        

!       OPEN (UNIT=19,FILE='rhw7_31.dat',STATUS='OLD')  ! W7 @ day 31
!       OPEN (UNIT=19,FILE='co110modfor03bg.dat',STATUS='OLD')
!       OPEN (UNIT=19,FILE='co21.dat',STATUS='OLD')
!       OPEN (UNIT=19,
!     &    FILE='/users/mazzali/explo/w7.HYDRO',STATUS='OLD')   ! Wdd1 @ 1020s
!        OPEN (UNIT=19,FILE='co110n-5.dat',STATUS='OLD')        ! CO6138, 10d
        OPEN (UNIT=41,FILE=datapath(1:len_trim(datapath))//'/filters.dat',STATUS='OLD')          ! Filter fns.
!       OPEN (UNIT=70,FILE='tbib.dat',STATUS='OLD')             ! old T(r)
!
!
!      output files
!
        OPEN (UNIT=3,FILE='stst.dat',STATUS='NEW')              ! (pbib) checks
        OPEN (UNIT=10,FILE='diagn.dat',STATUS='NEW')            ! casc. diagn.
        OPEN (UNIT=20,FILE='yhea.dat',STATUS='NEW')             ! Ioniz. strat.
        OPEN (UNIT=21,FILE='lhea.dat',STATUS='NEW')             ! Wrk line list
        OPEN (UNIT=29,FILE='spcp.dat',STATUS='NEW')             ! spec. (MC,f)
        OPEN (UNIT=30,FILE='sica.dat',STATUS='NEW')             ! MC spectrum
        OPEN (UNIT=31,FILE='eica.dat',STATUS='NEW')             ! Emitted lines
        OPEN (UNIT=32,FILE='spct.dat',STATUS='NEW')             ! sp for plot
        OPEN (UNIT=37,FILE='atmd.oud',STATUS='NEW')             ! at. models
        OPEN (UNIT=42,FILE='ptfn.dat',STATUS='NEW')             ! partn fns
        OPEN (UNIT=43,FILE='taul.dat',STATUS='NEW')             ! Sob. opt dep
        OPEN (UNIT=46,FILE='sbib.dat',STATUS='NEW')             ! lst sp lines
!       SHChange: Open additional data files 
!        IF (WHIST==1) OPEN (UNIT=50,FILE='/scratch/shaching'//workdir(1:len_trim(workdir))//'/AdditionalData/hist.dat',STATUS='REPLACE')
        
        if (workdir(len_trim(workdir):len_trim(workdir))=='/') then 
          wdlast = workdir(1:(len_trim(workdir)-1))
          WRITE(*,'(A)') 'Interesting: snipped workdir'
        ELSE
          WDLAST = workdir(1:(len_trim(workdir)))
        END IF
        wdlast = wdlast((scan(wdlast,'/',.true.)+1):len(wdlast))
        IF (WHIST==1) OPEN (UNIT=50,FILE='/scratch/shaching/hist_'//wdlast(1:len_trim(wdlast))//'.dat',STATUS='REPLACE')

!        OPEN (UNIT=51,FILE='./AdditionalData/encint.dat',STATUS='REPLACE')
!        OPEN (UNIT=54,FILE='./AdditionalData/passint_wlsel.dat',STATUS='REPLACE')
!        OPEN (UNIT=55,FILE='./AdditionalData/passstayint_wlsel.dat',STATUS='REPLACE')
!        OPEN (UNIT=56,FILE='./AdditionalData/encint_wlsel.dat',STATUS='REPLACE')
!        OPEN (UNIT=57,FILE='./AdditionalData/em_wlsel.dat',STATUS='REPLACE')
!        OPEN (UNIT=58,FILE='./AdditionalData/passreabs_wlsel.dat',STATUS='REPLACE')
!        OPEN (UNIT=61,FILE='./AdditionalData/emprocint_wlsel.dat',STATUS='REPLACE')

!        OPEN (UNIT=52,FILE='./AdditionalData/imm_wlsel.dat',STATUS='REPLACE')
!        OPEN (UNIT=53,FILE='./AdditionalData/imm_bl_wlsel.dat',STATUS='REPLACE')
!        OPEN (UNIT=62,FILE='./AdditionalData/encext_wlsel.dat',STATUS='REPLACE')
!        OPEN (UNIT=64,FILE='./AdditionalData/passext_wlsel.dat',STATUS='REPLACE')
!        OPEN (UNIT=66,FILE='./AdditionalData/passstayext_wlsel.dat',STATUS='REPLACE')
!        OPEN (UNIT=67,FILE='./AdditionalData/passstayext_bl_wlsel.dat',STATUS='REPLACE')
!        OPEN (UNIT=68,FILE='./AdditionalData/emprocext_wlsel.dat',STATUS='REPLACE')
!        OPEN (UNIT=69,FILE='./AdditionalData/emprocext_bl_wlsel.dat',STATUS='REPLACE')

!
        RETURN
        END

