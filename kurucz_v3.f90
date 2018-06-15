program kurucz
 !
 implicit none
 integer npoints,ipoint,info,i
 integer :: ib,ie,nfiles
 integer,parameter :: nlevels_max=200000
 !
 real(8) temp,freql,freqr,gns,halfwidth,meanmass,partfunc,dfreq,dfreq_,acoef,abscoef,&
         tranfreq,cmcoef,beta,ln2,pi,dpwcoef,absthresh,emcoef,exp_coef,gf,loggf,intband
 integer(4) :: iso,iso_,nlevels,j,ilevel
 real(8) :: energy_l,energy_u,j_l,j_u, j0, pf, energy, lambda,x1,x2
 real(8) ::  de,x0,xp,xm,elow=100000.

 real(8) :: thresh = 1e-16
 real(8),allocatable :: freq(:),intens(:),energies(:)
 real(8),allocatable :: jrot(:),energy_(:)
 !
 character(60) intfilename(200),enrfilename
 character(30) proftype,specttype
 character(2) isochar
 character(5) label_u,label_l
 character(5),allocatable :: label_(:)
 logical :: new_u, new_l
 !
 real(8),parameter :: planck=6.6260693d-27,avogno=6.0221415d+23,vellgt=2.99792458d+10,boltz=1.380658d-16

! initialize variables 
!
halfwidth = 1e-2 ; meanmass = 1e-1 ; absthresh = 1e-8
!
!read number of intensity files
 read*,nfiles

if (nfiles>200) then 
  print('("Too many files (>200):",i8)'),nfiles
  stop 'Too many files'
endif

!read intensities filenames
do i = 1,nfiles
 read*,intfilename(i)
 !print('(1x,50a)'),intfilename(i)
enddo
!
read*,enrfilename
!
!read the isotopologue number 
 read*,iso
 !
 if (all(iso/=(/51,28,29,30,18,12,13/))) then 
   print('(" Illegal iso-number = ",i4,", 12,13,28,29,30, and 18 are only allowed")'), iso
   stop 'Illegal iso-number'
 endif
 !
!read nuclear spin statistics weights
 read*,gns

!read temperature, partition function
 read*,temp,partfunc

!read frequency range and number of points
 read*,freql,freqr,npoints

!read type of profiling
 read*,proftype

!read type of spectra
 read*,specttype

!read gaussian half-width or molecular mean mass
 select case (proftype(1:5))
    case ('gauss'); read*,halfwidth
    case ('doppl'); read*,meanmass
    case ('rect ','bin  '); 
       read*,halfwidth
       halfwidth = min((freqr-freql)/real(npoints,8),halfwidth)
    case ('stick','HITRA','hitra'); read*,absthresh
 end select

!read the threshold
 read*,thresh

!compute constants
 ln2=log(2.0d0)
 pi=acos(-1.0d0)
 beta=planck*vellgt/(boltz*temp)
 cmcoef=1.0d0/(8.0d0*pi*vellgt)
 emcoef=planck*vellgt/(4.0d0*pi)
 dpwcoef=sqrt(2.0*ln2*boltz*avogno)/vellgt
 dpwcoef = dpwcoef*sqrt(temp/meanmass)


!set up freq. grid
 allocate(freq(npoints),stat=info); if (info/=0) stop 'error: freq is out of memory'
 allocate(intens(npoints),stat=info); if (info/=0) stop 'error: intens is out of memory'
 dfreq=(freqr-freql)/real(npoints-1,8)
 forall(ipoint=1:npoints) freq(ipoint)=freql+real(ipoint-1,8)*dfreq


 !
 ! prepare computing the partition function
 if (partfunc<=0) then
    !
    if (trim(enrfilename)=='none'.or.trim(enrfilename)=='NONE') then 
       !
       nlevels = 0
       allocate(energy_(nlevels_max),label_(nlevels_max),stat=info)
       if (info/=0) stop 'error: energies_, label_ is out of memory'
       !
       partfunc = 0 
       !
       do i = 1,nfiles
         !
         open(unit=1,file=trim(intfilename(i)))
         do
            !
            if (nlevels+1>nlevels_max) then 
              write(6,"('nlevel>nlevels_max, ',2i9)") nlevels+1,nlevels_max
              !
              print('(f16.8,2x,a5)'),(energy_(ilevel),label_(ilevel),ilevel=1,nlevels-1)
              !
              stop 'nlevel>nlevels_max'
            endif 
            !
            read(1,"(F10.4,F7.3,F5.1,F10.3,F5.1,F11.3,4x,a5,3x,a5,3x,i2)",end=30) tranfreq,loggf,J_l,energy_l,J_u,energy_u,label_l,label_u,iso_
            !
            if (iso/=iso_) cycle
            !
            gf = exp(abs(loggf))
            energy_u = abs(energy_u)
            energy_l = abs(energy_l)
            !
            new_u = .true.
            new_l = .true.
            !
            cycle_levels: do ilevel = 1,nlevels
              !
              if ( new_u .and. abs( energy_u-energy_(ilevel) )<=1e-3 .and.trim(label_(ilevel))==trim(label_u)) then 
                !
                new_u = .false.
                !
                if ( trim(label_(ilevel) )==trim(label_u) .and. abs( energy_u-energy_(ilevel) )>1e-3 ) then 
                  write(6,"('upper state labels are the same  but energies are different:',2a5,2f18.4)") trim(label_(ilevel)),trim(label_u),energy_(ilevel),energy_u
                  stop 'input state labels disagree with energies'
                endif
                !
                if ( trim(label_(ilevel) )/=trim(label_u) .and. abs( energy_u-energy_(ilevel) )<=1e-3 ) then 
                  write(6,"('upper state labels are different but energies are the same:',2a5,2f18.4)") trim(label_(ilevel)),trim(label_u),energy_(ilevel),energy_u
                  stop 'input state labels disagree with energies'
                endif
                !
              endif
              !
              if ( new_l .and. abs( energy_l-energy_(ilevel) )<=1e-3 .and.trim(label_(ilevel))==trim(label_l)) then 
                !
                new_l = .false.
                !
                if ( trim(label_(ilevel))==trim(label_l) .and. abs( energy_l-energy_(ilevel) )>1e-3 ) then 
                  write(6,"('lower state labels are the same  but energies are different:',2a5,2f18.4)") trim(label_(ilevel)),trim(label_l),energy_(ilevel),energy_l
                  stop 'input state labels disagree with energies'
                endif
                !
                if ( trim(label_(ilevel))/=trim(label_l) .and. abs( energy_l-energy_(ilevel) )<=1e-3 ) then 
                  write(6,"('lower state labels are different but energies are the same:',2a5,2f18.4)") trim(label_(ilevel)),trim(label_l),energy_(ilevel),energy_l
                  stop 'input state labels disagree with energies'
                endif
                !
              endif
              !
              if (.not.new_l.and..not.new_u) exit cycle_levels
              !
            enddo cycle_levels
            !
            if (new_u) then
                nlevels = nlevels + 1
                energy_(nlevels) = energy_u
                label_(nlevels) = label_u
                !
                partfunc=partfunc+(2.0d0*J_u+1.0d0)*gns*exp(-beta*energy_u)
                !
                write(101,'(f16.8,2x,f9.1,2x,a5)') energy_u,j_l,label_u
                !
            endif
            !
            if (new_l) then
                nlevels = nlevels + 1
                energy_(nlevels) = energy_l
                label_(nlevels) = label_l
                !
                partfunc=partfunc+(2.0d0*J_l+1.0d0)*gns*exp(-beta*energy_l)
                !
                write(101,'(f16.8,2x,f9.1,2x,a5)') energy_l,J_l,label_l
                !
            endif
            !   
            cycle
         30  exit
         enddo
         !
         close(1)
         !
       enddo
       !
       deallocate(energy_,label_)
       !
       print('("Nlevels=",i8,"pf = ",f12.4)'),nlevels,partfunc
       !
    else
      !
      partfunc=0.0 ; j0=0
      pf = 0
      !count number of lines (levels) in the Energy files
      !
      write(isochar,"(i2)") iso
      !
      !enrfilename = "sioxx-"//trim(isochar)//".energies"
      !
      open(unit=1,file=trim(enrfilename))
      i = 0
      do
         !
         read(1,*,end=10) temp
         !
         i = i + 1
         !
         cycle
         !
      10  exit
      end do
      !
      rewind(1)
      !
      nlevels = i
      !
      allocate(energies(nlevels),Jrot(nlevels),stat=info)
      if (info/=0) stop 'error: sym,energies,Jrot is out of memory'
      !
      !
      ! start reading the energies 
      !
      do i=1,nlevels
         !
         read(1,*) energies(i),jrot(i)
         !
         j = jrot(i)
         energy = energies(i)
         !
         ! partition function
         !
         partfunc=partfunc+real((2*j+1),8)*gns*exp(-beta*energy)
         !
         j0 = j
         !
      end do
      close(1)
      !
      print('(1x,a,1x,es16.8)'),'# partition function value is',partfunc
      !
      deallocate(energies,Jrot)
      !
   endif 
   !
 endif
 !
 intens=0.0
 intband = 0.0
 !
 !start loop over all transitions
 !
 do i = 1,nfiles
   !
   open(unit=1,file=trim(intfilename(i)))
   do
      !   read new line from intensities file
      !
      !read(1,"(i5,i5,f6.1,f6.1,f6.1,2x,1a,i3,f10.2,f14.4,f10.2)",end=20) v_u,v_l,omega_u,omega_l,J_l,parity,ibranch,tranfreq,acoef,energy_l
      !
      !    0    4   3.5   3.5   4.5  +  1   3294.80        0.0340   6720.86
      !
      read(1,"(F10.4,F7.3,F5.1,F10.3,F5.1,F11.3,4x,a5,3x,a5,3x,i2)",end=20) lambda,loggf,J_l,energy_l,J_u,energy_u,label_l,label_u,iso_
      !
      if (iso/=iso_) cycle
      !
      X1=-1.510
      X2=-.001
      GF=EXP((LOGGF+X1+X2)*2.30258509299405E0)
      !
      !gf = exp((loggf))
      energy_u = abs(energy_u)
      energy_l = abs(energy_l)
      !
      !gf = (1.4992 x 10^4)(A/s-1)(lambda/nm)^2
      !
      ! A = gf/(1.4992 x 10^-14)/ (lambda/nm)^2
      !
      acoef = gf/(1.4992e-14)/lambda**2
      !
      ! correct for the statistical weight 
      !
      !acoef = acoef/(2.0d0*j_u+1.0d0)/gns
      !
      ! change to cm-1
      !
      tranfreq = energy_u-energy_l
      !
      !
!     ++++++++++^^^^^^^+++++^^^^^^^^^^+++++^^^^^^^^^^^++++^++^+^^^+^^+^+++^^
! 433.0318 -3.524 19.5-10563.271 20.5 -33649.772 106X02F2   A02F1   13
! wl(nm)   log gf  J    E(cm-1)   J'   E'(cm-1) code  V      V'     iso
!                                                   label   label'

      !
      if (tranfreq<freql.or.tranfreq>freqr) cycle 
      !
      !   half width for Doppler profiling
      if (proftype(1:5)=='doppl') halfwidth=dpwcoef*tranfreq
      !
      x0 = sqrt(ln2)/halfwidth*dfreq*0.5d0
      !
      !   if transition frequency is out of selected range
      if (tranfreq>freqr+10.0*halfwidth.or.tranfreq<freql-10.0*halfwidth) cycle
      !
      if (tranfreq<epsilon(1.0d0)) cycle
      !
      select case (trim(specttype))
        !
      case default
        !
        print('(a,2x,a)'),'Illegal key:',trim(specttype)
        stop 'Illegal specttype-key'
        !
      case ('absorption','ABSORPTION','absorpt','ABSORPT')
        !
        ! compute absorption coefficient [cm/molecule]
        !
        abscoef=cmcoef*acoef*gns*(2.0d0*j_u+1.0d0)*exp(-beta*energy_l)*(1.0d0-exp(-beta*tranfreq))/(tranfreq**2*partfunc)
        !
      case ('emission','EMISSION','EMISS','emiss')
        !
        ! emission coefficient [Ergs/mol/Sr]
        !
        abscoef=emcoef*acoef*gns*(2.0d0*j_u+1.0d0)*exp(-beta*energy_u)*tranfreq/(partfunc)
        !
      end select
      !
      intband = intband + abscoef
      !
      ! if only stick spectrum needed
      if (proftype(1:5)=='stick') then
        !
        if (abscoef>absthresh) then 
          print('((1x,2es16.8,2x,a5,2x,a5,f6.1,2x,f6.1,f14.4,f12.4))'), &
                   tranfreq,abscoef,label_u,label_l,j_u,j_l,acoef,energy_l
        endif
        cycle
      end if
      !
      !   if only stick spectrum needed
      if (proftype(1:6)=='hitran'.or.proftype(1:6)=='HITRAN') then
         if (abscoef>absthresh) & 
             print('(I2,I1,F12.6,E10.3,E10.3,F5.4,F5.4,F10.4,F4.2,F8.6,a7,8x,a7,8x,F5.1,8x,F5.1,8x,6I1,6I2,1x,F7.1,F7.1)'),&
             iso_,1,tranfreq,abscoef,acoef,0.0,0.0,energy_l,0.0,0.0,label_u,label_l,j_u,j_l,0,0,0,0,0,0,0,0,0,0,0,0,gns*(2.0d0*j_u+1.0d0),gns*(2.0d0*j_l+1.0d0)
         cycle
      end if
      !
      !   normalize absorption coefficient
      !
      !abscoef=abscoef*sqrt(ln2/pi)/halfwidth
      !
      if (abscoef<thresh) cycle 
      !
      !exp_coef = -ln2/halfwidth**2
      !
      ib =  max(nint( ( tranfreq-halfwidth*10.0-freql)/dfreq )+1,1)
      ie =  min(nint( ( tranfreq+halfwidth*10.0-freql)/dfreq )+1,npoints)
      !   normalize absorption coefficient
      !
      !read gaussian half-width or molecular mean mass
      select case (trim(proftype(1:5)))
          !
      case ('gauss','doppl')
          !
          !abscoef=abscoef*sqrt(ln2/pi)/halfwidth
          !
          !$omp parallel do private(ipoint,dfreq_,xp,xm,de) shared(intens) schedule(dynamic)
          do ipoint=ib,ie
             !
             dfreq_=freq(ipoint)-tranfreq
             !if (abs(dfreq_)>halfwidth*10.0) cycle
             !
             xp = sqrt(ln2)/halfwidth*(dfreq_)+x0
             xm = sqrt(ln2)/halfwidth*(dfreq_)-x0
             !
             de = erf(xp)-erf(xm)
             !
             intens(ipoint)=intens(ipoint)+abscoef*0.5d0/dfreq*de
             !
          enddo
          !$omp end parallel do 
          !
          !print('(2(1x,es16.8),i4,i4,f11.4,f11.4,f11.4,es16.8,i6,i6)'),freq(1),intens(1),ji,jf,energyi,energyf,tranfreq,linestr,ileveli,ilevelf
          !
          !ji_ = ji ; jf_ = jf
          !
      case ('rect');
        !
        !$omp parallel do private(ipoint,dfreq_) shared(intens) schedule(dynamic)
        do ipoint=ib,ie
           dfreq_=freq(ipoint)-tranfreq
           !if (abs(dfreq_)>halfwidth*10.0) cycle
           intens(ipoint)=max(intens(ipoint),abscoef)
        enddo
        !$omp end parallel do
        !
        !print('(2(1x,es16.8),i4,i4,f11.4,f11.4,f11.4,es16.8,i6,i6)'),freq(1),intens(1),ji,jf,energyi,energyf,tranfreq,linestr,ileveli,ilevelf
        !
      case ('bin');
        !
        !$omp parallel do private(ipoint,dfreq_) shared(intens) schedule(dynamic)
        do ipoint=ib,ie
           dfreq_=freq(ipoint)-tranfreq
           !if (abs(dfreq_)>halfwidth*10.0) cycle
           !
           intens(ipoint)=intens(ipoint)+abscoef/dfreq
           !
        enddo
        !$omp end parallel do
        !
        !print('(2(1x,es16.8),i4,i4,f11.4,f11.4,f11.4,es16.8,i6,i6)'),freq(1),intens(1),ji,jf,energyi,energyf,tranfreq,linestr,ileveli,ilevelf
          !
      end select





      !
      cycle
   20  exit
   enddo
   !
 enddo 

 close(1)

 
 if (proftype(1:5)=='doppl'.or.proftype(1:5)=='gauss'.or.proftype(1:4)=='rect'.or.proftype(1:3)=='bin') then 
   !
   !print out the generated intensities convoluted with a given profile
   !
   print('(2(1x,es16.8))'),(freq(ipoint),intens(ipoint),ipoint=1,npoints)
   !
   print('("Total intensity  (sum):",es16.8," (int):",es16.8)'), intband,sum(intens)*dfreq
   !
 else
   !
   print('("Total intensity = ",es16.8)'), intband
   !
 endif 


end program kurucz
