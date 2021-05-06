! This program dumps an entire model at the specified range and spacing
program dumpmodel
  use types
  use util
  use constants, only : R_E, PI, VERSION
  use ngo_dens_model_adapter, only : fngo=>funcPlasmaParams, ngoStateData, &
       ngoStateDataP, ngosetup=>setup
  use ngo_3d_dens_model_adapter, only : fngo3d=>funcPlasmaParams, ngo3dStateData, &
       ngo3dStateDataP, ngo3dsetup=>setup
  use gcpm_dens_model_adapter, only : fgcpm=>funcPlasmaParams, &
       gcpmStateData, gcpmStateDataP
  use simple_3d_model_adapter, only : fsimple=>funcPlasmaParams, simpleStateData, &
       simpleStateDataP
  use interp_dens_model_adapter, only : finterp=>funcPlasmaParams, &
       interpStateData, interpStateDataP, interpsetup=>setup
  use scattered_interp_dens_model_adapter, only : &
       fscatteredinterp=>funcPlasmaParams, &
       scatteredinterpStateData, scatteredinterpStateDataP, &
       scatteredinterpsetup=>setup
  use AT64ThCh_adapter, only : fAT64ThCh=>funcPlasmaParams, &
       AT64ThChStateData, AT64ThChStateDataP
  USE ISO_FORTRAN_ENV ! for OUTPUT_UNIT definition, from 2003 standard
  implicit none
  
  character(len=10000) :: filename, interp_interpfile, ngo_configfile
  integer,parameter :: outfile=50
  character (len=10000) :: buffer
  real(kind=DP), allocatable :: x(:), y(:), z(:) 
  integer :: nx,ny,nz, ind, nspec, i,j,k
  character, allocatable :: data(:)
  real(kind=DP) :: minx,maxx,miny,maxy,minz,maxz, tmpinput, delx,dely,delz
  integer :: modelnum, sz, foundopt
  real(kind=DP), allocatable :: f(:,:,:,:)

  type(ngoStateData),target :: ngo_state_data
  type(ngoStateDataP) :: ngo_state_datap

  type(ngo3dStateData),target :: ngo3d_state_data
  type(ngo3dStateDataP) :: ngo3d_state_datap

  type(gcpmStateData),target :: gcpm_state_data
  type(gcpmStateDataP) :: gcpm_state_datap

  type(simpleStateData),target :: simple_state_data
  type(simpleStateDataP),target :: simple_state_datap

  type(interpStateData),target :: interp_state_data
  type(interpStateDataP) :: interp_state_datap

  type(scatteredinterpStateData),target :: scattered_interp_state_data
  type(scatteredinterpStateDataP) :: scattered_interp_state_datap

  type(AT64ThChStateData),target :: AT64ThCh_state_data
  type(AT64ThChStateDataP) :: AT64ThCh_state_datap

  real(kind=DP), allocatable :: qs(:), Ns(:), ms(:), nus(:)
  real(kind=DP) :: B0(3)

  integer(kind=DP) :: mag_coords  ! If 1, use magnetic dipole coordinates (cartesian). Else, SM.
  real(kind=DP) :: x_mag(3)
  modelnum = 0
  ind=0

  print '(a,f5.2)', 'Dumpmodel version', VERSION

  if( iargc() == 0 ) then
     print *, 'Usage:'
     print *, '  program --param1=value1 --param2=value2 ...'
     print *, '  '
     print *, '--minx      minimum x coordinate'
     print *, '--maxx      maximum x coordinate'
     print *, '--miny      minimum y coordinate'
     print *, '--maxy      maximum y coordinate'
     print *, '--minz      minimum z coordinate'
     print *, '--maxz      maximum z coordinate'
     print *, '--nx        number of points in x direction'
     print *, '--ny        number of points in y direction'
     print *, '--nz        number of points in z direction'
     print *, '--filename  output filename'
     print *, '--modelnum  (1) Ngo model'
     print *, '            (2) GCPM ionosphere model'
     print *, '            (3) Interpolated model (gridded data)'
     print *, '            (4) Interpolated model (scattered data)'
     
     print *, ' Ngo Parameters (required if model 1 is chosen):'
     print *, '   --ngo_configfile     newray input filename'
     print *, '   --yearday            year and day, e.g., 1999098'
     print *, '   --milliseconds_day   milliseconds of day'
     print *, '   --use_tsyganenko     (1=use, 0=do not use)'
     print *, '   --use_igrf           (1=use, 0=do not use)'
     print *, '   --tsyganenko_Pdyn    between 0.5 and 10 nPa'
     print *, '   --tsyganenko_Dst     between -100 and +20 in nT'
     print *, '   --tsyganenko_ByIMF   between -10 and +10 nT'
     print *, '   --tsyganenko_BzIMF   between -10 and +10 nT'
     print *, '   --tsyganenko_W1      TS04_s W1'
     print *, '   --tsyganenko_W2      TS04_s W2'
     print *, '   --tsyganenko_W3      TS04_s W3'
     print *, '   --tsyganenko_W4      TS04_s W4'
     print *, '   --tsyganenko_W5      TS04_s W5'
     print *, '   --tsyganenko_W6      TS04_s W6'
     ! GCPM parameters
     print *, ' GCPM Parameters (required if model 2 is chosen):'
     print *, '   --gcpm_kp            kp index'
     print *, '   --yearday            year and day, e.g., 1999098'
     print *, '   --milliseconds_day   milliseconds of day'
     print *, '   --use_tsyganenko     (1=use, 0=do not use)'
     print *, '   --use_igrf           (1=use, 0=do not use)'
     print *, '   --tsyganenko_Pdyn    between 0.5 and 10 nPa'
     print *, '   --tsyganenko_Dst     between -100 and +20 in nT'
     print *, '   --tsyganenko_ByIMF   between -10 and +10 nT'
     print *, '   --tsyganenko_BzIMF   between -10 and +10 nT'
     print *, '   --tsyganenko_W1      TS04_s W1'
     print *, '   --tsyganenko_W2      TS04_s W2'
     print *, '   --tsyganenko_W3      TS04_s W3'
     print *, '   --tsyganenko_W4      TS04_s W4'
     print *, '   --tsyganenko_W5      TS04_s W5'
     print *, '   --tsyganenko_W6      TS04_s W6'
     ! Interpolated parameters
     print *, ' Interp parameters (required if model 3 is chosen):'
     print *, '   --interp_interpfile  grid filename'
     print *, '   --yearday            year and day, e.g., 1999098'
     print *, '   --milliseconds_day   milliseconds of day'
     print *, '   --use_tsyganenko     (1=use, 0=do not use)'
     print *, '   --use_igrf           (1=use, 0=do not use)'
     print *, '   --tsyganenko_Pdyn    between 0.5 and 10 nPa'
     print *, '   --tsyganenko_Dst     between -100 and +20 in nT'
     print *, '   --tsyganenko_ByIMF   between -10 and +10 nT'
     print *, '   --tsyganenko_BzIMF   between -10 and +10 nT'
     print *, '   --tsyganenko_W1      TS04_s W1'
     print *, '   --tsyganenko_W2      TS04_s W2'
     print *, '   --tsyganenko_W3      TS04_s W3'
     print *, '   --tsyganenko_W4      TS04_s W4'
     print *, '   --tsyganenko_W5      TS04_s W5'
     print *, '   --tsyganenko_W6      TS04_s W6'
     ! Scattered interpolator parameters
     print *, ' Scattered interp parameters (required if model 4 is chosen):'
     print *, '   --interp_interpfile  data filename'
     print *, '   --yearday            year and day, e.g., 1999098'
     print *, '   --milliseconds_day   milliseconds of day'
     print *, '   --use_tsyganenko     (1=use, 0=do not use)'
     print *, '   --use_igrf           (1=use, 0=do not use)'
     print *, '   --tsyganenko_Pdyn    between 0.5 and 10 nPa'
     print *, '   --tsyganenko_Dst     between -100 and +20 in nT'
     print *, '   --tsyganenko_ByIMF   between -10 and +10 nT'
     print *, '   --tsyganenko_BzIMF   between -10 and +10 nT'
     print *, '   --tsyganenko_W1      TS04_s W1'
     print *, '   --tsyganenko_W2      TS04_s W2'
     print *, '   --tsyganenko_W3      TS04_s W3'
     print *, '   --tsyganenko_W4      TS04_s W4'
     print *, '   --tsyganenko_W5      TS04_s W5'
     print *, '   --tsyganenko_W6      TS04_s W6'
     print *, '   --scattered_interp_window_scale'
     print *, '                        window radius scale factor above'
     print *, '                        maximum sample spacing'
     print *, '   --scattered_interp_order' 
     print *, '                        monomial order'
     print *, '   --scattered_interp_exact'
     print *, '                        exact(1) or inexact(0)'
     print *, '   --scattered_interp_local_window_scale'
     print *, '                        scale factor above minimum average'
     print *, '                        sample spacing'
     
     stop
  end if

  ! Read the arguments
  call getopt_named( 'minx', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) minx
  end if
  call getopt_named( 'maxx', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) maxx
  end if
  call getopt_named( 'miny', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) miny
  end if
  call getopt_named( 'maxy', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) maxy
  end if
  call getopt_named( 'minz', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) minz
  end if
  call getopt_named( 'maxz', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) maxz
  end if
  call getopt_named( 'nx', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) tmpinput
     nx = floor(tmpinput)
  end if
  call getopt_named( 'ny', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) tmpinput
     ny = floor(tmpinput)
  end if
  call getopt_named( 'nz', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) tmpinput
     nz = floor(tmpinput)
  end if

  call getopt_named( 'filename', filename, foundopt )

  call getopt_named( 'modelnum', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) tmpinput
     modelnum = floor(tmpinput)
  end if

  allocate(x(nx))
  allocate(y(ny))
  allocate(z(nz))

  ! equivalent to linspace
  delx = (maxx-minx)/(nx-1.0_DP)
  dely = (maxy-miny)/(ny-1.0_DP)
  delz = (maxz-minz)/(nz-1.0_DP)
  if( nx /= 1 ) then
     x = (/ (ind, ind=0,nx-1) /)*delx + minx
  else
     x = (/ minx /)
  end if
  if( ny /= 1 ) then
     y = (/ (ind, ind=0,ny-1) /)*dely + miny
  else
     y = (/ miny /)
  end if
  if( nz /= 1 ) then
     z = (/ (ind, ind=0,nz-1) /)*delz + minz
  else
     z = (/ minz /)
  end if

  print *, 'Common parameters'
  print *, '     minx: ', minx
  print *, '     maxx: ', maxx
  print *, '     miny: ', miny
  print *, '     maxy: ', maxy
  print *, '     minz: ', minz
  print *, '     maxz: ', maxz
  print *, '       nx: ', nx
  print *, '       ny: ', ny
  print *, '       nz: ', nz
  print *, ' filename: ', trim(filename)
  flush(OUTPUT_UNIT)
     
  if( modelnum == 1 ) then
     !!!!!!!!!!!!!!!!!!!!!!! Ngo setup
     ! The Ngo model is the old raytracer plasmasphere model
     nspec=4
     
     ! we need to set up the model paramaters and marshall the setup data
     ! Read the arguments
     ! configuration file
     call getopt_named( 'ngo_configfile', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,'(a)') ngo_configfile
     end if
     ! yearday
     call getopt_named( 'yearday', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        ngo_state_data%itime(1) = floor(tmpinput)
     end if
     ! milliseconds_day
     call getopt_named( 'milliseconds_day', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        ngo_state_data%itime(2) = floor(tmpinput)
     end if
     ! use_tsyganenko
     call getopt_named( 'use_tsyganenko', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        ngo_state_data%use_tsyganenko = floor(tmpinput)
     end if
     ! use_igrf
     call getopt_named( 'use_igrf', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        ngo_state_data%use_igrf = floor(tmpinput)
     end if
     ! tsyganenko_Pdyn
     call getopt_named( 'tsyganenko_Pdyn', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) ngo_state_data%Pdyn
     end if
     ! tsyganenko_Dst
     call getopt_named( 'tsyganenko_Dst', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) ngo_state_data%Dst
     end if
     ! tsyganenko_ByIMF
     call getopt_named( 'tsyganenko_ByIMF', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) ngo_state_data%ByIMF
     end if
     ! tsyganenko_BzIMF
     call getopt_named( 'tsyganenko_BzIMF', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) ngo_state_data%BzIMF
     end if
     ! tsyganenko_W1
     call getopt_named( 'tsyganenko_W1', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) ngo_state_data%W1
     end if
     ! tsyganenko_W2
     call getopt_named( 'tsyganenko_W2', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) ngo_state_data%W2
     end if
     ! tsyganenko_W3
     call getopt_named( 'tsyganenko_W3', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) ngo_state_data%W3
     end if
     ! tsyganenko_W4
     call getopt_named( 'tsyganenko_W4', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) ngo_state_data%W4
     end if
     ! tsyganenko_W5
     call getopt_named( 'tsyganenko_W5', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) ngo_state_data%W5
     end if
     ! tsyganenko_W6
     call getopt_named( 'tsyganenko_W6', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) ngo_state_data%W6
     end if

     ! Marshall our data to the callback
     ! associate a pointer to the state data provided by the user
     ngo_state_dataP%p => ngo_state_data
     ! marshall the data pointer to our function
     sz = size(transfer(ngo_state_datap, data))
     allocate(data(sz))
     data = transfer(ngo_state_dataP, data)
     
     ! Call the setup routine to open a source file
     call ngosetup(ngo_state_data, ngo_configfile)

     print *, 'Ngo Model parameters:'
     print *, '   ngo_configfile:   ', trim(ngo_configfile)
     print *, '   yearday:          ', ngo_state_data%itime(1)
     print *, '   milliseconds_day: ', ngo_state_data%itime(2)
     print *, '   use_tsyganenko:   ', ngo_state_data%use_tsyganenko
     print *, '   use_igrf:         ', ngo_state_data%use_igrf
     print *, '   tsyganenko_Pdyn:  ', ngo_state_data%Pdyn
     print *, '   tsyganenko_Dst:   ', ngo_state_data%Dst
     print *, '   tsyganenko_ByIMF: ', ngo_state_data%ByIMF
     print *, '   tsyganenko_BzIMF: ', ngo_state_data%BzIMF
     print *, '   tsyganenko_W1:    ', ngo_state_data%W1
     print *, '   tsyganenko_W2:    ', ngo_state_data%W2
     print *, '   tsyganenko_W3:    ', ngo_state_data%W3
     print *, '   tsyganenko_W4:    ', ngo_state_data%W4
     print *, '   tsyganenko_W5:    ', ngo_state_data%W5
     print *, '   tsyganenko_W6:    ', ngo_state_data%W6
     flush(OUTPUT_UNIT)

     ! Allocate space for the data
     ! number of species, times 4 (charge, number density, mass, 
     ! collision frequency), plus 3 (3 components of B0)
     allocate(f(nspec*4+3, nx, ny, nz))
     
     ! Grab the data
     do k=1,nz
        flush(OUTPUT_UNIT)
        do j=1,ny
           do i=1,nx
              call fngo((/x(i),y(j),z(k)/), qs, Ns, ms, nus, B0, data)
              f(:,i,j,k) = (/qs, Ns, Ms, nus, B0/)
           end do
        end do
     end do
  elseif( modelnum == 2 ) then
     !!!!!!!!!!!!!!!!!!!!!!! GCPM SETUP
     ! GCPM is a complete plasmasphere model.
     nspec=4
     
     ! we need to set up the model paramaters and marshall the setup data
     ! Read the arguments
     call getopt_named( 'gcpm_kp', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) gcpm_state_data%akp
     end if
     ! yearday
     call getopt_named( 'yearday', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        gcpm_state_data%itime(1) = floor(tmpinput)
     end if
     ! milliseconds_day
     call getopt_named( 'milliseconds_day', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        gcpm_state_data%itime(2) = floor(tmpinput)
     end if
     ! use_tsyganenko
     call getopt_named( 'use_tsyganenko', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        gcpm_state_data%use_tsyganenko = floor(tmpinput)
     end if
     ! use_igrf
     call getopt_named( 'use_igrf', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        gcpm_state_data%use_igrf = floor(tmpinput)
     end if
     ! tsyganenko_Pdyn
     call getopt_named( 'tsyganenko_Pdyn', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) gcpm_state_data%Pdyn
     end if
     ! tsyganenko_Dst
     call getopt_named( 'tsyganenko_Dst', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) gcpm_state_data%Dst
     end if
     ! tsyganenko_ByIMF
     call getopt_named( 'tsyganenko_ByIMF', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) gcpm_state_data%ByIMF
     end if
     ! tsyganenko_BzIMF
     call getopt_named( 'tsyganenko_BzIMF', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) gcpm_state_data%BzIMF
     end if
     ! tsyganenko_W1
     call getopt_named( 'tsyganenko_W1', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) gcpm_state_data%W1
     end if
     ! tsyganenko_W2
     call getopt_named( 'tsyganenko_W2', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) gcpm_state_data%W2
     end if
     ! tsyganenko_W3
     call getopt_named( 'tsyganenko_W3', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) gcpm_state_data%W3
     end if
     ! tsyganenko_W4
     call getopt_named( 'tsyganenko_W4', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) gcpm_state_data%W4
     end if
     ! tsyganenko_W5
     call getopt_named( 'tsyganenko_W5', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) gcpm_state_data%W5
     end if
     ! tsyganenko_W6
     call getopt_named( 'tsyganenko_W6', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) gcpm_state_data%W6
     end if
     ! Fixed MLT:
     call getopt_named( 'MLT', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) gcpm_state_data%MLT
     end if
     call getopt_named( 'fixed_MLT', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        gcpm_state_data%fixed_MLT = floor(tmpinput)
     end if

     ! Marshall our data to the callback
     ! associate a pointer to the state data provided by the user
     gcpm_state_dataP%p => gcpm_state_data
     ! marshall the data pointer to our function
     sz = size(transfer(gcpm_state_datap, data))
     allocate(data(sz))
     data = transfer(gcpm_state_dataP, data)
     
     print *, 'GCPM model parameters:'
     print *, '   gcpm_kp:          ', gcpm_state_data%akp
     print *, '   yearday:          ', gcpm_state_data%itime(1)
     print *, '   milliseconds_day: ', gcpm_state_data%itime(2)
     print *, '   use_tsyganenko:   ', gcpm_state_data%use_tsyganenko
     print *, '   use_igrf:         ', gcpm_state_data%use_igrf
     print *, '   tsyganenko_Pdyn:  ', gcpm_state_data%Pdyn
     print *, '   tsyganenko_Dst:   ', gcpm_state_data%Dst
     print *, '   tsyganenko_ByIMF: ', gcpm_state_data%ByIMF
     print *, '   tsyganenko_BzIMF: ', gcpm_state_data%BzIMF
     print *, '   tsyganenko_W1:    ', gcpm_state_data%W1
     print *, '   tsyganenko_W2:    ', gcpm_state_data%W2
     print *, '   tsyganenko_W3:    ', gcpm_state_data%W3
     print *, '   tsyganenko_W4:    ', gcpm_state_data%W4
     print *, '   tsyganenko_W5:    ', gcpm_state_data%W5
     print *, '   tsyganenko_W6:    ', gcpm_state_data%W6
     print *, '   fixed_MLT:        ', gcpm_state_data%fixed_MLT
     print *, '   MLT:              ', gcpm_state_data%MLT
     flush(OUTPUT_UNIT)

     ! Allocate space for the data
     ! number of species, times 4 (charge, number density, mass, 
     ! collision frequency), plus 3 (3 components of B0)
     allocate(f(nspec*4+3, nx, ny, nz))
     
     ! Grab the data
     do k=1,nz
        do j=1,ny
           do i=1,nx
              ! print *, 'i=',i, '/', nx, 'j=',j, '/', ny, 'k=',k, '/', nz
              flush(OUTPUT_UNIT)
              call fgcpm((/x(i),y(j),z(k)/), qs, Ns, ms, nus, B0, data)
              f(:,i,j,k) = (/qs, Ns, Ms, nus, B0/)
           end do
        end do
     end do
  elseif( modelnum == 3 ) then
     !!!!!!!!!!!!!!!!!!!!!!! INTERPOLATED SETUP
     ! This is an interpolated plasmasphere model, which works from a
     ! gridded source file (e.g. like that generated by
     ! gcpm_dens_model_buildgrid)

     ! we need to set up the model paramaters and marshall the setup data
     ! Read the arguments
     call getopt_named( 'interp_interpfile', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,'(a)') interp_interpfile
     end if
     ! yearday
     call getopt_named( 'yearday', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        interp_state_data%itime(1) = floor(tmpinput)
     end if
     ! milliseconds_day
     call getopt_named( 'milliseconds_day', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        interp_state_data%itime(2) = floor(tmpinput)
     end if
     ! use_tsyganenko
     call getopt_named( 'use_tsyganenko', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        interp_state_data%use_tsyganenko = floor(tmpinput)
     end if
     ! use_igrf
     call getopt_named( 'use_igrf', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        interp_state_data%use_igrf = floor(tmpinput)
     end if
     ! tsyganenko_Pdyn
     call getopt_named( 'tsyganenko_Pdyn', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) interp_state_data%Pdyn
     end if
     ! tsyganenko_Dst
     call getopt_named( 'tsyganenko_Dst', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) interp_state_data%Dst
     end if
     ! tsyganenko_ByIMF
     call getopt_named( 'tsyganenko_ByIMF', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) interp_state_data%ByIMF
     end if
     ! tsyganenko_BzIMF
     call getopt_named( 'tsyganenko_BzIMF', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) interp_state_data%BzIMF
     end if
     ! tsyganenko_W1
     call getopt_named( 'tsyganenko_W1', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) interp_state_data%W1
     end if
     ! tsyganenko_W2
     call getopt_named( 'tsyganenko_W2', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) interp_state_data%W2
     end if
     ! tsyganenko_W3
     call getopt_named( 'tsyganenko_W3', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) interp_state_data%W3
     end if
     ! tsyganenko_W4
     call getopt_named( 'tsyganenko_W4', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) interp_state_data%W4
     end if
     ! tsyganenko_W5
     call getopt_named( 'tsyganenko_W5', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) interp_state_data%W5
     end if
     ! tsyganenko_W6
     call getopt_named( 'tsyganenko_W6', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) interp_state_data%W6
     end if

     ! Marshall our data to the callback
     ! associate a pointer to the state data provided by the user
     interp_state_dataP%p => interp_state_data
     ! marshall the data pointer to our function
     sz = size(transfer(interp_state_datap, data))
     allocate(data(sz))
     data = transfer(interp_state_dataP, data)
     
     print *, 'Interpolated model parameters:'
     print *, '   interp_interpfile:', trim(interp_interpfile)
     print *, '   yearday:          ', interp_state_data%itime(1)
     print *, '   milliseconds_day: ', interp_state_data%itime(2)
     print *, '   use_tsyganenko:   ', interp_state_data%use_tsyganenko
     print *, '   use_igrf:         ', interp_state_data%use_igrf
     print *, '   tsyganenko_Pdyn:  ', interp_state_data%Pdyn
     print *, '   tsyganenko_Dst:   ', interp_state_data%Dst
     print *, '   tsyganenko_ByIMF: ', interp_state_data%ByIMF
     print *, '   tsyganenko_BzIMF: ', interp_state_data%BzIMF
     print *, '   tsyganenko_W1:    ', interp_state_data%W1
     print *, '   tsyganenko_W2:    ', interp_state_data%W2
     print *, '   tsyganenko_W3:    ', interp_state_data%W3
     print *, '   tsyganenko_W4:    ', interp_state_data%W4
     print *, '   tsyganenko_W5:    ', interp_state_data%W5
     print *, '   tsyganenko_W6:    ', interp_state_data%W6

     ! Additional model setup
     print *, 'Reading input file'
     flush(OUTPUT_UNIT)
     call interpsetup(interp_state_data, interp_interpfile)
     print *, 'Done'
     flush(OUTPUT_UNIT)

     ! Allocate space for the data
     ! number of species, times 4 (charge, number density, mass, 
     ! collision frequency), plus 3 (3 components of B0)
     nspec = interp_state_data%nspec
     allocate(f(interp_state_data%nspec*4+3, nx, ny, nz))
     
     ! Grab the data
     do k=1,nz
        flush(OUTPUT_UNIT)
        do j=1,ny
           do i=1,nx
              call finterp((/x(i),y(j),z(k)/), qs, Ns, ms, nus, B0, data)
              f(:,i,j,k) = (/qs, Ns, Ms, nus, B0/)
           end do
        end do
     end do
  elseif( modelnum == 4 ) then
     !!!!!!!!!!!!!!!!!!!!!!! SCATTERED INTERPOLATED SETUP
     ! This is an interpolated plasmasphere model, which works from a
     ! scattered source file (e.g. like that generated by
     ! gcpm_dens_model_buildgrid_random)
     
     ! we need to set up the model paramaters and marshall the setup data
     ! Read the arguments
     call getopt_named( 'interp_interpfile', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,'(a)') interp_interpfile
     end if
     ! yearday
     call getopt_named( 'yearday', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        scattered_interp_state_data%itime(1) = floor(tmpinput)
     end if
     ! milliseconds_day
     call getopt_named( 'milliseconds_day', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        scattered_interp_state_data%itime(2) = floor(tmpinput)
     end if
     ! use_tsyganenko
     call getopt_named( 'use_tsyganenko', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        scattered_interp_state_data%use_tsyganenko = floor(tmpinput)
     end if
     ! use_igrf
     call getopt_named( 'use_igrf', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        scattered_interp_state_data%use_igrf = floor(tmpinput)
     end if
     ! tsyganenko_Pdyn
     call getopt_named( 'tsyganenko_Pdyn', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) scattered_interp_state_data%Pdyn
     end if
     ! tsyganenko_Dst
     call getopt_named( 'tsyganenko_Dst', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) scattered_interp_state_data%Dst
     end if
     ! tsyganenko_ByIMF
     call getopt_named( 'tsyganenko_ByIMF', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) scattered_interp_state_data%ByIMF
     end if
     ! tsyganenko_BzIMF
     call getopt_named( 'tsyganenko_BzIMF', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) scattered_interp_state_data%BzIMF
     end if
     ! tsyganenko_W1
     call getopt_named( 'tsyganenko_W1', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) scattered_interp_state_data%W1
     end if
     ! tsyganenko_W2
     call getopt_named( 'tsyganenko_W2', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) scattered_interp_state_data%W2
     end if
     ! tsyganenko_W3
     call getopt_named( 'tsyganenko_W3', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) scattered_interp_state_data%W3
     end if
     ! tsyganenko_W4
     call getopt_named( 'tsyganenko_W4', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) scattered_interp_state_data%W4
     end if
     ! tsyganenko_W5
     call getopt_named( 'tsyganenko_W5', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) scattered_interp_state_data%W5
     end if
     ! tsyganenko_W6
     call getopt_named( 'tsyganenko_W6', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) scattered_interp_state_data%W6
     end if
     ! scattered_interp_radius
     call getopt_named( 'scattered_interp_window_scale', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) scattered_interp_state_data%window_scale
     end if
     ! scattered_interp_order
     call getopt_named( 'scattered_interp_order', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        scattered_interp_state_data%order = floor(tmpinput)
     end if
     ! scattered_interp_exact
     call getopt_named( 'scattered_interp_exact', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        scattered_interp_state_data%exact = floor(tmpinput)
     end if
     ! scattered_interp_scaled
     scattered_interp_state_data%scaled = 0
     ! scattered_interp_radius_scalefactor
     call getopt_named( 'scattered_interp_local_window_scale', &
          buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) scattered_interp_state_data%local_window_scale
     end if
     
     
     ! Marshall our data to the callback
     ! associate a pointer to the state data provided by the user
     scattered_interp_state_dataP%p => scattered_interp_state_data
     ! marshall the data pointer to our function
     sz = size(transfer(scattered_interp_state_datap, data))
     allocate(data(sz))
     data = transfer(scattered_interp_state_dataP, data)
     
     print *, 'Scattered interpolator model parameters:'
     print *, '   interp_interpfile:', trim(interp_interpfile)
     print *, '   yearday:          ', scattered_interp_state_data%itime(1)
     print *, '   milliseconds_day: ', scattered_interp_state_data%itime(2)
     print *, '   use_tsyganenko:   ', &
          scattered_interp_state_data%use_tsyganenko
     print *, '   use_igrf:         ', scattered_interp_state_data%use_igrf
     print *, '   tsyganenko_Pdyn:  ', scattered_interp_state_data%Pdyn
     print *, '   tsyganenko_Dst:   ', scattered_interp_state_data%Dst
     print *, '   tsyganenko_ByIMF: ', scattered_interp_state_data%ByIMF
     print *, '   tsyganenko_BzIMF: ', scattered_interp_state_data%BzIMF
     print *, '   tsyganenko_W1:    ', scattered_interp_state_data%W1
     print *, '   tsyganenko_W2:    ', scattered_interp_state_data%W2
     print *, '   tsyganenko_W3:    ', scattered_interp_state_data%W3
     print *, '   tsyganenko_W4:    ', scattered_interp_state_data%W4
     print *, '   tsyganenko_W5:    ', scattered_interp_state_data%W5
     print *, '   tsyganenko_W6:    ', scattered_interp_state_data%W6
     print *, '   scattered_interp_window_scale: ', &
          scattered_interp_state_data%window_scale
     print *, '   scattered_interp_order: ', &
          scattered_interp_state_data%order
     print *, '   scattered_interp_exact: ', &
          scattered_interp_state_data%exact
     print *, '   scattered_interp_local_window_scale: ', &
          scattered_interp_state_data%local_window_scale

     ! Additional model setup
     print *, 'Reading input file'
     flush(OUTPUT_UNIT)
     call scatteredinterpsetup(&
          scattered_interp_state_data, interp_interpfile)
     print *, 'Done'
     flush(OUTPUT_UNIT)

     ! Allocate space for the data
     ! number of species, times 4 (charge, number density, mass, 
     ! collision frequency), plus 3 (3 components of B0)
     allocate(f(scattered_interp_state_data%nspec*4+3, nx, ny, nz))
     allocate(Ns(scattered_interp_state_data%nspec))
     nspec = scattered_interp_state_data%nspec
     
     ! Grab the data
     do k=1,nz
        flush(OUTPUT_UNIT)
        do j=1,ny
           do i=1,nx
              call fscatteredinterp((/x(i),y(j),z(k)/), &
                   qs, Ns, ms, nus, B0, data)
              f(:,i,j,k) = (/qs, Ns, Ms, nus, B0/)
           end do
        end do
     end do
  elseif( modelnum == 5 ) then
     !!!!!!!!!!!!!!!!!!!!!!! Ngo setup
     ! The Ngo model is the old raytracer plasmasphere model
     nspec=4
     
     ! we need to set up the model paramaters and marshall the setup data
     ! Read the arguments
     ! configuration file
     call getopt_named( 'ngo_configfile', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,'(a)') ngo_configfile
     end if
     ! yearday
     call getopt_named( 'yearday', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        ngo3d_state_data%itime(1) = floor(tmpinput)
     end if
     ! milliseconds_day
     call getopt_named( 'milliseconds_day', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        ngo3d_state_data%itime(2) = floor(tmpinput)
     end if
     ! use_tsyganenko
     call getopt_named( 'use_tsyganenko', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        ngo3d_state_data%use_tsyganenko = floor(tmpinput)
     end if
     ! use_igrf
     call getopt_named( 'use_igrf', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        ngo3d_state_data%use_igrf = floor(tmpinput)
     end if
     call getopt_named( 'kp', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        ngo3d_state_data%kp = floor(tmpinput)
     end if
     ! tsyganenko_Pdyn
     call getopt_named( 'tsyganenko_Pdyn', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) ngo3d_state_data%Pdyn
     end if
     ! tsyganenko_Dst
     call getopt_named( 'tsyganenko_Dst', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) ngo3d_state_data%Dst
     end if
     ! tsyganenko_ByIMF
     call getopt_named( 'tsyganenko_ByIMF', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) ngo3d_state_data%ByIMF
     end if
     ! tsyganenko_BzIMF
     call getopt_named( 'tsyganenko_BzIMF', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) ngo3d_state_data%BzIMF
     end if
     ! tsyganenko_W1
     call getopt_named( 'tsyganenko_W1', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) ngo3d_state_data%W1
     end if
     ! tsyganenko_W2
     call getopt_named( 'tsyganenko_W2', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) ngo3d_state_data%W2
     end if
     ! tsyganenko_W3
     call getopt_named( 'tsyganenko_W3', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) ngo3d_state_data%W3
     end if
     ! tsyganenko_W4
     call getopt_named( 'tsyganenko_W4', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) ngo3d_state_data%W4
     end if
     ! tsyganenko_W5
     call getopt_named( 'tsyganenko_W5', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) ngo3d_state_data%W5
     end if
     ! tsyganenko_W6
     call getopt_named( 'tsyganenko_W6', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) ngo3d_state_data%W6
     end if
     ! Fixed MLT:
     call getopt_named( 'MLT', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) ngo3d_state_data%MLT
     end if
     call getopt_named( 'fixed_MLT', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        ngo3d_state_data%fixed_MLT = floor(tmpinput)
     end if


     ! Marshall our data to the callback
     ! associate a pointer to the state data provided by the user
     ngo3d_state_dataP%p => ngo3d_state_data
     ! marshall the data pointer to our function
     sz = size(transfer(ngo3d_state_datap, data))
     allocate(data(sz))
     data = transfer(ngo3d_state_dataP, data)
     
     ! Call the setup routine to open a source file
     call ngo3dsetup(ngo3d_state_data, ngo_configfile)

     print *, 'Ngo Model parameters:'
     print *, '   ngo_configfile:   ', trim(ngo_configfile)
     print *, '   yearday:          ', ngo3d_state_data%itime(1)
     print *, '   milliseconds_day: ', ngo3d_state_data%itime(2)
     print *, '   use_tsyganenko:   ', ngo3d_state_data%use_tsyganenko
     print *, '   use_igrf:         ', ngo3d_state_data%use_igrf
     print *, '   kp:           ', ngo3d_state_data%kp
     print *, '   tsyganenko_Pdyn:  ', ngo3d_state_data%Pdyn
     print *, '   tsyganenko_Dst:   ', ngo3d_state_data%Dst
     print *, '   tsyganenko_ByIMF: ', ngo3d_state_data%ByIMF
     print *, '   tsyganenko_BzIMF: ', ngo3d_state_data%BzIMF
     print *, '   tsyganenko_W1:    ', ngo3d_state_data%W1
     print *, '   tsyganenko_W2:    ', ngo3d_state_data%W2
     print *, '   tsyganenko_W3:    ', ngo3d_state_data%W3
     print *, '   tsyganenko_W4:    ', ngo3d_state_data%W4
     print *, '   tsyganenko_W5:    ', ngo3d_state_data%W5
     print *, '   tsyganenko_W6:    ', ngo3d_state_data%W6
     print *, '   fixed_MLT:        ', ngo3d_state_data%fixed_MLT
     print *, '   MLT:              ', ngo3d_state_data%MLT
     flush(OUTPUT_UNIT)


     ! Allocate space for the data
     ! number of species, times 4 (charge, number density, mass, 
     ! collision frequency), plus 3 (3 components of B0)
     allocate(f(nspec*4+3, nx, ny, nz))
     
     ! Grab the data
     do k=1,nz
        flush(OUTPUT_UNIT)
        do j=1,ny
           do i=1,nx
              call fngo3d((/x(i),y(j),z(k)/), qs, Ns, ms, nus, B0, data)
              f(:,i,j,k) = (/qs, Ns, Ms, nus, B0/)
           end do
        end do
     end do




  elseif( modelnum == 6 ) then
     !!!!!!!!!!!!!!!!!!!!!!! Simple 3d plasmasphere setup
     ! Austin's simplified plasmasphere model
     
     ! we need to set up the model paramaters and marshall the setup data
     ! Read the arguments
     ! configuration file
     call getopt_named( 'ngo_configfile', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,'(a)') ngo_configfile
     end if
     ! yearday
     call getopt_named( 'yearday', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        simple_state_data%itime(1) = floor(tmpinput)
     end if
     ! milliseconds_day
     call getopt_named( 'milliseconds_day', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        simple_state_data%itime(2) = floor(tmpinput)
     end if
     ! use_tsyganenko
     call getopt_named( 'use_tsyganenko', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        simple_state_data%use_tsyganenko = floor(tmpinput)
     end if
     ! use_igrf
     call getopt_named( 'use_igrf', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        simple_state_data%use_igrf = floor(tmpinput)
     end if
     call getopt_named( 'kp', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        simple_state_data%kp = tmpinput
     end if
     ! tsyganenko_Pdyn
     call getopt_named( 'tsyganenko_Pdyn', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) simple_state_data%Pdyn
     end if
     ! tsyganenko_Dst
     call getopt_named( 'tsyganenko_Dst', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) simple_state_data%Dst
     end if
     ! tsyganenko_ByIMF
     call getopt_named( 'tsyganenko_ByIMF', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) simple_state_data%ByIMF
     end if
     ! tsyganenko_BzIMF
     call getopt_named( 'tsyganenko_BzIMF', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) simple_state_data%BzIMF
     end if
     ! tsyganenko_W1
     call getopt_named( 'tsyganenko_W1', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) simple_state_data%W1
     end if
     ! tsyganenko_W2
     call getopt_named( 'tsyganenko_W2', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) simple_state_data%W2
     end if
     ! tsyganenko_W3
     call getopt_named( 'tsyganenko_W3', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) simple_state_data%W3
     end if
     ! tsyganenko_W4
     call getopt_named( 'tsyganenko_W4', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) simple_state_data%W4
     end if
     ! tsyganenko_W5
     call getopt_named( 'tsyganenko_W5', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) simple_state_data%W5
     end if
     ! tsyganenko_W6
     call getopt_named( 'tsyganenko_W6', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) simple_state_data%W6
     end if
     ! Fixed MLT:
     call getopt_named( 'MLT', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) simple_state_data%MLT
     end if
     call getopt_named( 'fixed_MLT', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        simple_state_data%fixed_MLT = floor(tmpinput)
     end if
    call getopt_named( 'mag_coords', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) mag_coords

        ! simple_state_data%fixed_MLT = floor(tmpinput)
     end if


     ! Marshall our data to the callback
     ! associate a pointer to the state data provided by the user
     simple_state_dataP%p => simple_state_data
     ! marshall the data pointer to our function
     sz = size(transfer(simple_state_dataP, data))
     allocate(data(sz))
     data = transfer(simple_state_dataP, data)


     
     ! ! Call the setup routine to open a source file
     ! call simple3dsetup("/shared/users/asousa/software/raytracer_DP/bin/iri_dayside_600km.dat","/shared/users/asousa/software/raytracer_DP/bin/iri_nightside_600km.dat")

     print *, 'Simplified 3d Model parameters:'
     ! print *, '   ngo_configfile:   ', trim(ngo_configfile)
     print *, '   yearday:          ', simple_state_data%itime(1)
     print *, '   milliseconds_day: ', simple_state_data%itime(2)
     print *, '   use_tsyganenko:   ', simple_state_data%use_tsyganenko
     print *, '   use_igrf:         ', simple_state_data%use_igrf
     print *, '   kp:               ', simple_state_data%kp
     print *, '   tsyganenko_Pdyn:  ', simple_state_data%Pdyn
     print *, '   tsyganenko_Dst:   ', simple_state_data%Dst
     print *, '   tsyganenko_ByIMF: ', simple_state_data%ByIMF
     print *, '   tsyganenko_BzIMF: ', simple_state_data%BzIMF
     print *, '   tsyganenko_W1:    ', simple_state_data%W1
     print *, '   tsyganenko_W2:    ', simple_state_data%W2
     print *, '   tsyganenko_W3:    ', simple_state_data%W3
     print *, '   tsyganenko_W4:    ', simple_state_data%W4
     print *, '   tsyganenko_W5:    ', simple_state_data%W5
     print *, '   tsyganenko_W6:    ', simple_state_data%W6
     print *, '   mag_coords:       ', mag_coords
     flush(OUTPUT_UNIT)



     ! Allocate space for the data
     ! number of species, times 4 (charge, number density, mass, 
     ! collision frequency), plus 3 (3 components of B0)
     nspec = 4 ! always four species
     allocate(f(nspec*4+3, nx, ny, nz))
     
     ! Grab the data
     do k=1,nz
        flush(OUTPUT_UNIT)
        do j=1,ny
           do i=1,nx
              if (mag_coords==1) then
                ! geomag coordinates
                call SM_TO_MAG_D(simple_state_data%itime, (/x(i),y(j),z(k)/), x_mag)
                call fsimple(x_mag, qs, Ns, ms, nus, B0, data)
              else
                ! SM coordinates
                call fsimple((/x(i),y(j),z(k)/), qs, Ns, ms, nus, B0, data)
              end if
              f(:,i,j,k) = (/qs, Ns, Ms, nus, B0/)
           end do
        end do
     end do
  
  elseif( modelnum == 7 ) then
     !!! diffusive eq.
     call getopt_named( 'ngo_configfile', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,'(a)') ngo_configfile
     end if
     ! yearday
     call getopt_named( 'yearday', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        AT64ThCh_state_data%itime(1) = floor(tmpinput)
     end if
     ! milliseconds_day
     call getopt_named( 'milliseconds_day', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        AT64ThCh_state_data%itime(2) = floor(tmpinput)
     end if
     ! use_tsyganenko
     call getopt_named( 'use_tsyganenko', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        AT64ThCh_state_data%use_tsyganenko = floor(tmpinput)
     end if
     ! use_igrf
     call getopt_named( 'use_igrf', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        AT64ThCh_state_data%use_igrf = floor(tmpinput)
     end if
     call getopt_named( 'gcpm_kp', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        AT64ThCh_state_data%gcpm_kp = floor(tmpinput)
     end if
     ! tsyganenko_Pdyn
     call getopt_named( 'tsyganenko_Pdyn', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) AT64ThCh_state_data%Pdyn
     end if
     ! tsyganenko_Dst
     call getopt_named( 'tsyganenko_Dst', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) AT64ThCh_state_data%Dst
     end if
     ! tsyganenko_ByIMF
     call getopt_named( 'tsyganenko_ByIMF', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) AT64ThCh_state_data%ByIMF
     end if
     ! tsyganenko_BzIMF
     call getopt_named( 'tsyganenko_BzIMF', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) AT64ThCh_state_data%BzIMF
     end if
     ! tsyganenko_W1
     call getopt_named( 'tsyganenko_W1', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) AT64ThCh_state_data%W1
     end if
     ! tsyganenko_W2
     call getopt_named( 'tsyganenko_W2', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) AT64ThCh_state_data%W2
     end if
     ! tsyganenko_W3
     call getopt_named( 'tsyganenko_W3', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) AT64ThCh_state_data%W3
     end if
     ! tsyganenko_W4
     call getopt_named( 'tsyganenko_W4', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) AT64ThCh_state_data%W4
     end if
     ! tsyganenko_W5
     call getopt_named( 'tsyganenko_W5', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) AT64ThCh_state_data%W5
     end if
     ! tsyganenko_W6
     call getopt_named( 'tsyganenko_W6', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) AT64ThCh_state_data%W6
     end if
    call getopt_named( 'mag_coords', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) mag_coords

     end if


     ! Marshall our data to the callback
     ! associate a pointer to the state data provided by the user
     AT64ThCh_state_dataP%p => AT64ThCh_state_data
     ! marshall the data pointer to our function
     sz = size(transfer(AT64ThCh_state_dataP, data))
     allocate(data(sz))
     data = transfer(AT64ThCh_state_dataP, data)

     flush(OUTPUT_UNIT)



     ! Allocate space for the data
     ! number of species, times 4 (charge, number density, mass, 
     ! collision frequency), plus 3 (3 components of B0)
     nspec = 3 
     allocate(f(nspec*4+3, nx, ny, nz))
     
     ! Grab the data
     do k=1,nz
        flush(OUTPUT_UNIT)
        do j=1,ny
           do i=1,nx
              if (mag_coords==1) then
                ! geomag coordinates
                call SM_TO_MAG_D(AT64ThCh_state_data%itime, (/x(i),y(j),z(k)/), x_mag)
                call fAT64ThCh(x_mag, qs, Ns, ms, nus, B0, data)
              else
                ! SM coordinates
                call fAT64ThCh((/x(i),y(j),z(k)/), qs, Ns, ms, nus, B0, data)
              end if
              f(:,i,j,k) = (/qs, Ns, Ms, nus, B0/)
           end do
        end do
     end do


  end if

  open(unit=outfile, file=filename, status="replace")
  ! print out the headers
  write(outfile, '(5i10)'), nspec, nx,ny,nz
  write(outfile, '(6es24.15e3)'), minx,maxx, miny,maxy, minz,maxz
  ! Print out the data
  write(outfile, '(es24.15e3)'), f
  flush(outfile)
  ! close
  close(unit=outfile)
  
end program dumpmodel
