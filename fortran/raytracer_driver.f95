program raytracer_driver
  use types
  use util
  use constants, only : R_E, PI, VERSION
  use ngo_dens_model_adapter, only : fngo=>funcPlasmaParams, ngoStateData, &
       ngoStateDataP, ngosetup=>setup
  use ngo_3d_dens_model_adapter, only : fngo3d=>funcPlasmaParams, ngo3dStateData, &
       ngo3dStateDataP, ngo3dsetup=>setup

  use simple_3d_model_adapter, only : fsimple=>funcPlasmaParams, simpleStateData, &
       simpleStateDataP

  use gcpm_dens_model_adapter, only : fgcpm=>funcPlasmaParams, &
       gcpmStateData, gcpmStateDataP
  use interp_dens_model_adapter, only : finterp=>funcPlasmaParams, &
       interpStateData, interpStateDataP, interpsetup=>setup
  


  use raytracer, only : raytracer_run, raytracer_stopconditions

  ! Temporary file for me to safely fuck with the gradients
  ! use raytracer_spherical, only : raytracer_run, raytracer_stopconditions

  use scattered_interp_dens_model_adapter, only : &
       fscatteredinterp=>funcPlasmaParams, &
       scatteredinterpStateData, scatteredinterpStateDataP, &
       scatteredinterpsetup=>setup
  USE ISO_FORTRAN_ENV ! for OUTPUT_UNIT definition, from 2003 standard
  implicit none

  real(kind=DP) :: pos0(3), w, dir0(3), dt0, dtmax, maxerr, tmax
  integer :: fixedstep, root

  real(kind=DP), allocatable :: pos(:,:), time(:), vprel(:,:), vgrel(:,:), &
                                n(:,:), B0(:,:), &
                                qs(:,:), ms(:,:), Ns(:,:), nus(:,:)

  character (len=10000) :: buffer, inputraysfile, outputfile
  integer :: stopcond, i, maxsteps, j

  real(kind=DP) :: delSP, delDP, minalt
  
  character(len=10000) :: interp_interpfile, ngo_configfile
  integer,parameter :: infile=51, outfile=50
  character, allocatable :: data(:)
  real(kind=DP) :: tmpinput
  integer :: modelnum, sz, status, raynum, foundopt, outputper

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

  modelnum = 0

  ! print '(a,f5.2)', 'Raytracer version', VERSION

  if( iargc() == 0 ) then
     ! print *, 'Usage:'
     ! print *, '  program --param1=value1 --param2=value2 ...'
     ! print *, '  '
     ! print *, '--dt0           initial timestep in seconds'
     ! print *, '--dtmax         maximum timestep in seconds'
     ! print *, '--tmax          maximum time for simulation'
     ! print *, '--root          root number (2=whistler at VLF frequencies)'
     ! print *, '--fixedstep     fixed timesteps (1) or adaptive (0)'
     ! print *, '--maxerr        maximum error for adaptive timestepping'
     ! print *, '--maxsteps      maximum number of timesteps'
     ! print *, '--minalt        minimum altitude'
     ! print *, '--inputraysfile ray input filename'
     ! print *, '--outputfile    output filename'
     ! print *, '--outputper     output fields only every "outputper" steps'
     ! print *, '--modelnum      (1) Ngo model'
     ! print *, '                (2) GCPM ionosphere model'
     ! print *, '                (3) Interpolated model (gridded data)'
     ! print *, '                (4) Interpolated model (scattered data)'
     
     ! Ngo parameters
     ! print *, ' Ngo Parameters (required if model 1 is chosen):'
     ! print *, '   --ngo_configfile     newray input filename'
     ! print *, '   --yearday            year and day, e.g., 1999098'
     ! print *, '   --milliseconds_day   milliseconds of day'
     ! print *, '   --use_tsyganenko     (1=use, 0=do not use)'
     ! print *, '   --use_igrf           (1=use, 0=do not use)'
     ! print *, '   --tsyganenko_Pdyn    between 0.5 and 10 nPa'
     ! print *, '   --tsyganenko_Dst     between -100 and +20 in nT'
     ! print *, '   --tsyganenko_ByIMF   between -10 and +10 nT'
     ! print *, '   --tsyganenko_BzIMF   between -10 and +10 nT'
     ! print *, '   --tsyganenko_W1      TS04_s W1'
     ! print *, '   --tsyganenko_W2      TS04_s W2'
     ! print *, '   --tsyganenko_W3      TS04_s W3'
     ! print *, '   --tsyganenko_W4      TS04_s W4'
     ! print *, '   --tsyganenko_W5      TS04_s W5'
     ! print *, '   --tsyganenko_W6      TS04_s W6'
     ! GCPM parameters
     ! print *, ' GCPM Parameters (required if model 2 is chosen):'
     ! print *, '   --gcpm_kp            kp index'
     ! print *, '   --yearday            year and day, e.g., 1999098'
     ! print *, '   --milliseconds_day   milliseconds of day'
     ! print *, '   --use_tsyganenko     (1=use, 0=do not use)'
     ! print *, '   --use_igrf           (1=use, 0=do not use)'
     ! print *, '   --tsyganenko_Pdyn    between 0.5 and 10 nPa'
     ! print *, '   --tsyganenko_Dst     between -100 and +20 in nT'
     ! print *, '   --tsyganenko_ByIMF   between -10 and +10 nT'
     ! print *, '   --tsyganenko_BzIMF   between -10 and +10 nT'
     ! print *, '   --tsyganenko_W1      TS04_s W1'
     ! print *, '   --tsyganenko_W2      TS04_s W2'
     ! print *, '   --tsyganenko_W3      TS04_s W3'
     ! print *, '   --tsyganenko_W4      TS04_s W4'
     ! print *, '   --tsyganenko_W5      TS04_s W5'
     ! print *, '   --tsyganenko_W6      TS04_s W6'
     ! Interpolated parameters
     ! print *, ' Interp parameters (required if model 3 is chosen):'
     ! print *, '   --interp_interpfile  grid filename'
     ! print *, '   --yearday            year and day, e.g., 1999098'
     ! print *, '   --milliseconds_day   milliseconds of day'
     ! print *, '   --use_tsyganenko     (1=use, 0=do not use)'
     ! print *, '   --use_igrf           (1=use, 0=do not use)'
     ! print *, '   --tsyganenko_Pdyn    between 0.5 and 10 nPa'
     ! print *, '   --tsyganenko_Dst     between -100 and +20 in nT'
     ! print *, '   --tsyganenko_ByIMF   between -10 and +10 nT'
     ! print *, '   --tsyganenko_BzIMF   between -10 and +10 nT'
     ! print *, '   --tsyganenko_W1      TS04_s W1'
     ! print *, '   --tsyganenko_W2      TS04_s W2'
     ! print *, '   --tsyganenko_W3      TS04_s W3'
     ! print *, '   --tsyganenko_W4      TS04_s W4'
     ! print *, '   --tsyganenko_W5      TS04_s W5'
     ! print *, '   --tsyganenko_W6      TS04_s W6'
     ! Scattered interpolator parameters
     ! print *, ' Scattered interp parameters (required if model 4 is chosen):'
     ! print *, '   --interp_interpfile  data filename'
     ! print *, '   --yearday            year and day, e.g., 1999098'
     ! print *, '   --milliseconds_day   milliseconds of day'
     ! print *, '   --use_tsyganenko     (1=use, 0=do not use)'
     ! print *, '   --use_igrf           (1=use, 0=do not use)'
     ! print *, '   --tsyganenko_Pdyn    between 0.5 and 10 nPa'
     ! print *, '   --tsyganenko_Dst     between -100 and +20 in nT'
     ! print *, '   --tsyganenko_ByIMF   between -10 and +10 nT'
     ! print *, '   --tsyganenko_BzIMF   between -10 and +10 nT'
     ! print *, '   --tsyganenko_W1      TS04_s W1'
     ! print *, '   --tsyganenko_W2      TS04_s W2'
     ! print *, '   --tsyganenko_W3      TS04_s W3'
     ! print *, '   --tsyganenko_W4      TS04_s W4'
     ! print *, '   --tsyganenko_W5      TS04_s W5'
     ! print *, '   --tsyganenko_W6      TS04_s W6'
     ! print *, '   --scattered_interp_window_scale'
     ! print *, '                        window radius scale factor above'
     ! print *, '                        maximum sample spacing'
     ! print *, '   --scattered_interp_order' 
     ! print *, '                        monomial order'
     ! print *, '   --scattered_interp_exact'
     ! print *, '                        exact(1) or inexact(0)'
     ! print *, '   --scattered_interp_local_window_scale'
     ! print *, '                        scale factor above minimum average'
     ! print *, '                        sample spacing'

     flush(OUTPUT_UNIT)

     stop
  end if

  ! Initial dt
  call getopt_named( 'dt0', buffer, foundopt )
  if( foundopt == 1 ) then
     read(buffer,*) dt0
  end if
  call getopt_named( 'dtmax', buffer, foundopt )
  if( foundopt == 1 ) then
     read(buffer,*) dtmax
  end if
  call getopt_named( 'tmax', buffer, foundopt )
  if( foundopt == 1 ) then
     read(buffer,*) tmax
  end if
  call getopt_named( 'root', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) tmpinput
     root = floor(tmpinput)
  end if
  call getopt_named( 'fixedstep', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) tmpinput
     fixedstep = floor(tmpinput)
  end if
  call getopt_named( 'maxerr', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) maxerr
  end if
  call getopt_named( 'maxsteps', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) maxsteps
  end if
  call getopt_named( 'minalt', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) minalt
  end if
  call getopt_named( 'inputraysfile', inputraysfile, foundopt )
  call getopt_named( 'outputfile', outputfile, foundopt )
  call getopt_named( 'modelnum', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) tmpinput
     modelnum = floor(tmpinput)
  end if
  
  outputper = 1
  call getopt_named( 'outputper', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) tmpinput
     outputper = floor(tmpinput)
  end if
  
  ! print *, 'Common parameters'
  ! print *, '          dt0: ', dt0
  ! print *, '        dtmax: ', dtmax
  ! print *, '         tmax: ', tmax
  ! print *, '         root: ', root
  ! print *, '    fixedstep: ', fixedstep
  ! print *, '       maxerr: ', maxerr
  ! print *, '     maxsteps: ', maxsteps
  ! print *, '       minalt: ', minalt
  ! print *, 'inputraysfile: ', trim(inputraysfile)
  ! print *, '   outputfile: ', trim(outputfile)
  ! print *, '    outputper: ', outputper
  ! print *, '     modelnum: ', modelnum
  flush(OUTPUT_UNIT)

  ! Such a large value is necessary because many of the
  ! plasmasphere/magnetosphere models use only real precision 7
  ! significant digits in a single.  From numerical tests, 1e-3 is too
  ! large.  1e-6 is inaccurate.  1e-4 seems about a good usable
  ! number.  For double-precision models, obviously something higher
  ! can be used (1e-8 or thereabouts).
  delSP = 1.0e-4_DP
  delDP = 1.0e-6_DP
  ! del = 1.0e-6_DP
  ! del = 1.0e-8_DP

  if( modelnum == 1 ) then
     !!!!!!!!!!!!!!!!!!!!!!! Ngo setup
     ! The Ngo model is the old raytracer plasmasphere model
     
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

     ! print *, 'Ngo Model parameters:'
     ! print *, '   ngo_configfile:   ', trim(ngo_configfile)
     ! print *, '   yearday:          ', ngo_state_data%itime(1)
     ! print *, '   milliseconds_day: ', ngo_state_data%itime(2)
     ! print *, '   use_tsyganenko:   ', ngo_state_data%use_tsyganenko
     ! print *, '   use_igrf:         ', ngo_state_data%use_igrf
     ! print *, '   tsyganenko_Pdyn:  ', ngo_state_data%Pdyn
     ! print *, '   tsyganenko_Dst:   ', ngo_state_data%Dst
     ! print *, '   tsyganenko_ByIMF: ', ngo_state_data%ByIMF
     ! print *, '   tsyganenko_BzIMF: ', ngo_state_data%BzIMF
     ! print *, '   tsyganenko_W1:    ', ngo_state_data%W1
     ! print *, '   tsyganenko_W2:    ', ngo_state_data%W2
     ! print *, '   tsyganenko_W3:    ', ngo_state_data%W3
     ! print *, '   tsyganenko_W4:    ', ngo_state_data%W4
     ! print *, '   tsyganenko_W5:    ', ngo_state_data%W5
     ! print *, '   tsyganenko_W6:    ', ngo_state_data%W6
     flush(OUTPUT_UNIT)

  elseif( modelnum == 2 ) then
     !!!!!!!!!!!!!!!!!!!!!!! GCPM SETUP
     ! GCPM is a complete plasmasphere model.
     
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
     
     ! print *, 'GCPM model parameters:'
     ! print *, '   gcpm_kp:          ', gcpm_state_data%akp
     ! print *, '   yearday:          ', gcpm_state_data%itime(1)
     ! print *, '   milliseconds_day: ', gcpm_state_data%itime(2)
     ! print *, '   use_tsyganenko:   ', gcpm_state_data%use_tsyganenko
     ! print *, '   use_igrf:         ', gcpm_state_data%use_igrf
     ! print *, '   tsyganenko_Pdyn:  ', gcpm_state_data%Pdyn
     ! print *, '   tsyganenko_Dst:   ', gcpm_state_data%Dst
     ! print *, '   tsyganenko_ByIMF: ', gcpm_state_data%ByIMF
     ! print *, '   tsyganenko_BzIMF: ', gcpm_state_data%BzIMF
     ! print *, '   tsyganenko_W1:    ', gcpm_state_data%W1
     ! print *, '   tsyganenko_W2:    ', gcpm_state_data%W2
     ! print *, '   tsyganenko_W3:    ', gcpm_state_data%W3
     ! print *, '   tsyganenko_W4:    ', gcpm_state_data%W4
     ! print *, '   tsyganenko_W5:    ', gcpm_state_data%W5
     ! print *, '   tsyganenko_W6:    ', gcpm_state_data%W6
     ! print *, '   fixed_MLT:        ', gcpm_state_data%fixed_MLT
     ! print *, '   MLT:              ', gcpm_state_data%MLT


     flush(OUTPUT_UNIT)

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

     print *, 'Interpolatod model parameters:'
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
     flush(OUTPUT_UNIT)

     ! Additional model setup
     print *, 'Reading input file'
     flush(OUTPUT_UNIT)
     call interpsetup(interp_state_data, interp_interpfile)
     print *, 'Done'
     flush(OUTPUT_UNIT)

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

  elseif( modelnum == 5 ) then
     !!!!!!!!!!!!!!!!!!!!!!! Ngo setup
     ! The Ngo model is the old raytracer plasmasphere model
     
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
        ngo3d_state_data%kp = tmpinput
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
     print *, '   ngo_kp:           ', ngo3d_state_data%kp
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
     flush(OUTPUT_UNIT)

  elseif( modelnum == 6 ) then
     !!!!!!!!!!!!!!!!!!!!!!! Ngo setup
     ! The Ngo model is the old raytracer plasmasphere model
     
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

     ! Marshall our data to the callback
     ! associate a pointer to the state data provided by the user
     simple_state_dataP%p => simple_state_data
     ! marshall the data pointer to our function
     sz = size(transfer(simple_state_dataP, data))
     allocate(data(sz))
     data = transfer(simple_state_dataP, data)
     
     ! ! Call the setup routine to open a source file
     ! call ngo3dsetup(ngo3d_state_data, ngo_configfile)

     print *, 'Ngo Model parameters:'
     ! print *, '   ngo_configfile:   ', trim(ngo_configfile)
     print *, '   yearday:          ', simple_state_data%itime(1)
     print *, '   milliseconds_day: ', simple_state_data%itime(2)
     print *, '   use_tsyganenko:   ', simple_state_data%use_tsyganenko
     print *, '   use_igrf:         ', simple_state_data%use_igrf
     print *, '   ngo_kp:           ', simple_state_data%kp
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
     flush(OUTPUT_UNIT)




  end if

  open(unit=infile, file=inputraysfile, status="old")
  open(unit=outfile, file=outputfile, status="replace")
  print *, ''
  raynum = 1
  do
     ! Read in the next ray starting point and direction
     read(infile, *,iostat=status), pos0, dir0, w
     if( status /= 0 ) then
        ! end of file or error, abort the loop
        exit
     end if
     print *, 'ray ', raynum, ', pos0=', pos0, ', dir0=', dir0, ', w=', w
     flush(OUTPUT_UNIT)
     if( modelnum == 1 ) then
        call raytracer_run( &
             pos,time,vprel,vgrel,n,&
             B0, qs, ms, Ns, nus, stopcond, &
             pos0, dir0, w, dt0, dtmax, maxerr, maxsteps, minalt, root, tmax, &
             fixedstep, delSP, fngo, data, raytracer_stopconditions)
     elseif( modelnum == 2 ) then
        call raytracer_run( &
             pos,time,vprel,vgrel,n, &
             B0, qs, ms, Ns, nus, stopcond, &
             pos0, dir0, w, dt0, dtmax, maxerr, maxsteps, minalt, root, tmax, &
             fixedstep, delSP, fgcpm, data, raytracer_stopconditions)
     elseif( modelnum == 3 ) then
        call raytracer_run( &
             pos,time,vprel,vgrel,n, &
             B0, qs, ms, Ns, nus, stopcond, &
             pos0, dir0, w, dt0, dtmax, maxerr, maxsteps, minalt, root, tmax, &
             fixedstep, delDP, finterp, data, raytracer_stopconditions)
     elseif( modelnum == 4 ) then
        call raytracer_run( &
             pos,time,vprel,vgrel,n, &
             B0, qs, ms, Ns, nus, stopcond, &
             pos0, dir0, w, dt0, dtmax, maxerr, maxsteps, minalt, root, tmax, &
             fixedstep, delDP, fscatteredinterp, data, raytracer_stopconditions)
     elseif( modelnum == 5 ) then
        call raytracer_run( &
             pos,time,vprel,vgrel,n,&
             B0, qs, ms, Ns, nus, stopcond, &
             pos0, dir0, w, dt0, dtmax, maxerr, maxsteps, minalt, root, tmax, &
             fixedstep, delDP, fngo3d, data, raytracer_stopconditions)
     elseif( modelnum == 6 ) then
        call raytracer_run( &
             pos,time,vprel,vgrel,n,&
             B0, qs, ms, Ns, nus, stopcond, &
             pos0, dir0, w, dt0, dtmax, maxerr, maxsteps, minalt, root, tmax, &
             fixedstep, delDP, fsimple, data, raytracer_stopconditions)
     end if
     ! Write the data to the output file
     do i=1,size(time,1),outputper
        write(outfile, &
             fmt='(i10, i10, 17es24.15e3, i10)', &
             advance='no'), &
             raynum, stopcond, &
             time(i), pos(:,i), vprel(:,i), vgrel(:,i), n(:,i), &
             B0(:,i), w, size(qs,1)
        do j=1,size(qs,1)
           write(outfile, fmt='(es24.15e3)',  advance='no'), qs(j,i)
        end do
        do j=1,size(qs,1)
           write(outfile, fmt='(es24.15e3)',  advance='no'), ms(j,i)
        end do
        do j=1,size(qs,1)
           write(outfile, fmt='(es24.15e3)',  advance='no'), Ns(j,i)
        end do
        do j=1,size(qs,1)
           write(outfile, fmt='(es24.15e3)',  advance='no'), nus(j,i)
        end do
        write(outfile, fmt='(a)'), ''
     end do
        
     deallocate(pos)
     deallocate(time)
     deallocate(vprel)
     deallocate(vgrel)
     deallocate(n)
     deallocate(B0)
     deallocate(qs)
     deallocate(ms)
     deallocate(Ns)
     deallocate(nus)

     raynum = raynum+1
  end do
  flush(outfile)
  flush(infile)
  
  close(unit=outfile)
  close(unit=infile)



end program raytracer_driver

