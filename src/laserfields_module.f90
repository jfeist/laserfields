! Copyright (c) 2012, Johannes Feist
! licensed under the MIT open source license, see LICENSE file

module laserfields_module
  use nrtype
  implicit none
  save

  real(dp), parameter :: gaussian_time_cutoff_fwhm = 3.5d0
  integer, parameter :: minimum_steps_per_laser_period = 75

  !> Datatype describing a single laserfield, should usually be created through make_laserfield routine
  type :: laserfield
     !> The shape of the envelope, can be 'gaussianF', 'gaussianI', 'sin2', 'sin_exp', 'linear', 'linear2', 'readin'

     !> - gaussianF, gaussianI: Gaussian pulses with form exp(-z t<sup>2</sup>) (see laserfield%duration doc for meaning of z)
     !> - sin2: sin<sup>2</sup> envelope
     !> - sin_exp: sin<sup>form_exponent</sup> envelope
     !> - linear: linear rampon, then constant amplitude, then linear rampoff
     !> - linear2: same as linear, but with sin<sup>2</sup> rampon/off
     !> - readin: read numerical data from an input file and interpolate
     character(len=30) :: form = ''
     !> Only for form=='readin': file to read numerical data from
     character(len=200) :: datafile = ''
     !> If .true., describe the vector potential A(t) with the given envelope - E(t) will contain derivative terms.
     !> This is often useful to make sure that A(-infty) = A(infty), a requirement for propagating fields.
     logical :: is_vecpot = .false.
     ! these are the input laserfields.in "convenient" units - W/cm^2, nm, as
     !> peak intensity in W/cm<sup>2</sup>
     real(dp) :: intensity_Wcm2 = 0.d0
     !> carrier wavelength in nm
     real(dp) :: lambda_nm = 0.d0
     !> peak time of envelope in as
     real(dp) :: peak_time_as = 0.d0
     !> duration in as - depending on the value of form, this has different meanings!

     !> - gaussianF: FWHM of the field envelope
     !> - gaussianI: FWHM of the intensity envelope
     !> - sin2, sin_exp: total duration of pulse
     !> - linear, linear2: time during which pulse has peak amplitude
     !> - readin: ignored
     real(dp) :: duration_as = 0.d0
     !> Only relevant for form linear and linear2: duration of rampon/rampoff at beginning and end of pulse, in as
     real(dp) :: rampon_as = 0.d0
     !> Carrier-envelope phase (CEP), in multiples of &pi;.<br> phase_pi=0 corresponds to a pulse with sin(w (t-t<sub>peak</sub>)) oscillation.
     real(dp) :: phase_pi = 0.d0
     !> Only relevant for form sin_exp: exponent of envelope
     real(dp) :: form_exponent = 1.d0
     !> linear temporal chirp rate, in units of omega_0/as

     !> &omega; = &omega;<sub>0</sub> (1 + linear_chirp_rate_w0 (t-t<sub>peak</sub>))
     !> e.g. for a value of 1.d-3, the frequency will change by &omega; over 1000 as.
     !> Since the carrier wave should not have negative frequencies, you must take care that the chirp is not too large.
     real(dp) :: linear_chirp_rate_w0as = 0.d0

     ! these are "derived" parameters that should not be entered directly. they are all in atomic units
     !> <b>Derived</b>: peak electric field strength in a.u.
     real(dp) :: E0 = 0.d0
     !> <b>Derived</b>: carrier angular frequency in a.u.
     real(dp) :: omega = 0.d0
     !> <b>Derived</b>: carrier period in a.u.
     real(dp) :: TX = 0.d0
     !> <b>Derived</b>: peak time of envelope in a.u.
     real(dp) :: peak_time = 0.d0
     !> <b>Derived</b>: duration in a.u.
     real(dp) :: duration = 0.d0
     !> <b>Derived</b>: rampon/rampoff in a.u.
     real(dp) :: rampon = 0.d0

     !> Numerical arrays to save the fields for interpolation for read-in fields.
     !> also used for numerical integration of E(t) to obtain A(t) and Z(t) for is_vecpot=.false.
     !> and numerical integration of A(t) to obtain Z(t) for is_vecpot=.true.
     real(dp), dimension(:), allocatable :: tt, ee, aa, zz
  end type laserfield

  !> Global array saving the laserfields read from parameter files and added with add_laserfield
  type(laserfield), dimension(100) :: all_laserfields
  !> Number of laserfields in the global array all_laserfields
  integer :: n_laserfields = 0

  !> Make a new type(laserfield), either from parameters or reading from a datafile
  interface make_laserfield
     module procedure make_laserfield_params
     module procedure make_laserfield_datafile
  end interface make_laserfield
  private :: make_laserfield_params, make_laserfield_datafile
  !> Return the electric field E(t). if called with a type(laserfield) argument,
  !> gets E(t) for just that one field, otherwise gets sum of all fields in all_laserfields
  interface get_EL
     module procedure laserfields_get_EL
     module procedure lf_get_EL
  end interface
  private :: laserfields_get_EL, lf_get_EL
  !> Return the vector potential A(t). if called with a type(laserfield) argument,
  !> gets A(t) for just that one field, otherwise gets sum of all fields in all_laserfields
  interface get_AL
     module procedure laserfields_get_AL
     module procedure lf_get_AL
  end interface
  private :: laserfields_get_AL, lf_get_AL
  !> Return the free-space displacement Z(t) of an electron. if called with a type(laserfield) argument,
  !> gets Z(t) for just that one field, otherwise gets sum of all fields in all_laserfields
  interface get_ZL
     module procedure laserfields_get_ZL
     module procedure lf_get_ZL
  end interface
  private :: laserfields_get_ZL, lf_get_ZL
  !> Return the fourier transform E(&omega;). if called with a type(laserfield) argument,
  !> gets E(&omega;) for just that one field, otherwise gets sum of all fields in all_laserfields
  interface get_EL_fourier_transform
     module procedure laserfields_get_EL_fourier_transform
     module procedure lf_get_EL_fourier_transform
  end interface
  private :: laserfields_get_EL_fourier_transform, lf_get_EL_fourier_transform
  !> Return the fourier transform A(&omega;). if called with a type(laserfield) argument,
  !> gets A(&omega;) for just that one field, otherwise gets sum of all fields in all_laserfields
  interface get_AL_fourier_transform
     module procedure laserfields_get_AL_fourier_transform
     module procedure lf_get_AL_fourier_transform
  end interface
  private :: laserfields_get_AL_fourier_transform, lf_get_AL_fourier_transform
  !> Return the fourier transform E(&omega;) as a string that can be used as a function in gnuplot.
  !> if called with a type(laserfield) argument,
  !> gets E(&omega;) for just that one field, otherwise gets sum of all fields in all_laserfields
  interface get_EL_fourier_transform_string
     module procedure laserfields_get_EL_fourier_transform_string
     module procedure lf_get_EL_fourier_transform_string
  end interface
  private :: laserfields_get_EL_fourier_transform_string, lf_get_EL_fourier_transform_string

contains
  !---------------------------------------------------------------------------
  !> add a type(laserfield) to the global list all_laserfields
  subroutine add_laserfield(lf)
    type(laserfield), intent(in) :: lf
    if (n_laserfields == size(all_laserfields)) STOP 'ERROR: added too many laser fields!'
    n_laserfields = n_laserfields + 1
    all_laserfields(n_laserfields) = lf
  end subroutine add_laserfield
  !---------------------------------------------------------------------------
  !> make a new type(laserfield). for the meaning of the parameter see the documentation of laserfield
  type(laserfield) function make_laserfield_params(form,intensity_Wcm2,lambda_nm,peak_time_as,&
       duration_as,rampon_as,form_exponent,phase_pi,is_vecpot,linear_chirp_rate_w0as) result(lf)
    character(len=*), intent(in) :: form
    real(dp), intent(in) :: intensity_Wcm2
    real(dp), intent(in) :: lambda_nm
    real(dp), intent(in) :: peak_time_as
    real(dp), intent(in) :: duration_as
    real(dp), intent(in), optional :: rampon_as, phase_pi, form_exponent, linear_chirp_rate_w0as
    logical, intent(in), optional :: is_vecpot

    lf%form = form
    lf%intensity_Wcm2 = intensity_Wcm2
    lf%lambda_nm = lambda_nm
    lf%peak_time_as = peak_time_as
    lf%duration_as = duration_as
    if (present(form_exponent)) lf%form_exponent = form_exponent
    if (present(rampon_as)) lf%rampon_as = rampon_as
    if (present(phase_pi)) lf%phase_pi = phase_pi
    if (present(is_vecpot)) lf%is_vecpot = is_vecpot
    if (present(linear_chirp_rate_w0as)) lf%linear_chirp_rate_w0as = linear_chirp_rate_w0as

    call laserfield_set_dependent(lf)
  end function make_laserfield_params
  !---------------------------------------------------------------------------
  type(laserfield) function make_laserfield_datafile(datafile,is_vecpot) result(lf)
    character(len=*), intent(in) :: datafile
    logical, intent(in), optional :: is_vecpot

    lf%form = 'readin'
    lf%datafile = datafile
    if (present(is_vecpot)) lf%is_vecpot = is_vecpot
    call laserfield_set_dependent(lf)
  end function make_laserfield_datafile
  !---------------------------------------------------------------------------
  subroutine laserfield_set_dependent(lf)
    use atomic_units
    type(laserfield), intent(inout) :: lf

    ! make sure we know how to deal with this laser field
    select case (lf%form)
    case('gaussianF','gaussianI','linear','linear2','sin2','sin_exp','readin')
       continue
    case default
       write(0,'(2A)') 'ERROR! laser field form unknown, form = ', trim(lf%form)
       STOP 516
    end select

    ! for laser fields that we read from file, all the parameters apart from is_vecpot
    ! are undefined at first, and we have to infer them from the file that we read in
    if (lf%form == 'readin') then
       call read_laserfield_from_file(lf)
    end if

    lf%E0 = sqrt(lf%intensity_Wcm2 * au_wcm2toEL2)
    lf%TX = lf%lambda_nm * au_nm / au_c
    if (lf%TX/=0) then
       lf%omega = TWOPI / lf%TX
    else ! special provision: if lambda=0, we take a laser field without any oscillation, so also omega=0
       lf%omega = 0
    end if
    lf%peak_time = lf%peak_time_as * au_as
    lf%duration = lf%duration_as * au_as
    lf%rampon = lf%rampon_as * au_as

    ! handle exceptions and warnings
    if (lf%is_vecpot) then
       select case (lf%form)
       case('linear','sin_exp')
          write(0,*) 'ERROR: Selected laser field cannot be given by A-envelope: ', trim(lf%form)
          write(0,*) '       -> E-field would be discontinuous!'
          STOP 311
       end select
    end if

    select case (lf%form)
    case('linear','linear2')
       if (lf%linear_chirp_rate_w0as /= 0.d0) then
          write(0,'(3A)') 'WARNING: be careful with chirped laser pulses with ', trim(lf%form), &
               & ' envelopes! Check that E(t) and A(t) are really as you expect them!'
       end if
    end select

    if (lf%omega==0 .and. lf%linear_chirp_rate_w0as/=0) then
       write(0,'(A)') 'WARNING: laserfield has lambda=0, i.e. no oscillation, but chirp/=0. chirp will be ignored!'
    end if

    if (lf%omega==0 .and. lf%is_vecpot) then
       write(0,'(2A)') 'WARNING: laserfield has lambda=0, i.e. no oscillation, but is_vecpot is true. ', &
            & 'The peak intensity is currently not treated correctly!'
    end if

  end subroutine laserfield_set_dependent
  !---------------------------------------------------------------------------
  subroutine read_laserfield_from_file(lf)
    use laserfields_miscfuncs
    use atomic_units
    type(laserfield), intent(inout) :: lf
    integer ::  ii, io_error, unit, npoints
    character(len=200) :: tmpstr
    real(dp), dimension(:), allocatable :: tmptt, tmpvals
    real(dp) :: dt, lastzerocrossing, starttime, endtime

    write(6,*) '# Reading laserfield from file: ', trim(lf%datafile)

    unit = get_unused_unit()
    open(unit,file=trim(lf%datafile),status='old',action='read',iostat=io_error)
    if (io_error /= 0) then
       write(0,*) 'ERROR opening file for laserfield, filename =', trim(lf%datafile)
       STOP 511
    end if

    npoints=0
    ! count number of valid lines in file
    do while (io_error==0) ! read lines until an IO error ocurs, i.e. (normally) EOF
       read(unit,'(a200)',iostat=io_error) tmpstr
       tmpstr = adjustl(tmpstr)
       ! at a comment, empty line or end of file, cycle the loop
       if (tmpstr(1:1)=='#' .or. tmpstr=='' .or. io_error/=0) cycle
       ! otherwise, the line has data
       npoints=npoints+1
    end do

    if (npoints == 0) then
       write(0,*) 'ERROR: No data points found in file for laserfield. Stop.'
       STOP 512
    end if

    write(6,*) '# Number of data points found:', npoints
    allocate(tmptt(npoints), tmpvals(npoints))

    ! now read in
    ii = 0
    rewind(unit, iostat=io_error)
    do while (io_error==0)
       read(unit,'(a200)',iostat=io_error) tmpstr
       tmpstr = adjustl(tmpstr)
       if (tmpstr(1:1)=='#' .or. tmpstr=='' .or. io_error/=0) cycle

       ii = ii + 1
       read(tmpstr,*) tmptt(ii), tmpvals(ii)
    end do
    close(unit)

    !**************
    ! analyze the data we have read to guess some information about the field
    starttime = tmptt(1)
    endtime = tmptt(npoints)

    lf%TX = huge(1.d0)
    dt = huge(1.d0)
    lastzerocrossing = starttime
    do ii = 2, npoints
       if (tmptt(ii)<=tmptt(ii-1)) then
          write(0,*) 'ERROR: times in data file for laser field must be monotonously increasing! filename =', trim(lf%datafile)
          STOP 514
       end if

       dt = min(dt,tmptt(ii)-tmptt(ii-1))
       if (sign(1.d0,tmpvals(ii))/=sign(1.d0,tmpvals(ii-1))) then
          ! we have crossed a zero, estimate TX with this
          lf%TX = min(lf%TX,2*(tmptt(ii)-lastzerocrossing))
          lastzerocrossing = tmptt(ii)
       end if
    end do

    ! ensure small enough dt to allow for good interpolation and numeric integration/differentiation
    ! note that TX will usually be much smaller than the "real" TX already, so this should
    ! be small enough by far
    dt = min(dt,lf%TX/200)

    !**************
    ! transform the field with arbitrary times to a grid with equal spacing for the time steps
    ! watch out, we are reusing npoints
    npoints = nint((endtime - starttime) / dt) + 1  ! the "+ 1" is for the endpoint
    ! this is the actual dt we use to have exactly npoints from starttime to endtime
    dt = (endtime - starttime) / (npoints-1)
    allocate(lf%tt(npoints), lf%EE(npoints), lf%AA(npoints), lf%ZZ(npoints))
    do ii = 1, npoints
       lf%tt(ii) = starttime + (ii-1) * dt
    end do

    ! Interpolate to get laser field on regularly spaced points
    if (lf%is_vecpot) then
       lf%AA(:) = interpolate(tmptt,tmpvals,lf%tt,degree=6)
       call laserfield_numerical_derivation_E_from_A(lf)
    else
       lf%EE(:) = interpolate(tmptt,tmpvals,lf%tt,degree=6)
       call laserfield_numerical_integration_A_from_E(lf)
    end if
    call laserfield_numerical_integration_Z_from_A(lf)

    ! set guessed parameters
    lf%E0 = maxval(abs(lf%EE))
    lf%intensity_Wcm2 = lf%E0**2 / au_Wcm2toEL2
    lf%omega = TWOPI / lf%TX
    lf%lambda_nm = lf%TX / au_nm * au_c
    lf%peak_time = (starttime + endtime) / 2
    lf%peak_time_as = lf%peak_time / au_as
    lf%duration = endtime - starttime
    lf%duration_as = lf%duration / au_as

  end subroutine read_laserfield_from_file
  !---------------------------------------------------------------------------
  subroutine laserfield_numerical_derivation_E_from_A(lf)
    type(laserfield), intent(inout) :: lf
    real(dp) :: dt, zeit
    integer :: ii

    ! we already have A(t) for the laser field and want to calculate lf%EE = -dA/dt

    ! use a smaller dt than the grid on which A(t) is given
    dt = (lf%tt(2) - lf%tt(1))/50

    do ii = 1, size(lf%tt)
       zeit = lf%tt(ii)
       ! do not forget the minus sign!
       lf%EE(ii) = -(get_AL(lf,zeit+dt) - get_AL(lf,zeit-dt)) / (2*dt)
    end do

  end subroutine laserfield_numerical_derivation_E_from_A
  !---------------------------------------------------------------------------
  subroutine laserfield_numerical_integration_A_from_E(lf)
    type(laserfield), intent(inout) :: lf
    integer :: ii, jj, nsteps, simpson_fac
    real(dp) :: EL, zeit
    ! we have a laser field where we can get the electric field and
    ! want to numerically integrate that to get the vector potential A(t) = -Int E(t) dt
    ! we use the composite Simpson rule

    ! nsteps HAS to be ODD for Simpson rule, so that the number of intervals is even
    nsteps = 21

    ! initialize first value to zero
    lf%AA(1) = 0.d0
    do ii = 2, size(lf%AA)
       ! for each step we want to take, we do nsteps integration steps with the trapezoid rule
       zeit = lf%tt(ii-1)
       EL = get_EL(lf,zeit)

       ! do integration from tt(ii-1) to tt(ii)
       ! Int_a^b f(t) dt ~ (f(a) + 4*f(a+h) + 2*f(a+2*h) + 4*f(a+3*h) + ... + 4*f(b-h) + f(b)) * h/3
       ! where h = (b-a)/(nsteps-1)
       lf%AA(ii) = -EL
       simpson_fac = 4
       do jj = 2, nsteps-1
          zeit = (lf%tt(ii-1) * (nsteps-jj) + lf%tt(ii) * (jj-1)) / (nsteps-1)
          ! jj = 1      -> zeit = lf%tt(ii-1) * (nsteps-1) / (nsteps-1) -> lf%tt(ii-1) -> ok
          ! jj = nsteps -> zeit = lf%tt(ii)   * (nsteps-1) / (nsteps-1) -> lf%tt(ii)   -> ok
          EL = get_EL(lf,zeit)
          ! do not forget the minus sign!
          lf%AA(ii) = lf%AA(ii) - simpson_fac * EL
          ! 6 - 2 = 4,  6 - 4 = 2 -> alternating between 4 and 2
          simpson_fac = 6 - simpson_fac
       end do
       ! last step
       zeit = lf%tt(ii)
       EL = get_EL(lf,zeit)
       lf%AA(ii) = lf%AA(ii) - EL
       ! int = sum(...) * h/3
       lf%AA(ii) = lf%AA(ii) * (lf%tt(ii)-lf%tt(ii-1))/(3*(nsteps-1))
       ! A(t+dt) = A(t) - Int_t^{t+dt} E(t) dt
       lf%AA(ii) = lf%AA(ii) + lf%AA(ii-1)
       !write(60,'(99g25.12)') lf%tt(ii), lf%AA(ii)
    end do
  end subroutine laserfield_numerical_integration_A_from_E
  !---------------------------------------------------------------------------
  subroutine laserfield_numerical_integration_Z_from_A(lf)
    type(laserfield), intent(inout) :: lf
    integer :: ii, jj, nsteps, simpson_fac
    real(dp) :: AL, zeit
    ! we have a laser field where we can get the vector potential
    ! want to numerically integrate that to get the classical position Z(t) = -Int A(t) dt
    ! we use the composite Simpson rule

    ! nsteps HAS to be ODD for Simpson rule, so that the number of intervals is even
    nsteps = 21

    ! initialize first value to zero
    lf%ZZ(1) = 0.d0
    do ii = 2, size(lf%ZZ)
       ! for each step we want to take, we do nsteps integration steps with the trapezoid rule
       zeit = lf%tt(ii-1)
       AL = get_AL(lf,zeit)

       ! do integration from tt(ii-1) to tt(ii)
       ! Int_a^b f(t) dt ~ (f(a) + 4*f(a+h) + 2*f(a+2*h) + 4*f(a+3*h) + ... + 4*f(b-h) + f(b)) * h/3
       ! where h = (b-a)/(nsteps-1)
       lf%ZZ(ii) = -AL
       simpson_fac = 4
       do jj = 2, nsteps-1
          zeit = (lf%tt(ii-1) * (nsteps-jj) + lf%tt(ii) * (jj-1)) / (nsteps-1)
          ! jj = 1      -> zeit = lf%tt(ii-1) * (nsteps-1) / (nsteps-1) -> lf%tt(ii-1) -> ok
          ! jj = nsteps -> zeit = lf%tt(ii)   * (nsteps-1) / (nsteps-1) -> lf%tt(ii)   -> ok
          AL = get_AL(lf,zeit)
          ! do not forget the minus sign!
          lf%ZZ(ii) = lf%ZZ(ii) - simpson_fac * AL
          ! 6 - 2 = 4,  6 - 4 = 2 -> alternating between 4 and 2
          simpson_fac = 6 - simpson_fac
       end do
       ! last step
       zeit = lf%tt(ii)
       AL = get_AL(lf,zeit)
       lf%ZZ(ii) = lf%ZZ(ii) - AL
       ! int = sum(...) * h/3
       lf%ZZ(ii) = lf%ZZ(ii) * (lf%tt(ii)-lf%tt(ii-1))/(3*(nsteps-1))
       ! Z(t+dt) = Z(t) - Int_t^{t+dt} A(t) dt
       lf%ZZ(ii) = lf%ZZ(ii) + lf%ZZ(ii-1)
       !write(60,'(99g25.12)') lf%tt(ii), lf%ZZ(ii)
    end do
  end subroutine laserfield_numerical_integration_Z_from_A
  !---------------------------------------------------------------------------
  subroutine lf_setup_AA_interpolation(lf)
    type(laserfield), intent(inout) :: lf
    integer :: ii, npoints
    real(dp) :: starttime, endtime, dt
    ! we have a laser field where we can get the electric field and
    ! want to numerically integrate that to get the vector potential
    ! we use the trapezoid rule

    starttime = lf_get_starttime(lf)
    endtime   = lf_get_endtime  (lf)
    ! ensure small enough dt to allow for good interpolation
    dt = lf%TX/250
    if (dt == 0.d0) then
       ! this is a field without oscillation - interpolate envelope with 500 points
       dt = (endtime-starttime) / 500
    end if

    ! find number of points closest to wanted dt
    npoints = nint((endtime-starttime) / dt + 1)  ! the "+ 1" is for the endpoint
    allocate(lf%tt(npoints), lf%AA(npoints))
    do ii = 1, npoints
       ! set the times through linear interpolation from starttime to endtime, without reference to dt
       ! this avoids numerical addition problems if starttime is very large and dt is small
       ! the actual time step is not exactly equal to dt, as "(endtime-starttime) / dt + 1" is not precisely integer
       lf%tt(ii) = (starttime * (npoints-ii) + endtime * (ii-1)) / (npoints-1)
    end do

    call laserfield_numerical_integration_A_from_E(lf)

  end subroutine lf_setup_AA_interpolation
  !---------------------------------------------------------------------------
  subroutine lf_setup_ZZ_interpolation(lf)
    type(laserfield), intent(inout) :: lf
    integer :: ii, npoints
    real(dp) :: starttime, endtime, dt
    ! we have a laser field where we can get the vector potential and
    ! want to numerically integrate that to get the classical position for the acceleration gauge / Kramers frame hamiltonian
    ! we use the trapezoid rule

    starttime = lf_get_starttime(lf)
    endtime   = lf_get_endtime  (lf)
    ! ensure small enough dt to allow for good interpolation
    dt = lf%TX/250
    if (dt == 0.d0) then
       ! this is a field without oscillation - interpolate envelope with 500 points
       dt = (endtime-starttime) / 500
    end if

    ! find number of points closest to wanted dt
    npoints = nint((endtime-starttime) / dt + 1)  ! the "+ 1" is for the endpoint
    if (.not.allocated(lf%tt)) then
       allocate(lf%tt(npoints))
    else if (size(lf%tt)/=npoints) then
       write(0,*) 'ERROR: size(lf%tt)/=npoints in lf_setup_ZZ_interpolation!'
       STOP 356
    end if

    allocate(lf%ZZ(npoints))
    do ii = 1, npoints
       ! set the times through linear interpolation from starttime to endtime, without reference to dt
       ! this avoids numerical addition problems if starttime is very large and dt is small
       ! the actual time step is not exactly equal to dt, as "(endtime-starttime) / dt + 1" is not precisely integer
       lf%tt(ii) = (starttime * (npoints-ii) + endtime * (ii-1)) / (npoints-1)
    end do

    call laserfield_numerical_integration_Z_from_A(lf)

  end subroutine lf_setup_ZZ_interpolation
  !---------------------------------------------------------------------------
  ! this calculates the envelope of the laser field
  ! as well as the first derivative (by time) of the envelope
  ! the derivative is needed for calculating E(t) = -dA/dt if the laserfield describes A(t) (lf%is_vecpot is .true.)
  subroutine lf_get_envelope(lf,zeit,env,envpr)
    type(laserfield), intent(in) :: lf
    real(dp), intent(in) :: zeit
    real(dp), intent(out) :: env, envpr
    real(dp) :: trel, ttmp

    trel = zeit - lf%peak_time

    select case (lf%form)
    case('gaussianF') ! laser pulse with Gaussian envelope, lf%duration is FWHM of electric field envelope
       env   = exp(-trel**2*log(16.d0)/lf%duration**2)
       envpr = -2*trel*log(16.d0)/lf%duration**2 * env
    case('gaussianI') ! laser pulse with Gaussian envelope, lf%duration is FWHM of intensity envelope
       env   = exp(-trel**2*log( 4.d0)/lf%duration**2)
       envpr = -2*trel*log( 4.d0)/lf%duration**2 * env
    case('linear')
       ! "linear" is a field with constant intensity for lf%duration, and linear rampon/rampoff of duration lf%rampon at the beginning and end
       if (abs(trel) < lf%duration/2) then
          env   = 1.d0
          envpr = 0.d0
       else if (abs(trel) < lf%duration/2 + lf%rampon) then
          env   = 1.d0 - (abs(trel)-lf%duration/2) / lf%rampon
          envpr = -sign(1.d0,trel) / lf%rampon
       else
          env   = 0.d0
          envpr = 0.d0
       end if
    case('linear2')
       ! same as linear, but with sin2 rampon/rampoff
       if (abs(trel) < lf%duration/2) then
          env   = 1.d0
          envpr = 0.d0
       else if (abs(trel) < lf%duration/2 + lf%rampon) then
          ttmp  = 1.d0 - (abs(trel)-lf%duration/2) / lf%rampon
          env   = sin(PIO2*ttmp)**2
          envpr = -sign(1.d0,trel) * sin(PI*trel) * PIO2/lf%rampon
       else
          env   = 0.d0
          envpr = 0.d0
       end if
    case ('sin2')
       if (abs(trel) < lf%duration/2) then
          env   =  cos(PI*trel/lf%duration)**2
          !        2 sin(wt) cos(wt) = sin(2wt)
          envpr = -sin(TWOPI*trel/lf%duration) * PI/lf%duration
       else
          env   = 0.d0
          envpr = 0.d0
       end if
    case ('sin_exp')
       if (abs(trel) < lf%duration/2) then
          env   =  cos(PI*trel/lf%duration)**lf%form_exponent
          envpr = -sin(PI*trel/lf%duration) * lf%form_exponent * &
               &   cos(PI*trel/lf%duration)**(lf%form_exponent-1) * PI/lf%duration
       else
          env   = 0.d0
          envpr = 0.d0
       end if
    case default
       write(0,'(2A)') 'ERROR! laser field form unknown, form = ', trim(lf%form)
       STOP 516
    end select

    env   = lf%E0 * env
    envpr = lf%E0 * envpr

  end subroutine lf_get_envelope
  !---------------------------------------------------------------------------
  complex(dpc) function lf_envelope_fourier(lf,omega) result(val)
    use atomic_units
    ! return the fourier transform of the envelope of the laser field
    ! we write the whole pulse as
    ! f(t) = (env(t) exp(IU*(phi0 + w0*tp + chirp*tp**2)) + c.c. ) / (2*IU), where tp = t-tpeak
    ! for the fourier transform we include the chirp term exp(i chirp (t-tpeak)**2) in the envelope.
    ! the fourier transform of the envelope is then a complex function.
    ! however, for unchirped pulses, the result will be purely real!
    type(laserfield), intent(in) :: lf
    real(dp), intent(in) :: omega
    real(dp) :: chirp
    complex(dpc) :: z

    ! important - use omega, not lf%omega as the argument here!!

    ! instantaneous frequency:
    ! w(t) = lf%omega * (1.d0 + lf%linear_chirp_rate_w0as/au_as * (zeit-lf%peak_time))
    chirp =  lf%omega * lf%linear_chirp_rate_w0as / au_as

    ! for the various calculations, see chirped_fourier.nb in the mathematica subversion directory
    select case (lf%form)
    case('gaussianF')
       ! F[exp(-z*t**2)] = exp(-w**2/4*z)/sqrt(2*z) (for real(z)>0)
       z = log(16.d0)/lf%duration**2 - IU * chirp
       val = exp(-omega**2/(4*z)) / sqrt(2*z)
    case('gaussianI')
       z = log( 4.d0)/lf%duration**2 - IU * chirp
       val = exp(-omega**2/(4*z)) / sqrt(2*z)
    case('linear')
       if (chirp /= 0.d0) then
          write(0,'(A)') 'ERROR! fourier transform of "linear" envelope with chirp not implemented!'
          STOP 676
       end if
       val = sqrt(8.d0/PI) * sin(omega*lf%rampon/2) * sin(omega*(lf%rampon+lf%duration)/2) / (lf%rampon * omega**2)
    case('linear2')
       if (chirp /= 0.d0) then
          write(0,'(A)') 'ERROR! linear2 fourier transform with chirp not implemented!'
          STOP 676
       end if
       val = sqrt(2*PI**3) * cos(omega*lf%rampon/2) * sin(omega*(lf%rampon+lf%duration)/2) / (PI**2*omega - lf%rampon**2*omega**3)
    case ('sin2')
       if (chirp == 0.d0) then ! the expression with chirp can not be evaluated with chirp == 0, so we take this as a special case
          val = sqrt(8.d0*PI**3) * sin(omega*lf%duration/2)/(TWOPI**2*omega - omega**3*lf%duration**2)
       else
          ! now we use that cos(pi*t/T)**2 * exp(i*c*t**2) can be written as 0.5 exp(i*c*t**2) + 0.25 exp(i*c*t**2 - 2*i*pi*t/T) + 0.25 exp(i*c*t**2 + 2*i*pi*t/T)
          ! the integral of exp(IU*(a*t+b*t**2)) from t=-T/2 to t=T/2 can be calculated analytically and is implemented in the function below
          ! the arguments are a={-omega, -2*pi/T-omega, 2*pi/T-omega} and b=chirp
          val =    expiatbt2_intT(                   - omega, chirp, lf%duration)/2  &
               & + expiatbt2_intT(-TWOPI/lf%duration - omega, chirp, lf%duration)/4  &
               & + expiatbt2_intT( TWOPI/lf%duration - omega, chirp, lf%duration)/4
       end if
    case ('sin_exp')
       STOP 'ERROR! sin_exp fourier transform too complicated to implement'
    case default
       STOP 'ERROR! unknown laser field form'
    end select

    val = lf%E0 * val

  end function lf_envelope_fourier
  !---------------------------------------------------------------------------
  ! returns the result of the integral Int(exp(IU*(a*t+b*t**2)),{t,-T/2,T/2}) / sqrt(TWOPI)
  complex(dpc) function expiatbt2_intT(a,b,T) result(res)
    use faddeeva, ONLY : erf
    real(dp), intent(in) :: a, b, T
    res = erf((a-b*T)/sqrt(4*IU*b)) - erf((a+b*T)/sqrt(4*IU*b))
    res = res * (-IU/sqrt(8*IU*b)) * exp(-IU*a**2/(4*b))
  end function expiatbt2_intT
  !---------------------------------------------------------------------------
  character(1000) function lf_envelope_fourier_string(lf,omegastr) result(val)
    use atomic_units
    use laserfields_miscfuncs
    ! return the fourier transform of the envelope of the laser field
    ! we write the whole pulse as
    ! f(t) = (env(t) exp(IU*(phi0 + w0*tp + chirp*tp**2)) + c.c. ) / (2*IU), where tp = t-tpeak
    ! for the fourier transform of the envelope, we include the chirp term
    ! exp(i chirp (t-tpeak)**2) in the envelope.
    ! the fourier transform is then a complex function.
    ! however, for unchirped pulses, the result will be purely real!
    type(laserfield), intent(in) :: lf
    character(*), intent(in) :: omegastr
    character(1000) :: tmpstr
    real(dp) :: chirp
    complex(dpc) :: z

    ! important - use omega, not lf%omega as the argument here!!

    ! instantaneous frequency:
    ! w(t) = lf%omega * (1.d0 + lf%linear_chirp_rate_w0as/au_as * (zeit-lf%peak_time))
    chirp =  lf%omega * lf%linear_chirp_rate_w0as / au_as

    ! for the various calculations, see chirped_fourier.nb in the mathematica subversion directory

    select case (lf%form)
    case('gaussianF')
       ! F[exp(-z*t**2)] = exp(-w**2/4*z)/sqrt(2*z) (for real(z)>0)
       z = log(16.d0)/lf%duration**2 - IU * chirp
       write(tmpstr,'(999A)') gnuplotstring(1/sqrt(2*z)), ' * exp(', gnuplotstring(-1/(4*z)),' * ('//omegastr//')**2)'
    case('gaussianI')
       z = log( 4.d0)/lf%duration**2 - IU * chirp
       write(tmpstr,'(999A)') gnuplotstring(1/sqrt(2*z)), ' * exp(', gnuplotstring(-1/(4*z)),' * ('//omegastr//')**2)'
    case('linear')
       if (chirp /= 0.d0) then
          write(0,'(A)') 'ERROR! fourier transform of "linear" envelope with chirp not implemented!'
          STOP 676
       end if
       write(tmpstr,'(SP,999(ES15.8,A))') 8.d0/(PI*lf%rampon), ' * sin(', lf%rampon/2,'*('//omegastr//')) * sin(', &
            & (lf%rampon+lf%duration)/2,'*('//omegastr//')) / ('//omegastr//')**2'
    case('linear2')
       if (chirp /= 0.d0) then
          write(0,'(A)') 'ERROR! linear2 fourier transform with chirp not implemented!'
          STOP 676
       end if
       write(tmpstr,'(SP,999(ES15.8,A))') sqrt(2*PI**3), ' * cos(', lf%rampon/2,'*('//omegastr//')) * sin(', &
            & (lf%rampon+lf%duration)/2,'*('//omegastr//')) / (',PI**2,'*('//omegastr//') - ',lf%rampon**2,'*('//omegastr//')**3)'
    case ('sin2')
       if (chirp == 0.d0) then ! the expression with chirp can not be evaluated with chirp == 0, so we take this as a special case
          write(tmpstr,'(SP,999(ES15.8,A))') sqrt(8*PI**3), ' * sin(',lf%duration/2,'*('//omegastr//')) / (', &
               & TWOPI**2,'*('//omegastr//') - ',lf%duration**2,'*('//omegastr//')**3)'
       else
          write(0,'(A)') 'ERROR! sin2 fourier transform function string with chirp not implemented!'
          STOP 676
       end if
    case ('sin_exp')
       write(0,'(2A)') 'ERROR! sin_exp fourier transform too complicated to implement'
       STOP 517
    case default
       write(0,'(2A)') 'ERROR! laser field form unknown, form = ', trim(lf%form)
       STOP 518
    end select

    ! multiply with sqrt(2.) to get normalisation such that integral
    ! over omega from 0 to infinity gives the total power in the laser field
    ! (same as integral over t from -infinity to infinity)
    ! as the field E(t) is real, we know that E(omega) is symmetric about 0
    write(val,'(SP,A,ES15.8,3A)') '(',lf%E0,' * ',trim(tmpstr),')'

  end function lf_envelope_fourier_string
  !---------------------------------------------------------------------------
  subroutine lf_get_omega(lf,zeit,omega,env)
    ! Calculate current frequency (which changes with time when a chirp rate is specified)
    ! returns with an error if omega is 0 or smaller
    ! optional argument env can be given to suppress this error when env (laser field envelope) is smaller than epsilon
    use atomic_units
    type(laserfield), intent(in) :: lf
    real(dp), intent(in) :: zeit
    real(dp), intent(out) :: omega
    real(dp), intent(inout), optional :: env

    ! The chirp rate is given in omega_0/as
    omega = lf%omega * (1.d0 + lf%linear_chirp_rate_w0as/au_as * (zeit-lf%peak_time))
    ! if omega is larger than or equal to zero, everything is okay and we can return
    if (omega >= 0.d0) return

    ! omega < 0.d0 -> check envelope if present
    if (present(env)) then
       if (env <= epsilon(1.d0)) then
          ! envelope <= epsilon, recover by setting envelope to ZERO
          env = 0.d0
          ! now we can return. we do not write error message to prevent flooding e.g. fourierlaserfield output
          return
       end if
    end if

    ! omega < 0.d0 and no envelope or envelope too large
    write(0,'(a)') ' ERROR: Field chirp is too large, omega(t) <= zero!'
    write(0,'(a,9g15.5)') ' t, omega = ', zeit, omega
    if (present(env)) write(0,'(a,9g15.5)') ' envelope > epsilon, env, epsilon =', env, epsilon(1.d0)
    write(0,*) 'STOPPING!'
    STOP 519

  end subroutine lf_get_omega
  !---------------------------------------------------------------------------
  real(dp) function lf_get_EL(lf,zeit,env_out) result(EL)
    use laserfields_miscfuncs
    type(laserfield), intent(in) :: lf
    real(dp), intent(in) :: zeit
    real(dp), intent(out), optional :: env_out
    real(dp) :: env, envpr, osc, oscpr
    real(dp) :: omega

    if (lf%form == 'readin') then
       if (zeit < lf%tt(1) .or. zeit > lf%tt(size(lf%tt))) then
          EL = 0.d0
       else
          ! interpolate a value and cycle to next laser field
          EL = interpolate(lf%tt,lf%EE,zeit,degree=6)
       end if
       ! set the envelope part to zero, as the readin files can not make this distinction
       env = 0.d0
    else
       call lf_get_envelope(lf,zeit,env,envpr)

       if (lf%omega==0) then ! No oscillation enclosed by envelope
          if (lf%is_vecpot) then
             EL = -envpr
          else
             EL = env
          end if
       else
          call lf_get_omega(lf,zeit,omega,env)
          osc   = sin(omega * (zeit-lf%peak_time) + PI*lf%phase_pi)
          ! d(w(t)*(t-peak))/dt = w(t) + w'(t)*(t-peak) = omega + lf%omega * lf%chirp * (t-peak) = 2 * omega - lf%omega
          oscpr = (2.d0*omega-lf%omega)*cos(omega*(zeit-lf%peak_time)+PI*lf%phase_pi)

          if (lf%is_vecpot) then
             ! Divide out derivative of oscillation to ensure peak amplitude of E0
             ! NOTE: this does not actually give a peak amplitude of E0 for chirped linear/linear2 pulses
             EL = -(env * oscpr + envpr * osc) / lf%omega
          else ! describes electric field directly
             EL = env * osc
          end if
       end if
    end if

    if (present(env_out)) env_out = env

  end function lf_get_EL
  !---------------------------------------------------------------------------
  real(dp) function lf_get_AL(lf,zeit,env_out) result(AL)
    use laserfields_miscfuncs
    type(laserfield), intent(inout) :: lf
    real(dp), intent(in) :: zeit
    real(dp), intent(out), optional :: env_out
    real(dp) :: env, envpr, osc
    real(dp) :: omega

    if (lf%form == 'readin' .or. (.not.lf%is_vecpot)) then
       ! setup lf%AA/lf%tt if it does not already exist
       if (.not.allocated(lf%AA)) call lf_setup_AA_interpolation(lf)

       if      (zeit < lf%tt(1))           then
          AL = 0.d0
       else if (zeit > lf%tt(size(lf%tt))) then
          ! AL (should be zero after pulse for propagating pulse!) does not change after end of pulse
          AL = lf%AA(size(lf%AA))
       else
          ! interpolate a value from the numerical array lf%AA
          AL = interpolate(lf%tt,lf%AA,zeit,degree=6)
       end if

       ! set the envelope part to zero, as the readin files can not make this distinction
       env = 0.d0
    else
       call lf_get_envelope(lf,zeit,env,envpr)

       if (lf%omega==0) then ! No oscillation enclosed by envelope
          AL = env
       else
          call lf_get_omega(lf,zeit,omega,env)
          osc = sin(omega * (zeit-lf%peak_time) + PI*lf%phase_pi)
          ! Divide out derivative of oscillation to ensure peak amplitude of E0 for electric field
          AL = env*osc / lf%omega
       end if
    end if
    if (present(env_out)) env_out = env
  end function lf_get_AL
  !---------------------------------------------------------------------------
  real(dp) function lf_get_ZL(lf,zeit) result(ZL)
    use laserfields_miscfuncs
    type(laserfield), intent(inout) :: lf
    real(dp), intent(in) :: zeit
    integer :: np
    if (.not.allocated(lf%ZZ)) call lf_setup_ZZ_interpolation(lf)
    np = size(lf%tt)
    if      (zeit < lf%tt(1))  then
       ZL = 0.d0
    else if (zeit > lf%tt(np)) then
       ! ZL changes linearly after end of pulse (usually AL(endtime) should be zero, so ZL stays constant)
       ZL = lf%ZZ(np) - (zeit - lf%tt(np)) * get_AL(lf,lf%tt(np))
    else
       ! interpolate a value from the numerical array lf%ZZ
       ZL = interpolate(lf%tt,lf%ZZ,zeit,degree=6)
    end if
  end function lf_get_ZL
  !---------------------------------------------------------------------------
  real(dp) function laserfields_get_EL(zeit,envarr,partarr) result(EL)
    real(dp), intent(in) :: zeit
    real(dp), dimension(n_laserfields), intent(out), optional :: envarr, partarr
    integer :: i_field
    real(dp) :: env, part

    EL = 0.d0
    if (present(envarr)) envarr = 0.d0
    if (present(partarr)) partarr = 0.d0

    do i_field = 1, n_laserfields
       part = get_EL(all_laserfields(i_field),zeit,env)
       EL = EL + part
       if (present(envarr)) envarr(i_field) = env
       if (present(partarr)) partarr(i_field) = part
    end do

  end function laserfields_get_EL
  !---------------------------------------------------------------------------
  real(dp) function laserfields_get_AL(zeit,envarr,partarr) result(AL)
    real(dp), intent(in) :: zeit
    real(dp), dimension(n_laserfields), intent(out), optional :: envarr, partarr
    integer :: i_field
    real(dp) :: env, part

    AL = 0.d0
    do i_field = 1, n_laserfields
       part = get_AL(all_laserfields(i_field),zeit,env)
       AL = AL + part
       if (present(envarr)) envarr(i_field) = env
       if (present(partarr)) partarr(i_field) = part
    end do
  end function laserfields_get_AL
  !---------------------------------------------------------------------------
  real(dp) function laserfields_get_ZL(zeit) result(ZL)
    real(dp), intent(in) :: zeit
    integer :: i_field
    ZL = 0.d0
    do i_field = 1, n_laserfields
       ZL = ZL + get_ZL(all_laserfields(i_field),zeit)
    end do
  end function laserfields_get_ZL
  !---------------------------------------------------------------------------
  real(dp) function lf_get_halftotaltime(lf) result(tt)
    type(laserfield), intent(in) :: lf
    select case (lf%form)
    case('gaussianF')
       ! for gaussian, take peak_time + fwhm*gaussian_time_cutoff_fwhm as cutoff
       tt = lf%duration * gaussian_time_cutoff_fwhm
    case('gaussianI')
       ! take into account that for gaussianI, duration is FWHM of intensity
       tt = lf%duration * gaussian_time_cutoff_fwhm * sqrt(2.d0)
    case('linear','linear2')
       tt = lf%duration * 0.5d0 + lf%rampon
    case ('sin2','sin_exp')
       tt = lf%duration * 0.5d0
    case default
       write(0,*) 'ERROR: laser field form unknown, form = ',lf%form
       STOP 521
    end select
  end function lf_get_halftotaltime
  !---------------------------------------------------------------------------
  real(dp) function lf_get_starttime(lf) result(tt)
    type(laserfield), intent(in) :: lf
    if (lf%form=='readin') then
       tt = lf%tt(1)
    else
       tt = lf%peak_time - lf_get_halftotaltime(lf)
    end if
  end function lf_get_starttime
  !---------------------------------------------------------------------------
  real(dp) function lf_get_endtime(lf) result(tt)
    type(laserfield), intent(in) :: lf
    if (lf%form=='readin') then
       tt = lf%tt(size(lf%tt))
    else
       tt = lf%peak_time + lf_get_halftotaltime(lf)
    end if
  end function lf_get_endtime
  !---------------------------------------------------------------------------
  real(dp) function laserfields_starttime() result(tt)
    integer :: ii
    tt = huge(1.d0)
    do ii = 1, n_laserfields
       tt = min(tt, lf_get_starttime(all_laserfields(ii)))
    end do
  end function laserfields_starttime
  !---------------------------------------------------------------------------
  real(dp) function laserfields_endtime() result(tt)
    integer :: ii
    tt = -huge(1.d0)
    do ii = 1, n_laserfields
       tt = max(tt, lf_get_endtime(all_laserfields(ii)))
    end do
  end function laserfields_endtime
  !---------------------------------------------------------------------------
  function laserfields_smallest_TX() result(TX)
    real(dp) :: TX
    integer :: ii
    TX = huge(1.d0)
    do ii = 1, n_laserfields
       if (all_laserfields(ii)%TX < TX) TX = all_laserfields(ii)%TX
    end do
  end function laserfields_smallest_TX
  !---------------------------------------------------------------------------
  !> This function gives the maximum timestep that resolves the current oscillation well.
  !> the global parameter minimum_steps_per_laser_period determines how large the timestep here is allowed to be.
  real(dp) function laserfields_largest_possible_dt(zeit,endzeit,global) result(dt)
    real(dp), intent(in) :: zeit
    real(dp), intent(in), optional :: endzeit
    real(dp) :: tstart, tend, omega, TX
    integer :: ii
    logical, intent(in), optional :: global

    ! if endzeit is specified, that is the time until which we want to propagate at most
    if (present(endzeit)) then
       if (endzeit <= zeit) then
          ! we have already passed endzeit, but still want to do an extra step - shouldn't normally happen!
          ! do at least a very small step - the smallest possible that still changes zeit
          dt = nearest(zeit,1.d0)-zeit
          return
       end if
       ! normally, have the time until endzeit be the maximum timestep
       dt = endzeit-zeit
    else
       ! otherwise, do steps of at most 100 atomic units
       dt = 100.d0
    end if

    ! if global is requested, we just get the global maximum timestep
    if (present(global)) then
       if (global) then
          dt = min(dt,minval(all_laserfields(1:n_laserfields)%TX)/minimum_steps_per_laser_period)
          return ! exit function
       end if
    end if

    ! otherwise, loop over laser fields and get the maximum timestep allowed at this time (doing at most minimum_steps_per_laser_period steps per period)
    do ii = 1, n_laserfields
       tstart = lf_get_starttime(all_laserfields(ii))
       tend   = lf_get_endtime  (all_laserfields(ii))

       ! this laser field is only relevant if we are before its endtime
       if (zeit < tend) then
          if (zeit >= tstart) then
             if (all_laserfields(ii)%omega==0) then
                ! for fields without oscillation, we want 500 steps for the whole field
                dt = min(dt,(all_laserfields(ii)%duration+2*all_laserfields(ii)%rampon)/500)
             else
                ! if it's currently active, dt has to be at most TX(zeit)/minimum_steps_per_laser_period, where TX depends on zeit for chirped pulses
                call lf_get_omega(all_laserfields(ii),zeit,omega)
                TX = TWOPI / omega
                dt = min(dt,TX/minimum_steps_per_laser_period)
             end if
          else if (zeit+dt > tstart) then
             ! if we would go beyond the starttime, decrease dt to go only there
             dt = tstart - zeit
          end if
       end if
    end do

  end function laserfields_largest_possible_dt
  !---------------------------------------------------------------------------
  subroutine write_laserfields(unit,timestep)
    use atomic_units
    integer, intent(in) :: unit
    real(dp), intent(in), optional :: timestep
    real(dp) :: tt, endtime
    real(dp) :: EL, AL, ZL
    real(dp), dimension(n_laserfields) :: env, part

    write(unit,'(A)') '# t  E(t)  A(t)  Z(t)  I(t)  [E_ii(t)] [envelope_ii(t)]'

    tt = laserfields_starttime()
    endtime = laserfields_endtime()
    do
       EL = get_EL(tt, env, part)
       AL = get_AL(tt)
       ZL = get_ZL(tt)
       write(unit,'(9999(1x,g22.14e3))') tt, EL, AL, ZL, EL**2/au_Wcm2toEL2, part(:), env(:)
       if (tt >= endtime) exit

       if (present(timestep)) then
          tt = tt + timestep
       else
          tt = tt + laserfields_largest_possible_dt(tt,endtime)
       end if
    end do
  end subroutine write_laserfields
  !---------------------------------------------------------------------------
  real(dp) function laserfield_teff(lf,n_photon) result(teff)
    type(laserfield), intent(in) :: lf
    integer, intent(in) :: n_photon
    ! returns the "effective duration" of a laser field for n-photon processes
    ! the values for T_eff are calculated according to
    ! I_0^n * T_eff = \Int_0^T I(t)^n dt = \Int_0^T envelope(t)^(2n) dt
    ! calculations done in intensity_integrals.nb
    select case (lf%form)
    case('gaussianF')
       ! result for exp(-t^2/(r*T^2)) = T sqrt(r PI / 2n)
       ! for gaussian, r = 1/log(16)
       teff = lf%duration * sqrt(PIO2/(n_photon*log(16.d0)))
    case('gaussianI')
       ! for gaussianI, r = 1/log(4)
       teff = lf%duration * sqrt(PIO2/(n_photon*log(4.d0)))
    case('sin2')
       teff = lf%duration * gamma(0.5d0+n_photon*2) / (sqrt(PI)*gamma(1.d0+n_photon*2))
    case('sin_exp')
       teff = lf%duration * gamma(0.5d0+n_photon*lf%form_exponent) / (sqrt(PI)*gamma(1.d0+n_photon*lf%form_exponent))
    case('linear')
       teff = lf%duration + 2*lf%rampon / (1+2.d0*n_photon)
    case('linear2')
       ! flat top for lf%duration with a sin2 pulse of duration 2*lf%rampon as rampon/off
       teff = lf%duration + 2*lf%rampon * gamma(0.5d0+n_photon*2) / (sqrt(PI)*gamma(1.d0+n_photon*2))
    case default
       write(0,*) 'ERROR: effective duration not implemented for laser field form ', lf%form
       STOP
    end select
  end function laserfield_teff
  !---------------------------------------------------------------------------
  real(dp) function laserfield_int(lf,n_photon)
    use atomic_units
    type(laserfield), intent(in) :: lf
    integer, intent(in) :: n_photon
    laserfield_int = laserfield_teff(lf,n_photon) * (lf%intensity_Wcm2 * au_Wcm2)**n_photon
  end function laserfield_int
  !---------------------------------------------------------------------------
  real(dp) function tdcs_factor(lf,n_photon)
    use atomic_units
    real(dp), parameter :: csunit(2) = (/ 1.d24*au_eV/((au_cm)**2), 1.d55*au_eV/((au_cm)**4*au_as*1.d18) /)
    type(laserfield), intent(in) :: lf
    integer, intent(in) :: n_photon
    if (n_photon < 1 .or. n_photon > 2) STOP 'ERROR: Only n_photon = {1,2} allowed for tdcs_factor in laserfields module'
    tdcs_factor = csunit(n_photon) * lf%omega**n_photon / laserfield_int(lf,n_photon)
  end function tdcs_factor
  !-------------------------------------------------------------------
  real(dp) function CS_factor(lf,n_photon)
    use atomic_units
    real(dp), parameter :: csunit(2) = (/ 1.d21/(au_cm)**2, 1.d52/((au_cm)**4*au_as*1.d18) /)
    type(laserfield), intent(in) :: lf
    integer, intent(in) :: n_photon
    if (n_photon < 1 .or. n_photon > 2) STOP 'ERROR: Only n_photon = {1,2} allowed for CS_factor in laserfields module'
    CS_factor = csunit(n_photon) * lf%omega**n_photon / laserfield_int(lf,n_photon)
  end function CS_factor
  !-------------------------------------------------------------------------
  logical function laserfields_can_get_fourier() result(yes_we_can)
    integer :: i_field
    yes_we_can = .true.
    do i_field = 1, n_laserfields
       yes_we_can = yes_we_can .and. lf_can_get_fourier(all_laserfields(i_field))
    end do
  end function laserfields_can_get_fourier
  !-------------------------------------------------------------------------
  logical function lf_can_get_fourier(lf) result(yes_we_can)
    type(laserfield), intent(in) :: lf
    yes_we_can = .true.
    select case (lf%form)
    case('gaussianF','gaussianI','sin2')
       continue
    case('sin_exp','readin')
       write(0,'(9a)') 'WARNING: can not get analytical fourier transform for ', trim(lf%form), ' fields!'
       yes_we_can = .false.
    case('linear','linear2')
       if (lf%linear_chirp_rate_w0as /= 0.d0) then
          write(0,'(9a)') 'WARNING: can not get analytical fourier transform for chirped ', trim(lf%form), ' fields!'
          yes_we_can = .false.
       end if
    end select
  end function lf_can_get_fourier
  !-------------------------------------------------------------------------
  complex(dpc) function lf_get_EL_fourier_transform(lf,omega) result(ELFT)
    ! analytically determine the fourier transform of the defined laser fields
    ! determined as Int exp(-i*omega*t) E(t) dt
    use atomic_units
    type(laserfield), intent(in) :: lf
    real(dp), intent(in) :: omega

    if (.not.lf_can_get_fourier(lf)) then
       write(0,'(A,I2,A)') 'ERROR! can not get analytic fourier transform for this laser field! stopping!'
       STOP 1103
    end if

    if (lf%omega==0) then ! No oscillation enclosed by envelope
       ELFT = lf_envelope_fourier(lf,omega)
    else
       ! with tp = t-tpeak, the whole pulse is
       ! f(t) =  env(t) sin    (phi0 + w0*tp + chirp*tp**2)
       !      = (env(t) exp(IU*(phi0 + w0*tp + chirp*tp**2)) - c.c. ) / (2*IU)
       ! for the fourier transform, we include the chirp term exp(i chirp tp**2) in the envelope.
       ! this part is transformed in lf_envelope_fourier.
       ! exp(IU*phi0) is just a constant prefactor, and the linear phase w0*tp just gives a shift in frequency,
       ! F[f(t) exp(IU w0 t)](w) = F[f(t)](w-w0)
       ! complex conjugation of the transformed function gives complex conjugation + reversal of the argument in the transform, so
       ! F[conjg(f(t) exp(IU w0 t))](w) = conjg(F[f(t) exp(IU w0 t)](-w)) = conjg(F[f(t)](-w-w0))

       ELFT = (        lf_envelope_fourier(lf,  omega - lf%omega)  * exp( IU*PI*lf%phase_pi)             &
            &  - conjg(lf_envelope_fourier(lf, -omega - lf%omega)) * exp(-IU*PI*lf%phase_pi) ) / (2*IU)
    end if

    ! the fourier transform of the part was determined as if it was centered around t=0
    ! shift in time now -- just adds a phase exp(-IU*omega*peak_time), as F[f(t-a)] = exp(-IU*omega*a) F[f(t)]
    ELFT = ELFT * exp(-IU*omega*lf%peak_time)

    if (lf%is_vecpot) then
       ! if this laser field was defined as a vector potential, we need to multiply with -IU*omega to get the fourier transform of the electric field, E=-dA/dt
       ! F[-dA/dt] = -iw F[A]
       ELFT = -IU * omega * ELFT 
       ! in addition, we need to take into account that A0 = E0 / lf%omega if we have an oscillation
       if (lf%omega/=0) ELFT = ELFT / lf%omega
    end if

  end function lf_get_EL_fourier_transform
  !-------------------------------------------------------------------------
  character(4000) function lf_get_EL_fourier_transform_string(lf) result(Ff)
    ! analytically determine the fourier transform of the defined laser fields
    ! determined as Int exp(-i*omega*t) E(t) dt
    use atomic_units
    use laserfields_miscfuncs
    type(laserfield), intent(in) :: lf

    if (.not.lf_can_get_fourier(lf)) then
       write(0,'(A,I2,A)') 'ERROR! can not get analytic fourier transform for this laser field! stopping!'
       STOP 1103
    end if

    if (lf%omega==0) then ! No oscillation enclosed by envelope
       Ff = trim(lf_envelope_fourier_string(lf,'w'))
    else
       ! with tp = t-tpeak, the whole pulse is
       ! f(t) =  env(t) sin    (phi0 + w0*tp + chirp*tp**2)
       !      = (env(t) exp(IU*(phi0 + w0*tp + chirp*tp**2)) - c.c. ) / (2*IU)
       ! for the fourier transform, we include the chirp term exp(i chirp tp**2) in the envelope.
       ! this part is transformed in lf_envelope_fourier.
       ! exp(IU*phi0) is just a constant prefactor, and the linear phase w0*tp just gives a shift in frequency,
       ! F[f(t) exp(IU w0 t)](w) = F[f(t)](w-w0)
       ! complex conjugation of the transformed function gives complex conjugation + reversal of the argument in the transform, so
       ! F[conjg(f(t) exp(IU w0 t))](w) = conjg(F[f(t) exp(IU w0 t)](-w)) = conjg(F[f(t)](-w-w0))
       Ff   =         '('//trim(lf_envelope_fourier_string(lf, ' w - '//gnuplotstring(lf%omega))) &
            &                                                  //' * '//gnuplotstring(exp( IU*PI*lf%phase_pi))//  &
            &  '- conjg('//trim(lf_envelope_fourier_string(lf, '-w - '//gnuplotstring(lf%omega))) &
            &                                                 //') * '//gnuplotstring(exp(-IU*PI*lf%phase_pi))//') / {0,2}'
    end if

    ! the fourier transform of the part was determined as if it was centered around t=0
    ! shift in time now -- just adds a phase exp(-IU*omega*peak_time), as F[f(t-a)] = exp(-IU*omega*a) F[f(t)]
    Ff = trim(Ff) // ' * exp('//gnuplotstring(-IU*lf%peak_time)//'*w)'

    if (lf%is_vecpot) then
       ! if this laser field was defined as a vector potential, we need to multiply with -IU*omega to get the fourier transform of the electric field, E=-dA/dt
       ! F[-dA/dt] = -iw F[A]
       ! in addition, we need to take into account that A0 = E0 / lf%omega
       Ff = '('//trim(Ff)//' * '//gnuplotstring(-IU/lf%omega)//'*w)'
    end if
  end function lf_get_EL_fourier_transform_string
  !-------------------------------------------------------------------------
  complex(dpc) function lf_get_AL_fourier_transform(lf,omega) result(ALFT)
    ! analytically determine the fourier transform of the defined laser fields
    ! determined as Int exp(-i*w*t) A(t) dt
    ! use connection with electric field:
    ! A(t)   = Int exp(i w t)       A(w) dw
    ! -dA/dt = Int exp(i w t) (-iw) A(w) dw = E(t) = Int exp(i w t) E(omega) dw
    ! --> E(omega) = (-iw) A(omega)
    type(laserfield), intent(in) :: lf
    real(dp), intent(in) :: omega
    ! is undefined if called with omega=0.d0, so set to zero for that value
    ALFT = 0
    if (omega /= 0.d0) ALFT = get_EL_fourier_transform(lf,omega) / (-IU*omega)
  end function lf_get_AL_fourier_transform
  !-------------------------------------------------------------------------
  complex(dpc) function laserfields_get_EL_fourier_transform(omega) result(ELFT)
    real(dp), intent(in) :: omega
    integer :: i_field
    ELFT = 0
    do i_field = 1, n_laserfields
       ELFT = ELFT + get_EL_fourier_transform(all_laserfields(i_field),omega)
    end do
  end function laserfields_get_EL_fourier_transform
  !-------------------------------------------------------------------------
  character(4000) function laserfields_get_EL_fourier_transform_string() result(Ff)
    integer :: i_field
    Ff = get_EL_fourier_transform_string(all_laserfields(1))
    do i_field = 2, n_laserfields
       Ff = trim(Ff)//' + '//get_EL_fourier_transform_string(all_laserfields(i_field))
    end do
  end function laserfields_get_EL_fourier_transform_string
  !-------------------------------------------------------------------------
  complex(dpc) function laserfields_get_AL_fourier_transform(omega) result(ALFT)
    real(dp), intent(in) :: omega
    integer :: i_field
    ALFT = 0
    do i_field = 1, n_laserfields
       ALFT = ALFT + get_AL_fourier_transform(all_laserfields(i_field),omega)
    end do
  end function laserfields_get_AL_fourier_transform
  !-------------------------------------------------------------------------
end module laserfields_module
