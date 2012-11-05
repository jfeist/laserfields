! Copyright (c) 2012, Johannes Feist
! licensed under the MIT open source license, see LICENSE file

program plotlaserfourier
  use laserfields
  implicit none
  real(dp), parameter :: TWOPI = 6.283185307179586476925286766559005768394_dp
  integer :: ii, i_field, nsteps
  real(dp) :: omega, de, emax
  complex(dpc), dimension(:), allocatable :: lf_ft
  complex(dpc) :: alft
  character(len=200) :: parfile

  parfile = 'laserfields.in'
  if (command_argument_count() >= 1) then
     call get_command_argument(1,parfile)
  end if
  call laserfields_read_parameters(trim(parfile))

  if (.not.laserfields_can_get_fourier()) then
     write(0,*) 'ERROR: can not get analytical fourier transform for these laser fields. Stopping!'
     STOP 733
  end if

  allocate(lf_ft(0:n_laserfields))

  ! determine energy resolution by estimating fourier bandwidth
  de = 1.d0
  emax = 0.d0
  do i_field = 1, n_laserfields
     de = min(de,TWOPI/(100*all_laserfields(i_field)%duration))
     ! go up to twice the maximum frequency - for really broadband pulses, this might not actually be enough
     emax = max(emax,all_laserfields(i_field)%omega*2.01) ! use 2.01 here so we do not get omega==lf%omega as an argument, where some of the routines have problems
  end do
  nsteps = nint(emax/de)
  de = emax/(nsteps+1)

  do ii = 1, nsteps
     omega = ii * de
     do i_field = 1, n_laserfields
        lf_ft(i_field) = get_EL_fourier_transform(all_laserfields(i_field), omega)
     end do
     lf_ft(0) = sum(lf_ft(1:))
     alft = get_AL_fourier_transform(omega)

     ! we are not unrolling the phases here (i.e. removing phase jumps by 2 pi)
     ! neither are we shifting in time to get a flatter phase overall
     write(6,'(9999(1x,g22.14e3))') omega, lf_ft(0), alft, lf_ft(1:)
  end do

end program plotlaserfourier
