! Copyright (c) 2012, Johannes Feist
! licensed under the MIT open source license, see LICENSE file

program printlaserfourier
  use laserfields
  implicit none
  integer :: i_field
  character(len=200) :: parfile

  parfile = 'laserfields.in'
  if (command_argument_count() >= 1) then
     call get_command_argument(1,parfile)
  end if

  call laserfields_read_parameters(trim(parfile))

  if (.not.laserfields_can_get_fourier()) then
     write(0,*) 'error: can not get analytical fourier transform for these laser fields. stopping!'
     stop 733
  end if

  write(6,'(3a)') 'ffall(w) = ',trim(get_el_fourier_transform_string()),';'
  do i_field = 1, n_laserfields
     write(6,'(a,i2.2,3a)') 'ff_',i_field,'(w) = ',trim(get_el_fourier_transform_string(all_laserfields(i_field))),';'
  end do

end program printlaserfourier
