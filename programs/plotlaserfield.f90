! Copyright (c) 2012, Johannes Feist
! licensed under the MIT open source license, see LICENSE file

program plotlaserfield
  use laserfields
  implicit none
  character(len=200) :: parfile

  parfile = 'laserfields.in'
  if (command_argument_count() >= 1) then
     call get_command_argument(1,parfile)
  end if
  call laserfields_read_parameters(trim(parfile))
  call laserfields_write_parameters('laserfields_parameters.dat')

  call write_laserfields(6)

end program plotlaserfield
