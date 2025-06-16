! Copyright (c) 2012, Johannes Feist
! licensed under the MIT open source license, see LICENSE file

program plotlaserfield
  use laserfields
  implicit none
  character(len=200) :: parfile, arg
  ! Default value for timestep is -1, which indicates that it is not set.
  real(dp) :: timestep = -1.d0
  integer :: iarg

  parfile = 'laserfields.in'
  iarg = 1
  do while (iarg <= command_argument_count())
    call get_command_argument(iarg, arg)
    if (trim(arg) == '--dt') then
      iarg = iarg + 1
      call get_command_argument(iarg, arg)
      read(arg,*) timestep
      if (timestep <= 0) then
        write(6,*) 'Error: Time step must be positive.'
        stop
      end if
    else
      parfile = trim(arg)
    end if
    iarg = iarg + 1
  end do
  call laserfields_read_parameters(trim(parfile))
  call laserfields_write_parameters('laserfields_parameters.dat')

  if (timestep > 0) then
    call write_laserfields(6, timestep)
  else
    call write_laserfields(6)
  end if

end program plotlaserfield
