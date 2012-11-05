! Copyright (c) 2012, Johannes Feist
! licensed under the MIT open source license, see LICENSE file

program main
  use laserfields
  implicit none
  real(dp) :: tt, dt

  call laserfields_read_parameters('laserfields.in')

  !             make_laserfield(form,intensity_Wcm2,lambda_nm,peak_time_as,duration_as,rampon_as,form_exponent,phase_pi,is_vecpot,group,linear_chirp_rate_w0as)
  call add_laserfield(make_laserfield('sin2',1.d15,20.d0,100.d0,200.d0))
  call add_laserfield(make_laserfield('weirdfield.dat'))

  dt = laserfields_smallest_TX() / 50.
  tt = laserfields_starttime() - dt
  do while (tt<=laserfields_endtime())
     tt = tt + dt
     write(6,*) tt, get_EL(tt)
  end do

end program main
