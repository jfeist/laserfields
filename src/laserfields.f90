! Copyright (c) 2012, Johannes Feist
! licensed under the MIT open source license, see LICENSE file

module laserfields
  ! this is a trick for making use with intel fortran easier:
  ! we explicitly list the variables/functions we want to import, which
  ! causes the generated laserfields.mod file to contain the definitions themselves and
  ! not to depend on laserfields_module.mod and laserfields_paramfilehandling.mod
  ! in addition, this gives us the opportunity not to export all functions in these
  ! modules for "external" use
  use nrtype, ONLY : dp, dpc
  use laserfields_module, only : &
       & laserfield, all_laserfields, n_laserfields, get_EL, get_AL, get_ZL, &
       & get_EL_fourier_transform, get_AL_fourier_transform, get_EL_fourier_transform_string, &
       & add_laserfield, make_laserfield, &
       & lf_get_envelope, lf_envelope_fourier, lf_envelope_fourier_string, lf_get_omega, &
       & lf_get_starttime, lf_get_endtime, laserfields_starttime, laserfields_endtime, &
       & laserfields_smallest_TX, laserfields_largest_possible_dt, write_laserfields, &
       & laserfield_teff, laserfield_int, tdcs_factor, CS_factor, laserfields_can_get_fourier, &
       & lf_can_get_fourier
  use laserfields_paramfilehandling, only : laserfields_read_parameters, laserfields_write_parameters
end module laserfields
