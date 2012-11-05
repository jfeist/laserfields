! Copyright (c) 2012, Johannes Feist
! licensed under the MIT open source license, see LICENSE file

module atomic_units
  use nrtype
  ! these units are based on the latest CODATA values from physics.nist.gov/constants from 2012-10-28
  real(dp), parameter :: au_as   = 1.d0/24.18884326211233d0  ! attosecond in a.u.
  real(dp), parameter :: au_Wcm2toEL2 = 1.d0/3.5094452159384d16 ! from W/cm^2 to E^2 in a.u.
  real(dp), parameter :: au_Wcm2 = 1.d0/6.436409342904736d15  ! W/cm^2 in a.u.
  real(dp), parameter :: au_cm   = 1.d0/5.2917721092d-9      ! cm in a.u.
  real(dp), parameter :: au_nm   = 1.d0/5.2917721092d-2      ! nm in a.u.
  real(dp), parameter :: au_c    = 137.035999074             ! c (speed of light) in a.u. == 1/alpha
  real(dp), parameter :: au_eV   = 1.d0/27.21138505732838d0  ! eV in a.u.
end module atomic_units
