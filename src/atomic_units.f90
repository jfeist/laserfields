! Copyright (c) 2012, Johannes Feist
! licensed under the MIT open source license, see LICENSE file

!> conversion factors for various units to/from atomic units
!> based on the 2012-10-28 CODATA values from physics.nist.gov/constants
module atomic_units
  use nrtype
  !> t [a.u.] = t [attoseconds] * au_as
  real(dp), parameter :: au_as   = 1.d0/24.18884326211233d0

  !> E<sub>F</sub> [a.u.] = sqrt(I [W/cm<sup>2</sup>] * au_Wcm2toEL2),

  !> i.e. convert from intensity in W/cm<sup>2</sup> to peak electric field squared in a.u.
  real(dp), parameter :: au_Wcm2toEL2 = 1.d0/3.5094452159384d16

  !> I [a.u.] = I [W/cm<sup>2</sup>] * au_Wcm2
  real(dp), parameter :: au_Wcm2 = 1.d0/6.436409342904736d15

  !> x [a.u.] = x [cm] * au_cm
  real(dp), parameter :: au_cm   = 1.d0/5.2917721092d-9

  !> x [a.u.] = x [nm] * au_nm
  real(dp), parameter :: au_nm   = 1.d0/5.2917721092d-2

  !> speed of light in a.u. == 1/alpha
  real(dp), parameter :: au_c    = 137.035999074

  !> E [a.u.] = E [eV] * au_eV (in Hartree atomic units)
  real(dp), parameter :: au_eV   = 1.d0/27.21138505732838d0
end module atomic_units
