! Copyright (c) 2012, Johannes Feist
! licensed under the MIT open source license, see LICENSE file

module nrtype
  ! based on nrtype.f90 from numerical recipes, which is in the public domain (http://www.nr.com/public-domain.html)
  integer, parameter :: dp  = kind(1.0d0)
  integer, parameter :: dpc = kind((1.0d0,1.0d0))
  real(dp), parameter :: PI    = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: PIO2  = 1.57079632679489661923132169163975144209858_dp
  real(dp), parameter :: TWOPI = 6.283185307179586476925286766559005768394_dp
  complex(dpc), parameter :: IU = (0.d0,1.d0)
end module nrtype
