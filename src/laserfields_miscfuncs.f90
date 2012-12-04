! Copyright (c) 2012, Johannes Feist
! licensed under the MIT open source license, see LICENSE file

module laserfields_miscfuncs
  use nrtype
  implicit none

  !> 1-dimensional polynomial interpolation
  interface interpolate
     module procedure interpolate_single
     module procedure interpolate_array
  end interface
  private interpolate_array, interpolate_single

  !> convert real or complex numbers to strings that can be used in gnuplot
  interface gnuplotstring
     module procedure gnuplotstring_d
     module procedure gnuplotstring_z
  end interface
  private gnuplotstring_d, gnuplotstring_z
contains
  !----------------------------------------------------------------------
  !> obtain the phase (=argument) of a complex number
  real(dp) pure elemental function phase(c)
    complex(dpc), intent(in) :: c
    phase = atan2(aimag(c),real(c))
  end function phase
  !---------------------------------------------------------------------------
  !> return a complex number containing absolute value squared and phase
  !> as real and imaginary part: {|c|<sup>2</sup>,phase(c)}
  complex(dpc) pure elemental function abs2_phase(c)
    complex(dpc), intent(in) :: c
    abs2_phase = cmplx(abs(c)**2,phase(c),dpc)
  end function abs2_phase
  !---------------------------------------------------------------------------
  !> converts a double to a character string for gnuplot
  character(15) function gnuplotstring_d(val) result(ch)
    real(dp), intent(in) :: val
    write(ch,'(sp,es15.8)') val
  end function gnuplotstring_d
  !---------------------------------------------------------------------------
  !> converts a double complex to a character string for gnuplot
  character(33) function gnuplotstring_z(val) result(ch)
    complex(dpc), intent(in) :: val
    write(ch,'(sp,a,es15.8,a,es15.8,a)') '{',real(val),',',aimag(val),'}'
  end function gnuplotstring_z
  !---------------------------------------------------------------------------
  !> for a set of known values (x_old,f_old) calculate the value of f(x) at
  !> point x_new using polynomial interpolation
  real(dp) function interpolate_single(x_old,f_old,x_new,degree) result(f_new)
    real(dp), dimension(:), intent(in) :: x_old,f_old
    real(dp), intent(in) :: x_new
    !> order of the interpolating polynomial
    integer, intent(in) :: degree
    real(dp), dimension(1) :: x_new_arr, f_new_arr
    x_new_arr(1) = x_new
    f_new_arr = interpolate(x_old,f_old,x_new_arr,degree)
    f_new = f_new_arr(1)
  end function interpolate_single
  !---------------------------------------------------------------------------
  !> for a set of known values (x_old,f_old) calculate the values of f(x) at
  !> points x_new using polynomial interpolation
  function interpolate_array(x_old,f_old,x_new,degree) result(f_new)
    real(dp), dimension(:), intent(in) :: x_old, f_old, x_new
    real(dp), dimension(size(x_new)) :: f_new
    !> order of the interpolating polynomial
    integer, intent(in) :: degree
    integer :: ii, npoints_old, npoints_new, nearest_node, spoint, npoints_int

    if (size(x_old)/=size(f_old)) STOP 'ERROR: number of nodes and function values do not agree in interpolate!'

    npoints_new = size(x_new)
    npoints_old = size(x_old)
    ! number of points used for interpolation
    npoints_int=degree+1

    do ii = 1, npoints_new
       if (x_new(ii) < x_old(1) .or. x_new(ii) > x_old(npoints_old)) then
          STOP 'ERROR: extrapolation would be necessary in interpolate!'
       end if

       ! get nearest node smaller then x_new(ii)
       nearest_node = binary_search(x_old,x_new(ii))
       ! get starting point for interpolation
       spoint = min(max(nearest_node-(npoints_int-1)/2,1),npoints_old+1-npoints_int)
       ! call interpolation routine
       f_new(ii) = lagrange_interpol(x_old(spoint:spoint+degree),f_old(spoint:spoint+degree),npoints_int,x_new(ii))
    end do
  end function interpolate_array
  !-----------------------------------------------------------------------
  !> lagrange polynomial interpolation, using the barycentric form (see e.g. wikipedia)
  real(dp) pure function lagrange_interpol(xk,yk,n,x) result(y)
    !> size of xk and yk.
    integer, intent(in) :: n
    !> x-values of points to interpolate.
    real(dp), intent(in) :: xk(n)
    !> y-values of points to interpolate.
    real(dp), intent(in) :: yk(n)
    !> position at which to evaluate the interpolating polynomial.
    real(dp), intent(in) :: x
    real(dp) :: wk(n), nom, denom, diff
    integer :: ii, jj
    do jj = 1, n
       wk(jj) = 1
       do ii = 1, n
          if (ii/=jj) wk(jj) = wk(jj) / (xk(jj)-xk(ii))
       end do
    end do
    nom = 0
    denom = 0
    do jj = 1, N
       diff = abs(x-xk(jj))
       ! this guards against x=xk(jj)=0
       if (diff/=0) diff = 2*diff / (abs(x)+abs(xk(jj)))
       if (diff < 1.e-12) then
          y = yk(jj)
          return
       end if
       nom   = nom   + wk(jj)*yk(jj) / (x - xk(jj))
       denom = denom + wk(jj)        / (x - xk(jj))
    end do
    y = nom/denom
  end function lagrange_interpol
  !-----------------------------------------------------------------------
  !> do a binary search and return index il for which xs(il) <= x <= xs(il+1)
  pure integer function binary_search(xs,x) result(il)
    !> input array to search, has to be sorted (either in ascending or descending order)
    real(dp), intent(in) :: xs(:)
    !> value to search for
    real(dp), intent(in) :: x
    integer :: n, im, iu
    logical :: increasing

    n = size(xs)
    ! we do not check whether the sequence is actually monotonic, but we handle increasing and decreasing sequences
    increasing = xs(1) < xs(n)

    if (increasing) then
       if      (x < xs(1)) then; il = 1; return
       else if (x > xs(n)) then; il = n; return
       end if
    else
       if      (x > xs(1)) then; il = 1; return
       else if (x < xs(n)) then; il = n; return
       end if
    end if

    il=0   ! initial value for lower border
    iu=n+1 ! initial value for upper border
    do while (iu-il>1)
       im=(iu+il)/2
       if (increasing.eqv.(x > xs(im))) then
          il=im
       else
          iu=im
       end if
    end do
  end function binary_search
  !-----------------------------------------------------------------------
end module laserfields_miscfuncs
