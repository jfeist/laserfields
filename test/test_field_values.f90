program test_field_values
  use nrtype, only: dp, dpc, twopi
  use atomic_units, only: au_as, au_Wcm2toEL2, au_nm, au_c
  use laserfields
  implicit none

  integer, parameter :: nfields = 6, nsamples = 3
  real(dp), parameter :: atol = 1.d-11
  real(dp), parameter :: atol_small = 1.d-12

  type(laserfield) :: lfs(nfields)
  real(dp) :: sample_shifts(nsamples)
  real(dp) :: sample_times(nsamples)
  real(dp) :: expected_e(nfields, nsamples)
  real(dp) :: expected_a(nfields, nsamples)
  complex(dpc) :: expected_ep(nfields, nsamples)
  complex(dpc) :: expected_ap(nfields, nsamples)
  real(dp) :: t, t_before, t_after, min_start, max_end
  real(dp) :: intensity_wcm2, lambda_nm, phase_pi
  real(dp) :: e0_ref, omega0_ref, t0_ref, sigma_ref, t_ref, tflat_ref, tramp_ref
  real(dp) :: duration_gaussian_as, duration_sin_as, duration_flat_as, rampon_as
  real(dp) :: sum_e, sum_a
  complex(dpc) :: sum_ep, sum_ap, rec
  integer :: i, j, fails

  fails = 0
  sample_shifts = (/ -0.2d0, 0.d0, 0.3d0 /)

  e0_ref = 1.5d0
  omega0_ref = 0.12d0
  t0_ref = 500.d0
  sigma_ref = 100.d0
  t_ref = 800.d0
  tflat_ref = 400.d0
  tramp_ref = 150.d0
  phase_pi = 0.8d0

  intensity_wcm2 = e0_ref**2 / au_Wcm2toEL2
  lambda_nm = twopi * au_c / (omega0_ref * au_nm)
  duration_gaussian_as = sigma_ref * sqrt(log(16.d0)) / au_as
  duration_sin_as = t_ref / au_as
  duration_flat_as = tflat_ref / au_as
  rampon_as = tramp_ref / au_as

  lfs(1) = make_laserfield('gaussianI', intensity_wcm2, lambda_nm, t0_ref/au_as, duration_gaussian_as, &
  & phase_pi=phase_pi, is_vecpot=.true., linear_chirp_rate_w0as=0.d0)
  lfs(2) = make_laserfield('sin2', intensity_wcm2, lambda_nm, t0_ref/au_as, duration_sin_as, &
  & phase_pi=phase_pi, is_vecpot=.true., linear_chirp_rate_w0as=0.d0)
  lfs(3) = make_laserfield('sin_exp', intensity_wcm2, lambda_nm, t0_ref/au_as, duration_sin_as, &
  & form_exponent=4.d0, phase_pi=phase_pi, is_vecpot=.true., linear_chirp_rate_w0as=0.d0)
  lfs(4) = make_laserfield('sin_exp', intensity_wcm2, lambda_nm, t0_ref/au_as, duration_sin_as, &
  & form_exponent=7.d0, phase_pi=phase_pi, is_vecpot=.true., linear_chirp_rate_w0as=0.d0)
  lfs(5) = make_laserfield('linear', intensity_wcm2, lambda_nm, t0_ref/au_as, duration_flat_as, &
  & rampon_as=rampon_as, phase_pi=phase_pi, is_vecpot=.false., linear_chirp_rate_w0as=0.d0)
  lfs(6) = make_laserfield('linear2', intensity_wcm2, lambda_nm, t0_ref/au_as, duration_flat_as, &
  & rampon_as=rampon_as, phase_pi=phase_pi, is_vecpot=.true., linear_chirp_rate_w0as=0.d0)

  expected_e(1,:) = (/ -0.4733721103340174d0, 1.213525491562421d0, 0.43939711942241655d0 /)
  expected_e(2,:) = (/ -0.46657740126715147d0, 1.213525491562421d0, 0.4560190763122897d0 /)
  expected_e(3,:) = (/ -0.4696176702845257d0, 1.213525491562421d0, 0.44856301902567197d0 /)
  expected_e(4,:) = (/ -0.47415631450844775d0, 1.213525491562421d0, 0.4374727545356417d0 /)
  expected_e(5,:) = (/ 1.4265847744427302d0, 0.8816778784387098d0, -1.4265847744427274d0 /)
  expected_e(6,:) = (/ -0.4635254915624215d0, 1.213525491562421d0, 0.4635254915624302d0 /)

  expected_a(1,:) = (/ 11.823200448244128d0, 7.347315653655915d0, -11.742442578919206d0 /)
  expected_a(2,:) = (/ 11.868113281037923d0, 7.347315653655915d0, -11.843028666243946d0 /)
  expected_a(3,:) = (/ 11.848054069403913d0, 7.347315653655915d0, -11.798022564282457d0 /)
  expected_a(4,:) = (/ 11.818028803324419d0, 7.347315653655915d0, -11.730833895083071d0 /)
  expected_a(5,:) = 0.d0
  expected_a(6,:) = (/ 11.888206453689419d0, 7.347315653655915d0, -11.888206453689396d0 /)

  expected_ep(1,:) = (/ (-0.2366860551670087d0, 0.7073805747087495d0), (0.6067627457812105d0, 0.4408389392193549d0), (0.21969855971120827d0, -0.7075431243057029d0) /)
  expected_ep(2,:) = (/ (-0.23328870063357574d0, 0.7114637065144673d0), (0.6067627457812105d0, 0.4408389392193549d0), (0.22800953815614486d0, -0.7115150383005721d0) /)
  expected_ep(3,:) = (/ (-0.23480883514226286d0, 0.7096391697345191d0), (0.6067627457812105d0, 0.4408389392193549d0), (0.22428150951283599d0, -0.7097408968808412d0) /)
  expected_ep(4,:) = (/ (-0.23707815725422388d0, 0.7069101152176173d0), (0.6067627457812105d0, 0.4408389392193549d0), (0.21873637726782086d0, -0.7070857016210783d0) /)
  expected_ep(5,:) = (/ (0.7132923872213651d0, 0.23176274578121075d0), (0.4408389392193549d0, -0.6067627457812105d0), (-0.7132923872213637d0, -0.23176274578121514d0) /)
  expected_ep(6,:) = (/ (-0.23176274578121075d0, 0.7132923872213651d0), (0.6067627457812105d0, 0.4408389392193549d0), (0.2317627457812151d0, -0.7132923872213638d0) /)

  expected_ap(1,:) = (/ (5.911600224122064d0, 1.9207953490721235d0), (3.6736578268279576d0, -5.0563562148434205d0), (-5.871221289459603d0, -1.907675437887428d0) /)
  expected_ap(2,:) = (/ (5.934056640518961d0, 1.9280918810662835d0), (3.6736578268279576d0, -5.0563562148434205d0), (-5.921514333121973d0, -1.9240166383568338d0) /)
  expected_ap(3,:) = (/ (5.924027034701957d0, 1.924833064590886d0), (3.6736578268279576d0, -5.0563562148434205d0), (-5.899011282141228d0, -1.9167049538678569d0) /)
  expected_ap(4,:) = (/ (5.909014401662209d0, 1.9199551644239556d0), (3.6736578268279576d0, -5.0563562148434205d0), (-5.865416947541536d0, -1.9057894928745776d0) /)
  expected_ap(5,:) = (0.d0, 0.d0)
  expected_ap(6,:) = (/ (5.944103226844709d0, 1.9313562148434231d0), (3.6736578268279576d0, -5.0563562148434205d0), (-5.944103226844698d0, -1.9313562148434595d0) /)

  do i = 1, nfields
    sample_times = lfs(i)%peak_time + sample_shifts * lfs(i)%TX

    do j = 1, nsamples
      t = sample_times(j)
      call check_close_real('E field', get_EL(lfs(i), t), expected_e(i,j), atol, fails)
      call check_close_complex('E posfreq', get_EL_posfreq(lfs(i), t), expected_ep(i,j), atol, fails)

      rec = get_EL_posfreq(lfs(i), t) + conjg(get_EL_posfreq(lfs(i), t))
      call check_close_real('E reconstruction real part', real(rec, dp), get_EL(lfs(i), t), atol, fails)
      call check_small('E reconstruction imaginary part', abs(aimag(rec)), atol, fails)

      if (lfs(i)%is_vecpot) then
        call check_close_real('A field', get_AL(lfs(i), t), expected_a(i,j), atol, fails)
        call check_close_complex('A posfreq', get_AL_posfreq(lfs(i), t), expected_ap(i,j), atol, fails)

        rec = get_AL_posfreq(lfs(i), t) + conjg(get_AL_posfreq(lfs(i), t))
        call check_close_real('A reconstruction real part', real(rec, dp), get_AL(lfs(i), t), atol, fails)
        call check_small('A reconstruction imaginary part', abs(aimag(rec)), atol, fails)
      end if
    end do

    t_before = lf_get_starttime(lfs(i)) - 0.1d0 * lfs(i)%TX
    t_after = lf_get_endtime(lfs(i)) + 0.1d0 * lfs(i)%TX
    if (trim(lfs(i)%form) == 'gaussianI') then
      call check_small('Gaussian E before pulse', abs(get_EL(lfs(i), t_before)), atol_small, fails)
      call check_small('Gaussian E after pulse', abs(get_EL(lfs(i), t_after)), atol_small, fails)
      if (lfs(i)%is_vecpot) then
        call check_small('Gaussian A before pulse', abs(get_AL(lfs(i), t_before)), atol_small, fails)
        call check_small('Gaussian A after pulse', abs(get_AL(lfs(i), t_after)), atol_small, fails)
      end if
    else
      call check_small('Compact-support E before pulse', abs(get_EL(lfs(i), t_before)), atol, fails)
      call check_small('Compact-support E after pulse', abs(get_EL(lfs(i), t_after)), atol, fails)
      if (lfs(i)%is_vecpot) then
        call check_small('Compact-support A before pulse', abs(get_AL(lfs(i), t_before)), atol, fails)
        call check_small('Compact-support A after pulse', abs(get_AL(lfs(i), t_after)), atol, fails)
      end if
    end if
  end do

  n_laserfields = 0
  do i = 1, nfields
    call add_laserfield(lfs(i))
  end do

  t = 300.d0
  sum_e = 0.d0
  sum_ep = 0.d0
  do i = 1, nfields
    sum_e = sum_e + get_EL(lfs(i), t)
    sum_ep = sum_ep + get_EL_posfreq(lfs(i), t)
  end do
  call check_close_real('Collection E sum', get_EL(t), sum_e, atol, fails)
  call check_close_complex('Collection E posfreq sum', get_EL_posfreq(t), sum_ep, atol, fails)

  do j = 1, nsamples
    t = lfs(1)%peak_time + sample_shifts(j) * lfs(1)%TX
    rec = get_EL_posfreq(t) + conjg(get_EL_posfreq(t))
    call check_close_real('Collection E reconstruction real part', real(rec, dp), get_EL(t), atol, fails)
    call check_small('Collection E reconstruction imaginary part', abs(aimag(rec)), atol, fails)

  end do

  n_laserfields = 0
  do i = 1, 4
    call add_laserfield(lfs(i))
  end do
  call add_laserfield(lfs(6))

  t = 300.d0
  sum_a = 0.d0
  sum_ap = 0.d0
  do i = 1, 4
    sum_a = sum_a + get_AL(lfs(i), t)
    sum_ap = sum_ap + get_AL_posfreq(lfs(i), t)
  end do
  sum_a = sum_a + get_AL(lfs(6), t)
  sum_ap = sum_ap + get_AL_posfreq(lfs(6), t)
  call check_close_real('Collection A sum', get_AL(t), sum_a, atol, fails)
  call check_close_complex('Collection A posfreq sum', get_AL_posfreq(t), sum_ap, atol, fails)

  do j = 1, nsamples
    t = lfs(1)%peak_time + sample_shifts(j) * lfs(1)%TX
    rec = get_AL_posfreq(t) + conjg(get_AL_posfreq(t))
    call check_close_real('Collection A reconstruction real part', real(rec, dp), get_AL(t), atol, fails)
    call check_small('Collection A reconstruction imaginary part', abs(aimag(rec)), atol, fails)
  end do

  min_start = minval((/ (lf_get_starttime(lfs(i)), i=1,nfields) /))
  max_end = maxval((/ (lf_get_endtime(lfs(i)), i=1,nfields) /))
  call check_close_real('Collection start time', laserfields_starttime(), min_start, atol, fails)
  call check_close_real('Collection end time', laserfields_endtime(), max_end, atol, fails)

  if (fails > 0) then
    write(*,'(A,I0)') 'FAIL test_field_values: ', fails
    error stop 1
  end if

  write(*,'(A)') 'PASS test_field_values'

contains

  subroutine check_close_real(label, value, expected, tolerance, nfail)
    character(*), intent(in) :: label
    real(dp), intent(in) :: value, expected, tolerance
    integer, intent(inout) :: nfail
    if (abs(value-expected) > tolerance) then
      nfail = nfail + 1
      write(*,'(A,2X,A,2X,ES24.15,2X,ES24.15,2X,ES10.3)') 'FAIL', trim(label), value, expected, abs(value-expected)
    end if
  end subroutine check_close_real

  subroutine check_close_complex(label, value, expected, tolerance, nfail)
    character(*), intent(in) :: label
    complex(dpc), intent(in) :: value, expected
    real(dp), intent(in) :: tolerance
    integer, intent(inout) :: nfail
    real(dp) :: err
    err = abs(value-expected)
    if (err > tolerance) then
      nfail = nfail + 1
      write(*,'(A,2X,A,2X,2(ES24.15,1X),2X,2(ES24.15,1X),2X,ES10.3)') 'FAIL', trim(label), &
        & real(value,dp), aimag(value), real(expected,dp), aimag(expected), err
    end if
  end subroutine check_close_complex

  subroutine check_small(label, value, tolerance, nfail)
    character(*), intent(in) :: label
    real(dp), intent(in) :: value, tolerance
    integer, intent(inout) :: nfail
    if (value > tolerance) then
      nfail = nfail + 1
      write(*,'(A,2X,A,2X,ES24.15,2X,ES10.3)') 'FAIL', trim(label), value, tolerance
    end if
  end subroutine check_small

end program test_field_values