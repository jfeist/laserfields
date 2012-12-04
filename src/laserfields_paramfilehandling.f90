! Copyright (c) 2012, Johannes Feist
! licensed under the MIT open source license, see LICENSE file

module laserfields_paramfilehandling
  use nrtype
  implicit none
  private ! make everything private by default
  public :: laserfields_read_parameters, laserfields_write_parameters

  !> read laserfield parameters and put into global all_laser_fields array

  !> can read from either a file or a unit
  interface laserfields_read_parameters
     module procedure read_parameters_from_file
     module procedure read_parameters_from_unit
  end interface
  private :: read_parameters_from_file, read_parameters_from_unit

  !> write parameters of all laserfields in global all_laser_fields array

  !> can write to either a file or a unit
  interface laserfields_write_parameters
     module procedure write_parameters_to_file
     module procedure write_parameters_to_unit
  end interface
  private :: write_parameters_to_file, write_parameters_to_unit

  interface write_param
     module procedure write_int_param
     module procedure write_double_param
     module procedure write_logical_param
     module procedure write_char_param
  end interface
  interface read_param
     module procedure read_int_param
     module procedure read_double_param
     module procedure read_logical_param
     module procedure read_char_param
  end interface

contains
  !---------------------------------------------------------------------------
  subroutine write_parameters_to_file(filename)
    use laserfields_miscfuncs, only: get_unused_unit
    character(len=*), intent(in) :: filename
    integer :: unit
    unit = get_unused_unit()
    open(unit,file=filename,status='unknown',action='write')
    call write_parameters_to_unit(unit)
    close(unit)
  end subroutine write_parameters_to_file
  !---------------------------------------------------------------------------
  subroutine write_parameters_to_unit(unit)
    use laserfields_module
    integer, intent(in) :: unit
    type(laserfield) :: lf
    integer :: ii
    character(len=2) :: chii

    do ii = 1, n_laserfields
       write(unit,*) ''
       write(chii,'(i2)') ii
       call write_header(unit,'laser '//chii//' parameters')
       lf = all_laserfields(ii)
       call write_param(unit,'laserfield',lf%form)

       call write_param(unit,"datafile",lf%datafile)
       call write_param(unit,"is_vecpot",lf%is_vecpot)
       call write_param(unit,"intensity_Wcm2",lf%intensity_Wcm2)
       call write_param(unit,"lambda_nm",lf%lambda_nm)
       call write_param(unit,"peak_time_as",lf%peak_time_as)
       call write_param(unit,"duration_as",lf%duration_as)
       call write_param(unit,"rampon_as",lf%rampon_as)
       call write_param(unit,"phase_pi",lf%phase_pi)
       call write_param(unit,"form_exponent",lf%form_exponent)
       call write_param(unit,"linear_chirp_rate_w0as",lf%linear_chirp_rate_w0as)
       write(unit,*) ''

       call write_param(unit,"E0",lf%E0)
       call write_param(unit,"omega",lf%omega)
       call write_param(unit,"TX",lf%TX)
       call write_param(unit,"peak_time",lf%peak_time)
       call write_param(unit,"duration",lf%duration)
       call write_param(unit,"rampon",lf%rampon)
       call write_param(unit,'laserfield_end',1)
    end do
  end subroutine write_parameters_to_unit
  !---------------------------------------------------------------------------
  subroutine read_parameters_from_file(filename)
    use laserfields_miscfuncs, only: get_unused_unit
    character(len=*), intent(in) :: filename
    integer :: unit
    unit = get_unused_unit()
    open(unit,file=filename,status='old',action='read')
    call read_parameters_from_unit(unit)
    close(unit)
  end subroutine read_parameters_from_file
  !---------------------------------------------------------------------------
  subroutine read_parameters_from_unit(unit)
    integer, intent(in) :: unit
    integer :: status
    character(len=200) :: name, valuestr
    readparams : do
       status = get_next_param(unit,name,valuestr)
       if (status /= 0) exit readparams
       if (name == 'laserfield') then
          call read_laserfield_parameters(unit,valuestr)
       else
          write(0,'(4A)') 'ERROR: unknown parameter encountered, name = ', trim(adjustl(name)), &
               & ' value = ', trim(adjustl(valuestr))
          STOP 14
       end if
    end do readparams
  end subroutine read_parameters_from_unit
  !---------------------------------------------------------------------------
  subroutine read_laserfield_parameters(unit, form)
    use laserfields_module
    integer, intent(in) :: unit
    character(len=*), intent(in) :: form
    type(laserfield) :: lf
    integer :: status
    character(len=200) :: name, valuestr

    if (len_trim(form) > len(lf%form)) then
       write(0,*) 'laser field form too long in read_laserfield_parameters, form =', form
       STOP 951
    end if

    lf%form = form

    status = get_next_param(unit,name,valuestr)
    do while (name /= 'laserfield_end')
       if (name == "form") read(valuestr,*) lf%form
       if (name == "datafile") read(valuestr,*) lf%datafile
       if (name == "is_vecpot") read(valuestr,*) lf%is_vecpot
       if (name == "intensity_Wcm2") read(valuestr,*) lf%intensity_Wcm2
       if (name == "lambda_nm") read(valuestr,*) lf%lambda_nm
       if (name == "peak_time_as") read(valuestr,*) lf%peak_time_as
       if (name == "duration_as") read(valuestr,*) lf%duration_as
       if (name == "rampon_as") read(valuestr,*) lf%rampon_as
       if (name == "phase_pi") read(valuestr,*) lf%phase_pi
       if (name == "form_exponent") read(valuestr,*) lf%form_exponent
       if (name == "linear_chirp_rate_w0as") read(valuestr,*) lf%linear_chirp_rate_w0as
       status = get_next_param(unit,name,valuestr)
    end do

    call laserfield_set_dependent(lf)
    if (lf%E0 /= 0.d0) call add_laserfield(lf)

  end subroutine read_laserfield_parameters
  !---------------------------------------------------------------------------
  subroutine write_header(unit,name)
    integer, parameter :: headerlen = 60
    integer, intent(in) :: unit
    character(len=*), intent(in) :: name
    character(len=headerlen) :: header
    integer :: inlen, padl, padr
    inlen = len_trim(name)
    header(1:headerlen) = repeat('#',headerlen)
    if (inlen > 0) then
       padl = (headerlen - inlen)/2
       padr = headerlen - inlen - padl ! strlen = padl + inlen + padr
       if (padr < 0) write(0,'(a,a)') 'len(name) too large in write_header, name =',name
       header(padl:padl+inlen+1) = ' '//trim(name)//' '
    end if
    write(unit,'(a)') header
  end subroutine write_header
  !---------------------------------------------------------------------------
  function param_name_string(name)
    integer, parameter :: strlen = 30
    character(len=strlen) :: param_name_string
    character(len=*), intent(in) :: name
    integer :: inlen, padl, padr
    inlen = len_trim(name)
    padl = 3
    padr = strlen - inlen - padl ! strlen = padl + inlen + padr
    if (padr < 0) write(0,'(a,a)') 'len(name) too large in param_name_string, name =',name
    param_name_string(1:strlen) = repeat(' ',strlen)
    param_name_string(padl+1:padl+inlen) = trim(name)
  end function param_name_string
  !---------------------------------------------------------------------------
  subroutine write_int_param(unit,name,value)
    integer, intent(in) :: unit, value
    character(len=*), intent(in) :: name
    character(len=30) :: tmpch
    write(tmpch,'(i30)') value
    write(unit,'(3a)') param_name_string(name), '  =  ', trim(adjustl(tmpch))
  end subroutine write_int_param
  !---------------------------------------------------------------------------
  subroutine write_double_param(unit,name,value)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: value
    character(len=30) :: tmpch
    write(tmpch,'(g30.18e3)') value
    write(unit,'(3a)') param_name_string(name), '  =  ', trim(adjustl(tmpch))
  end subroutine write_double_param
  !---------------------------------------------------------------------------
  subroutine write_logical_param(unit,name,value)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: name
    logical, intent(in) :: value
    character(len=30) :: tmpch
    write(tmpch,'(l1)') value
    write(unit,'(3a)') param_name_string(name), '  =  ', trim(adjustl(tmpch))
  end subroutine write_logical_param
  !---------------------------------------------------------------------------
  subroutine write_char_param(unit,name,value)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: name, value
    write(unit,'(3a)') param_name_string(name), '  =  ', trim(adjustl(value))
  end subroutine write_char_param
  !---------------------------------------------------------------------------
  subroutine read_int_param(valstr,par)
    character(len=*), intent(in) :: valstr
    integer, intent(out) :: par
    read(valstr,*) par
  end subroutine read_int_param
  !---------------------------------------------------------------------------
  subroutine read_double_param(valstr,par)
    character(len=*), intent(in) :: valstr
    real(dp), intent(out) :: par
    read(valstr,*) par
  end subroutine read_double_param
  !---------------------------------------------------------------------------
  subroutine read_logical_param(valstr,par)
    character(len=*), intent(in) :: valstr
    logical, intent(out) :: par
    read(valstr,*) par
  end subroutine read_logical_param
  !---------------------------------------------------------------------------
  subroutine read_char_param(valstr,par)
    character(len=*), intent(in) :: valstr
    character(len=*), intent(out) :: par
    par = valstr
  end subroutine read_char_param
  !---------------------------------------------------------------------------
  subroutine get_paramstr(instr,name,valuestr)
    character(len=*), intent(in) :: instr
    character(len=*), intent(out) :: name, valuestr
    integer :: ii

    ii = 1
    do while (instr(ii:ii) /= '=')
       ii = ii + 1
       if (ii > len(instr)) then
          write(0,*) 'no equal sign (=) in parameter string!'
          stop 531
       end if
    end do

    name = trim(adjustl(instr(:ii-1)))
    valuestr = trim(adjustl(instr(ii+1:)))

  end subroutine get_paramstr
  !---------------------------------------------------------------------------
  function get_next_param(unit,name,valuestr)
    integer :: get_next_param
    integer, intent(in) :: unit
    character(len=*), intent(out) :: name, valuestr
    character(len=200) :: tmpstr

    get_next_param = 0

    tmpstr(1:1) = '#'
    do while (tmpstr(1:1) == '#' .or. len(trim(tmpstr))==0)
       read(unit,'(a200)',iostat=get_next_param) tmpstr
       ! if end of file occured
       if (get_next_param < 0) return
       ! if error occured
       if (get_next_param > 0) then
          inquire(unit=unit, name=tmpstr)
          write(0,'(a,i5,2a)') 'error in get_next_param while reading from unit ', unit, ' = file ',tmpstr
          stop 19
       end if
    end do
    call get_paramstr(tmpstr,name,valuestr)

  end function get_next_param
  !---------------------------------------------------------------------------
end module laserfields_paramfilehandling
