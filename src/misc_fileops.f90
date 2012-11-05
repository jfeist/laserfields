! Copyright (c) 2012, Johannes Feist
! licensed under the MIT open source license, see LICENSE file

module misc_fileops
  use nrtype
  implicit none
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
  function get_unused_unit() result(unit)
    integer :: unit
    logical :: opened
    ! find an unused unit to use
    unit = 16
    opened = .true.
    do while (opened)
       unit = unit + 1
       inquire(unit=unit,opened=opened)
    end do
  end function get_unused_unit
  !---------------------------------------------------------------------------
end module misc_fileops
