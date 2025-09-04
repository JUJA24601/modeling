program main
  implicit none
  integer(kind=4) :: emflag = 1
  character(len=20) :: inputfilename
  namelist /esorem/ emflag

  call get_command_argument(1,inputfilename)

  open(1,file=inputfilename)
  read(1,esorem)
  close(1)

  if(emflag.eq.0) then
    call esses
  elseif(emflag.eq.1) then
    call emses
  else if(emflag.eq.2) then
    call imses
  end if

  stop
end program main
