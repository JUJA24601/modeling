#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine extfld(ps)
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   E X T F L D
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   ............................................................
!
!-------------------- parameters and variables
  use oh_type
  use paramt
  use allcom
! #define MCW MPI_COMM_WORLD
#define MCW CTCA_subcomm
#define MSS MPI_STATUS_SIZE
  implicit none
!
  integer(kind=4) :: i,j,k, iw
  integer(kind=4) :: xl,xu, yl,yu, zl,zu
  integer(kind=4) :: ps
  real(kind=8) :: progrs
  real(kind=8) :: x,y,z
  real(kind=8) :: orge(3), orgb(3)


!-------------------- 
      xl = sdoms(1,1,sdid(ps)+1); xu = sdoms(2,1,sdid(ps)+1)
      yl = sdoms(1,2,sdid(ps)+1); yu = sdoms(2,2,sdid(ps)+1)
      zl = sdoms(1,3,sdid(ps)+1); zu = sdoms(2,3,sdid(ps)+1)


!-------------------- zero clear
      do k=-1,zu-zl+1
      do j=-1,yu-yl+1
      do i=-1,xu-xl+1
        eb(EX:EZ,i,j,k,nebfld*2+ps) = (/0.0d0,0.0d0,0.0d0/)
      end do
      end do
      end do


!-------------------- 
      progrs = cv*2.0d0*istep
      do iw=1,nwave
        if(twave(iw).eq.1) then
          orge(1) = orgwv(1,iw) + progrs*sin(angwv(1,iw))*cos(angwv(2,iw))
          orge(2) = orgwv(2,iw) + progrs*sin(angwv(1,iw))*sin(angwv(2,iw))
          orge(3) = orgwv(3,iw) + progrs*cos(angwv(1,iw))
          orgb(1) = orgwv(1,iw) + progrs*sin(angwv(1,iw))*cos(angwv(2,iw))
          orgb(2) = orgwv(2,iw) + progrs*sin(angwv(1,iw))*sin(angwv(2,iw))
          orgb(3) = orgwv(3,iw) + progrs*cos(angwv(1,iw))
!
          do k=-1,zu-zl+1
          do j=-1,yu-yl+1
          do i=-1,xu-xl+1
!           --------- E-components
            x = (i + xl + 0.5d0)*dr; y = (j + yl)*dr; z = (k + zl)*dr
            x = modulo(x,slx); y = modulo(y,sly); z = modulo(z,slz)
            eb(EX,i,j,k,nebfld*2+ps) = &
           &  ebamp(EX,iw)*sinfld(x,y,z,ldwv(1,iw),orge,angwv(:,iw))
            x = (i + xl)*dr; y = (j + yl + 0.5d0)*dr; z = (k + zl)*dr
            x = modulo(x,slx); y = modulo(y,sly); z = modulo(z,slz)
            eb(EY,i,j,k,nebfld*2+ps) = &
           &  ebamp(EY,iw)*sinfld(x,y,z,ldwv(1,iw),orge,angwv(:,iw))
            x = (i + xl)*dr; y = (j + yl)*dr; z = (k + zl + 0.5d0)*dr
            x = modulo(x,slx); y = modulo(y,sly); z = modulo(z,slz)
            eb(EZ,i,j,k,nebfld*2+ps) = &
           &  ebamp(EZ,iw)*sinfld(x,y,z,ldwv(1,iw),orge,angwv(:,iw))
!           --------- B-components
            x = (i + xl)*dr; y = (j + yl + 0.5d0)*dr; z = (k + zl + 0.5d0)*dr
            x = modulo(x,slx); y = modulo(y,sly); z = modulo(z,slz)
            eb(BX,i,j,k,nebfld*2+ps) = &
           &  ebamp(BX,iw)*sinfld(x,y,z,ldwv(1,iw),orgb,angwv(:,iw))
            x = (i + xl + 0.5d0)*dr; y = (j + yl)*dr; z = (k + zl + 0.5d0)*dr
            x = modulo(x,slx); y = modulo(y,sly); z = modulo(z,slz)
            eb(BY,i,j,k,nebfld*2+ps) = &
           &  ebamp(BY,iw)*sinfld(x,y,z,ldwv(1,iw),orgb,angwv(:,iw))
            x = (i + xl + 0.5d0)*dr; y = (j + yl + 0.5d0)*dr; z = (k + zl)*dr
            x = modulo(x,slx); y = modulo(y,sly); z = modulo(z,slz)
            eb(BZ,i,j,k,nebfld*2+ps) = &
           &  ebamp(BZ,iw)*sinfld(x,y,z,ldwv(1,iw),orgb,angwv(:,iw))
          end do
          end do
          end do
        else if(twave(iw).eq.2) then
          orge(1) = orgwv(1,iw) + progrs*sin(angwv(1,iw))*cos(angwv(2,iw))
          orge(2) = orgwv(2,iw) + progrs*sin(angwv(1,iw))*sin(angwv(2,iw))
          orge(3) = orgwv(3,iw) + progrs*cos(angwv(1,iw))
          orgb(1) = orgwv(1,iw) + progrs*sin(angwv(1,iw))*cos(angwv(2,iw))
          orgb(2) = orgwv(2,iw) + progrs*sin(angwv(1,iw))*sin(angwv(2,iw))
          orgb(3) = orgwv(3,iw) + progrs*cos(angwv(1,iw))
!
          do k=-1,zu-zl+1
          do j=-1,yu-yl+1
          do i=-1,xu-xl+1
!           --------- E-components
            x = (i + xl + 0.5d0)*dr; y = (j + yl)*dr; z = (k + zl)*dr
            x = modulo(x,slx); y = modulo(y,sly); z = modulo(z,slz)
            eb(EX,i,j,k,nebfld*2+ps) = &
           &  ebamp(EX,iw)*gaus1dA(x,y,z,ldwv(1,iw),orge,angwv(:,iw))
            x = (i + xl)*dr; y = (j + yl + 0.5d0)*dr; z = (k + zl)*dr
            x = modulo(x,slx); y = modulo(y,sly); z = modulo(z,slz)
            eb(EY,i,j,k,nebfld*2+ps) = &
           &  ebamp(EY,iw)*gaus1dA(x,y,z,ldwv(1,iw),orge,angwv(:,iw))
            x = (i + xl)*dr; y = (j + yl)*dr; z = (k + zl + 0.5d0)*dr
            x = modulo(x,slx); y = modulo(y,sly); z = modulo(z,slz)
            eb(EZ,i,j,k,nebfld*2+ps) = &
           &  ebamp(EZ,iw)*gaus1dA(x,y,z,ldwv(1,iw),orge,angwv(:,iw))
!           --------- B-components
            x = (i + xl)*dr; y = (j + yl + 0.5d0)*dr; z = (k + zl + 0.5d0)*dr
            x = modulo(x,slx); y = modulo(y,sly); z = modulo(z,slz)
            eb(BX,i,j,k,nebfld*2+ps) = &
           &  ebamp(BX,iw)*gaus1dA(x,y,z,ldwv(1,iw),orgb,angwv(:,iw))
            x = (i + xl + 0.5d0)*dr; y = (j + yl)*dr; z = (k + zl + 0.5d0)*dr
            x = modulo(x,slx); y = modulo(y,sly); z = modulo(z,slz)
            eb(BY,i,j,k,nebfld*2+ps) = &
           &  ebamp(BY,iw)*gaus1dA(x,y,z,ldwv(1,iw),orgb,angwv(:,iw))
            x = (i + xl + 0.5d0)*dr; y = (j + yl + 0.5d0)*dr; z = (k + zl)*dr
            x = modulo(x,slx); y = modulo(y,sly); z = modulo(z,slz)
            eb(BZ,i,j,k,nebfld*2+ps) = &
           &  ebamp(BZ,iw)*gaus1dA(x,y,z,ldwv(1,iw),orgb,angwv(:,iw))
          end do
          end do
          end do
        else if(twave(iw).eq.3) then
          orge(1) = orgwv(1,iw) + progrs*sin(angwv(1,iw))*cos(angwv(2,iw))
          orge(2) = orgwv(2,iw) + progrs*sin(angwv(1,iw))*sin(angwv(2,iw))
          orge(3) = orgwv(3,iw) + progrs*cos(angwv(1,iw))
          orgb(1) = orgwv(1,iw) + progrs*sin(angwv(1,iw))*cos(angwv(2,iw))
          orgb(2) = orgwv(2,iw) + progrs*sin(angwv(1,iw))*sin(angwv(2,iw))
          orgb(3) = orgwv(3,iw) + progrs*cos(angwv(1,iw))
!
          do k=-1,zu-zl+1
          do j=-1,yu-yl+1
          do i=-1,xu-xl+1
!           --------- E-components
            x = (i + xl + 0.5d0)*dr; y = (j + yl)*dr; z = (k + zl)*dr
            x = modulo(x,slx); y = modulo(y,sly); z = modulo(z,slz)
            eb(EX,i,j,k,nebfld*2+ps) = &
           &  ebamp(EX,iw)*gaus1dB(x,y,z,ldwv(1,iw),orge,angwv(:,iw))
            x = (i + xl)*dr; y = (j + yl + 0.5d0)*dr; z = (k + zl)*dr
            x = modulo(x,slx); y = modulo(y,sly); z = modulo(z,slz)
            eb(EY,i,j,k,nebfld*2+ps) = &
           &  ebamp(EY,iw)*gaus1dB(x,y,z,ldwv(1,iw),orge,angwv(:,iw))
            x = (i + xl)*dr; y = (j + yl)*dr; z = (k + zl + 0.5d0)*dr
            x = modulo(x,slx); y = modulo(y,sly); z = modulo(z,slz)
            eb(EZ,i,j,k,nebfld*2+ps) = &
           &  ebamp(EZ,iw)*gaus1dB(x,y,z,ldwv(1,iw),orge,angwv(:,iw))
!           --------- B-components
            x = (i + xl)*dr; y = (j + yl + 0.5d0)*dr; z = (k + zl + 0.5d0)*dr
            x = modulo(x,slx); y = modulo(y,sly); z = modulo(z,slz)
            eb(BX,i,j,k,nebfld*2+ps) = &
           &  ebamp(BX,iw)*gaus1dB(x,y,z,ldwv(1,iw),orgb,angwv(:,iw))
            x = (i + xl + 0.5d0)*dr; y = (j + yl)*dr; z = (k + zl + 0.5d0)*dr
            x = modulo(x,slx); y = modulo(y,sly); z = modulo(z,slz)
            eb(BY,i,j,k,nebfld*2+ps) = &
           &  ebamp(BY,iw)*gaus1dB(x,y,z,ldwv(1,iw),orgb,angwv(:,iw))
            x = (i + xl + 0.5d0)*dr; y = (j + yl + 0.5d0)*dr; z = (k + zl)*dr
            x = modulo(x,slx); y = modulo(y,sly); z = modulo(z,slz)
            eb(BZ,i,j,k,nebfld*2+ps) = &
           &  ebamp(BZ,iw)*gaus1dB(x,y,z,ldwv(1,iw),orgb,angwv(:,iw))
          end do
          end do
          end do
        else if(twave(iw).eq.4) then
          orge(1) = orgwv(1,iw) + progrs*sin(angwv(1,iw))*cos(angwv(2,iw))
          orge(2) = orgwv(2,iw) + progrs*sin(angwv(1,iw))*sin(angwv(2,iw))
          orge(3) = orgwv(3,iw) + progrs*cos(angwv(1,iw))
          orgb(1) = orgwv(1,iw) + progrs*sin(angwv(1,iw))*cos(angwv(2,iw))
          orgb(2) = orgwv(2,iw) + progrs*sin(angwv(1,iw))*sin(angwv(2,iw))
          orgb(3) = orgwv(3,iw) + progrs*cos(angwv(1,iw))
!
          do k=-1,zu-zl+1
          do j=-1,yu-yl+1
          do i=-1,xu-xl+1
!           --------- E-components
            x = (i + xl + 0.5d0)*dr; y = (j + yl)*dr; z = (k + zl)*dr
            x = modulo(x,slx); y = modulo(y,sly); z = modulo(z,slz)
            eb(EX,i,j,k,nebfld*2+ps) = &
           &  ebamp(EX,iw)*gssinfld(x,y,z,ldwv(1,iw),ldwv(2,iw),orge,angwv(:,iw))
            x = (i + xl)*dr; y = (j + yl + 0.5d0)*dr; z = (k + zl)*dr
            x = modulo(x,slx); y = modulo(y,sly); z = modulo(z,slz)
            eb(EY,i,j,k,nebfld*2+ps) = &
           &  ebamp(EY,iw)*gssinfld(x,y,z,ldwv(1,iw),ldwv(2,iw),orge,angwv(:,iw))
            x = (i + xl)*dr; y = (j + yl)*dr; z = (k + zl + 0.5d0)*dr
            x = modulo(x,slx); y = modulo(y,sly); z = modulo(z,slz)
            eb(EZ,i,j,k,nebfld*2+ps) = &
           &  ebamp(EZ,iw)*gssinfld(x,y,z,ldwv(1,iw),ldwv(2,iw),orge,angwv(:,iw))
!           --------- B-components
            x = (i + xl)*dr; y = (j + yl + 0.5d0)*dr; z = (k + zl + 0.5d0)*dr
            x = modulo(x,slx); y = modulo(y,sly); z = modulo(z,slz)
            eb(BX,i,j,k,nebfld*2+ps) = &
           &  ebamp(BX,iw)*gssinfld(x,y,z,ldwv(1,iw),ldwv(2,iw),orgb,angwv(:,iw))
            x = (i + xl + 0.5d0)*dr; y = (j + yl)*dr; z = (k + zl + 0.5d0)*dr
            x = modulo(x,slx); y = modulo(y,sly); z = modulo(z,slz)
            eb(BY,i,j,k,nebfld*2+ps) = &
           &  ebamp(BY,iw)*gssinfld(x,y,z,ldwv(1,iw),ldwv(2,iw),orgb,angwv(:,iw))
            x = (i + xl + 0.5d0)*dr; y = (j + yl + 0.5d0)*dr; z = (k + zl)*dr
            x = modulo(x,slx); y = modulo(y,sly); z = modulo(z,slz)
            eb(BZ,i,j,k,nebfld*2+ps) = &
           &  ebamp(BZ,iw)*gssinfld(x,y,z,ldwv(1,iw),ldwv(2,iw),orgb,angwv(:,iw))
          end do
          end do
          end do
        end if
      end do


  return


  contains


  function sinfld(x,y,z,ld,org,ang)
    implicit none
    real(kind=8) :: sinfld
    real(kind=8),intent(in) :: x,y,z,ld,org(:),ang(:)
    real(kind=8) :: displ

    displ = (x - org(1))*sin(ang(1))*cos(ang(2)) &
   &      + (y - org(2))*sin(ang(1))*sin(ang(2)) &
   &      + (z - org(3))*cos(ang(1))

    sinfld = sin(2.0d0*pi*displ/ld)

    return
  end function sinfld


  function gaus1dA(x,y,z,pw,org,ang)
    implicit none
    real(kind=8) :: gaus1dA
    real(kind=8),intent(in) :: x,y,z,pw,org(:),ang(:)
    real(kind=8) :: displ

    displ = (x - org(1))*sin(ang(1))*cos(ang(2)) &
   &      + (y - org(2))*sin(ang(1))*sin(ang(2)) &
   &      + (z - org(3))*cos(ang(1))

    gaus1dA = exp(-displ*displ/pw/pw)

    return
  end function gaus1dA


  function gaus1dB(x,y,z,pw,org,ang)
    implicit none
    real(kind=8) :: gaus1dB
    real(kind=8),intent(in) :: x,y,z,pw,org(:),ang(:)
    real(kind=8) :: tpw, ipw, displ

    displ = (x - org(1))*sin(ang(1))*cos(ang(2)) &
   &      + (y - org(2))*sin(ang(1))*sin(ang(2)) &
   &      + (z - org(3))*cos(ang(1))

    tpw = pw*0.152d0
    ipw = 1.0d0/pw
    displ = -displ*ipw

    if(displ.gt.0.0d0) then
      gaus1dB = ipw*displ*displ*displ*exp(-displ) &
     &         *(4.0d0 - displ)*tpw*0.0625d0/exp(-2.0d0)
    else
      gaus1dB = 0.0d0
    end if

    return
  end function gaus1dB


  function gssinfld(x,y,z,pw,ld,org,ang)
    implicit none
    real(kind=8) :: gssinfld
    real(kind=8),intent(in) :: x,y,z,pw,ld,org(:),ang(:)
    real(kind=8) :: displ

    displ = (x - org(1))*sin(ang(1))*cos(ang(2)) &
   &      + (y - org(2))*sin(ang(1))*sin(ang(2)) &
   &      + (z - org(3))*cos(ang(1))

    gssinfld = exp(-displ*displ/pw/pw)*sin(2.0d0*pi*displ/ld)

    return
  end function gssinfld


  end subroutine extfld
