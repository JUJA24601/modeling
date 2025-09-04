#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine vpushh(ps,ustep,dir,mode,order)
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   V P U S H H
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   ............................................................
!
!
!-------------------- parameter and common block
  use oh_type
  use paramt
  use allcom
  implicit none
!
  integer(kind=8) :: m, ns,ne, nee
  integer(kind=4) :: i,j,k
  integer(kind=4) :: i1,j1,k1
  integer(kind=4) :: is, ustep, itch
  integer(kind=4) :: xl,xu, yl,yu, zl,zu
  integer(kind=4) :: ps
  integer(kind=4) :: mode, order, dir
  real(kind=8) :: dxl,dxu, dyl,dyu, dzl,dzu
  real(kind=8) :: xlocal,ylocal,zlocal
  real(kind=8) :: x1,y1,z1, z2, xy1,xz1,yz1, xz2,yz2
  real(kind=8) :: gustep,qmp
  real(kind=8) :: v1,v2,v3,v4,v5,v6,v7,v8
  real(kind=8) :: eex,eey,eez, bbx,bby,bbz, detinv, vvt(3)
  real(kind=8) :: ewx,ewy,tew
  real(kind=8) :: displx,disply,displz,displ3i
  logical :: explic,forwar


!-------------------- 
      xl = sdoms(1,1,sdid(ps)+1); xu = sdoms(2,1,sdid(ps)+1)
      yl = sdoms(1,2,sdid(ps)+1); yu = sdoms(2,2,sdid(ps)+1)
      zl = sdoms(1,3,sdid(ps)+1); zu = sdoms(2,3,sdid(ps)+1)
      dxl = xl; dxu = xu; dyl = yl; dyu = yu; dzl = zl; dzu = zu


!-------------------- 
      if(mode.le.0) then
        explic = .true.
      else
        explic = .false.
      end if
      if(order.le.0) then
        forwar = .true.
      else
        forwar = .false.
      end if


!============================== species loop
      nee = pbase(ps)
ISL1: do is=1,nspec
        ns = nee + 1
        ne = nee + totalp(is,ps)
        nee = ne
!       ------------- 
        if(ewmodel.eq.1.or.ewmodel.eq.2) then
          tew = t - dt*nretard
          ewx = Ew(1,is,1)*cos(omegaw(1)*tew)
          ewy = Ew(2,is,1)*sin(omegaw(1)*tew)
        else
          ewx = 0.0d0
          ewy = 0.0d0
        end if
        gustep = ustep
        qmp = dir*qm(is)*gustep*0.5d0
!       ------------- inner loop
        do m=ns,ne
          if(pbuf(m)%nid.eq.-1) cycle
!
          xlocal = pbuf(m)%x - dxl
          ylocal = pbuf(m)%y - dyl
          zlocal = pbuf(m)%z - dzl
!
          i = int(xlocal)
          j = int(ylocal)
          k = int(zlocal)
!
          i1 = i + 1
          j1 = j + 1
          k1 = k + 1
!
          x1 = xlocal - i
          y1 = ylocal - j
          z1 = zlocal - k
!
          xy1 = x1*y1
          xz1 = x1*z1
          yz1 = y1*z1
          z2 = 1.0d0 - z1
          xz2 = x1*z2
          yz2 = y1*z2
          v3 = xy1*z1
          v2 = xz1 - v3
          v4 = yz1 - v3
          v1 = z1 - xz1 - v4
          v7 = xy1*z2
          v6 = xz2 - v7
          v8 = yz2 - v7
          v5 = z2 - xz2 - v8
!
!         ----------- eb-fields interpolation
          eex = wrk(EX,i ,j ,k1)*v1 + wrk(EX,i1,j ,k1)*v2 &
         &    + wrk(EX,i1,j1,k1)*v3 + wrk(EX,i ,j1,k1)*v4 &
         &    + wrk(EX,i ,j ,k )*v5 + wrk(EX,i1,j ,k )*v6 &
         &    + wrk(EX,i1,j1,k )*v7 + wrk(EX,i ,j1,k )*v8
          eey = wrk(EY,i ,j ,k1)*v1 + wrk(EY,i1,j ,k1)*v2 &
         &    + wrk(EY,i1,j1,k1)*v3 + wrk(EY,i ,j1,k1)*v4 &
         &    + wrk(EY,i ,j ,k )*v5 + wrk(EY,i1,j ,k )*v6 &
         &    + wrk(EY,i1,j1,k )*v7 + wrk(EY,i ,j1,k )*v8
          eez = wrk(EZ,i ,j ,k1)*v1 + wrk(EZ,i1,j ,k1)*v2 &
         &    + wrk(EZ,i1,j1,k1)*v3 + wrk(EZ,i ,j1,k1)*v4 &
         &    + wrk(EZ,i ,j ,k )*v5 + wrk(EZ,i1,j ,k )*v6 &
         &    + wrk(EZ,i1,j1,k )*v7 + wrk(EZ,i ,j1,k )*v8
          bbx = wrk(BX,i ,j ,k1)*v1 + wrk(BX,i1,j ,k1)*v2 &
         &    + wrk(BX,i1,j1,k1)*v3 + wrk(BX,i ,j1,k1)*v4 &
         &    + wrk(BX,i ,j ,k )*v5 + wrk(BX,i1,j ,k )*v6 &
         &    + wrk(BX,i1,j1,k )*v7 + wrk(BX,i ,j1,k )*v8
          bby = wrk(BY,i ,j ,k1)*v1 + wrk(BY,i1,j ,k1)*v2 &
         &    + wrk(BY,i1,j1,k1)*v3 + wrk(BY,i ,j1,k1)*v4 &
         &    + wrk(BY,i ,j ,k )*v5 + wrk(BY,i1,j ,k )*v6 &
         &    + wrk(BY,i1,j1,k )*v7 + wrk(BY,i ,j1,k )*v8
          bbz = wrk(BZ,i ,j ,k1)*v1 + wrk(BZ,i1,j ,k1)*v2 &
         &    + wrk(BZ,i1,j1,k1)*v3 + wrk(BZ,i ,j1,k1)*v4 &
         &    + wrk(BZ,i ,j ,k )*v5 + wrk(BZ,i1,j ,k )*v6 &
         &    + wrk(BZ,i1,j1,k )*v7 + wrk(BZ,i ,j1,k )*v8
!
!         ----------- 
          do itch=1,ntch
            displx = pbuf(m)%x - rtch(1,itch)
            disply = pbuf(m)%y - rtch(2,itch)
            displz = pbuf(m)%z - rtch(3,itch)
            displ3i = displx*displx + disply*disply + displz*displz
            displ3i = 1.0d0/((displ3i + r2cutoff(itch))*sqrt(displ3i))
            eex = eex + e1tch(itch)*displx*displ3i
            eey = eey + e1tch(itch)*disply*displ3i
            eez = eez + e1tch(itch)*displz*displ3i
          end do
!
!         ----------- charge-to-mass ratio is taken into consideration
          eex = (eex + ewx)*qmp
          eey = (eey + ewy)*qmp
          eez = eez*qmp
          bbx = bbx*qmp
          bby = bby*qmp
          bbz = bbz*qmp
!
!         ----------- update particle velocities
          vvt(1) = pbuf(m)%vx
          vvt(2) = pbuf(m)%vy
          vvt(3) = pbuf(m)%vz
          if(explic) then
            if(forwar) then
              vvt = eaccel(vvt(1),vvt(2),vvt(3),eex,eey,eez)
              vvt = gyroex(vvt(1),vvt(2),vvt(3),bbx,bby,bbz)
            else
              vvt = gyroex(vvt(1),vvt(2),vvt(3),bbx,bby,bbz)
              vvt = eaccel(vvt(1),vvt(2),vvt(3),eex,eey,eez)
            end if
          else
            detinv = 1.0d0/(1.0d0 + bbx*bbx + bby*bby + bbz*bbz)
            if(forwar) then
              vvt = eaccel(vvt(1),vvt(2),vvt(3),eex,eey,eez)
              vvt = gyroim(vvt(1),vvt(2),vvt(3),bbx,bby,bbz,detinv)
            else
              vvt = gyroim(vvt(1),vvt(2),vvt(3),bbx,bby,bbz,detinv)
              vvt = eaccel(vvt(1),vvt(2),vvt(3),eex,eey,eez)
            end if
          end if
          pbuf(m)%vx = vvt(1)
          pbuf(m)%vy = vvt(2)
          pbuf(m)%vz = vvt(3)
        end do
      end do ISL1

  return

  contains

    function eaccel(vvx,vvy,vvz,etx,ety,etz) result(vv)
      implicit none
      real(kind=8),intent(in) :: vvx,vvy,vvz
      real(kind=8),intent(in) :: etx,ety,etz
      real(kind=8) :: vv(3)

      vv(1) = vvx + etx
      vv(2) = vvy + ety
      vv(3) = vvz + etz

      return
    end function eaccel

    function gyroex(vvx,vvy,vvz,btx,bty,btz) result(vv)
      implicit none
      real(kind=8),intent(in) :: vvx,vvy,vvz
      real(kind=8),intent(in) :: btx,bty,btz
      real(kind=8) :: vv(3)

      vv(1) = vvx + vvy*btz - vvz*bty
      vv(2) = vvy + vvz*btx - vvx*btz
      vv(3) = vvz + vvx*bty - vvy*btx

      return
    end function gyroex

    function gyroim(vvx,vvy,vvz,btx,bty,btz,deti) result(vv)
      implicit none
      real(kind=8),intent(in) :: vvx,vvy,vvz
      real(kind=8),intent(in) :: btx,bty,btz
      real(kind=8),intent(in) :: deti
      real(kind=8) :: vv(3)

      vv(1) = (1.0d0 + btx*btx)*vvx &
     &      + (btx*bty + btz)*vvy &
     &      + (btz*btx - bty)*vvz
      vv(2) = (btx*bty - btz)*vvx &
     &      + (1.0d0 + bty*bty)*vvy &
     &      + (bty*btz + btx)*vvz
      vv(3) = (btz*btx + bty)*vvx &
     &      + (bty*btz - btx)*vvy &
     &      + (1.0d0 + btz*btz)*vvz

      vv(1) = vv(1)*deti
      vv(2) = vv(2)*deti
      vv(3) = vv(3)*deti

      return
    end function gyroim

  end subroutine vpushh
