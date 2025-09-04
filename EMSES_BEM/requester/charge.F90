#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine charge(ps,func)
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   C H A R G E
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .      this subroutine gives a charge distribution to      .
!   .      grid points from particle locations by the          .
!   .      "area-sharing" or "area-weighting" scheme.          .
!   ............................................................
!
!-------------------- parameter and common block
  use oh_type
  use paramt
  use allcom
  implicit none
!
  integer(kind=8) :: m, ns,ne, nee
  integer(kind=4) :: i,j,k, i1,j1,k1, is, nis
  integer(kind=4) :: ngx,ngy,ngz
  integer(kind=4) :: func
  integer(kind=4) :: ps
  real(kind=8) :: xl,xu, yl,yu, zl,zu
  real(kind=8) :: xlocal,ylocal,zlocal
  real(kind=8) :: x1,y1,z1, z2, xy1,xz1,yz1, xz2,yz2
  real(kind=8) :: v1,v2,v3,v4,v5,v6,v7,v8
!  character(len=8) :: finpfn


!-------------------- 
      xl = sdoms(1,1,sdid(ps)+1); xu = sdoms(2,1,sdid(ps)+1)
      yl = sdoms(1,2,sdid(ps)+1); yu = sdoms(2,2,sdid(ps)+1)
      zl = sdoms(1,3,sdid(ps)+1); zu = sdoms(2,3,sdid(ps)+1)
      ngx = sdoms(2,1,sdid(ps)+1) - sdoms(1,1,sdid(ps)+1)
      ngy = sdoms(2,2,sdid(ps)+1) - sdoms(1,2,sdid(ps)+1)
      ngz = sdoms(2,3,sdid(ps)+1) - sdoms(1,3,sdid(ps)+1)


!-------------------- for three dimensional system
!     --------------- zero clear rho
      if(ps.eq.1) then
        rho(1:1,:,:,:,1:2) = 0.0d0
        rhobk(1:2,:,:,:,1:2) = 0.0d0
      end if
      if(func.eq.1) &
     &  rhodg(:,:,:,:,ps) = 0.0d0


!-------------------- sum up charge component
!============================== species loop
      nee = pbase(ps)
 ISL: do is=1,nspec
!-------------------- inner loop
        ns = nee + 1
        ne = nee + totalp(is,ps)
        nee = ne
        if(func.eq.0) then
          do m=ns,ne
            if(pbuf(m)%nid.eq.-1) cycle
!
            xlocal = pbuf(m)%x - xl
            ylocal = pbuf(m)%y - yl
            zlocal = pbuf(m)%z - zl
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
!           --------- temporary value
            xy1 = x1 * y1
            xz1 = x1 * z1
            yz1 = y1 * z1
            z2 = 1.0d0 - z1
            xz2 = x1 * z2
            yz2 = y1 * z2
!
            v3 = xy1 * z1
            v2 = xz1 - v3
            v4 = yz1 - v3
            v1 = z1 - xz1 - v4
!
            v7 = xy1 * z2
            v6 = xz2 - v7
            v8 = yz2 - v7
            v5 = z2 - xz2 - v8
!           --------- distribute charge with weight
            if(pbuf(m)%preside.ge.0) then
              rho(1,i ,j ,k1,ps) = rho(1,i ,j ,k1,ps) + v1*q(is)
              rho(1,i1,j ,k1,ps) = rho(1,i1,j ,k1,ps) + v2*q(is)
              rho(1,i1,j1,k1,ps) = rho(1,i1,j1,k1,ps) + v3*q(is)
              rho(1,i ,j1,k1,ps) = rho(1,i ,j1,k1,ps) + v4*q(is)
              rho(1,i ,j ,k ,ps) = rho(1,i ,j ,k ,ps) + v5*q(is)
              rho(1,i1,j ,k ,ps) = rho(1,i1,j ,k ,ps) + v6*q(is)
              rho(1,i1,j1,k ,ps) = rho(1,i1,j1,k ,ps) + v7*q(is)
              rho(1,i ,j1,k ,ps) = rho(1,i ,j1,k ,ps) + v8*q(is)
            else if(pbuf(m)%preside.eq.IEBD.or. &
           &        pbuf(m)%preside.eq.IPCB) then
              rhobk(2,i ,j ,k1,ps) = rhobk(2,i ,j ,k1,ps) - v1*q(is)
              rhobk(2,i1,j ,k1,ps) = rhobk(2,i1,j ,k1,ps) - v2*q(is)
              rhobk(2,i1,j1,k1,ps) = rhobk(2,i1,j1,k1,ps) - v3*q(is)
              rhobk(2,i ,j1,k1,ps) = rhobk(2,i ,j1,k1,ps) - v4*q(is)
              rhobk(2,i ,j ,k ,ps) = rhobk(2,i ,j ,k ,ps) - v5*q(is)
              rhobk(2,i1,j ,k ,ps) = rhobk(2,i1,j ,k ,ps) - v6*q(is)
              rhobk(2,i1,j1,k ,ps) = rhobk(2,i1,j1,k ,ps) - v7*q(is)
              rhobk(2,i ,j1,k ,ps) = rhobk(2,i ,j1,k ,ps) - v8*q(is)
!              rhobk(1:2,i ,j ,k1,ps) = rhobk(1:2,i ,j ,k1,ps) - v1*q(is)
!              rhobk(1:2,i1,j ,k1,ps) = rhobk(1:2,i1,j ,k1,ps) - v2*q(is)
!              rhobk(1:2,i1,j1,k1,ps) = rhobk(1:2,i1,j1,k1,ps) - v3*q(is)
!              rhobk(1:2,i ,j1,k1,ps) = rhobk(1:2,i ,j1,k1,ps) - v4*q(is)
!              rhobk(1:2,i ,j ,k ,ps) = rhobk(1:2,i ,j ,k ,ps) - v5*q(is)
!              rhobk(1:2,i1,j ,k ,ps) = rhobk(1:2,i1,j ,k ,ps) - v6*q(is)
!              rhobk(1:2,i1,j1,k ,ps) = rhobk(1:2,i1,j1,k ,ps) - v7*q(is)
!              rhobk(1:2,i ,j1,k ,ps) = rhobk(1:2,i ,j1,k ,ps) - v8*q(is)
              nphgram(pbuf(m)%nid+1,is,ps) = nphgram(pbuf(m)%nid+1,is,ps) - 1
              pbuf(m)%nid = -1
              pbuf(m)%preside = 0
            else if(pbuf(m)%preside.eq.INCB) then
              rhobk(1,i ,j ,k1,ps) = rhobk(1,i ,j ,k1,ps) + v1*q(is)
              rhobk(1,i1,j ,k1,ps) = rhobk(1,i1,j ,k1,ps) + v2*q(is)
              rhobk(1,i1,j1,k1,ps) = rhobk(1,i1,j1,k1,ps) + v3*q(is)
              rhobk(1,i ,j1,k1,ps) = rhobk(1,i ,j1,k1,ps) + v4*q(is)
              rhobk(1,i ,j ,k ,ps) = rhobk(1,i ,j ,k ,ps) + v5*q(is)
              rhobk(1,i1,j ,k ,ps) = rhobk(1,i1,j ,k ,ps) + v6*q(is)
              rhobk(1,i1,j1,k ,ps) = rhobk(1,i1,j1,k ,ps) + v7*q(is)
              rhobk(1,i ,j1,k ,ps) = rhobk(1,i ,j1,k ,ps) + v8*q(is)
              nphgram(pbuf(m)%nid+1,is,ps) = nphgram(pbuf(m)%nid+1,is,ps) - 1
              pbuf(m)%nid = -1
              pbuf(m)%preside = 0
            else if(pbuf(m)%preside.eq.OEBD.or. &
           &        pbuf(m)%preside.eq.OPCB) then
              rhobk(2,i ,j ,k1,ps) = rhobk(2,i ,j ,k1,ps) + v1*q(is)
              rhobk(2,i1,j ,k1,ps) = rhobk(2,i1,j ,k1,ps) + v2*q(is)
              rhobk(2,i1,j1,k1,ps) = rhobk(2,i1,j1,k1,ps) + v3*q(is)
              rhobk(2,i ,j1,k1,ps) = rhobk(2,i ,j1,k1,ps) + v4*q(is)
              rhobk(2,i ,j ,k ,ps) = rhobk(2,i ,j ,k ,ps) + v5*q(is)
              rhobk(2,i1,j ,k ,ps) = rhobk(2,i1,j ,k ,ps) + v6*q(is)
              rhobk(2,i1,j1,k ,ps) = rhobk(2,i1,j1,k ,ps) + v7*q(is)
              rhobk(2,i ,j1,k ,ps) = rhobk(2,i ,j1,k ,ps) + v8*q(is)
!              rhobk(1:2,i ,j ,k1,ps) = rhobk(1:2,i ,j ,k1,ps) + v1*q(is)
!              rhobk(1:2,i1,j ,k1,ps) = rhobk(1:2,i1,j ,k1,ps) + v2*q(is)
!              rhobk(1:2,i1,j1,k1,ps) = rhobk(1:2,i1,j1,k1,ps) + v3*q(is)
!              rhobk(1:2,i ,j1,k1,ps) = rhobk(1:2,i ,j1,k1,ps) + v4*q(is)
!              rhobk(1:2,i ,j ,k ,ps) = rhobk(1:2,i ,j ,k ,ps) + v5*q(is)
!              rhobk(1:2,i1,j ,k ,ps) = rhobk(1:2,i1,j ,k ,ps) + v6*q(is)
!              rhobk(1:2,i1,j1,k ,ps) = rhobk(1:2,i1,j1,k ,ps) + v7*q(is)
!              rhobk(1:2,i ,j1,k ,ps) = rhobk(1:2,i ,j1,k ,ps) + v8*q(is)
              nphgram(pbuf(m)%nid+1,is,ps) = nphgram(pbuf(m)%nid+1,is,ps) - 1
              pbuf(m)%nid = -1
              pbuf(m)%preside = 0
            else if(pbuf(m)%preside.eq.ONCB) then
              rhobk(1,i ,j ,k1,ps) = rhobk(1,i ,j ,k1,ps) - v1*q(is)
              rhobk(1,i1,j ,k1,ps) = rhobk(1,i1,j ,k1,ps) - v2*q(is)
              rhobk(1,i1,j1,k1,ps) = rhobk(1,i1,j1,k1,ps) - v3*q(is)
              rhobk(1,i ,j1,k1,ps) = rhobk(1,i ,j1,k1,ps) - v4*q(is)
              rhobk(1,i ,j ,k ,ps) = rhobk(1,i ,j ,k ,ps) - v5*q(is)
              rhobk(1,i1,j ,k ,ps) = rhobk(1,i1,j ,k ,ps) - v6*q(is)
              rhobk(1,i1,j1,k ,ps) = rhobk(1,i1,j1,k ,ps) - v7*q(is)
              rhobk(1,i ,j1,k ,ps) = rhobk(1,i ,j1,k ,ps) - v8*q(is)
              nphgram(pbuf(m)%nid+1,is,ps) = nphgram(pbuf(m)%nid+1,is,ps) - 1
              pbuf(m)%nid = -1
              pbuf(m)%preside = 0
            end if
          end do
        else if(func.eq.1) then
          nis = nspec + is
          do m=ns,ne
            if(pbuf(m)%nid.eq.-1) cycle
!
            xlocal = pbuf(m)%x - xl
            ylocal = pbuf(m)%y - yl
            zlocal = pbuf(m)%z - zl
!
            i = int(xlocal)
            j = int(ylocal)
            k = int(zlocal)
!
            i1= i + 1
            j1= j + 1
            k1= k + 1
!
            x1 = xlocal - i
            y1 = ylocal - j
            z1 = zlocal - k
!           --------- temporary value
            xy1 = x1 * y1
            xz1 = x1 * z1
            yz1 = y1 * z1
            z2 = 1.0d0 - z1
            xz2 = x1 * z2
            yz2 = y1 * z2
!
            v3 = xy1 * z1
            v2 = xz1 - v3
            v4 = yz1 - v3
            v1 = z1 - xz1 - v4
!
            v7 = xy1 * z2
            v6 = xz2 - v7
            v8 = yz2 - v7
            v5 = z2 - xz2 - v8
!           --------- distribute charge with weight
            if(pbuf(m)%preside.ge.0) then
              rho(1,i ,j ,k1,ps) = rho(1,i ,j ,k1,ps) + v1*q(is)
              rho(1,i1,j ,k1,ps) = rho(1,i1,j ,k1,ps) + v2*q(is)
              rho(1,i1,j1,k1,ps) = rho(1,i1,j1,k1,ps) + v3*q(is)
              rho(1,i ,j1,k1,ps) = rho(1,i ,j1,k1,ps) + v4*q(is)
              rho(1,i ,j ,k ,ps) = rho(1,i ,j ,k ,ps) + v5*q(is)
              rho(1,i1,j ,k ,ps) = rho(1,i1,j ,k ,ps) + v6*q(is)
              rho(1,i1,j1,k ,ps) = rho(1,i1,j1,k ,ps) + v7*q(is)
              rho(1,i ,j1,k ,ps) = rho(1,i ,j1,k ,ps) + v8*q(is)
              rhodg(is,i ,j ,k1,ps) = rhodg(is,i ,j ,k1,ps) + v1*q(is)
              rhodg(is,i1,j ,k1,ps) = rhodg(is,i1,j ,k1,ps) + v2*q(is)
              rhodg(is,i1,j1,k1,ps) = rhodg(is,i1,j1,k1,ps) + v3*q(is)
              rhodg(is,i ,j1,k1,ps) = rhodg(is,i ,j1,k1,ps) + v4*q(is)
              rhodg(is,i ,j ,k ,ps) = rhodg(is,i ,j ,k ,ps) + v5*q(is)
              rhodg(is,i1,j ,k ,ps) = rhodg(is,i1,j ,k ,ps) + v6*q(is)
              rhodg(is,i1,j1,k ,ps) = rhodg(is,i1,j1,k ,ps) + v7*q(is)
              rhodg(is,i ,j1,k ,ps) = rhodg(is,i ,j1,k ,ps) + v8*q(is)
            else if(pbuf(m)%preside.eq.IEBD) then
              rhobk(2,i ,j ,k1,ps) = rhobk(2,i ,j ,k1,ps) - v1*q(is)
              rhobk(2,i1,j ,k1,ps) = rhobk(2,i1,j ,k1,ps) - v2*q(is)
              rhobk(2,i1,j1,k1,ps) = rhobk(2,i1,j1,k1,ps) - v3*q(is)
              rhobk(2,i ,j1,k1,ps) = rhobk(2,i ,j1,k1,ps) - v4*q(is)
              rhobk(2,i ,j ,k ,ps) = rhobk(2,i ,j ,k ,ps) - v5*q(is)
              rhobk(2,i1,j ,k ,ps) = rhobk(2,i1,j ,k ,ps) - v6*q(is)
              rhobk(2,i1,j1,k ,ps) = rhobk(2,i1,j1,k ,ps) - v7*q(is)
              rhobk(2,i ,j1,k ,ps) = rhobk(2,i ,j1,k ,ps) - v8*q(is)
!              rhobk(1:2,i ,j ,k1,ps) = rhobk(1:2,i ,j ,k1,ps) - v1*q(is)
!              rhobk(1:2,i1,j ,k1,ps) = rhobk(1:2,i1,j ,k1,ps) - v2*q(is)
!              rhobk(1:2,i1,j1,k1,ps) = rhobk(1:2,i1,j1,k1,ps) - v3*q(is)
!              rhobk(1:2,i ,j1,k1,ps) = rhobk(1:2,i ,j1,k1,ps) - v4*q(is)
!              rhobk(1:2,i ,j ,k ,ps) = rhobk(1:2,i ,j ,k ,ps) - v5*q(is)
!              rhobk(1:2,i1,j ,k ,ps) = rhobk(1:2,i1,j ,k ,ps) - v6*q(is)
!              rhobk(1:2,i1,j1,k ,ps) = rhobk(1:2,i1,j1,k ,ps) - v7*q(is)
!              rhobk(1:2,i ,j1,k ,ps) = rhobk(1:2,i ,j1,k ,ps) - v8*q(is)
              nphgram(pbuf(m)%nid+1,is,ps) = nphgram(pbuf(m)%nid+1,is,ps) - 1
              pbuf(m)%nid = -1
              pbuf(m)%preside = 0
            else if(pbuf(m)%preside.eq.INCB) then
              rhobk(1,i ,j ,k1,ps) = rhobk(1,i ,j ,k1,ps) + v1*q(is)
              rhobk(1,i1,j ,k1,ps) = rhobk(1,i1,j ,k1,ps) + v2*q(is)
              rhobk(1,i1,j1,k1,ps) = rhobk(1,i1,j1,k1,ps) + v3*q(is)
              rhobk(1,i ,j1,k1,ps) = rhobk(1,i ,j1,k1,ps) + v4*q(is)
              rhobk(1,i ,j ,k ,ps) = rhobk(1,i ,j ,k ,ps) + v5*q(is)
              rhobk(1,i1,j ,k ,ps) = rhobk(1,i1,j ,k ,ps) + v6*q(is)
              rhobk(1,i1,j1,k ,ps) = rhobk(1,i1,j1,k ,ps) + v7*q(is)
              rhobk(1,i ,j1,k ,ps) = rhobk(1,i ,j1,k ,ps) + v8*q(is)
              rhodg(nis,i ,j ,k1,ps) = rhodg(nis,i ,j ,k1,ps) + v1*q(is)
              rhodg(nis,i1,j ,k1,ps) = rhodg(nis,i1,j ,k1,ps) + v2*q(is)
              rhodg(nis,i1,j1,k1,ps) = rhodg(nis,i1,j1,k1,ps) + v3*q(is)
              rhodg(nis,i ,j1,k1,ps) = rhodg(nis,i ,j1,k1,ps) + v4*q(is)
              rhodg(nis,i ,j ,k ,ps) = rhodg(nis,i ,j ,k ,ps) + v5*q(is)
              rhodg(nis,i1,j ,k ,ps) = rhodg(nis,i1,j ,k ,ps) + v6*q(is)
              rhodg(nis,i1,j1,k ,ps) = rhodg(nis,i1,j1,k ,ps) + v7*q(is)
              rhodg(nis,i ,j1,k ,ps) = rhodg(nis,i ,j1,k ,ps) + v8*q(is)
              nphgram(pbuf(m)%nid+1,is,ps) = nphgram(pbuf(m)%nid+1,is,ps) - 1
              pbuf(m)%nid = -1
              pbuf(m)%preside = 0
            else if(pbuf(m)%preside.eq.IPCB) then
              rhobk(2,i ,j ,k1,ps) = rhobk(2,i ,j ,k1,ps) - v1*q(is)
              rhobk(2,i1,j ,k1,ps) = rhobk(2,i1,j ,k1,ps) - v2*q(is)
              rhobk(2,i1,j1,k1,ps) = rhobk(2,i1,j1,k1,ps) - v3*q(is)
              rhobk(2,i ,j1,k1,ps) = rhobk(2,i ,j1,k1,ps) - v4*q(is)
              rhobk(2,i ,j ,k ,ps) = rhobk(2,i ,j ,k ,ps) - v5*q(is)
              rhobk(2,i1,j ,k ,ps) = rhobk(2,i1,j ,k ,ps) - v6*q(is)
              rhobk(2,i1,j1,k ,ps) = rhobk(2,i1,j1,k ,ps) - v7*q(is)
              rhobk(2,i ,j1,k ,ps) = rhobk(2,i ,j1,k ,ps) - v8*q(is)
!              rhobk(1:2,i ,j ,k1,ps) = rhobk(1:2,i ,j ,k1,ps) - v1*q(is)
!              rhobk(1:2,i1,j ,k1,ps) = rhobk(1:2,i1,j ,k1,ps) - v2*q(is)
!              rhobk(1:2,i1,j1,k1,ps) = rhobk(1:2,i1,j1,k1,ps) - v3*q(is)
!              rhobk(1:2,i ,j1,k1,ps) = rhobk(1:2,i ,j1,k1,ps) - v4*q(is)
!              rhobk(1:2,i ,j ,k ,ps) = rhobk(1:2,i ,j ,k ,ps) - v5*q(is)
!              rhobk(1:2,i1,j ,k ,ps) = rhobk(1:2,i1,j ,k ,ps) - v6*q(is)
!              rhobk(1:2,i1,j1,k ,ps) = rhobk(1:2,i1,j1,k ,ps) - v7*q(is)
!              rhobk(1:2,i ,j1,k ,ps) = rhobk(1:2,i ,j1,k ,ps) - v8*q(is)
              rhodg(nis,i ,j ,k1,ps) = rhodg(nis,i ,j ,k1,ps) + v1*q(is)
              rhodg(nis,i1,j ,k1,ps) = rhodg(nis,i1,j ,k1,ps) + v2*q(is)
              rhodg(nis,i1,j1,k1,ps) = rhodg(nis,i1,j1,k1,ps) + v3*q(is)
              rhodg(nis,i ,j1,k1,ps) = rhodg(nis,i ,j1,k1,ps) + v4*q(is)
              rhodg(nis,i ,j ,k ,ps) = rhodg(nis,i ,j ,k ,ps) + v5*q(is)
              rhodg(nis,i1,j ,k ,ps) = rhodg(nis,i1,j ,k ,ps) + v6*q(is)
              rhodg(nis,i1,j1,k ,ps) = rhodg(nis,i1,j1,k ,ps) + v7*q(is)
              rhodg(nis,i ,j1,k ,ps) = rhodg(nis,i ,j1,k ,ps) + v8*q(is)
              nphgram(pbuf(m)%nid+1,is,ps) = nphgram(pbuf(m)%nid+1,is,ps) - 1
              pbuf(m)%nid = -1
              pbuf(m)%preside = 0
            else if(pbuf(m)%preside.eq.OEBD) then
              rhobk(2,i ,j ,k1,ps) = rhobk(2,i ,j ,k1,ps) + v1*q(is)
              rhobk(2,i1,j ,k1,ps) = rhobk(2,i1,j ,k1,ps) + v2*q(is)
              rhobk(2,i1,j1,k1,ps) = rhobk(2,i1,j1,k1,ps) + v3*q(is)
              rhobk(2,i ,j1,k1,ps) = rhobk(2,i ,j1,k1,ps) + v4*q(is)
              rhobk(2,i ,j ,k ,ps) = rhobk(2,i ,j ,k ,ps) + v5*q(is)
              rhobk(2,i1,j ,k ,ps) = rhobk(2,i1,j ,k ,ps) + v6*q(is)
              rhobk(2,i1,j1,k ,ps) = rhobk(2,i1,j1,k ,ps) + v7*q(is)
              rhobk(2,i ,j1,k ,ps) = rhobk(2,i ,j1,k ,ps) + v8*q(is)
!              rhobk(1:2,i ,j ,k1,ps) = rhobk(1:2,i ,j ,k1,ps) + v1*q(is)
!              rhobk(1:2,i1,j ,k1,ps) = rhobk(1:2,i1,j ,k1,ps) + v2*q(is)
!              rhobk(1:2,i1,j1,k1,ps) = rhobk(1:2,i1,j1,k1,ps) + v3*q(is)
!              rhobk(1:2,i ,j1,k1,ps) = rhobk(1:2,i ,j1,k1,ps) + v4*q(is)
!              rhobk(1:2,i ,j ,k ,ps) = rhobk(1:2,i ,j ,k ,ps) + v5*q(is)
!              rhobk(1:2,i1,j ,k ,ps) = rhobk(1:2,i1,j ,k ,ps) + v6*q(is)
!              rhobk(1:2,i1,j1,k ,ps) = rhobk(1:2,i1,j1,k ,ps) + v7*q(is)
!              rhobk(1:2,i ,j1,k ,ps) = rhobk(1:2,i ,j1,k ,ps) + v8*q(is)
              nphgram(pbuf(m)%nid+1,is,ps) = nphgram(pbuf(m)%nid+1,is,ps) - 1
              pbuf(m)%nid = -1
              pbuf(m)%preside = 0
            else if(pbuf(m)%preside.eq.ONCB) then
              rhobk(1,i ,j ,k1,ps) = rhobk(1,i ,j ,k1,ps) - v1*q(is)
              rhobk(1,i1,j ,k1,ps) = rhobk(1,i1,j ,k1,ps) - v2*q(is)
              rhobk(1,i1,j1,k1,ps) = rhobk(1,i1,j1,k1,ps) - v3*q(is)
              rhobk(1,i ,j1,k1,ps) = rhobk(1,i ,j1,k1,ps) - v4*q(is)
              rhobk(1,i ,j ,k ,ps) = rhobk(1,i ,j ,k ,ps) - v5*q(is)
              rhobk(1,i1,j ,k ,ps) = rhobk(1,i1,j ,k ,ps) - v6*q(is)
              rhobk(1,i1,j1,k ,ps) = rhobk(1,i1,j1,k ,ps) - v7*q(is)
              rhobk(1,i ,j1,k ,ps) = rhobk(1,i ,j1,k ,ps) - v8*q(is)
              rhodg(nis,i ,j ,k1,ps) = rhodg(nis,i ,j ,k1,ps) - v1*q(is)
              rhodg(nis,i1,j ,k1,ps) = rhodg(nis,i1,j ,k1,ps) - v2*q(is)
              rhodg(nis,i1,j1,k1,ps) = rhodg(nis,i1,j1,k1,ps) - v3*q(is)
              rhodg(nis,i ,j1,k1,ps) = rhodg(nis,i ,j1,k1,ps) - v4*q(is)
              rhodg(nis,i ,j ,k ,ps) = rhodg(nis,i ,j ,k ,ps) - v5*q(is)
              rhodg(nis,i1,j ,k ,ps) = rhodg(nis,i1,j ,k ,ps) - v6*q(is)
              rhodg(nis,i1,j1,k ,ps) = rhodg(nis,i1,j1,k ,ps) - v7*q(is)
              rhodg(nis,i ,j1,k ,ps) = rhodg(nis,i ,j1,k ,ps) - v8*q(is)
              nphgram(pbuf(m)%nid+1,is,ps) = nphgram(pbuf(m)%nid+1,is,ps) - 1
              pbuf(m)%nid = -1
              pbuf(m)%preside = 0
            else if(pbuf(m)%preside.eq.OPCB) then
              rhobk(2,i ,j ,k1,ps) = rhobk(2,i ,j ,k1,ps) + v1*q(is)
              rhobk(2,i1,j ,k1,ps) = rhobk(2,i1,j ,k1,ps) + v2*q(is)
              rhobk(2,i1,j1,k1,ps) = rhobk(2,i1,j1,k1,ps) + v3*q(is)
              rhobk(2,i ,j1,k1,ps) = rhobk(2,i ,j1,k1,ps) + v4*q(is)
              rhobk(2,i ,j ,k ,ps) = rhobk(2,i ,j ,k ,ps) + v5*q(is)
              rhobk(2,i1,j ,k ,ps) = rhobk(2,i1,j ,k ,ps) + v6*q(is)
              rhobk(2,i1,j1,k ,ps) = rhobk(2,i1,j1,k ,ps) + v7*q(is)
              rhobk(2,i ,j1,k ,ps) = rhobk(2,i ,j1,k ,ps) + v8*q(is)
!              rhobk(1:2,i ,j ,k1,ps) = rhobk(1:2,i ,j ,k1,ps) + v1*q(is)
!              rhobk(1:2,i1,j ,k1,ps) = rhobk(1:2,i1,j ,k1,ps) + v2*q(is)
!              rhobk(1:2,i1,j1,k1,ps) = rhobk(1:2,i1,j1,k1,ps) + v3*q(is)
!              rhobk(1:2,i ,j1,k1,ps) = rhobk(1:2,i ,j1,k1,ps) + v4*q(is)
!              rhobk(1:2,i ,j ,k ,ps) = rhobk(1:2,i ,j ,k ,ps) + v5*q(is)
!              rhobk(1:2,i1,j ,k ,ps) = rhobk(1:2,i1,j ,k ,ps) + v6*q(is)
!              rhobk(1:2,i1,j1,k ,ps) = rhobk(1:2,i1,j1,k ,ps) + v7*q(is)
!              rhobk(1:2,i ,j1,k ,ps) = rhobk(1:2,i ,j1,k ,ps) + v8*q(is)
              rhodg(nis,i ,j ,k1,ps) = rhodg(nis,i ,j ,k1,ps) - v1*q(is)
              rhodg(nis,i1,j ,k1,ps) = rhodg(nis,i1,j ,k1,ps) - v2*q(is)
              rhodg(nis,i1,j1,k1,ps) = rhodg(nis,i1,j1,k1,ps) - v3*q(is)
              rhodg(nis,i ,j1,k1,ps) = rhodg(nis,i ,j1,k1,ps) - v4*q(is)
              rhodg(nis,i ,j ,k ,ps) = rhodg(nis,i ,j ,k ,ps) - v5*q(is)
              rhodg(nis,i1,j ,k ,ps) = rhodg(nis,i1,j ,k ,ps) - v6*q(is)
              rhodg(nis,i1,j1,k ,ps) = rhodg(nis,i1,j1,k ,ps) - v7*q(is)
              rhodg(nis,i ,j1,k ,ps) = rhodg(nis,i ,j1,k ,ps) - v8*q(is)
              nphgram(pbuf(m)%nid+1,is,ps) = nphgram(pbuf(m)%nid+1,is,ps) - 1
              pbuf(m)%nid = -1
              pbuf(m)%preside = 0
            end if
          end do
        end if
!
      end do ISL
!============================== species loop end


!-------------------- masking of charge
!      if(iimsk.eq.1) call fsmask(5)


  return
  end subroutine



!
  subroutine add_boundary_charge(rhodat,ps,fid,cid,ncomp)
!
!   ____________________________________________________________
!
!                       S U B R O U T I N E
!              A D D _ B O U N D A R Y _ C H A R G E
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .      this subroutine gives a charge distribution to      .
!   .      grid points from particle locations by the          .
!   .      "area-sharing" or "area-weighting" scheme.          .
!   ............................................................
!
!
!-------------------- parameter and common block
  use oh_type
  use paramt
  use allcom
  implicit none
!
  integer(kind=4),intent(in) :: ps,fid,cid,ncomp
  real(kind=8),intent(inout) :: &
 &  rhodat(1:,fsizes(1,1,fid):,fsizes(1,2,fid):,fsizes(1,3,fid):)
  integer(kind=4) :: i,j,k
  integer(kind=4) :: xu,yu,zu
  integer(kind=4) :: sl,su, dl,du, nl,nu


!-------------------- 
      xu = sdoms(2,1,sdid(ps)+1) - sdoms(1,1,sdid(ps)+1)
      yu = sdoms(2,2,sdid(ps)+1) - sdoms(1,2,sdid(ps)+1)
      zu = sdoms(2,3,sdid(ps)+1) - sdoms(1,3,sdid(ps)+1)


!-------------------- 
      sl = ctypes(2,2,1,cid)	!=-1
      su = ctypes(2,1,1,cid)	!=+1
!
      nl = ctypes(3,2,1,cid)	!= 1
      nu = ctypes(3,1,1,cid)	!= 1
!
      dl = sl + nl	!= 0
      du = su - nu	!= 0


!-------------------- 
      if((bounds(1,3,sdid(ps)+1).eq.1).or. &
     &   (bounds(1,3,sdid(ps)+1).eq.2.and.nfbnd(3).eq.0.and.nz.eq.1)) then
        do k=0,nl-1		!(i.e., do k=0,0)
        do j=0,yu+(su+nu-sl)-1	!(i.e., do j=0,yu+2)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+2)
          rhodat(1:ncomp,sl+i,sl+j,dl+k) = rhodat(1:ncomp,sl+i,sl+j,dl+k) &
         &                               + rhodat(1:ncomp,sl+i,sl+j,sl+k)
        end do
        end do
        end do
      end if

!-------------------- 
      if((bounds(2,3,sdid(ps)+1).eq.1).or. &
     &   (bounds(2,3,sdid(ps)+1).eq.2.and.nfbnd(3).eq.0.and.nz.eq.1)) then
        do k=0,nu-1		!(i.e., do k=0,0)
        do j=0,yu+(su+nu-sl)-1	!(i.e., do j=0,yu+2)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+2)
          rhodat(1:ncomp,sl+i,sl+j,zu+du+k) = rhodat(1:ncomp,sl+i,sl+j,zu+du+k) &
         &                                  + rhodat(1:ncomp,sl+i,sl+j,zu+su+k)
        end do
        end do
        end do
      end if

!-------------------- 
      if((bounds(1,2,sdid(ps)+1).eq.1).or. &
     &   (bounds(1,2,sdid(ps)+1).eq.2.and.nfbnd(2).eq.0.and.ny.eq.1)) then
        do k=0,zu+(du+nu-dl)-1	!(i.e., do k=0,zu)
        do j=0,nl-1		!(i.e., do j=0,0)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+2)
          rhodat(1:ncomp,sl+i,dl+j,dl+k) = rhodat(1:ncomp,sl+i,dl+j,dl+k) &
         &                               + rhodat(1:ncomp,sl+i,sl+j,dl+k)
        end do
        end do
        end do
      end if

!-------------------- 
      if((bounds(2,2,sdid(ps)+1).eq.1).or. &
     &   (bounds(2,2,sdid(ps)+1).eq.2.and.nfbnd(2).eq.0.and.ny.eq.1)) then
        do k=0,zu+(du+nu-dl)-1	!(i.e., do k=0,zu)
        do j=0,nu-1		!(i.e., do j=0,0)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+2)
          rhodat(1:ncomp,sl+i,yu+du+j,dl+k) = rhodat(1:ncomp,sl+i,yu+du+j,dl+k) &
         &                                  + rhodat(1:ncomp,sl+i,yu+su+j,dl+k)
        end do
        end do
        end do
      end if

!-------------------- 
      if((bounds(1,1,sdid(ps)+1).eq.1).or. &
     &   (bounds(1,1,sdid(ps)+1).eq.2.and.nfbnd(1).eq.0.and.nx.eq.1)) then
        do k=0,zu+(du+nu-dl)-1	!(i.e., do k=0,zu)
        do j=0,yu+(du+nu-dl)-1	!(i.e., do j=0,yu)
        do i=0,nl-1		!(i.e., do i=0,0)
          rhodat(1:ncomp,dl+i,dl+j,dl+k) = rhodat(1:ncomp,dl+i,dl+j,dl+k) &
         &                               + rhodat(1:ncomp,sl+i,dl+j,dl+k)
        end do
        end do
        end do
      end if

!-------------------- 
      if((bounds(2,1,sdid(ps)+1).eq.1).or. &
     &   (bounds(2,1,sdid(ps)+1).eq.2.and.nfbnd(1).eq.0.and.nx.eq.1)) then
        do k=0,zu+(du+nu-dl)-1	!(i.e., do k=0,zu)
        do j=0,yu+(du+nu-dl)-1	!(i.e., do j=0,yu)
        do i=0,nu-1		!(i.e., do i=0,0)
          rhodat(1:ncomp,xu+du+i,dl+j,dl+k) = rhodat(1:ncomp,xu+du+i,dl+j,dl+k) &
       &                                    + rhodat(1:ncomp,xu+su+i,dl+j,dl+k)
        end do
        end do
        end do
      end if


  return
  end subroutine add_boundary_charge



!
  subroutine exchange_lchg_borders(rhodat,ps,fid,cid,ncomp)
!
!   ____________________________________________________________
!
!                       S U B R O U T I N E
!            E X C H A N G E _ L C H G _ B O R D E R S
!   ____________________________________________________________
!
!   ............................................................
!   ............................................................
!
!
!-------------------- parameter and common block
  use oh_type
  use paramt
  use allcom
  implicit none
!
  integer(kind=4),intent(in) :: ps,fid,cid,ncomp
  real(kind=8),intent(inout) :: &
 &  rhodat(1:,fsizes(1,1,fid):,fsizes(1,2,fid):,fsizes(1,3,fid):)
  integer(kind=4) :: i,j,k
  integer(kind=4) :: xu,yu,zu
  integer(kind=4) :: sl,su, dl,du, nl,nu


!-------------------- 
      xu = sdoms(2,1,sdid(ps)+1) - sdoms(1,1,sdid(ps)+1)
      yu = sdoms(2,2,sdid(ps)+1) - sdoms(1,2,sdid(ps)+1)
      zu = sdoms(2,3,sdid(ps)+1) - sdoms(1,3,sdid(ps)+1)


!-------------------- 
      sl = ctypes(2,2,1,cid)	!=-1
      su = ctypes(2,1,1,cid)	!=+1
!
      nl = ctypes(3,2,1,cid)	!= 1
      nu = ctypes(3,1,1,cid)	!= 1
!
      dl = sl + nl	!= 0
      du = su - nu	!= 0


!-------------------- 
      if(bounds(1,3,sdid(ps)+1).eq.2.and.nfbnd(3).eq.0.and.nz.eq.1) then
        do k=0,nl-1		!(i.e., do k=0,0)
        do j=0,yu+(su+nu-sl)-1	!(i.e., do j=0,yu+2)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+2)
          rhodat(1:ncomp,sl+i,sl+j,zu+su+k) = rhodat(1:ncomp,sl+i,sl+j,dl+k)
        end do
        end do
        end do
      end if

!-------------------- 
      if(bounds(2,3,sdid(ps)+1).eq.2.and.nfbnd(3).eq.0.and.nz.eq.1) then
        do k=0,nu-1		!(i.e., do k=0,0)
        do j=0,yu+(su+nu-sl)-1	!(i.e., do j=0,yu+2)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+2)
          rhodat(1:ncomp,sl+i,sl+j,sl+k) = rhodat(1:ncomp,sl+i,sl+j,zu+du+k)
        end do
        end do
        end do
      end if

!-------------------- 
      if(bounds(1,2,sdid(ps)+1).eq.2.and.nfbnd(2).eq.0.and.ny.eq.1) then
        do k=0,zu+(su+nu-sl)-1	!(i.e., do k=0,zu+2)
        do j=0,nl-1		!(i.e., do j=0,0)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+2)
          rhodat(1:ncomp,sl+i,yu+su+j,sl+k) = rhodat(1:ncomp,sl+i,dl+j,sl+k)
        end do
        end do
        end do
      end if

!-------------------- 
      if(bounds(2,2,sdid(ps)+1).eq.2.and.nfbnd(2).eq.0.and.ny.eq.1) then
        do k=0,zu+(su+nu-sl)-1	!(i.e., do k=0,zu+2)
        do j=0,nu-1		!(i.e., do j=0,0)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+2)
          rhodat(1:ncomp,sl+i,sl+j,sl+k) = rhodat(1:ncomp,sl+i,yu+du+j,sl+k)
        end do
        end do
        end do
      end if

!-------------------- 
      if(bounds(1,1,sdid(ps)+1).eq.2.and.nfbnd(1).eq.0.and.nx.eq.1) then
        do k=0,zu+(su+nu-sl)-1	!(i.e., do k=0,zu+2)
        do j=0,yu+(su+nu-sl)-1	!(i.e., do j=0,yu+2)
        do i=0,nl-1		!(i.e., do i=0,0)
          rhodat(1:ncomp,xu+su+i,sl+j,sl+k) = rhodat(1:ncomp,dl+i,sl+j,sl+k)
        end do
        end do
        end do
      end if

!-------------------- 
      if(bounds(2,1,sdid(ps)+1).eq.2.and.nfbnd(1).eq.0.and.nx.eq.1) then
        do k=0,zu+(su+nu-sl)-1	!(i.e., do k=0,zu+2)
        do j=0,yu+(su+nu-sl)-1	!(i.e., do j=0,yu+2)
        do i=0,nu-1		!(i.e., do i=0,0)
          rhodat(1:ncomp,sl+i,sl+j,sl+k) = rhodat(1:ncomp,xu+du+i,sl+j,sl+k)
        end do
        end do
        end do
      end if


  return
  end subroutine exchange_lchg_borders



!
  subroutine add_background_charge(ps)
!
!   ____________________________________________________________
!
!                       S U B R O U T I N E
!            A D D _ B A C K G R O U N D _ C H A R G E
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .      this subroutine gives a charge distribution to      .
!   .      grid points from particle locations by the          .
!   .      "area-sharing" or "area-weighting" scheme.          .
!   ............................................................
!
!
!-------------------- parameter and common block
  use oh_type
  use paramt
  use allcom
  implicit none
!
  integer(kind=4) :: i,j,k
  integer(kind=4) :: xl,xu, yl,yu, zl,zu
  integer(kind=4) :: ps


!-------------------- 
      xl = sdoms(1,1,sdid(ps)+1); xu = sdoms(2,1,sdid(ps)+1)
      yl = sdoms(1,2,sdid(ps)+1); yu = sdoms(2,2,sdid(ps)+1)
      zl = sdoms(1,3,sdid(ps)+1); zu = sdoms(2,3,sdid(ps)+1)


!-------------------- 
      if(nflag_testp.ne.1) then
        do k=0,zu-zl
        do j=0,yu-yl
        do i=0,xu-xl
          rho(1,i,j,k,ps) = rho(1,i,j,k,ps) + rhobk(1,i,j,k,3)
        end do
        end do
        end do
      else
        do k=0,zu-zl
        do j=0,yu-yl
        do i=0,xu-xl
          rho(1,i,j,k,ps) = rhobk(1,i,j,k,3)
        end do
        end do
        end do
      end if


  return
  end subroutine add_background_charge
