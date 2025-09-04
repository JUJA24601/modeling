#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine velcty(ps)
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   V E L C T Y
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .      this subroutine updates praticle velocity.          .
!   .      buneman-boris method is used.                       .
!   ............................................................
!
!-------------------- parameter and common blocks
  use oh_type
  use paramt
  use allcom
  implicit none
!
  integer(kind=8) :: m, ns,ne, nee
  integer(kind=4) :: i,j,k, i1,j1,k1, is
  integer(kind=4) :: ps
  real(kind=8) :: xl,yl,zl
  real(kind=8) :: xlocal,ylocal,zlocal
  real(kind=8) :: x1,y1,z1, z2, xy1,xz1,yz1, xz2,yz2
  real(kind=8) :: v1,v2,v3,v4,v5,v6,v7,v8
  real(kind=8) :: eex,eey,eez, bbx,bby,bbz, boris, vxt,vyt,vzt


!-------------------- 
      xl = sdoms(1,1,sdid(ps)+1)
      yl = sdoms(1,2,sdid(ps)+1)
      zl = sdoms(1,3,sdid(ps)+1)


!============================== species loop
      nee = pbase(ps)
 SPL: do is=1,nspec
        ns = nee + 1
        ne = nee + totalp(is,ps)
        nee = ne
        do m=ns,ne
          if(pbuf(m)%nid.eq.-1) cycle
!
          xlocal = pbuf(m)%x - xl
          ylocal = pbuf(m)%y - yl
          zlocal = pbuf(m)%z - zl
!
          x1 = dint(xlocal)
          y1 = dint(ylocal)
          z1 = dint(zlocal)
!
          i = x1
          j = y1
          k = z1
!
          i1 = i + 1
          j1 = j + 1
          k1 = k + 1
!
          x1 = xlocal - x1
          y1 = ylocal - y1
          z1 = zlocal - z1
!         ----------- temporary value
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
!         ----------- e-field interpolation
          eex = peb(EX,i ,j ,k1,ps)*v1 + peb(EX,i1,j ,k1,ps)*v2 &
         &    + peb(EX,i1,j1,k1,ps)*v3 + peb(EX,i ,j1,k1,ps)*v4 &
         &    + peb(EX,i ,j ,k ,ps)*v5 + peb(EX,i1,j ,k ,ps)*v6 &
         &    + peb(EX,i1,j1,k ,ps)*v7 + peb(EX,i ,j1,k ,ps)*v8
          eey = peb(EY,i ,j ,k1,ps)*v1 + peb(EY,i1,j ,k1,ps)*v2 &
         &    + peb(EY,i1,j1,k1,ps)*v3 + peb(EY,i ,j1,k1,ps)*v4 &
         &    + peb(EY,i ,j ,k ,ps)*v5 + peb(EY,i1,j ,k ,ps)*v6 &
         &    + peb(EY,i1,j1,k ,ps)*v7 + peb(EY,i ,j1,k ,ps)*v8
          eez = peb(EZ,i ,j ,k1,ps)*v1 + peb(EZ,i1,j ,k1,ps)*v2 &
         &    + peb(EZ,i1,j1,k1,ps)*v3 + peb(EZ,i ,j1,k1,ps)*v4 &
         &    + peb(EZ,i ,j ,k ,ps)*v5 + peb(EZ,i1,j ,k ,ps)*v6 &
         &    + peb(EZ,i1,j1,k ,ps)*v7 + peb(EZ,i ,j1,k ,ps)*v8
          bbx = peb(BX,i ,j ,k1,ps)*v1 + peb(BX,i1,j ,k1,ps)*v2 &
         &    + peb(BX,i1,j1,k1,ps)*v3 + peb(BX,i ,j1,k1,ps)*v4 &
         &    + peb(BX,i ,j ,k ,ps)*v5 + peb(BX,i1,j ,k ,ps)*v6 &
         &    + peb(BX,i1,j1,k ,ps)*v7 + peb(BX,i ,j1,k ,ps)*v8
          bby = peb(BY,i ,j ,k1,ps)*v1 + peb(BY,i1,j ,k1,ps)*v2 &
         &    + peb(BY,i1,j1,k1,ps)*v3 + peb(BY,i ,j1,k1,ps)*v4 &
         &    + peb(BY,i ,j ,k ,ps)*v5 + peb(BY,i1,j ,k ,ps)*v6 &
         &    + peb(BY,i1,j1,k ,ps)*v7 + peb(BY,i ,j1,k ,ps)*v8
          bbz = peb(BZ,i ,j ,k1,ps)*v1 + peb(BZ,i1,j ,k1,ps)*v2 &
         &    + peb(BZ,i1,j1,k1,ps)*v3 + peb(BZ,i ,j1,k1,ps)*v4 &
         &    + peb(BZ,i ,j ,k ,ps)*v5 + peb(BZ,i1,j ,k ,ps)*v6 &
         &    + peb(BZ,i1,j1,k ,ps)*v7 + peb(BZ,i ,j1,k ,ps)*v8
!         ----------- charge-to-mass ratio is taken into consideration
          eex = eex*qm(is)
          eey = eey*qm(is)
          eez = eez*qm(is)
          bbx = bbx*qm(is)
          bby = bby*qm(is)
          bbz = bbz*qm(is)
!         ----------- Buneman-Boris method
          boris=2.0d0/(1.0d0+bbx*bbx+bby*bby+bbz*bbz)
!
          pbuf(m)%vx = pbuf(m)%vx + eex
          pbuf(m)%vy = pbuf(m)%vy + eey
          pbuf(m)%vz = pbuf(m)%vz + eez
!
          vxt = pbuf(m)%vx + pbuf(m)%vy*bbz - pbuf(m)%vz*bby
          vyt = pbuf(m)%vy + pbuf(m)%vz*bbx - pbuf(m)%vx*bbz
          vzt = pbuf(m)%vz + pbuf(m)%vx*bby - pbuf(m)%vy*bbx
!
          pbuf(m)%vx = pbuf(m)%vx + boris*(vyt*bbz - vzt*bby)
          pbuf(m)%vy = pbuf(m)%vy + boris*(vzt*bbx - vxt*bbz)
          pbuf(m)%vz = pbuf(m)%vz + boris*(vxt*bby - vyt*bbx)
!
          pbuf(m)%vx = pbuf(m)%vx + eex
          pbuf(m)%vy = pbuf(m)%vy + eey
          pbuf(m)%vz = pbuf(m)%vz + eez
        end do
!
      end do SPL
!============================== species loop


  return
  end subroutine
