#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine medium
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   M E D I U M
!   ____________________________________________________________
!
!-------------------- parameter and common blocks
  use oh_type
  use paramt
  use allcom
  implicit none
!
  interface rfold
    integer(kind=4) function rfold_i(r,nr,dir)
      integer(kind=4),intent(in) :: r, nr
      integer(kind=4),intent(in) :: dir
    end function rfold_i
!
    real(kind=4) function rfold_r(r,nr,dir)
      real(kind=4),intent(in) :: r, nr
      integer(kind=4),intent(in) :: dir
    end function rfold_r
!
    real(kind=8) function rfold_d(r,nr,dir)
      real(kind=8),intent(in) :: r, nr
      integer(kind=4),intent(in) :: dir
    end function rfold_d
  end interface rfold
!
  type pmls
    integer(kind=4) :: nxl,nxu, nyl,nyu, nzl,nzu
  end type pmls
!
  integer(kind=4) :: i,j,k
  integer(kind=4) :: xl,xu, yl,yu, zl,zu
  integer(kind=4) :: ipc
  real(kind=8) :: x,y,z, x1,y1,z1
  real(kind=8) :: xlbd,xubd,ylbd,yubd,zlbd,zubd
  real(kind=8) :: cylaln, cylrad, cylaxs1,cylaxs2, cyledg1,cyledg2
  real(kind=8) :: sphrad, sphcnt1,sphcnt2,sphcnt3
  real(kind=8) :: radsq,radsqin
  real(kind=8) :: smax0(2,3)
  real(kind=8) :: sigmaxe,sigmaxm, sigmaye,sigmaym, sigmaze,sigmazm
  real(kind=8) :: mediatr, z0
  type(pmls) :: pml


!--------------------
      xl = sdoms(1,1,sdid(1)+1); xu = sdoms(2,1,sdid(1)+1)
      yl = sdoms(1,2,sdid(1)+1); yu = sdoms(2,2,sdid(1)+1)
      zl = sdoms(1,3,sdid(1)+1); zu = sdoms(2,3,sdid(1)+1)


!-------------------- 
      do ipc=npc,1,-1
        xlbd = xlpc(ipc); xubd = xupc(ipc)
        ylbd = ylpc(ipc); yubd = yupc(ipc)
        zlbd = zlpc(ipc); zubd = zupc(ipc)
!
        cylaln = cylinder(ipc)%align
        cylrad = cylinder(ipc)%radius
        cylaxs1 = cylinder(ipc)%axis(1)
        cylaxs2 = cylinder(ipc)%axis(2)
        cyledg1 = cylinder(ipc)%edge(1)
        cyledg2 = cylinder(ipc)%edge(2)
!
        sphrad = sphere(ipc)%radius
        sphcnt1 = sphere(ipc)%center(1)
        sphcnt2 = sphere(ipc)%center(2)
        sphcnt3 = sphere(ipc)%center(3)
!
!GEOTYPE0 or 1
        if(geotype(ipc).eq.0.or.geotype(ipc).eq.1) then
          do k=-1,zu-zl+1
          do j=-1,yu-yl+1
          do i=-1,xu-xl+1
            x  = (i + xl)*dr
            x1 = (i + 1 + xl)*dr
            y  = (j + yl)*dr
            y1 = (j + 1 + yl)*dr
            z  = (k + zl)*dr
            z1 = (k + 1 + zl)*dr
            if(x .ge.xlbd.and.x .le.xubd.and. &
           &   x1.ge.xlbd.and.x1.le.xubd.and. &
           &   y .ge.ylbd.and.y .le.yubd.and. &
           &   z .ge.zlbd.and.z .le.zubd) then
              mp(EX,i,j,k,1) = 1.0d0
            end if
            if(x .ge.xlbd.and.x .le.xubd.and. &
           &   y .ge.ylbd.and.y .le.yubd.and. &
           &   y1.ge.ylbd.and.y1.le.yubd.and. &
           &   z .ge.zlbd.and.z .le.zubd) then
              mp(EY,i,j,k,1) = 1.0d0
            end if
            if(x .ge.xlbd.and.x .le.xubd.and. &
           &   y .ge.ylbd.and.y .le.yubd.and. &
           &   z .ge.zlbd.and.z .le.zubd.and. &
           &   z1.ge.zlbd.and.z1.le.zubd) then
              mp(EZ,i,j,k,1) = 1.0d0
            end if
            if(x .ge.xlbd.and.x .le.xubd.and. &
           &   y .ge.ylbd.and.y .le.yubd.and. &
           &   y1.ge.ylbd.and.y1.le.yubd.and. &
           &   z .ge.zlbd.and.z .le.zubd.and. &
           &   z1.ge.zlbd.and.z1.le.zubd) then
              mp(BX,i,j,k,1) = 1.0d0
            end if
            if(x .ge.xlbd.and.x .le.xubd.and. &
           &   x1.ge.xlbd.and.x1.le.xubd.and. &
           &   y .ge.ylbd.and.y .le.yubd.and. &
           &   z .ge.zlbd.and.z .le.zubd.and. &
           &   z1.ge.zlbd.and.z1.le.zubd) then
              mp(BY,i,j,k,1) = 1.0d0
            end if
            if(x .ge.xlbd.and.x .le.xubd.and. &
           &   x1.ge.xlbd.and.x1.le.xubd.and. &
           &   y .ge.ylbd.and.y .le.yubd.and. &
           &   y1.ge.ylbd.and.y1.le.yubd.and. &
           &   z .ge.zlbd.and.z .le.zubd) then
              mp(BZ,i,j,k,1) = 1.0d0
            end if
!
            x  = rfold(i+xl,nx,1)*dr
            x1 = rfold(i+1+xl,nx,1)*dr
            y  = rfold(j+yl,ny,2)*dr
            y1 = rfold(j+1+yl,ny,2)*dr
            z  = rfold(k+zl,nz,3)*dr
            z1 = rfold(k+1+zl,nz,3)*dr
            if(x .ge.xlbd.and.x .le.xubd.and. &
           &   x1.ge.xlbd.and.x1.le.xubd.and. &
           &   y .ge.ylbd.and.y .le.yubd.and. &
           &   z .ge.zlbd.and.z .le.zubd) then
              mp(EX,i,j,k,1) = 1.0d0
            end if
            if(x .ge.xlbd.and.x .le.xubd.and. &
           &   y .ge.ylbd.and.y .le.yubd.and. &
           &   y1.ge.ylbd.and.y1.le.yubd.and. &
           &   z .ge.zlbd.and.z .le.zubd) then
              mp(EY,i,j,k,1) = 1.0d0
            end if
            if(x .ge.xlbd.and.x .le.xubd.and. &
           &   y .ge.ylbd.and.y .le.yubd.and. &
           &   z .ge.zlbd.and.z .le.zubd.and. &
           &   z1.ge.zlbd.and.z1.le.zubd) then
              mp(EZ,i,j,k,1) = 1.0d0
            end if
            if(x .ge.xlbd.and.x .le.xubd.and. &
           &   y .ge.ylbd.and.y .le.yubd.and. &
           &   y1.ge.ylbd.and.y1.le.yubd.and. &
           &   z .ge.zlbd.and.z .le.zubd.and. &
           &   z1.ge.zlbd.and.z1.le.zubd) then
              mp(BX,i,j,k,1) = 1.0d0
            end if
            if(x .ge.xlbd.and.x .le.xubd.and. &
           &   x1.ge.xlbd.and.x1.le.xubd.and. &
           &   y .ge.ylbd.and.y .le.yubd.and. &
           &   z .ge.zlbd.and.z .le.zubd.and. &
           &   z1.ge.zlbd.and.z1.le.zubd) then
              mp(BY,i,j,k,1) = 1.0d0
            end if
            if(x .ge.xlbd.and.x .le.xubd.and. &
           &   x1.ge.xlbd.and.x1.le.xubd.and. &
           &   y .ge.ylbd.and.y .le.yubd.and. &
           &   y1.ge.ylbd.and.y1.le.yubd.and. &
           &   z .ge.zlbd.and.z .le.zubd) then
              mp(BZ,i,j,k,1) = 1.0d0
            end if
          end do
          end do
          end do
!
!GEOTYPE2
        else if(geotype(ipc).eq.2) then
          radsq = cylrad**2
          radsqin = (cylrad - 1.0d0)**2
!GEOTYPE2-1
          if(cylaln.eq.1) then
            do k=-1,zu-zl+1
            do j=-1,yu-yl+1
            do i=-1,xu-xl+1
              x  = (i + xl)*dr
              x1 = (i + 1 + xl)*dr
              y  = (j + yl)*dr
              y1 = (j + 1 + yl)*dr
              z  = (k + zl)*dr
              z1 = (k + 1 + zl)*dr
              if(x .ge.cyledg1.and.x .le.cyledg2.and. &
             &   x1.ge.cyledg1.and.x1.le.cyledg2.and. &
             &   (y -cylaxs1)**2+(z -cylaxs2)**2.le.radsq) then
                mp(EX,i,j,k,1) = 1.0d0
              end if
              if(x .ge.cyledg1.and.x .le.cyledg2.and. &
             &   (y -cylaxs1)**2+(z -cylaxs2)**2.le.radsq.and. &
             &   (y1-cylaxs1)**2+(z -cylaxs2)**2.le.radsq) then
                mp(EY,i,j,k,1) = 1.0d0
              end if
              if(x .ge.cyledg1.and.x .le.cyledg2.and. &
             &   (y -cylaxs1)**2+(z -cylaxs2)**2.le.radsq.and. &
             &   (y -cylaxs1)**2+(z1-cylaxs2)**2.le.radsq) then
                mp(EZ,i,j,k,1) = 1.0d0
              end if
              if(x .ge.cyledg1.and.x .le.cyledg2.and. &
             &   (y -cylaxs1)**2+(z -cylaxs2)**2.le.radsq.and. &
             &   (y1-cylaxs1)**2+(z -cylaxs2)**2.le.radsq.and. &
             &   (y -cylaxs1)**2+(z1-cylaxs2)**2.le.radsq.and. &
             &   (y1-cylaxs1)**2+(z1-cylaxs2)**2.le.radsq) then
                mp(BX,i,j,k,1) = 1.0d0
              end if
              if(x .ge.cyledg1.and.x .le.cyledg2.and. &
             &   x1.ge.cyledg1.and.x1.le.cyledg2.and. &
             &   (y -cylaxs1)**2+(z -cylaxs2)**2.le.radsq.and. &
             &   (y -cylaxs1)**2+(z1-cylaxs2)**2.le.radsq) then
                mp(BY,i,j,k,1) = 1.0d0
              end if
              if(x .ge.cyledg1.and.x .le.cyledg2.and. &
             &   x1.ge.cyledg1.and.x1.le.cyledg2.and. &
             &   (y -cylaxs1)**2+(z -cylaxs2)**2.le.radsq.and. &
             &   (y1-cylaxs1)**2+(z -cylaxs2)**2.le.radsq) then
                mp(BZ,i,j,k,1) = 1.0d0
              end if
!
              x  = rfold(i+xl,nx,1)*dr
              x1 = rfold(i+1+xl,nx,1)*dr
              y  = rfold(j+yl,ny,2)*dr
              y1 = rfold(j+1+yl,ny,2)*dr
              z  = rfold(k+zl,nz,3)*dr
              z1 = rfold(k+1+zl,nz,3)*dr
              if(x .ge.cyledg1.and.x .le.cyledg2.and. &
             &   x1.ge.cyledg1.and.x1.le.cyledg2.and. &
             &   (y -cylaxs1)**2+(z -cylaxs2)**2.le.radsq) then
                mp(EX,i,j,k,1) = 1.0d0
              end if
              if(x .ge.cyledg1.and.x .le.cyledg2.and. &
             &   (y -cylaxs1)**2+(z -cylaxs2)**2.le.radsq.and. &
             &   (y1-cylaxs1)**2+(z -cylaxs2)**2.le.radsq) then
                mp(EY,i,j,k,1) = 1.0d0
              end if
              if(x .ge.cyledg1.and.x .le.cyledg2.and. &
             &   (y -cylaxs1)**2+(z -cylaxs2)**2.le.radsq.and. &
             &   (y -cylaxs1)**2+(z1-cylaxs2)**2.le.radsq) then
                mp(EZ,i,j,k,1) = 1.0d0
              end if
              if(x .ge.cyledg1.and.x .le.cyledg2.and. &
             &   (y -cylaxs1)**2+(z -cylaxs2)**2.le.radsq.and. &
             &   (y1-cylaxs1)**2+(z -cylaxs2)**2.le.radsq.and. &
             &   (y -cylaxs1)**2+(z1-cylaxs2)**2.le.radsq.and. &
             &   (y1-cylaxs1)**2+(z1-cylaxs2)**2.le.radsq) then
                mp(BX,i,j,k,1) = 1.0d0
              end if
              if(x .ge.cyledg1.and.x .le.cyledg2.and. &
             &   x1.ge.cyledg1.and.x1.le.cyledg2.and. &
             &   (y -cylaxs1)**2+(z -cylaxs2)**2.le.radsq.and. &
             &   (y -cylaxs1)**2+(z1-cylaxs2)**2.le.radsq) then
                mp(BY,i,j,k,1) = 1.0d0
              end if
              if(x .ge.cyledg1.and.x .le.cyledg2.and. &
             &   x1.ge.cyledg1.and.x1.le.cyledg2.and. &
             &   (y -cylaxs1)**2+(z -cylaxs2)**2.le.radsq.and. &
             &   (y1-cylaxs1)**2+(z -cylaxs2)**2.le.radsq) then
                mp(BZ,i,j,k,1) = 1.0d0
              end if
            end do
            end do
            end do
!GEOTYPE2-2
          else if(cylaln.eq.2) then
            do k=-1,zu-zl+1
            do j=-1,yu-yl+1
            do i=-1,xu-xl+1
              x  = (i + xl)*dr
              x1 = (i + 1 + xl)*dr
              y  = (j + yl)*dr
              y1 = (j + 1 + yl)*dr
              z  = (k + zl)*dr
              z1 = (k + 1 + zl)*dr
              if(y .ge.cyledg1.and.y .le.cyledg2.and. &
             &   (z -cylaxs1)**2+(x -cylaxs2)**2.le.radsq.and. &
             &   (z -cylaxs1)**2+(x1-cylaxs2)**2.le.radsq) then
                mp(EX,i,j,k,1) = 1.0d0
              end if
              if(y .ge.cyledg1.and.y .le.cyledg2.and. &
             &   y1.ge.cyledg1.and.y1.le.cyledg2.and. &
             &   (z -cylaxs1)**2+(x -cylaxs2)**2.le.radsq) then
                mp(EY,i,j,k,1) = 1.0d0
              end if
              if(y .ge.cyledg1.and.y .le.cyledg2.and. &
             &   (z -cylaxs1)**2+(x -cylaxs2)**2.le.radsq.and. &
             &   (z1-cylaxs1)**2+(x -cylaxs2)**2.le.radsq) then
                mp(EZ,i,j,k,1) = 1.0d0
              end if
              if(y .ge.cyledg1.and.y .le.cyledg2.and. &
             &   y1.ge.cyledg1.and.y1.le.cyledg2.and. &
             &   (z -cylaxs1)**2+(x -cylaxs2)**2.le.radsq.and. &
             &   (z1-cylaxs1)**2+(x -cylaxs2)**2.le.radsq) then
                mp(BX,i,j,k,1) = 1.0d0
              end if
              if(y .ge.cyledg1.and.y .le.cyledg2.and. &
             &   (z -cylaxs1)**2+(x -cylaxs2)**2.le.radsq.and. &
             &   (z -cylaxs1)**2+(x1-cylaxs2)**2.le.radsq.and. &
             &   (z1-cylaxs1)**2+(x -cylaxs2)**2.le.radsq.and. &
             &   (z1-cylaxs1)**2+(x1-cylaxs2)**2.le.radsq) then
                mp(BY,i,j,k,1) = 1.0d0
              end if
              if(y .ge.cyledg1.and.y .le.cyledg2.and. &
             &   y1.ge.cyledg1.and.y1.le.cyledg2.and. &
             &   (z -cylaxs1)**2+(x -cylaxs2)**2.le.radsq.and. &
             &   (z -cylaxs1)**2+(x1-cylaxs2)**2.le.radsq) then
                mp(BZ,i,j,k,1) = 1.0d0
              end if
!
              x  = rfold(i+xl,nx,1)*dr
              x1 = rfold(i+1+xl,nx,1)*dr
              y  = rfold(j+yl,ny,2)*dr
              y1 = rfold(j+1+yl,ny,2)*dr
              z  = rfold(k+zl,nz,3)*dr
              z1 = rfold(k+1+zl,nz,3)*dr
              if(y .ge.cyledg1.and.y .le.cyledg2.and. &
             &   (z -cylaxs1)**2+(x -cylaxs2)**2.le.radsq.and. &
             &   (z -cylaxs1)**2+(x1-cylaxs2)**2.le.radsq) then
                mp(EX,i,j,k,1) = 1.0d0
              end if
              if(y .ge.cyledg1.and.y .le.cyledg2.and. &
             &   y1.ge.cyledg1.and.y1.le.cyledg2.and. &
             &   (z -cylaxs1)**2+(x -cylaxs2)**2.le.radsq) then
                mp(EY,i,j,k,1) = 1.0d0
              end if
              if(y .ge.cyledg1.and.y .le.cyledg2.and. &
             &   (z -cylaxs1)**2+(x -cylaxs2)**2.le.radsq.and. &
             &   (z1-cylaxs1)**2+(x -cylaxs2)**2.le.radsq) then
                mp(EZ,i,j,k,1) = 1.0d0
              end if
              if(y .ge.cyledg1.and.y .le.cyledg2.and. &
             &   y1.ge.cyledg1.and.y1.le.cyledg2.and. &
             &   (z -cylaxs1)**2+(x -cylaxs2)**2.le.radsq.and. &
             &   (z1-cylaxs1)**2+(x -cylaxs2)**2.le.radsq) then
                mp(BX,i,j,k,1) = 1.0d0
              end if
              if(y .ge.cyledg1.and.y .le.cyledg2.and. &
             &   (z -cylaxs1)**2+(x -cylaxs2)**2.le.radsq.and. &
             &   (z -cylaxs1)**2+(x1-cylaxs2)**2.le.radsq.and. &
             &   (z1-cylaxs1)**2+(x -cylaxs2)**2.le.radsq.and. &
             &   (z1-cylaxs1)**2+(x1-cylaxs2)**2.le.radsq) then
                mp(BY,i,j,k,1) = 1.0d0
              end if
              if(y .ge.cyledg1.and.y .le.cyledg2.and. &
             &   y1.ge.cyledg1.and.y1.le.cyledg2.and. &
             &   (z -cylaxs1)**2+(x -cylaxs2)**2.le.radsq.and. &
             &   (z -cylaxs1)**2+(x1-cylaxs2)**2.le.radsq) then
                mp(BZ,i,j,k,1) = 1.0d0
              end if
            end do
            end do
            end do
!GEOTYPE2-3
          else if(cylaln.eq.3) then
            do k=-1,zu-zl+1
            do j=-1,yu-yl+1
            do i=-1,xu-xl+1
              x  = (i + xl)*dr
              x1 = (i + 1 + xl)*dr
              y  = (j + yl)*dr
              y1 = (j + 1 + yl)*dr
              z  = (k + zl)*dr
              z1 = (k + 1 + zl)*dr
              if(z .ge.cyledg1.and.z .le.cyledg2.and. &
             &   (x -cylaxs1)**2+(y -cylaxs2)**2.le.radsq.and. &
             &   (x1-cylaxs1)**2+(y -cylaxs2)**2.le.radsq) then
                mp(EX,i,j,k,1) = 1.0d0
              end if
              if(z .ge.cyledg1.and.z .le.cyledg2.and. &
             &   (x -cylaxs1)**2+(y -cylaxs2)**2.le.radsq.and. &
             &   (x -cylaxs1)**2+(y1-cylaxs2)**2.le.radsq) then
                mp(EY,i,j,k,1) = 1.0d0
              end if
              if(z .ge.cyledg1.and.z .le.cyledg2.and. &
             &   z1.ge.cyledg1.and.z1.le.cyledg2.and. &
             &   (x -cylaxs1)**2+(y -cylaxs2)**2.le.radsq) then
                mp(EZ,i,j,k,1) = 1.0d0
              end if
              if(z .ge.cyledg1.and.z .le.cyledg2.and. &
             &   z1.ge.cyledg1.and.z1.le.cyledg2.and. &
             &   (x -cylaxs1)**2+(y -cylaxs2)**2.le.radsq.and. &
             &   (x -cylaxs1)**2+(y1-cylaxs2)**2.le.radsq) then
                mp(BX,i,j,k,1) = 1.0d0
              end if
              if(z .ge.cyledg1.and.z .le.cyledg2.and. &
             &   z1.ge.cyledg1.and.z1.le.cyledg2.and. &
             &   (x -cylaxs1)**2+(y -cylaxs2)**2.le.radsq.and. &
             &   (x1-cylaxs1)**2+(y -cylaxs2)**2.le.radsq) then
                mp(BY,i,j,k,1) = 1.0d0
              end if
              if(z .ge.cyledg1.and.z .le.cyledg2.and. &
             &   (x -cylaxs1)**2+(y -cylaxs2)**2.le.radsq.and. &
             &   (x1-cylaxs1)**2+(y -cylaxs2)**2.le.radsq.and. &
             &   (x -cylaxs1)**2+(y1-cylaxs2)**2.le.radsq.and. &
             &   (x1-cylaxs1)**2+(y1-cylaxs2)**2.le.radsq) then
                mp(BZ,i,j,k,1) = 1.0d0
              end if
!
              x  = rfold(i+xl,nx,1)*dr
              x1 = rfold(i+1+xl,nx,1)*dr
              y  = rfold(j+yl,ny,2)*dr
              y1 = rfold(j+1+yl,ny,2)*dr
              z  = rfold(k+zl,nz,3)*dr
              z1 = rfold(k+1+zl,nz,3)*dr
              if(z .ge.cyledg1.and.z .le.cyledg2.and. &
             &   (x -cylaxs1)**2+(y -cylaxs2)**2.le.radsq.and. &
             &   (x1-cylaxs1)**2+(y -cylaxs2)**2.le.radsq) then
                mp(EX,i,j,k,1) = 1.0d0
              end if
              if(z .ge.cyledg1.and.z .le.cyledg2.and. &
             &   (x -cylaxs1)**2+(y -cylaxs2)**2.le.radsq.and. &
             &   (x -cylaxs1)**2+(y1-cylaxs2)**2.le.radsq) then
                mp(EY,i,j,k,1) = 1.0d0
              end if
              if(z .ge.cyledg1.and.z .le.cyledg2.and. &
             &   z1.ge.cyledg1.and.z1.le.cyledg2.and. &
             &   (x -cylaxs1)**2+(y -cylaxs2)**2.le.radsq) then
                mp(EZ,i,j,k,1) = 1.0d0
              end if
              if(z .ge.cyledg1.and.z .le.cyledg2.and. &
             &   z1.ge.cyledg1.and.z1.le.cyledg2.and. &
             &   (x -cylaxs1)**2+(y -cylaxs2)**2.le.radsq.and. &
             &   (x -cylaxs1)**2+(y1-cylaxs2)**2.le.radsq) then
                mp(BX,i,j,k,1) = 1.0d0
              end if
              if(z .ge.cyledg1.and.z .le.cyledg2.and. &
             &   z1.ge.cyledg1.and.z1.le.cyledg2.and. &
             &   (x -cylaxs1)**2+(y -cylaxs2)**2.le.radsq.and. &
             &   (x1-cylaxs1)**2+(y -cylaxs2)**2.le.radsq) then
                mp(BY,i,j,k,1) = 1.0d0
              end if
              if(z .ge.cyledg1.and.z .le.cyledg2.and. &
             &   (x -cylaxs1)**2+(y -cylaxs2)**2.le.radsq.and. &
             &   (x1-cylaxs1)**2+(y -cylaxs2)**2.le.radsq.and. &
             &   (x -cylaxs1)**2+(y1-cylaxs2)**2.le.radsq.and. &
             &   (x1-cylaxs1)**2+(y1-cylaxs2)**2.le.radsq) then
                mp(BZ,i,j,k,1) = 1.0d0
              end if
            end do
            end do
            end do
          end if
!
!GEOTYPE3
        else if(geotype(ipc).eq.3) then
          radsq = sphrad**2
          radsqin = (sphrad - 1.0d0)**2
!
          do k=-1,zu-zl+1
          do j=-1,yu-yl+1
          do i=-1,xu-xl+1
            x  = (i + xl)*dr
            x1 = (i + 1 + xl)*dr
            y  = (j + yl)*dr
            y1 = (j + 1 + yl)*dr
            z  = (k + zl)*dr
            z1 = (k + 1 + zl)*dr
            if((x -sphcnt1)**2+(y -sphcnt2)**2+(z -sphcnt3)**2.le.radsq.and. &
           &   (x1-sphcnt1)**2+(y -sphcnt2)**2+(z -sphcnt3)**2.le.radsq) then
              mp(EX,i,j,k,1) = 1.0d0
            end if
            if((x -sphcnt1)**2+(y -sphcnt2)**2+(z -sphcnt3)**2.le.radsq.and. &
           &   (x -sphcnt1)**2+(y1-sphcnt2)**2+(z -sphcnt3)**2.le.radsq) then
              mp(EY,i,j,k,1) = 1.0d0
            end if
            if((x -sphcnt1)**2+(y -sphcnt2)**2+(z -sphcnt3)**2.le.radsq.and. &
           &   (x -sphcnt1)**2+(y -sphcnt2)**2+(z1-sphcnt3)**2.le.radsq) then
              mp(EZ,i,j,k,1) = 1.0d0
            end if
            if((x -sphcnt1)**2+(y -sphcnt2)**2+(z -sphcnt3)**2.le.radsq.and. &
           &   (x -sphcnt1)**2+(y1-sphcnt2)**2+(z -sphcnt3)**2.le.radsq.and. &
           &   (x -sphcnt1)**2+(y -sphcnt2)**2+(z1-sphcnt3)**2.le.radsq.and. &
           &   (x -sphcnt1)**2+(y1-sphcnt2)**2+(z1-sphcnt3)**2.le.radsq) then
              mp(BX,i,j,k,1) = 1.0d0
            end if
            if((x -sphcnt1)**2+(y -sphcnt2)**2+(z -sphcnt3)**2.le.radsq.and. &
           &   (x1-sphcnt1)**2+(y -sphcnt2)**2+(z -sphcnt3)**2.le.radsq.and. &
           &   (x -sphcnt1)**2+(y -sphcnt2)**2+(z1-sphcnt3)**2.le.radsq.and. &
           &   (x1-sphcnt1)**2+(y -sphcnt2)**2+(z1-sphcnt3)**2.le.radsq) then
              mp(BY,i,j,k,1) = 1.0d0
            end if
            if((x -sphcnt1)**2+(y -sphcnt2)**2+(z -sphcnt3)**2.le.radsq.and. &
           &   (x1-sphcnt1)**2+(y -sphcnt2)**2+(z -sphcnt3)**2.le.radsq.and. &
           &   (x -sphcnt1)**2+(y1-sphcnt2)**2+(z -sphcnt3)**2.le.radsq.and. &
           &   (x1-sphcnt1)**2+(y1-sphcnt2)**2+(z -sphcnt3)**2.le.radsq) then
              mp(BZ,i,j,k,1) = 1.0d0
            end if
!
            x  = rfold(i+xl,nx,1)*dr
            x1 = rfold(i+1+xl,nx,1)*dr
            y  = rfold(j+yl,ny,2)*dr
            y1 = rfold(j+1+yl,ny,2)*dr
            z  = rfold(k+zl,nz,3)*dr
            z1 = rfold(k+1+zl,nz,3)*dr
            if((x -sphcnt1)**2+(y -sphcnt2)**2+(z -sphcnt3)**2.le.radsq.and. &
           &   (x1-sphcnt1)**2+(y -sphcnt2)**2+(z -sphcnt3)**2.le.radsq) then
              mp(EX,i,j,k,1) = 1.0d0
            end if
            if((x -sphcnt1)**2+(y -sphcnt2)**2+(z -sphcnt3)**2.le.radsq.and. &
           &   (x -sphcnt1)**2+(y1-sphcnt2)**2+(z -sphcnt3)**2.le.radsq) then
              mp(EY,i,j,k,1) = 1.0d0
            end if
            if((x -sphcnt1)**2+(y -sphcnt2)**2+(z -sphcnt3)**2.le.radsq.and. &
           &   (x -sphcnt1)**2+(y -sphcnt2)**2+(z1-sphcnt3)**2.le.radsq) then
              mp(EZ,i,j,k,1) = 1.0d0
            end if
            if((x -sphcnt1)**2+(y -sphcnt2)**2+(z -sphcnt3)**2.le.radsq.and. &
           &   (x -sphcnt1)**2+(y1-sphcnt2)**2+(z -sphcnt3)**2.le.radsq.and. &
           &   (x -sphcnt1)**2+(y -sphcnt2)**2+(z1-sphcnt3)**2.le.radsq.and. &
           &   (x -sphcnt1)**2+(y1-sphcnt2)**2+(z1-sphcnt3)**2.le.radsq) then
              mp(BX,i,j,k,1) = 1.0d0
            end if
            if((x -sphcnt1)**2+(y -sphcnt2)**2+(z -sphcnt3)**2.le.radsq.and. &
           &   (x1-sphcnt1)**2+(y -sphcnt2)**2+(z -sphcnt3)**2.le.radsq.and. &
           &   (x -sphcnt1)**2+(y -sphcnt2)**2+(z1-sphcnt3)**2.le.radsq.and. &
           &   (x1-sphcnt1)**2+(y -sphcnt2)**2+(z1-sphcnt3)**2.le.radsq) then
              mp(BY,i,j,k,1) = 1.0d0
            end if
            if((x -sphcnt1)**2+(y -sphcnt2)**2+(z -sphcnt3)**2.le.radsq.and. &
           &   (x1-sphcnt1)**2+(y -sphcnt2)**2+(z -sphcnt3)**2.le.radsq.and. &
           &   (x -sphcnt1)**2+(y1-sphcnt2)**2+(z -sphcnt3)**2.le.radsq.and. &
           &   (x1-sphcnt1)**2+(y1-sphcnt2)**2+(z -sphcnt3)**2.le.radsq) then
              mp(BZ,i,j,k,1) = 1.0d0
            end if
          end do
          end do
          end do
!
        end if
      end do


!--------------------
!      do k=-1,zu-zl+1
!      do j=-1,yu-yl+1
!      do i=-1,xu-xl+1
      do k=0,zu-zl
      do j=0,yu-yl
      do i=0,xu-xl
        x  = (i + xl)*dr
        x1 = (i + 1 + xl)*dr
        y  = (j + yl)*dr
        y1 = (j + 1 + yl)*dr
        z  = (k + zl)*dr
        z1 = (k + 1 + zl)*dr
        if(mp(EX,i,j,k,1).eq.0.0d0) then
!          x = (i + xl + 0.5d0)*dr
!          y = (j + yl)*dr
!          z = (k + zl)*dr
!          print*, "  ",x,y,z,"mp(EX) = 0", myid
!          if(medi(1,i,j,k).ne.1.or.medi(1,i+1,j,k).ne.1) then
!            print*, "  ",x,y,z,"illegal EX0", myid,medi(1,i,j,k),medi(1,i+1,j,k)
!          end if

        end if
        if(mp(EY,i,j,k,1).eq.0.0d0) then
          x = (i + xl)*dr
          y = (j + yl + 0.5d0)*dr
          z = (k + zl)*dr
!          print*, "  ",x,y,z,"mp(EY) = 0", myid
!          if(medi(1,i,j,k).ne.1.or.medi(1,i,j+1,k).ne.1) then
!            print*, "  ",x,y,z,"illegal EY0", myid,medi(1,i,j,k),medi(1,i,j+1,k)
!          end if
        end if
        if(mp(EZ,i,j,k,1).eq.0.0d0) then
          x = (i + xl)*dr
          y = (j + yl)*dr
          z = (k + zl + 0.5d0)*dr
!          print*, "  ",x,y,z,"mp(EZ) = 0", myid
!          if(medi(1,i,j,k).ne.1.or.medi(1,i,j,k+1).ne.1) then
!            print*, "  ",x,y,z,"illegal EZ0", myid,medi(1,i,j,k),medi(1,i,j,k+1)
!          end if
        end if
        if(mp(BX,i,j,k,1).eq.0.0d0) then
          x = (i + xl)*dr
          y = (j + yl + 0.5d0)*dr
          z = (k + zl + 0.5d0)*dr
!          print*, "  ",x,y,z,"mp(BX) = 0", myid
        end if
        if(mp(BY,i,j,k,1).eq.0.0d0) then
          x = (i + xl + 0.5d0)*dr
          y = (j + yl)*dr
          z = (k + zl + 0.5d0)*dr
!          print*, "  ",x,y,z,"mp(BY) = 0", myid
        end if
        if(mp(BZ,i,j,k,1).eq.0.0d0) then
          x = (i + xl + 0.5d0)*dr
          y = (j + yl + 0.5d0)*dr
          z = (k + zl)*dr
!          print*, "  ",x,y,z,"mp(BZ) = 0", myid
        end if
      end do
      end do
      end do


!--------------------
      z0 = 1.0d0/cv
!
      pml%nxl = max(pmlxl,0)
      pml%nxu = max(nx-pmlxu,0)
      pml%nyl = max(pmlyl,0)
      pml%nyu = max(ny-pmlyu,0)
      pml%nzl = max(pmlzl,0)
      pml%nzu = max(nz-pmlzu,0)
!
      if(pml%nxl.gt.0) then
        smax0(1,1) = -log(1.0d-6)*(pmlord+1)*0.5d0/z0/(pml%nxl*dr)
      else
        smax0(1,1) = 0.0d0
      end if
      if(pml%nxu.gt.0) then
        smax0(2,1) = -log(1.0d-6)*(pmlord+1)*0.5d0/z0/(pml%nxu*dr)
      else
        smax0(2,1) = 0.0d0
      end if
      if(pml%nyl.gt.0) then
        smax0(1,2) = -log(1.0d-6)*(pmlord+1)*0.5d0/z0/(pml%nyl*dr)
      else
        smax0(1,2) = 0.0d0
      end if
      if(pml%nyu.gt.0) then
        smax0(2,2) = -log(1.0d-6)*(pmlord+1)*0.5d0/z0/(pml%nyu*dr)
      else
        smax0(2,2) = 0.0d0
      end if
      if(pml%nzl.gt.0) then
        smax0(1,3) = -log(1.0d-6)*(pmlord+1)*0.5d0/z0/(pml%nzl*dr)
      else
        smax0(1,3) = 0.0d0
      end if
      if(pml%nzu.gt.0) then
        smax0(2,3) = -log(1.0d-6)*(pmlord+1)*0.5d0/z0/(pml%nzu*dr)
      else
        smax0(2,3) = 0.0d0
      end if
!
      do i=-1,nx+1
        if(pml%nxl.gt.0.and.i.lt.pmlxl) then
          sigmaxe = min(real(pml%nxl-i,8)/pml%nxl,1.0d0)**pmlord*smax0(1,1)
          sigmaxm = min((pml%nxl-i-0.5d0)/pml%nxl,1.0d0)**pmlord*smax0(1,1)
        else if(pml%nxu.gt.0.and.i.ge.pmlxu) then
          sigmaxe = min(real(i-nx+pml%nxu,8)/pml%nxu,1.0d0)**pmlord*smax0(2,1)
          sigmaxm = min((i-nx+pml%nxu+0.5d0)/pml%nxu,1.0d0)**pmlord*smax0(2,1)
        else
          sigmaxe = 0.0d0
          sigmaxm = 0.0d0
        end if
        mediatr = 0.5d0*sigmaxe*dt
        pmlc(1,i,1) = (1.0d0-mediatr)/(1.0d0+mediatr)
        pmlc(2,i,1) = cv*cv*dt/dr/(1.0d0+mediatr)
        mediatr = 0.25d0*sigmaxm*dt
        pmlc(1,i,2) = (1.0d0-mediatr)/(1.0d0+mediatr)
        pmlc(2,i,2) = 0.5d0*dt/dr/(1.0d0+mediatr)
      end do
!
      do j=-1,ny+1
        if(pml%nyl.gt.0.and.j.lt.pmlyl) then
          sigmaye = min(real(pml%nyl-j,8)/pml%nyl,1.0d0)**pmlord*smax0(1,2)
          sigmaym = min((pml%nyl-j-0.5d0)/pml%nyl,1.0d0)**pmlord*smax0(1,2)
       else if(pml%nyu.gt.0.and.j.ge.pmlyu) then
          sigmaye = min(real(j-ny+pml%nyu,8)/pml%nyu,1.0d0)**pmlord*smax0(2,2)
          sigmaym = min((j-ny+pml%nyu+0.5d0)/pml%nyu,1.0d0)**pmlord*smax0(2,2)
        else
          sigmaye = 0.0d0
          sigmaym = 0.0d0
        end if
        mediatr = 0.5d0*sigmaye*dt
        pmlc(3,j,1) = (1.0d0-mediatr)/(1.0d0+mediatr)
        pmlc(4,j,1) = cv*cv*dt/dr/(1.0d0+mediatr)
        mediatr = 0.25d0*sigmaym*dt
        pmlc(3,j,2) = (1.0d0-mediatr)/(1.0d0+mediatr)
        pmlc(4,j,2) = 0.5d0*dt/dr/(1.0d0+mediatr)
      end do
!
      do k=-1,nz+1
        if(pml%nzl.gt.0.and.k.lt.pmlzl) then
          sigmaze = min(real(pml%nzl-k,8)/pml%nzl,1.0d0)**pmlord*smax0(1,3)
          sigmazm = min((pml%nzl-k-0.5d0)/pml%nzl,1.0d0)**pmlord*smax0(1,3)
        else if(pml%nzu.gt.0.and.k.ge.pmlzu) then
          sigmaze = min(real(k-nz+pml%nzu,8)/pml%nzu,1.0d0)**pmlord*smax0(2,3)
          sigmazm = min((k-nz+pml%nzu+0.5d0)/pml%nzu,1.0d0)**pmlord*smax0(2,3)
        else
          sigmaze = 0.0d0
          sigmazm = 0.0d0
        end if
        mediatr = 0.5d0*sigmaze*dt
        pmlc(5,k,1) = (1.0d0-mediatr)/(1.0d0+mediatr)
        pmlc(6,k,1) = cv*cv*dt/dr/(1.0d0+mediatr)
        mediatr = 0.25d0*sigmazm*dt
        pmlc(5,k,2) = (1.0d0-mediatr)/(1.0d0+mediatr)
        pmlc(6,k,2) = 0.5d0*dt/dr/(1.0d0+mediatr)
      end do


    return
  end subroutine medium
