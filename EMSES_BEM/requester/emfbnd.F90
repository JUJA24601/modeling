#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine boundary_emfld(psd)
!
!   ____________________________________________________________
!
!                       S U B R O U T I N E
!                   B O U N D A R Y _ E M F L D
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   ............................................................
!
!-------------------- parameter and common blocks
  use oh_type
  use paramt
  use allcom
  implicit none
!
  integer(kind=4) :: i,j,k
  integer(kind=4) :: xl,xu, yl,yu, zl,zu
  integer(kind=4) :: ps, psd


!-------------------- test particle simulation
      if(juncan.ge.1000) return


!-------------------- 
      ps = modulo(psd-1,2) + 1
      xl = sdoms(1,1,sdid(ps)+1); xu = sdoms(2,1,sdid(ps)+1)
      yl = sdoms(1,2,sdid(ps)+1); yu = sdoms(2,2,sdid(ps)+1)
      zl = sdoms(1,3,sdid(ps)+1); zu = sdoms(2,3,sdid(ps)+1)


!-------------------- 
      if(bounds(1,3,sdid(ps)+1).eq.2.and.nfbnd(3).eq.0.and.nz.eq.1) then
        do j=-1,yu-yl+1
        do i=-1,xu-xl+1
!         ----------- zl-1-zl <- zu-1-zl
          eb(EX,i,j,-1,psd) = +eb(EX,i,j,zu-1-zl,psd)
          eb(EY,i,j,-1,psd) = +eb(EY,i,j,zu-1-zl,psd)
          eb(EZ,i,j,-1,psd) = +eb(EZ,i,j,zu-1-zl,psd)
          eb(BX,i,j,-1,psd) = +eb(BX,i,j,zu-1-zl,psd)
          eb(BY,i,j,-1,psd) = +eb(BY,i,j,zu-1-zl,psd)
          eb(BZ,i,j,-1,psd) = +eb(BZ,i,j,zu-1-zl,psd)
        end do
        end do
      else if(bounds(1,3,sdid(ps)+1).eq.2) then
        do j=-1,yu-yl+1
        do i=-1,xu-xl+1
          if(-1.ge.zl-1.and.-1.le.zu+1) then
            eb(EX,i,j,-1  -zl,psd) = -eb(EX,i,j,+1  -zl,psd)
            eb(EY,i,j,-1  -zl,psd) = -eb(EY,i,j,+1  -zl,psd)
            eb(BZ,i,j,-1  -zl,psd) = -eb(BZ,i,j,+1  -zl,psd)
!
            eb(EZ,i,j,-1  -zl,psd) = +eb(EZ,i,j, 0  -zl,psd)
            eb(BX,i,j,-1  -zl,psd) = +eb(BX,i,j, 0  -zl,psd)
            eb(BY,i,j,-1  -zl,psd) = +eb(BY,i,j, 0  -zl,psd)
          end if
!
          if( 0.ge.zl-1.and. 0.le.zu+1) then
            eb(EX,i,j, 0  -zl,psd) = 0.0d0
            eb(EY,i,j, 0  -zl,psd) = 0.0d0
            eb(BZ,i,j, 0  -zl,psd) = 0.0d0
          end if
        end do
        end do
      end if

!-------------------- 
      if(bounds(2,3,sdid(ps)+1).eq.2.and.nfbnd(3).eq.0.and.nz.eq.1) then
        do j=-1,yu-yl+1
        do i=-1,xu-xl+1
!         ----------- zu-zl <- zl-zl
          eb(EX,i,j,zu-zl,psd) = +eb(EX,i,j,0,psd)
          eb(EY,i,j,zu-zl,psd) = +eb(EY,i,j,0,psd)
          eb(EZ,i,j,zu-zl,psd) = +eb(EZ,i,j,0,psd)
          eb(BX,i,j,zu-zl,psd) = +eb(BX,i,j,0,psd)
          eb(BY,i,j,zu-zl,psd) = +eb(BY,i,j,0,psd)
          eb(BZ,i,j,zu-zl,psd) = +eb(BZ,i,j,0,psd)
!         ----------- zu+1-zl <- zl+1-zl
          eb(EX,i,j,zu+1-zl,psd) = +eb(EX,i,j,+1,psd)
          eb(EY,i,j,zu+1-zl,psd) = +eb(EY,i,j,+1,psd)
          eb(EZ,i,j,zu+1-zl,psd) = +eb(EZ,i,j,+1,psd)
          eb(BX,i,j,zu+1-zl,psd) = +eb(BX,i,j,+1,psd)
          eb(BY,i,j,zu+1-zl,psd) = +eb(BY,i,j,+1,psd)
          eb(BZ,i,j,zu+1-zl,psd) = +eb(BZ,i,j,+1,psd)
        end do
        end do
      else if(bounds(2,3,sdid(ps)+1).eq.2) then
        do j=-1,yu-yl+1
        do i=-1,xu-xl+1
          if(nz.ge.zl-1.and.nz.le.zu+1) then
            eb(EX,i,j,nz  -zl,psd) = 0.0d0
            eb(EY,i,j,nz  -zl,psd) = 0.0d0
            eb(BZ,i,j,nz  -zl,psd) = 0.0d0
!
            eb(EZ,i,j,nz  -zl,psd) = +eb(EZ,i,j,nz-1-zl,psd)
            eb(BX,i,j,nz  -zl,psd) = +eb(BX,i,j,nz-1-zl,psd)
            eb(BY,i,j,nz  -zl,psd) = +eb(BY,i,j,nz-1-zl,psd)
          end if
!
          if(nz+1.ge.zl-1.and.nz+1.le.zu+1) then
            eb(EX,i,j,nz+1-zl,psd) = -eb(EX,i,j,nz-1-zl,psd)
            eb(EY,i,j,nz+1-zl,psd) = -eb(EY,i,j,nz-1-zl,psd)
            eb(BZ,i,j,nz+1-zl,psd) = -eb(BZ,i,j,nz-1-zl,psd)
!
            eb(EZ,i,j,nz+1-zl,psd) = +eb(EZ,i,j,nz-2-zl,psd)
            eb(BX,i,j,nz+1-zl,psd) = +eb(BX,i,j,nz-2-zl,psd)
            eb(BY,i,j,nz+1-zl,psd) = +eb(BY,i,j,nz-2-zl,psd)
          end if
        end do
        end do
      end if

!-------------------- 
      if(bounds(1,2,sdid(ps)+1).eq.2.and.nfbnd(2).eq.0.and.ny.eq.1) then
        do k=-1,zu-zl+1
        do i=-1,xu-xl+1
!         ----------- yl-1-yl <- yu-1-yl
          eb(EX,i,-1,k,psd) = +eb(EX,i,yu-1-yl,k,psd)
          eb(EY,i,-1,k,psd) = +eb(EY,i,yu-1-yl,k,psd)
          eb(EZ,i,-1,k,psd) = +eb(EZ,i,yu-1-yl,k,psd)
          eb(BX,i,-1,k,psd) = +eb(BX,i,yu-1-yl,k,psd)
          eb(BY,i,-1,k,psd) = +eb(BY,i,yu-1-yl,k,psd)
          eb(BZ,i,-1,k,psd) = +eb(BZ,i,yu-1-yl,k,psd)
        end do
        end do
      else if(bounds(1,2,sdid(ps)+1).eq.2) then
        do k=-1,zu-zl+1
        do i=-1,xu-xl+1
          if(-1.ge.yl-1.and.-1.le.yu+1) then
            eb(EX,i,-1  -yl,k,psd) = -eb(EX,i,+1  -yl,k,psd)
            eb(EZ,i,-1  -yl,k,psd) = -eb(EZ,i,+1  -yl,k,psd)
            eb(BY,i,-1  -yl,k,psd) = -eb(BY,i,+1  -yl,k,psd)
!
            eb(EY,i,-1  -yl,k,psd) = +eb(EY,i, 0  -yl,k,psd)
            eb(BX,i,-1  -yl,k,psd) = +eb(BX,i, 0  -yl,k,psd)
            eb(BZ,i,-1  -yl,k,psd) = +eb(BZ,i, 0  -yl,k,psd)
          end if
!
          if( 0.ge.yl-1.and. 0.le.yu+1) then
            eb(EX,i, 0  -yl,k,psd) = 0.0d0
            eb(EZ,i, 0  -yl,k,psd) = 0.0d0
            eb(BY,i, 0  -yl,k,psd) = 0.0d0
          end if
        end do
        end do
      end if

!-------------------- 
      if(bounds(2,2,sdid(ps)+1).eq.2.and.nfbnd(2).eq.0.and.ny.eq.1) then
        do k=-1,zu-zl+1
        do i=-1,xu-xl+1
!         ----------- yu-yl <- yl-yl
          eb(EX,i,yu-yl,k,psd) = +eb(EX,i,0,k,psd)
          eb(EY,i,yu-yl,k,psd) = +eb(EY,i,0,k,psd)
          eb(EZ,i,yu-yl,k,psd) = +eb(EZ,i,0,k,psd)
          eb(BX,i,yu-yl,k,psd) = +eb(BX,i,0,k,psd)
          eb(BY,i,yu-yl,k,psd) = +eb(BY,i,0,k,psd)
          eb(BZ,i,yu-yl,k,psd) = +eb(BZ,i,0,k,psd)
!         ----------- yu+1-yl <- yl+1-yl
          eb(EX,i,yu+1-yl,k,psd) = +eb(EX,i,+1,k,psd)
          eb(EY,i,yu+1-yl,k,psd) = +eb(EY,i,+1,k,psd)
          eb(EZ,i,yu+1-yl,k,psd) = +eb(EZ,i,+1,k,psd)
          eb(BX,i,yu+1-yl,k,psd) = +eb(BX,i,+1,k,psd)
          eb(BY,i,yu+1-yl,k,psd) = +eb(BY,i,+1,k,psd)
          eb(BZ,i,yu+1-yl,k,psd) = +eb(BZ,i,+1,k,psd)
        end do
        end do
      else if(bounds(2,2,sdid(ps)+1).eq.2) then
        do k=-1,zu-zl+1
        do i=-1,xu-xl+1
          if(ny.ge.yl-1.and.ny.le.yu+1) then
            eb(EX,i,ny  -yl,k,psd) = 0.0d0
            eb(EZ,i,ny  -yl,k,psd) = 0.0d0
            eb(BY,i,ny  -yl,k,psd) = 0.0d0
!
            eb(EY,i,ny  -yl,k,psd) = +eb(EY,i,ny-1-yl,k,psd)
            eb(BX,i,ny  -yl,k,psd) = +eb(BX,i,ny-1-yl,k,psd)
            eb(BZ,i,ny  -yl,k,psd) = +eb(BZ,i,ny-1-yl,k,psd)
          end if
!
          if(ny+1.ge.yl-1.and.ny+1.le.yu+1) then
            eb(EX,i,ny+1-yl,k,psd) = -eb(EX,i,ny-1-yl,k,psd)
            eb(EZ,i,ny+1-yl,k,psd) = -eb(EZ,i,ny-1-yl,k,psd)
            eb(BY,i,ny+1-yl,k,psd) = -eb(BY,i,ny-1-yl,k,psd)
!
            eb(EY,i,ny+1-yl,k,psd) = +eb(EY,i,ny-2-yl,k,psd)
            eb(BX,i,ny+1-yl,k,psd) = +eb(BX,i,ny-2-yl,k,psd)
            eb(BZ,i,ny+1-yl,k,psd) = +eb(BZ,i,ny-2-yl,k,psd)
          end if
        end do
        end do
      end if

!-------------------- 
      if(bounds(1,1,sdid(ps)+1).eq.2.and.nfbnd(1).eq.0.and.nx.eq.1) then
        do k=-1,zu-zl+1
        do j=-1,yu-yl+1
!         ----------- xl-1-xl <- xu-1-xl
          eb(EX,-1,j,k,psd) = +eb(EX,xu-1-xl,j,k,psd)
          eb(EY,-1,j,k,psd) = +eb(EY,xu-1-xl,j,k,psd)
          eb(EZ,-1,j,k,psd) = +eb(EZ,xu-1-xl,j,k,psd)
          eb(BX,-1,j,k,psd) = +eb(BX,xu-1-xl,j,k,psd)
          eb(BY,-1,j,k,psd) = +eb(BY,xu-1-xl,j,k,psd)
          eb(BZ,-1,j,k,psd) = +eb(BZ,xu-1-xl,j,k,psd)
        end do
        end do
      else if(bounds(1,1,sdid(ps)+1).eq.2) then
        do k=-1,zu-zl+1
        do j=-1,yu-yl+1
          if(-1.ge.xl-1.and.-1.le.xu+1) then
            eb(EY,-1  -xl,j,k,psd) = -eb(EY,+1  -xl,j,k,psd)
            eb(EZ,-1  -xl,j,k,psd) = -eb(EZ,+1  -xl,j,k,psd)
            eb(BX,-1  -xl,j,k,psd) = -eb(BX,+1  -xl,j,k,psd)
!
            eb(EX,-1  -xl,j,k,psd) = +eb(EX, 0  -xl,j,k,psd)
            eb(BY,-1  -xl,j,k,psd) = +eb(BY, 0  -xl,j,k,psd)
            eb(BZ,-1  -xl,j,k,psd) = +eb(BZ, 0  -xl,j,k,psd)
          end if
!
          if( 0.ge.xl-1.and. 0.le.xu+1) then
            eb(EY, 0  -xl,j,k,psd) = 0.0d0
            eb(EZ, 0  -xl,j,k,psd) = 0.0d0
            eb(BX, 0  -xl,j,k,psd) = 0.0d0
          end if
        end do
        end do
      end if

!-------------------- 
      if(bounds(2,1,sdid(ps)+1).eq.2.and.nfbnd(1).eq.0.and.nx.eq.1) then
        do k=-1,zu-zl+1
        do j=-1,yu-yl+1
!         ----------- xu-xl <- xl-xl
          eb(EX,xu-xl,j,k,psd) = +eb(EX,0,j,k,psd)
          eb(EY,xu-xl,j,k,psd) = +eb(EY,0,j,k,psd)
          eb(EZ,xu-xl,j,k,psd) = +eb(EZ,0,j,k,psd)
          eb(BX,xu-xl,j,k,psd) = +eb(BX,0,j,k,psd)
          eb(BY,xu-xl,j,k,psd) = +eb(BY,0,j,k,psd)
          eb(BZ,xu-xl,j,k,psd) = +eb(BZ,0,j,k,psd)
!         ----------- xu+1-xl <- xl+1-xl
          eb(EX,xu+1-xl,j,k,psd) = +eb(EX,+1,j,k,psd)
          eb(EY,xu+1-xl,j,k,psd) = +eb(EY,+1,j,k,psd)
          eb(EZ,xu+1-xl,j,k,psd) = +eb(EZ,+1,j,k,psd)
          eb(BX,xu+1-xl,j,k,psd) = +eb(BX,+1,j,k,psd)
          eb(BY,xu+1-xl,j,k,psd) = +eb(BY,+1,j,k,psd)
          eb(BZ,xu+1-xl,j,k,psd) = +eb(BZ,+1,j,k,psd)
        end do
        end do
      else if(bounds(2,1,sdid(ps)+1).eq.2) then
        do k=-1,zu-zl+1
        do j=-1,yu-yl+1
          if(nx.ge.xl-1.and.nx.le.xu+1) then
            eb(EY,nx  -xl,j,k,psd) = 0.0d0
            eb(EZ,nx  -xl,j,k,psd) = 0.0d0
            eb(BX,nx  -xl,j,k,psd) = 0.0d0
!
            eb(EX,nx  -xl,j,k,psd) = +eb(EX,nx-1-xl,j,k,psd)
            eb(BY,nx  -xl,j,k,psd) = +eb(BY,nx-1-xl,j,k,psd)
            eb(BZ,nx  -xl,j,k,psd) = +eb(BZ,nx-1-xl,j,k,psd)
          end if
!
          if(nx+1.ge.xl-1.and.nx+1.le.xu+1) then
            eb(EY,nx+1-xl,j,k,psd) = -eb(EY,nx-1-xl,j,k,psd)
            eb(EZ,nx+1-xl,j,k,psd) = -eb(EZ,nx-1-xl,j,k,psd)
            eb(BX,nx+1-xl,j,k,psd) = -eb(BX,nx-1-xl,j,k,psd)
!
            eb(EX,nx+1-xl,j,k,psd) = +eb(EX,nx-2-xl,j,k,psd)
            eb(BY,nx+1-xl,j,k,psd) = +eb(BY,nx-2-xl,j,k,psd)
            eb(BZ,nx+1-xl,j,k,psd) = +eb(BZ,nx-2-xl,j,k,psd)
          end if
        end do
        end do
      end if


  return
  end subroutine boundary_emfld



!
  subroutine replace_boundary_field(psd)
!
!   ____________________________________________________________
!
!                       S U B R O U T I N E
!              A D D _ B O U N D A R Y _ C U R R E N T
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
  integer(kind=4) :: i,j,k
  integer(kind=4) :: xu,yu,zu
  integer(kind=4) :: sl,su, dl,du, nl,nu
  integer(kind=4) :: ps, psd


!-------------------- 
      ps = modulo(psd-1,2) + 1
      xu = sdoms(2,1,sdid(ps)+1) - sdoms(1,1,sdid(ps)+1)
      yu = sdoms(2,2,sdid(ps)+1) - sdoms(1,2,sdid(ps)+1)
      zu = sdoms(2,3,sdid(ps)+1) - sdoms(1,3,sdid(ps)+1)


!-------------------- 
      sl = ctypes(2,2,1,CAJ)	!=-2
      su = ctypes(2,1,1,CAJ)	!=+2
!
      nl = ctypes(3,2,1,CAJ)	!= 1
      nu = ctypes(3,1,1,CAJ)	!= 2
!
      dl = sl + nl		!=-1
      du = su - nu		!= 0

!-------------------- 
      if(bounds(1,3,sdid(ps)+1).eq.1) then
        do k=0,nl-1		!(i.e., do k=0,0)
        do j=0,yu+(su+nu-sl)-1	!(i.e., do j=0,yu+5)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+5)
          eb(EX,sl+i,sl+j,dl+k,psd) = eb(EX,sl+i,sl+j,sl+k,psd)
          eb(EY,sl+i,sl+j,dl+k,psd) = eb(EY,sl+i,sl+j,sl+k,psd)
          eb(EZ,sl+i,sl+j,dl+k,psd) = eb(EZ,sl+i,sl+j,sl+k,psd)
        end do
        end do
        end do
      end if

!-------------------- 
      if(bounds(2,3,sdid(ps)+1).eq.1) then
        do k=0,nu-1		!(i.e., do k=0,1)
        do j=0,yu+(su+nu-sl)-1	!(i.e., do j=0,yu+5)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+5)
          eb(EX,sl+i,sl+j,zu+du+k,psd) = eb(EX,sl+i,sl+j,zu+su+k,psd)
          eb(EY,sl+i,sl+j,zu+du+k,psd) = eb(EY,sl+i,sl+j,zu+su+k,psd)
          eb(EZ,sl+i,sl+j,zu+du+k,psd) = eb(EZ,sl+i,sl+j,zu+su+k,psd)
        end do
        end do
        end do
      end if

!-------------------- 
      if(bounds(1,2,sdid(ps)+1).eq.1) then
        do k=0,zu+(du+nu-dl)-1	!(i.e., do k=0,zu+2)
        do j=0,nl-1		!(i.e., do j=0,0)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+5)
          eb(EX,sl+i,dl+j,dl+k,psd) = eb(EX,sl+i,sl+j,dl+k,psd)
          eb(EY,sl+i,dl+j,dl+k,psd) = eb(EY,sl+i,sl+j,dl+k,psd)
          eb(EZ,sl+i,dl+j,dl+k,psd) = eb(EZ,sl+i,sl+j,dl+k,psd)
        end do
        end do
        end do
      end if

!-------------------- 
      if(bounds(2,2,sdid(ps)+1).eq.1) then
        do k=0,zu+(du+nu-dl)-1	!(i.e., do k=0,zu+2)
        do j=0,nu-1		!(i.e., do j=0,1)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+5)
          eb(EX,sl+i,yu+du+j,dl+k,psd) = eb(EX,sl+i,yu+su+j,dl+k,psd)
          eb(EY,sl+i,yu+du+j,dl+k,psd) = eb(EY,sl+i,yu+su+j,dl+k,psd)
          eb(EZ,sl+i,yu+du+j,dl+k,psd) = eb(EZ,sl+i,yu+su+j,dl+k,psd)
        end do
        end do
        end do
      end if

!-------------------- 
      if(bounds(1,1,sdid(ps)+1).eq.1) then
        do k=0,zu+(du+nu-dl)-1	!(i.e., do k=0,zu+2)
        do j=0,yu+(du+nu-dl)-1	!(i.e., do j=0,yu+2)
        do i=0,nl-1		!(i.e., do i=0,0)
          eb(EX,dl+i,dl+j,dl+k,psd) = eb(EX,sl+i,dl+j,dl+k,psd)
          eb(EY,dl+i,dl+j,dl+k,psd) = eb(EY,sl+i,dl+j,dl+k,psd)
          eb(EZ,dl+i,dl+j,dl+k,psd) = eb(EZ,sl+i,dl+j,dl+k,psd)
        end do
        end do
        end do
      end if

!-------------------- 
      if(bounds(2,1,sdid(ps)+1).eq.1) then
        do k=0,zu+(du+nu-dl)-1	!(i.e., do k=0,zu+2)
        do j=0,yu+(du+nu-dl)-1	!(i.e., do j=0,yu+2)
        do i=0,nu-1		!(i.e., do i=0,1)
          eb(EX,xu+du+i,dl+j,dl+k,psd) = eb(EX,xu+su+i,dl+j,dl+k,psd)
          eb(EY,xu+du+i,dl+j,dl+k,psd) = eb(EY,xu+su+i,dl+j,dl+k,psd)
          eb(EZ,xu+du+i,dl+j,dl+k,psd) = eb(EZ,xu+su+i,dl+j,dl+k,psd)
        end do
        end do
        end do
      end if


  return
  end subroutine replace_boundary_field
