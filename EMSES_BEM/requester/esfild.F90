#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
#include "oh_stats.h"
!
  subroutine esfld1
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   E S F I L D
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .   this subroutine gives a solution of electrostatic      .
!   .   by solving the Poisson's equation.                     .
!   ............................................................
!
!-------------------- parameter and common blocks
  use oh_type
  use paramt
  use allcom
  use interf
! #define MCW MPI_COMM_WORLD
#define MCW CTCA_subcomm
  implicit none
!
  integer(kind=4) :: i,j,k, i1,j1,k1, ii,jj,kk, ipc, icap
  integer(kind=4) :: xl,xu, yl,yu, zl,zu
!  integer(kind=4) :: siu,sju,sku
  integer(kind=4) :: xlower,xupper, ylower,yupper, zlower,zupper
  integer(kind=4) :: ierr
  real(kind=8) :: xlocal,ylocal,zlocal
  real(kind=8) :: x1,y1,z1, z2, xy1,xz1,yz1, xz2,yz2
  real(kind=8) :: v1,v2,v3,v4,v5,v6,v7,v8


!-------------------- 
      xl = sdoms(1,1,sdid(1)+1); xu = sdoms(2,1,sdid(1)+1)
      yl = sdoms(1,2,sdid(1)+1); yu = sdoms(2,2,sdid(1)+1)
      zl = sdoms(1,3,sdid(1)+1); zu = sdoms(2,3,sdid(1)+1)
!      if(pftmode.eq.3.and.myid.ne.snode-1) then
!        siu = sxu - sxl - 1
!      else
!        siu = sxu - sxl
!      end if
!      if(pftmode.eq.2.and.myid.ne.snode-1) then
!        sju = syu - syl - 1
!      else
!        sju = syu - syl
!      end if
!      if((pftmode.eq.1.or.pftmode.eq.4).and.myid.ne.snode-1) then
!        sku = szu - szl - 1
!      else
!        sku = szu - szl
!      end if


!-------------------- zero clear
      phi(:,:,:,:,:) = 0.0d0
      poi(:,:,:,:,:) = 0.0d0


!-------------------- Poisson Solver
      if(emmode.and.istep.gt.0) then
        rho(1,:,:,:,3:4) = rhobk(2,:,:,:,1:2)
        call poisson(3,rho(1,:,:,:,3:4))
      else
        call poisson(3,rho(1,:,:,:,1:2))
      end if


!-------------------- 
!      ngref(1) = 0
!      phiref(1) = 0.0d0
!      if(mtd_vbnd(1).ne.1.or.mtd_vbnd(2).ne.1.or.mtd_vbnd(3).ne.1) then
!        xlower = 0
!        if(bared(2,1,myid+1).eq.1) then
!          xupper = xu - xl
!        else
!          xupper = xu - xl - 1
!        end if
!
!        ylower = 0
!        if(bared(2,2,myid+1).eq.1) then
!          yupper = yu - yl
!        else
!          yupper = yu - yl - 1
!        end if
!
!        zlower = 0
!        if(bared(2,3,myid+1).eq.1) then
!          zupper = zu - zl
!        else
!          zupper = zu - zl - 1
!        end if
!
!        if(iphiref(1,1).eq.1.and.bared(1,1,myid+1).eq.1) then
!          do k=zlower,zupper
!          do j=ylower,yupper
!            ngref(1) = ngref(1) + 1
!            phiref(1) = phiref(1) + phi(1,0 ,j,k,1)
!          end do
!          end do
!          xlower = 1
!        end if
!
!        if(iphiref(2,1).eq.1.and.bared(2,1,myid+1).eq.1) then
!          do k=zlower,zupper
!          do j=ylower,yupper
!            ngref(1) = ngref(1) + 1
!            phiref(1) = phiref(1) + phi(1,xu-xl,j,k,1)
!          end do
!          end do
!          xupper = xu - xl - 1
!        end if
!
!        if(iphiref(1,2).eq.1.and.bared(1,2,myid+1).eq.1) then
!          do k=zlower,zupper
!          do i=xlower,xupper
!            ngref(1) = ngref(1) + 1
!            phiref(1) = phiref(1) + phi(1,i,0 ,k,1)
!          end do
!          end do
!          ylower = 1
!        end if
!
!        if(iphiref(2,2).eq.1.and.bared(2,2,myid+1).eq.1) then
!          do k=zlower,zupper
!          do i=xlower,xupper
!            ngref(1) = ngref(1) + 1
!            phiref(1) = phiref(1) + phi(1,i,yu-yl,k,1)
!          end do
!          end do
!          yupper = yu - yl - 1
!        end if
!
!        if(iphiref(1,3).eq.1.and.bared(1,3,myid+1).eq.1) then
!          do j=ylower,yupper
!          do i=xlower,xupper
!            ngref(1) = ngref(1) + 1
!            phiref(1) = phiref(1) + phi(1,i,j,0 ,1)
!          end do
!          end do
!        end if
!
!        if(iphiref(2,3).eq.1.and.bared(2,3,myid+1).eq.1) then
!          do j=ylower,yupper
!          do i=xlower,xupper
!            ngref(1) = ngref(1) + 1
!            phiref(1) = phiref(1) + phi(1,i,j,zu-zl,1)
!          end do
!          end do
!        end if
!
!        call MPI_Reduce(ngref(1),ngref(2),1,MPI_INTEGER,MPI_SUM,0,MCW,ierr)
!        call MPI_Reduce(phiref(1),phiref(2),1,MPI_REAL8,MPI_SUM,0,MCW,ierr)
!
!        if(myid.eq.0.and.ngref(2).ne.0) then
!          phiref(2) = phiref(2)/ngref(2)
!        else
!          phiref(2) = 0.0d0
!        end if
!      end if
!
!      do ipc=1,npprb
!        pprb(ipc,1) = 0.0d0
!        do k=dpprb(1,3,ipc),dpprb(2,3,ipc)
!        do j=dpprb(1,2,ipc),dpprb(2,2,ipc)
!        do i=dpprb(1,1,ipc),dpprb(2,1,ipc)
!          if(i.ge.xl.and.i.lt.xu.and. &
!             j.ge.yl.and.j.lt.yu.and. &
!             k.ge.zl.and.k.lt.zu) then
!            pprb(ipc,1) = pprb(ipc,1) + phi(1,i-xl,j-yl,k-zl,1)
!          end if
!        end do
!        end do
!        end do
!      end do
!
!      call MPI_Reduce(pprb(1:npprb,1),pprb(1:npprb,2),npprb,MPI_REAL8,MPI_SUM,0,MCW,ierr)
!      if(myid.eq.0) then
!        where(ngpprb(1:npprb).gt.0)
!          pprb(1:npprb,2) = pprb(1:npprb,2)/ngpprb(1:npprb)
!        elsewhere
!          pprb(1:npprb,2) = 0.0d0
!        end where
!      end if


  return
  end subroutine esfld1



!
  subroutine esfld2(ps)
!   ____________________________________________________________
!
!               S U B R O U T I N E   E S F I L D
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .   this subroutine gives a solution of electrostatic      .
!   .   by solving the Poisson's equation.                     .
!   ............................................................
!
!-------------------- parameter and common blocks
  use oh_type
  use paramt
  use allcom
  implicit none
!
  integer(kind=4) :: i,j,k, i1,j1,k1
  integer(kind=4) :: ipc, icap
  integer(kind=4) :: xl,yl,zl, xu,yu,zu
  integer(kind=4) :: ps
  real(kind=8) :: xlocal,ylocal,zlocal
  real(kind=8) :: x1,y1,z1, x2,y2,z2, etmp
  real(kind=8) :: vfactor, disp, r1, r2


!-------------------- 
      xl = sdoms(1,1,sdid(ps)+1)
      yl = sdoms(1,2,sdid(ps)+1)
      zl = sdoms(1,3,sdid(ps)+1)
      xu = sdoms(2,1,sdid(ps)+1) - sdoms(1,1,sdid(ps)+1)
      yu = sdoms(2,2,sdid(ps)+1) - sdoms(1,2,sdid(ps)+1)
      zu = sdoms(2,3,sdid(ps)+1) - sdoms(1,3,sdid(ps)+1)


!-------------------- smoothing potentials inside bodies
      if(sfecrrct.eq.2) then
        do ipc=1,npc
          if(boom(ipc)%align.eq.1) then
            vfactor = 0.5d0/log(2.0d0*boom(ipc)%hlength/boom(ipc)%rradius)
            do icap=nscpmx(ipc),nmxcpmx(ipc)-1
              disp = bdygrid(1,icap) - boom(ipc)%origin(1)
              r1 = sqrt((disp - boom(ipc)%hlength)*(disp - boom(ipc)%hlength) &
             &          + boom(ipc)%eradius*boom(ipc)%eradius)
              r2 = sqrt((disp + boom(ipc)%hlength)*(disp + boom(ipc)%hlength) &
             &          + boom(ipc)%eradius*boom(ipc)%eradius)
              i = bdygrid(1,icap) - xl
              j = bdygrid(2,icap) - yl
              k = bdygrid(3,icap) - zl
              if(i.ge.-1.and.i.le.xu+1.and. &
             &   j.ge.-1.and.j.le.yu+1.and. &
             &   k.ge.-1.and.k.le.zu+1) then
                phi(1,i,j,k,ps) = selfp(ipc)*vfactor &
             &            *log((+boom(ipc)%hlength - disp + r1) &
             &                /(-boom(ipc)%hlength - disp + r2))
              end if
            end do
          else if(boom(ipc)%align.eq.2) then
            vfactor = 0.5d0/log(2.0d0*boom(ipc)%hlength/boom(ipc)%rradius)
            do icap=nscpmx(ipc),nmxcpmx(ipc)-1
              disp = bdygrid(2,icap) - boom(ipc)%origin(2)
              r1 = sqrt((disp - boom(ipc)%hlength)*(disp - boom(ipc)%hlength) &
             &          + boom(ipc)%eradius*boom(ipc)%eradius)
              r2 = sqrt((disp + boom(ipc)%hlength)*(disp + boom(ipc)%hlength) &
             &          + boom(ipc)%eradius*boom(ipc)%eradius)
              i = bdygrid(1,icap) - xl
              j = bdygrid(2,icap) - yl
              k = bdygrid(3,icap) - zl
              if(i.ge.-1.and.i.le.xu+1.and. &
             &   j.ge.-1.and.j.le.yu+1.and. &
             &   k.ge.-1.and.k.le.zu+1) then
                phi(1,i,j,k,ps) = selfp(ipc)*vfactor &
             &            *log((+boom(ipc)%hlength - disp + r1) &
             &                /(-boom(ipc)%hlength - disp + r2))
              end if
            end do
          else if(boom(ipc)%align.eq.3) then
            vfactor = 0.5d0/log(2.0d0*boom(ipc)%hlength/boom(ipc)%rradius)
            do icap=nscpmx(ipc),nmxcpmx(ipc)-1
              disp = bdygrid(3,icap) - boom(ipc)%origin(3)
              r1 = sqrt((disp - boom(ipc)%hlength)*(disp - boom(ipc)%hlength) &
             &          + boom(ipc)%eradius*boom(ipc)%eradius)
              r2 = sqrt((disp + boom(ipc)%hlength)*(disp + boom(ipc)%hlength) &
             &          + boom(ipc)%eradius*boom(ipc)%eradius)
              i = bdygrid(1,icap) - xl
              j = bdygrid(2,icap) - yl
              k = bdygrid(3,icap) - zl
              if(i.ge.-1.and.i.le.xu+1.and. &
             &   j.ge.-1.and.j.le.yu+1.and. &
             &   k.ge.-1.and.k.le.zu+1) then
                phi(1,i,j,k,ps) = selfp(ipc)*vfactor &
             &            *log((+boom(ipc)%hlength - disp + r1) &
             &                /(-boom(ipc)%hlength - disp + r2))
              end if
            end do
          else
            do icap=nscpmx(ipc),nmxcpmx(ipc)-1
              i = bdygrid(1,icap) - xl
              j = bdygrid(2,icap) - yl
              k = bdygrid(3,icap) - zl
              if(i.ge.-1.and.i.le.xu+1.and. &
             &   j.ge.-1.and.j.le.yu+1.and. &
             &   k.ge.-1.and.k.le.zu+1) then
                phi(1,i,j,k,ps) = selfp(ipc)
              end if
            end do
          end if
        end do
      end if


!-------------------- zero clear
      eb(EX:EZ,:,:,:,ps) = 0.0d0


!-------------------- electrostatic field
      do k=-1,zu
      do j=-1,yu
      do i=-1,xu
        if(j.ne.-1.and.k.ne.-1) then
          eb(EX,i,j,k,ps) = &
       &       + ( phi(1,i,j,k,ps) - phi(1,i+1,j,k,ps) )*mltstp
        end if
        if(k.ne.-1.and.i.ne.-1) then
          eb(EY,i,j,k,ps) = &
       &       + ( phi(1,i,j,k,ps) - phi(1,i,j+1,k,ps) )*mltstp
        end if
        if(i.ne.-1.and.j.ne.-1) then
          eb(EZ,i,j,k,ps) = &
       &       + ( phi(1,i,j,k,ps) - phi(1,i,j,k+1,ps) )*mltstp
        end if
      end do
      end do
      end do


!-------------------- smoothing potentials inside bodies
      if(sfecrrct.eq.2) then
        do ipc=1,npc
          do icap=nmxcpmx(ipc),nmycpmx(ipc)-1
            xlocal = bdygrid(1,icap) - xl
            ylocal = bdygrid(2,icap) - yl
            zlocal = bdygrid(3,icap) - zl
!
            if(xlocal.ge.-1.and.xlocal.le.xu.and. &
           &   ylocal.ge. 0.and.ylocal.le.yu.and. &
           &   zlocal.ge. 0.and.zlocal.le.zu) then
              i = floor(xlocal)
              j = floor(ylocal)
              k = floor(zlocal)
!
              i1 = i + 1
              x1 = xlocal - i
              x2 = 1.0d0 - x1
!
              if(phi(1,i ,j,k,ps).eq.selfp(ipc)) then
                eb(EX,i,j,k,ps) = eb(EX,i,j,k,ps)/x2
              end if
!
              if(phi(1,i1,j,k,ps).eq.selfp(ipc)) then
                eb(EX,i,j,k,ps) = eb(EX,i,j,k,ps)/x1
              end if
            end if
          end do
!
          do icap=nmycpmx(ipc),nmzcpmx(ipc)-1
            xlocal = bdygrid(1,icap) - xl
            ylocal = bdygrid(2,icap) - yl
            zlocal = bdygrid(3,icap) - zl
!
            if(xlocal.ge. 0.and.xlocal.le.xu.and. &
           &   ylocal.ge.-1.and.ylocal.le.yu.and. &
           &   zlocal.ge. 0.and.zlocal.le.zu) then
              i = floor(xlocal)
              j = floor(ylocal)
              k = floor(zlocal)
!
              j1 = j + 1
              y1 = ylocal - j
              y2 = 1.0d0 - y1
!
              if(phi(1,i,j ,k,ps).eq.selfp(ipc)) then
                eb(EY,i,j,k,ps) = eb(EY,i,j,k,ps)/y2
              end if
!
              if(phi(1,i,j1,k,ps).eq.selfp(ipc)) then
                eb(EY,i,j,k,ps) = eb(EY,i,j,k,ps)/y1
              end if
            end if
          end do
!
          do icap=nmzcpmx(ipc),nmxycpmx(ipc)-1
            xlocal = bdygrid(1,icap) - xl
            ylocal = bdygrid(2,icap) - yl
            zlocal = bdygrid(3,icap) - zl
!
            if(xlocal.ge. 0.and.xlocal.le.xu.and. &
           &   ylocal.ge. 0.and.ylocal.le.yu.and. &
           &   zlocal.ge.-1.and.zlocal.le.zu) then
              i = floor(xlocal)
              j = floor(ylocal)
              k = floor(zlocal)
!
              k1 = k + 1
              z1 = zlocal - k
              z2 = 1.0d0 - z1
!
              if(phi(1,i,j,k ,ps).eq.selfp(ipc)) then
                eb(EZ,i,j,k,ps) = eb(EZ,i,j,k,ps)/z2
              end if
!
              if(phi(1,i,j,k1,ps).eq.selfp(ipc)) then
                eb(EZ,i,j,k,ps) = eb(EZ,i,j,k,ps)/z1
              end if
            end if
          end do
        end do
      end if


!--------------------- low pass filter
!      if(nfltrx.ge.1) then
!        call fltrf3(eb(EX,:,:,:,ps),work, &
!     &            ix,iy,iz,nx,ny,nz,nfltrx,nfltry,nfltrz,ixy)
!      end if
!      if(nfltry.ge.1) then
!        call fltrf3(eb(EY,:,:,:,ps),work, &
!     &            ix,iy,iz,nx,ny,nz,nfltrx,nfltry,nfltrz,ixy)
!      end if
!      if(nfltrz.ge.1) then
!        call fltrf3(eb(EZ,:,:,:,ps),work, &
!     &            ix,iy,iz,nx,ny,nz,nfltrx,nfltry,nfltrz,ixy)
!      end if
!--------------------- masking of longitudinal e-field
!      call fsmask(6)
!--------------------- boundary treatment of pe field
!      call fbound(3) 


!--------------------- 
!      if(mode_dipole.ne.0) then
!        do ig=1,ngap
!          kp(ig)=k_gap(ig)+nz0
!        end do
!#include "vical.fnc"
!      end if


  return
  end subroutine esfld2



!
  subroutine esfld3(ps)
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   E S F I L D
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .   this subroutine gives a solution of electrostatic      .
!   .   by solving the Poisson's equation.                     .
!   ............................................................
!
!-------------------- parameter and common blocks
  use oh_type
  use paramt
  use allcom
  implicit none
!
  integer(kind=4) :: i,j,k
  integer(kind=4) :: xu,yu,zu
  integer(kind=4) :: ps


!-------------------- 
      xu = sdoms(2,1,sdid(ps)+1) - sdoms(1,1,sdid(ps)+1)
      yu = sdoms(2,2,sdid(ps)+1) - sdoms(1,2,sdid(ps)+1)
      zu = sdoms(2,3,sdid(ps)+1) - sdoms(1,3,sdid(ps)+1)


!-------------------- electrostatic field
      do k=-1,zu
      do j=-1,yu
      do i=-1,xu
        if(j.ne.-1.and.k.ne.-1) then
          eb(EX,i,j,k,ps) = eb(EX,i,j,k,ps) &
     &         + ( phi(1,i,j,k,ps) - phi(1,i+1,j,k,ps) )*mltstp
        end if
        if(k.ne.-1.and.i.ne.-1) then
          eb(EY,i,j,k,ps) = eb(EY,i,j,k,ps) &
     &         + ( phi(1,i,j,k,ps) - phi(1,i,j+1,k,ps) )*mltstp
        end if
        if(i.ne.-1.and.j.ne.-1) then
          eb(EZ,i,j,k,ps) = eb(EZ,i,j,k,ps) &
     &         + ( phi(1,i,j,k,ps) - phi(1,i,j,k+1,ps) )*mltstp
        end if
      end do
      end do
      end do


!--------------------- low pass filter
!      if(nfltrx.ge.1) then
!        call fltrf3(eb(EX,:,:,:,ps),work, &
!     &            ix,iy,iz,nx,ny,nz,nfltrx,nfltry,nfltrz,ixy)
!      end if
!      if(nfltry.ge.1) then
!        call fltrf3(eb(EY,:,:,:,ps),work, &
!     &            ix,iy,iz,nx,ny,nz,nfltrx,nfltry,nfltrz,ixy)
!      end if
!      if(nfltrz.ge.1) then
!        call fltrf3(eb(EZ,:,:,:,ps),work, &
!     &            ix,iy,iz,nx,ny,nz,nfltrx,nfltry,nfltrz,ixy)
!      end if
!--------------------- masking of longitudinal e-field
!      call fsmask(6)
!--------------------- boundary treatment of pe field
!      call fbound(3) 


!--------------------- 
!      if(mode_dipole.ne.0) then
!        do ig=1,ngap
!          kp(ig)=k_gap(ig)+nz0
!        end do
!#include "vical.fnc"
!      end if


  return
  end subroutine esfld3
