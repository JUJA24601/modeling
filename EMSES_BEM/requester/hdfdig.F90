#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine hdfdig(func,idiag)
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   H D F D I G
!   ____________________________________________________________
!
!
!-------------------- parameter and common blocks
  use oh_type
  use paramt
  use allcom
  use hdf
  use interf
! #define MCW MPI_COMM_WORLD
#define MCW CTCA_subcomm
  implicit none
!
  integer,parameter :: PRECI=4
  integer(kind=8) :: m, nsP,nsS,neP,neS, neeP,neeS, nns,nne, nout
  integer(kind=8) :: lnoutp(ispec),noutp(ispec),poffset(ispec)
  integer(kind=4) :: func, idiag
  integer(kind=4) :: is, h5ind
  integer(kind=4) :: istat1,istat2
  integer(kind=4) :: i,j,k, ieb, ipc
  integer(kind=4) :: xl,xu, yl,yu, zl,zu
  integer(kind=4) :: xlower,xupper, ylower,yupper, zlower,zupper
  integer(kind=4) :: ierr, iid
  integer(kind=HSIZE_T) :: h5tab=0
  real(kind=8) :: cfbj
  real(kind=8) :: pout(1)
  real(kind=8) :: elemchg
  real(kind=8) :: massele
  real(kind=8) :: massion
  real(kind=8) :: epsiln0
  real(kind=8) :: boltfct
  real(kind=8) :: coe_v
  real(kind=8) :: coe_x
  real(kind=8) :: coe_t
  real(kind=8) :: coe_w
  real(kind=8) :: coe_ep
  real(kind=8) :: coe_qm
  real(kind=8) :: coe_j
  real(kind=8) :: coe_ph
  real(kind=8) :: coe_e
  real(kind=8) :: coe_m
  real(kind=8) :: expema
  character(len=30) :: plhdfn
  character(len=10) :: grname
  character(len=8) :: dsname

!-------------------- 
      xl = sdoms(1,1,sdid(1)+1); xu = sdoms(2,1,sdid(1)+1)
      yl = sdoms(1,2,sdid(1)+1); yu = sdoms(2,2,sdid(1)+1)
      zl = sdoms(1,3,sdid(1)+1); zu = sdoms(2,3,sdid(1)+1)


!-------------------- 
      if(idiag.gt.0.and.emarlxt.gt.0) then
        expema = exp(-1.0d0/emarlxt)
      else
        expema = 0.0d0
      end if


!--------------------
      if(func.eq.1) then
        if(idiag.lt.hdfdigstart.or.idiag.gt.hdfdigend) return
        elemchg = 1.602177d-19
        massele = 9.109382d-31
        massion = massele*massratio
        epsiln0 = 8.854188d-12
        boltfct = 1.380651E-23
        if(deltax.gt.0.0d0) then
          coe_v = cvsim/2.99792458d8
          coe_x = 1.0d0/deltax/1.0d-2
          coe_t = coe_x/coe_v
          coe_w = 1.0d0/coe_t
          coe_ep = 1.0d0/epsiln0
          coe_qm = massele/elemchg
          coe_j = coe_ep/coe_qm*coe_v*coe_v*coe_v/coe_x/coe_x
          coe_ph = coe_v*coe_v/coe_qm
          coe_e = coe_v*coe_v/coe_qm/coe_x
          coe_m = coe_v/coe_qm/coe_x
        else
          coe_v = 1.0d0
          coe_x = 1.0d0
          coe_t = 1.0d0
          coe_w = 1.0d0
          coe_ep = 1.0d0
          coe_qm = 1.0d0
          coe_j = 1.0d0
          coe_ph = 1.0d0
          coe_e = 1.0d0
          coe_m = 1.0d0
        end if
      end if


!-------------------- 
      if(func.eq.1) then
        if((ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0).or.daverg.ge.1) then
          phi(:,:,:,:,:) = 0.0d0
          poi(:,:,:,:,:) = 0.0d0
!
          call poisson(3,rho(1,:,:,:,1:2))
!
          do k=zl-1,zu+1
          do j=yl-1,yu+1
          do i=xl-1,xu+1
            phi(1,i-xl,j-yl,k-zl,1) = phi(1,i-xl,j-yl,k-zl,1) &
           &                        - (i - prefcrd(1))*e0x &
           &                        - (j - prefcrd(2))*e0y &
           &                        - (k - prefcrd(3))*e0z
          end do
          end do
          end do
        end if
!       call poichk(1)
!
!       ------------- 
        ngref(1) = 0
        phiref(1) = 0.0d0
        if(mtd_vbnd(1).ne.1.or.mtd_vbnd(2).ne.1.or.mtd_vbnd(3).ne.1) then
          xlower = 0
          if(bared(2,1,myid+1).eq.1) then
            xupper = xu - xl
          else
            xupper = xu - xl - 1
          end if
!
          ylower = 0
          if(bared(2,2,myid+1).eq.1) then
            yupper = yu - yl
          else
            yupper = yu - yl - 1
          end if
!
          zlower = 0
          if(bared(2,3,myid+1).eq.1) then
            zupper = zu - zl
          else
            zupper = zu - zl - 1
          end if
!
          if(iphiref(1,1).eq.1.and.bared(1,1,myid+1).eq.1) then
            do k=zlower,zupper
            do j=ylower,yupper
              ngref(1) = ngref(1) + 1
              phiref(1) = phiref(1) + phi(1,0 ,j,k,1)
            end do
            end do
            xlower = 1
          end if
!
          if(iphiref(2,1).eq.1.and.bared(2,1,myid+1).eq.1) then
            do k=zlower,zupper
            do j=ylower,yupper
              ngref(1) = ngref(1) + 1
              phiref(1) = phiref(1) + phi(1,xu-xl,j,k,1)
            end do
            end do
            xupper = xu - xl - 1
          end if
!
          if(iphiref(1,2).eq.1.and.bared(1,2,myid+1).eq.1) then
            do k=zlower,zupper
            do i=xlower,xupper
              ngref(1) = ngref(1) + 1
              phiref(1) = phiref(1) + phi(1,i,0 ,k,1)
            end do
            end do
            ylower = 1
          end if
!
          if(iphiref(2,2).eq.1.and.bared(2,2,myid+1).eq.1) then
            do k=zlower,zupper
            do i=xlower,xupper
              ngref(1) = ngref(1) + 1
              phiref(1) = phiref(1) + phi(1,i,yu-yl,k,1)
            end do
            end do
            yupper = yu - yl - 1
          end if
!
          if(iphiref(1,3).eq.1.and.bared(1,3,myid+1).eq.1) then
            do j=ylower,yupper
            do i=xlower,xupper
              ngref(1) = ngref(1) + 1
              phiref(1) = phiref(1) + phi(1,i,j,0 ,1)
            end do
            end do
          end if
!
          if(iphiref(2,3).eq.1.and.bared(2,3,myid+1).eq.1) then
            do j=ylower,yupper
            do i=xlower,xupper
              ngref(1) = ngref(1) + 1
              phiref(1) = phiref(1) + phi(1,i,j,zu-zl,1)
            end do
            end do
          end if
!
          call MPI_Reduce(ngref(1),ngref(2),1,MPI_INTEGER,MPI_SUM,0,MCW,ierr)
          call MPI_Reduce(phiref(1),phiref(2),1,MPI_REAL8,MPI_SUM,0,MCW,ierr)
!
          if(myid.eq.0.and.ngref(2).ne.0) then
            phiref(2) = phiref(2)/ngref(2)
          else
            phiref(2) = 0.0d0
          end if
        end if
!
        do ipc=1,npprb
          pprb(ipc,1) = 0.0d0
          do k=dpprb(1,3,ipc),dpprb(2,3,ipc)
          do j=dpprb(1,2,ipc),dpprb(2,2,ipc)
          do i=dpprb(1,1,ipc),dpprb(2,1,ipc)
            if(i.ge.xl.and.i.lt.xu.and. &
               j.ge.yl.and.j.lt.yu.and. &
               k.ge.zl.and.k.lt.zu) then
              pprb(ipc,1) = pprb(ipc,1) + phi(1,i-xl,j-yl,k-zl,1)
            end if
          end do
          end do
          end do
        end do
!
        call MPI_Reduce(pprb(1:npprb,1),pprb(1:npprb,2),npprb,MPI_REAL8,MPI_SUM,0,MCW,ierr)
        if(myid.eq.0) then
          where(ngpprb(1:npprb).gt.0)
            pprb(1:npprb,2) = pprb(1:npprb,2)/ngpprb(1:npprb)
          elsewhere
            pprb(1:npprb,2) = 0.0d0
          end where
        end if
!
!       ------------- 
        if(ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) then
          if(daverg.eq.0) then
            ebav(EX:BZ,:,:,:,1:nebfld+1) = eb(EX:BZ,:,:,:,1:(nebfld+1)*2:2)
            rhoav(1:nspec,:,:,:) = rhodg(1:nspec,:,:,:,1)
            rhoav(nspec+1:nspec+nspec,:,:,:) = &
           &  rhodg(nspec+1:nspec+nspec,:,:,:,1)
            rhoav(nspec*2+1,:,:,:) = rho(1,:,:,:,1)
            rhoav(nspec*2+2,:,:,:) = rhobk(1,:,:,:,3)           
!            rhoav(nspec*2+2,:,:,:) = rhobk(2,:,:,:,3)
            phiav(1,:,:,:) = phi(1,:,:,:,1)
          else if(daverg.eq.1) then
            h5fcount = h5fcount + 1
            ebav(EX:BZ,:,:,:,1:nebfld+1) = ebav(EX:BZ,:,:,:,1:nebfld+1) &
           &                             + eb(EX:BZ,:,:,:,1:(nebfld+1)*2:2)
            rhoav(1:nspec,:,:,:) = rhoav(1:nspec,:,:,:) &
           &                     + rhodg(1:nspec,:,:,:,1)
            rhoav(nspec+1:nspec+nspec,:,:,:) = &
           &  rhoav(nspec+1:nspec+nspec,:,:,:) &
           &  + rhodg(nspec+1:nspec+nspec,:,:,:,1)
            rhoav(nspec*2+1,:,:,:) = rhoav(nspec*2+1,:,:,:) &
           &                       + rho(1,:,:,:,1)
            rhoav(nspec*2+2,:,:,:) = rhoav(nspec*2+2,:,:,:) &
           &                       + rhobk(1,:,:,:,3)
!           &                       + rhobk(2,:,:,:,3)
            phiav(1,:,:,:) = phiav(1,:,:,:) &
           &               + phi(1,:,:,:,1)
!
            ebav(:,:,:,:,:) = ebav(:,:,:,:,:)/h5fcount
            rhoav(:,:,:,:) = rhoav(:,:,:,:)/h5fcount
            phiav(:,:,:,:) = phiav(:,:,:,:)/h5fcount
          else if(daverg.ge.2) then
            ebav(EX:BZ,:,:,:,1:nebfld+1) = &
           &  ebav(EX:BZ,:,:,:,1:nebfld+1)*expema &
           &  + eb(EX:BZ,:,:,:,1:(nebfld+1)*2:2)*(1.0d0-expema)
            rhoav(1:nspec,:,:,:) = rhoav(1:nspec,:,:,:)*expema &
           &                     + rhodg(1:nspec,:,:,:,1)*(1.0d0-expema)
            rhoav(nspec+1:nspec+nspec,:,:,:) = &
           &  rhoav(nspec+1:nspec+nspec,:,:,:)*expema &
           &  + rhodg(nspec+1:nspec+nspec,:,:,:,1)*(1.0d0-expema)
            rhoav(nspec*2+1,:,:,:) = rhoav(nspec*2+1,:,:,:)*expema &
           &                       + rho(1,:,:,:,1)*(1.0d0-expema)
            rhoav(nspec*2+2,:,:,:) = rhoav(nspec*2+2,:,:,:)*expema &
           &                       + rhobk(1,:,:,:,3)*(1.0d0-expema)
!           &                       + rhobk(2,:,:,:,3)*(1.0d0-expema)
            phiav(1,:,:,:) = phiav(1,:,:,:)*expema &
           &               + phi(1,:,:,:,1)*(1.0d0-expema)
          end if
        else
          if(daverg.eq.1) then
            h5fcount = h5fcount + 1
            ebav(EX:BZ,:,:,:,1:nebfld+1) = ebav(EX:BZ,:,:,:,1:nebfld+1) &
           &                             + eb(EX:BZ,:,:,:,1:(nebfld+1)*2:2)
            rhoav(1:nspec,:,:,:) = rhoav(1:nspec,:,:,:) &
           &                     + rhodg(1:nspec,:,:,:,1)
            rhoav(nspec+1:nspec+nspec,:,:,:) = &
           &  rhoav(nspec+1:nspec+nspec,:,:,:) &
           &  + rhodg(nspec+1:nspec+nspec,:,:,:,1)
            rhoav(nspec*2+1,:,:,:) = rhoav(nspec*2+1,:,:,:) &
           &                       + rho(1,:,:,:,1)
            rhoav(nspec*2+2,:,:,:) = rhoav(nspec*2+2,:,:,:) &
           &                       + rhobk(1,:,:,:,3)
!           &                       + rhobk(2,:,:,:,3)
            phiav(1,:,:,:) = phiav(1,:,:,:) &
           &               + phi(1,:,:,:,1)
          else if(daverg.ge.2) then
            ebav(EX:BZ,:,:,:,1:nebfld+1) = &
           &  ebav(EX:BZ,:,:,:,1:nebfld+1)*expema &
           &  + eb(EX:BZ,:,:,:,1:(nebfld+1)*2:2)*(1.0d0-expema)
            rhoav(1:nspec,:,:,:) = rhoav(1:nspec,:,:,:)*expema &
           &                     + rhodg(1:nspec,:,:,:,1)*(1.0d0-expema)
            rhoav(nspec+1:nspec+nspec,:,:,:) = &
           &  rhoav(nspec+1:nspec+nspec,:,:,:)*expema &
           &  + rhodg(nspec+1:nspec+nspec,:,:,:,1)*(1.0d0-expema)
            rhoav(nspec*2+1,:,:,:) = rhoav(nspec*2+1,:,:,:)*expema &
           &                       + rho(1,:,:,:,1)*(1.0d0-expema)
            rhoav(nspec*2+2,:,:,:) = rhoav(nspec*2+2,:,:,:)*expema &
           &                       + rhobk(1,:,:,:,3)*(1.0d0-expema)
!           &                       + rhobk(2,:,:,:,3)*(1.0d0-expema)
            phiav(1,:,:,:) = phiav(1,:,:,:)*expema &
           &               + phi(1,:,:,:,1)*(1.0d0-expema)
          end if
        end if
!
        if(idiag.le.1) then
          rhoav(nspec*2+3,:,:,:) = rhoav(nspec*2+3,:,:,:) &
         &  + sum(rhodg(nspec+1:nspec+nspec,:,:,:,1),dim=1)
        else
          rhoav(nspec*2+3,:,:,:) = rhoav(nspec*2+3,:,:,:)*expema &
         &  + sum(rhodg(nspec+1:nspec+nspec,:,:,:,1),dim=1)*(1.0d0-expema)
        end if
      end if
!
      if(func.eq.1.and.(ijxyz(JX).ge.1.or.ijxyz(JY).ge.1.or.ijxyz(JZ).ge.1)) then
        if(ijdiag.ne.0.and.mod(idiag,ijdiag).eq.0) then
          if(daverg.eq.0) then
            ajav(1:nspec*3,:,:,:) = ajdg(1:nspec*3,:,:,:,1)
            ajav(nspec*3+1:nspec*3+3,:,:,:) = aj(1:3,:,:,:,1)
            ajav(nspec*3+4:nspec*3+6,:,:,:) = aj(1:3,:,:,:,3)
          else if(daverg.eq.1) then
            h5jcount = h5jcount + 1
            ajav(1:nspec*3,:,:,:) = ajav(1:nspec*3,:,:,:) &
           &                      + ajdg(1:nspec*3,:,:,:,1)
            ajav(nspec*3+1:nspec*3+3,:,:,:) = &
           &  ajav(nspec*3+1:nspec*3+3,:,:,:) + aj(1:3,:,:,:,1)
            ajav(nspec*3+4:nspec*3+6,:,:,:) = &
           &  ajav(nspec*3+4:nspec*3+6,:,:,:) + aj(1:3,:,:,:,3)
!
            ajav(:,:,:,:) = ajav(:,:,:,:)/h5jcount
          else if(daverg.ge.2) then
            ajav(1:nspec*3,:,:,:) = ajav(1:nspec*3,:,:,:)*expema &
           &                      + ajdg(1:nspec*3,:,:,:,1)*(1.0d0-expema)
            ajav(nspec*3+1:nspec*3+3,:,:,:) = &
           &  ajav(nspec*3+1:nspec*3+3,:,:,:)*expema &
           &  + aj(1:3,:,:,:,1)*(1.0d0-expema)
            ajav(nspec*3+4:nspec*3+6,:,:,:) = &
           &  ajav(nspec*3+4:nspec*3+6,:,:,:)*expema &
           &  + aj(1:3,:,:,:,3)*(1.0d0-expema)
          end if
        else
          if(daverg.eq.1) then
            h5jcount = h5jcount + 1
            ajav(1:nspec*3,:,:,:) = ajav(1:nspec*3,:,:,:) &
           &                      + ajdg(1:nspec*3,:,:,:,1)
            ajav(nspec*3+1:nspec*3+3,:,:,:) = &
           &  ajav(nspec*3+1:nspec*3+3,:,:,:) + aj(1:3,:,:,:,1)
            ajav(nspec*3+4:nspec*3+6,:,:,:) = &
           &  ajav(nspec*3+4:nspec*3+6,:,:,:) + aj(1:3,:,:,:,3)
          else if(daverg.ge.2) then
            ajav(1:nspec*3,:,:,:) = ajav(1:nspec*3,:,:,:)*expema &
           &                      + ajdg(1:nspec*3,:,:,:,1)*(1.0d0-expema)
            ajav(nspec*3+1:nspec*3+3,:,:,:) = &
           &  ajav(nspec*3+1:nspec*3+3,:,:,:)*expema &
           &  + aj(1:3,:,:,:,1)*(1.0d0-expema)
            ajav(nspec*3+4:nspec*3+6,:,:,:) = &
           &  ajav(nspec*3+4:nspec*3+6,:,:,:)*expema &
           &  + aj(1:3,:,:,:,3)*(1.0d0-expema)
          end if
        end if
      end if


!-------------------- particle pushing fields
      do k=0,zu-zl
      do j=0,yu-yl
      do i=0,xu-xl
        wrk(EX:BZ,i,j,k) = (/e0x,e0y,e0z,b0x,b0y,b0z/)
      end do
      end do
      end do
!
      do ieb=1,nebfld+1
!       ------------- pe-fields
        do k=0,zu-zl
        do j=0,yu-yl
        do i=0,xu-xl
          wrk(EX,i,j,k) = wrk(EX,i,j,k) &
         &              + 0.5d0*(ebav(EX,i,j,k,ieb) + ebav(EX,i-1,j,k,ieb))
          wrk(EY,i,j,k) = wrk(EY,i,j,k) &
         &              + 0.5d0*(ebav(EY,i,j,k,ieb) + ebav(EY,i,j-1,k,ieb))
          wrk(EZ,i,j,k) = wrk(EZ,i,j,k) &
         &              + 0.5d0*(ebav(EZ,i,j,k,ieb) + ebav(EZ,i,j,k-1,ieb))
        end do
        end do
        end do
!
!       ------------- pb-fields
        do k=0,zu-zl
        do j=0,yu-yl
        do i=0,xu-xl
          wrk(BX,i,j,k) = wrk(BX,i,j,k) &
         &              + 0.25d0*(ebav(BX,i,j,k  ,ieb) + ebav(BX,i,j-1,k  ,ieb) &
         &                      + ebav(BX,i,j,k-1,ieb) + ebav(BX,i,j-1,k-1,ieb))
          wrk(BY,i,j,k) = wrk(BY,i,j,k) &
         &              + 0.25d0*(ebav(BY,i,j,k  ,ieb) + ebav(BY,i-1,j,k  ,ieb) &
         &                      + ebav(BY,i,j,k-1,ieb) + ebav(BY,i-1,j,k-1,ieb))
          wrk(BZ,i,j,k) = wrk(BZ,i,j,k) &
         &              + 0.25d0*(ebav(BZ,i,j  ,k,ieb) + ebav(BZ,i-1,j  ,k,ieb) &
         &                      + ebav(BZ,i,j-1,k,ieb) + ebav(BZ,i-1,j-1,k,ieb))
        end do
        end do
        end do
      end do


      h5ind = 0
!-------------------- fields
!     --------------- electromagnetic field:
      if(func.eq.1) write(dsname,'(i4.4)') lfdiag
      if(ifxyz(EX).ge.1) then
        h5ind = h5ind + 1
        if(func.eq.0) then
          write(grname,'(a)') 'ex'
          write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
          call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
          if(id_sd(1,h5ind).le.0) then
            write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : ex create h5 file'
          end if
        else if(func.eq.1.and.ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) then
          call wrt3mp(id_sd(2,h5ind),dsname,myid,nodes(1:3),idimfd(:,FEB),h5tab, &
         &            real(sum(ebav(EX,0:nxsd,0:nysd,0:nzsd,1:nebfld+1)/rene/coe_e,4),PRECI), &
         &            istat1,istat2,int(odoms,HSIZE_T))
        else if(func.eq.2) then
          call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
          if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
            write(*,*) h5ind, istats(:,h5ind), ' : ex close h5 file'
          end if
        end if
      end if
!
      if(ifxyz(EY).ge.1) then
        h5ind = h5ind + 1
        if(func.eq.0) then
          write(grname,'(a)') 'ey'
          write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
          call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
          if(id_sd(1,h5ind).le.0) then
            write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : ey create h5 file'
          end if
        else if(func.eq.1.and.ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) then
          call wrt3mp(id_sd(2,h5ind),dsname,myid,nodes(1:3),idimfd(:,FEB),h5tab, &
         &            real(sum(ebav(EY,0:nxsd,0:nysd,0:nzsd,1:nebfld+1)/rene/coe_e,4),PRECI), &
         &            istat1,istat2,int(odoms,HSIZE_T))
        else if(func.eq.2) then
          call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
          if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
            write(*,*) h5ind, istats(:,h5ind), ' : ey close h5 file'
          end if
        end if
      end if
!
      if(ifxyz(EZ).ge.1) then
        h5ind = h5ind + 1
        if(func.eq.0) then
          write(grname,'(a)') 'ez'
          write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
          call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
          if(id_sd(1,h5ind).le.0) then
            write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : ez create h5 file'
          end if
        else if(func.eq.1.and.ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) then
          call wrt3mp(id_sd(2,h5ind),dsname,myid,nodes(1:3),idimfd(:,FEB),h5tab, &
         &            real(sum(ebav(EZ,0:nxsd,0:nysd,0:nzsd,1:nebfld+1)/rene/coe_e,4),PRECI), &
         &            istat1,istat2,int(odoms,HSIZE_T))
        else if(func.eq.2) then
          call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
          if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
            write(*,*) h5ind, istats(:,h5ind), ' : ez close h5 file'
          end if
        end if
      end if
!
      if(ifxyz(BX).ge.1) then
        h5ind = h5ind + 1
        if(func.eq.0) then
          write(grname,'(a)') 'bx'
          write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
          call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
          if(id_sd(1,h5ind).le.0) then
            write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : bx create h5 file'
          end if
        else if(func.eq.1.and.ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) then
          call wrt3mp(id_sd(2,h5ind),dsname,myid,nodes(1:3),idimfd(:,FEB),h5tab, &
         &            real(sum(ebav(BX,0:nxsd,0:nysd,0:nzsd,1:nebfld+1)/renb/coe_m,4),PRECI), &
         &            istat1,istat2,int(odoms,HSIZE_T))
        else if(func.eq.2) then
          call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
          if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
            write(*,*) h5ind, istats(:,h5ind), ' : bx close h5 file'
          end if
        end if
      end if
!
      if(ifxyz(BY).ge.1) then
        h5ind = h5ind + 1
        if(func.eq.0) then
          write(grname,'(a)') 'by'
          write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
          call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
          if(id_sd(1,h5ind).le.0) then
            write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : by create h5 file'
          end if
        else if(func.eq.1.and.ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) then
          call wrt3mp(id_sd(2,h5ind),dsname,myid,nodes(1:3),idimfd(:,FEB),h5tab, &
         &            real(sum(ebav(BY,0:nxsd,0:nysd,0:nzsd,1:nebfld+1)/renb/coe_m,4),PRECI), &
         &            istat1,istat2,int(odoms,HSIZE_T))
        else if(func.eq.2) then
          call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
          if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
            write(*,*) h5ind, istats(:,h5ind), ' : by close h5 file'
          end if
        end if
      end if
!
      if(ifxyz(BZ).ge.1) then
        h5ind = h5ind + 1
        if(func.eq.0) then
          write(grname,'(a)') 'bz'
          write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
          call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
          if(id_sd(1,h5ind).le.0) then
            write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : bz create h5 file'
          end if
        else if(func.eq.1.and.ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) then
          call wrt3mp(id_sd(2,h5ind),dsname,myid,nodes(1:3),idimfd(:,FEB),h5tab, &
         &            real(sum(ebav(BZ,0:nxsd,0:nysd,0:nzsd,1:nebfld+1)/renb/coe_m,4),PRECI), &
         &            istat1,istat2,int(odoms,HSIZE_T))
        else if(func.eq.2) then
          call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
          if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
            write(*,*) h5ind, istats(:,h5ind), ' : bz close h5 file'
          end if
        end if
      end if
!
      if(ifxyz(EX).ge.2) then
        h5ind = h5ind + 1
        if(func.eq.0) then
          write(grname,'(a)') 'epx'
          write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
          call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
          if(id_sd(1,h5ind).le.0) then
            write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : epx create h5 file'
          end if
        else if(func.eq.1.and.ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) then
          call wrt3mp(id_sd(2,h5ind),dsname,myid,nodes(1:3),idimfd(:,FEB),h5tab, &
         &            real(wrk(EX,0:nxsd,0:nysd,0:nzsd)/rene/coe_e,PRECI), &
         &            istat1,istat2,int(odoms,HSIZE_T))
        else if(func.eq.2) then
          call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
          if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
            write(*,*) h5ind, istats(:,h5ind), ' : epx close h5 file'
          end if
        end if
      end if
!
      if(ifxyz(EY).ge.2) then
        h5ind = h5ind + 1
        if(func.eq.0) then
          write(grname,'(a)') 'epy'
          write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
          call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
          if(id_sd(1,h5ind).le.0) then
            write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : epy create h5 file'
          end if
        else if(func.eq.1.and.ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) then
          call wrt3mp(id_sd(2,h5ind),dsname,myid,nodes(1:3),idimfd(:,FEB),h5tab, &
         &            real(wrk(EY,0:nxsd,0:nysd,0:nzsd)/rene/coe_e,PRECI), &
         &            istat1,istat2,int(odoms,HSIZE_T))
        else if(func.eq.2) then
          call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
          if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
            write(*,*) h5ind, istats(:,h5ind), ' : epy close h5 file'
          end if
        end if
      end if
!
      if(ifxyz(EZ).ge.2) then
        h5ind = h5ind + 1
        if(func.eq.0) then
          write(grname,'(a)') 'epz'
          write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
          call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
          if(id_sd(1,h5ind).le.0) then
            write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : epz create h5 file'
          end if
        else if(func.eq.1.and.ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) then
          call wrt3mp(id_sd(2,h5ind),dsname,myid,nodes(1:3),idimfd(:,FEB),h5tab, &
         &            real(wrk(EZ,0:nxsd,0:nysd,0:nzsd)/rene/coe_e,PRECI), &
         &            istat1,istat2,int(odoms,HSIZE_T))
        else if(func.eq.2) then
          call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
          if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
            write(*,*) h5ind, istats(:,h5ind), ' : epz close h5 file'
          end if
        end if
      end if
!
      if(ifxyz(BX).ge.2) then
        h5ind = h5ind + 1
        if(func.eq.0) then
          write(grname,'(a)') 'bpx'
          write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
          call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
          if(id_sd(1,h5ind).le.0) then
            write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : bpx create h5 file'
          end if
        else if(func.eq.1.and.ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) then
          call wrt3mp(id_sd(2,h5ind),dsname,myid,nodes(1:3),idimfd(:,FEB),h5tab, &
         &            real(wrk(BX,0:nxsd,0:nysd,0:nzsd)/rene/coe_e,PRECI), &
         &            istat1,istat2,int(odoms,HSIZE_T))
        else if(func.eq.2) then
          call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
          if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
            write(*,*) h5ind, istats(:,h5ind), ' : bpx close h5 file'
          end if
        end if
      end if
!
      if(ifxyz(BY).ge.2) then
        h5ind = h5ind + 1
        if(func.eq.0) then
          write(grname,'(a)') 'bpy'
          write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
          call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
          if(id_sd(1,h5ind).le.0) then
            write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : bpy create h5 file'
          end if
        else if(func.eq.1.and.ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) then
          call wrt3mp(id_sd(2,h5ind),dsname,myid,nodes(1:3),idimfd(:,FEB),h5tab, &
         &            real(wrk(BY,0:nxsd,0:nysd,0:nzsd)/rene/coe_e,PRECI), &
         &            istat1,istat2,int(odoms,HSIZE_T))
        else if(func.eq.2) then
          call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
          if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
            write(*,*) h5ind, istats(:,h5ind), ' : bpy close h5 file'
          end if
        end if
      end if
!
      if(ifxyz(BZ).ge.2) then
        h5ind = h5ind + 1
        if(func.eq.0) then
          write(grname,'(a)') 'bpz'
          write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
          call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
          if(id_sd(1,h5ind).le.0) then
            write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : bpz create h5 file'
          end if
        else if(func.eq.1.and.ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) then
          call wrt3mp(id_sd(2,h5ind),dsname,myid,nodes(1:3),idimfd(:,FEB),h5tab, &
         &            real(wrk(BZ,0:nxsd,0:nysd,0:nzsd)/rene/coe_e,PRECI), &
         &            istat1,istat2,int(odoms,HSIZE_T))
        else if(func.eq.2) then
          call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
          if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
            write(*,*) h5ind, istats(:,h5ind), ' : bpz close h5 file'
          end if
        end if
      end if
!
      do ieb=1,nebfld+1
        if(ifxyz(EX).ge.3) then
          h5ind = h5ind + 1
          if(func.eq.0) then
            write(grname,'(a,i1.1,a)') 'e',ieb,'x'
            write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
            call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
            if(id_sd(1,h5ind).le.0) then
              write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : ex create h5 file'
            end if
          else if(func.eq.1.and.ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) then
            call wrt3mp(id_sd(2,h5ind),dsname,myid,nodes(1:3),idimfd(:,FEB),h5tab, &
           &            real(ebav(EX,0:nxsd,0:nysd,0:nzsd,ieb)/rene/coe_e,PRECI), &
           &            istat1,istat2,int(odoms,HSIZE_T))
          else if(func.eq.2) then
            call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
            if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
              write(*,*) h5ind, istats(:,h5ind), ' : ex close h5 file'
            end if
          end if
        end if
!
        if(ifxyz(EY).ge.3) then
          h5ind = h5ind + 1
          if(func.eq.0) then
            write(grname,'(a,i1.1,a)') 'e',ieb,'y'
            write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
            call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
            if(id_sd(1,h5ind).le.0) then
              write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : ey create h5 file'
            end if
          else if(func.eq.1.and.ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) then
            call wrt3mp(id_sd(2,h5ind),dsname,myid,nodes(1:3),idimfd(:,FEB),h5tab, &
           &            real(ebav(EY,0:nxsd,0:nysd,0:nzsd,ieb)/rene/coe_e,PRECI), &
           &            istat1,istat2,int(odoms,HSIZE_T))
          else if(func.eq.2) then
            call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
            if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
              write(*,*) h5ind, istats(:,h5ind), ' : ey close h5 file'
            end if
          end if
        end if
!
        if(ifxyz(EZ).ge.3) then
          h5ind = h5ind + 1
          if(func.eq.0) then
            write(grname,'(a,i1.1,a)') 'e',ieb,'z'
            write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
            call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
            if(id_sd(1,h5ind).le.0) then
              write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : ez create h5 file'
            end if
          else if(func.eq.1.and.ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) then
            call wrt3mp(id_sd(2,h5ind),dsname,myid,nodes(1:3),idimfd(:,FEB),h5tab, &
           &            real(ebav(EZ,0:nxsd,0:nysd,0:nzsd,ieb)/rene/coe_e,PRECI), &
           &            istat1,istat2,int(odoms,HSIZE_T))
          else if(func.eq.2) then
            call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
            if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
              write(*,*) h5ind, istats(:,h5ind), ' : ez close h5 file'
            end if
          end if
        end if
!
        if(ifxyz(BX).ge.3) then
          h5ind = h5ind + 1
          if(func.eq.0) then
            write(grname,'(a,i1.1,a)') 'b',ieb,'x'
            write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
            call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
            if(id_sd(1,h5ind).le.0) then
              write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : bx create h5 file'
            end if
          else if(func.eq.1.and.ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) then
            call wrt3mp(id_sd(2,h5ind),dsname,myid,nodes(1:3),idimfd(:,FEB),h5tab, &
           &            real(ebav(BX,0:nxsd,0:nysd,0:nzsd,ieb)/renb/coe_m,PRECI), &
           &            istat1,istat2,int(odoms,HSIZE_T))
          else if(func.eq.2) then
            call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
            if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
              write(*,*) h5ind, istats(:,h5ind), ' : bx close h5 file'
            end if
          end if
        end if
!
        if(ifxyz(BY).ge.3) then
          h5ind = h5ind + 1
          if(func.eq.0) then
            write(grname,'(a,i1.1,a)') 'b',ieb,'y'
            write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
            call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
            if(id_sd(1,h5ind).le.0) then
              write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : by create h5 file'
            end if
          else if(func.eq.1.and.ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) then
            call wrt3mp(id_sd(2,h5ind),dsname,myid,nodes(1:3),idimfd(:,FEB),h5tab, &
           &            real(ebav(BY,0:nxsd,0:nysd,0:nzsd,ieb)/renb/coe_m,PRECI), &
           &            istat1,istat2,int(odoms,HSIZE_T))
          else if(func.eq.2) then
            call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
            if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
              write(*,*) h5ind, istats(:,h5ind), ' : by close h5 file'
            end if
          end if
        end if
!
        if(ifxyz(BZ).ge.3) then
          h5ind = h5ind + 1
          if(func.eq.0) then
            write(grname,'(a,i1.1,a)') 'b',ieb,'z'
            write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
            call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
            if(id_sd(1,h5ind).le.0) then
              write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : bz create h5 file'
            end if
          else if(func.eq.1.and.ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) then
            call wrt3mp(id_sd(2,h5ind),dsname,myid,nodes(1:3),idimfd(:,FEB),h5tab, &
           &            real(ebav(BZ,0:nxsd,0:nysd,0:nzsd,ieb)/renb/coe_m,PRECI), &
           &            istat1,istat2,int(odoms,HSIZE_T))
          else if(func.eq.2) then
            call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
            if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
              write(*,*) h5ind, istats(:,h5ind), ' : bz close h5 file'
            end if
          end if
        end if
      end do

!     --------------- charge density + potential:
      if(func.eq.1) write(dsname,'(i4.4)') lfdiag
      if(ifxyz(7).ge.1) then
        h5ind = h5ind + 1
        if(func.eq.0) then
          write(grname,'(a)') 'rho'
          write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
          call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
          if(id_sd(1,h5ind).le.0) then
            write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : rho create h5 file'
          end if
        else if(func.eq.1.and.ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) then
          call wrt3mp(id_sd(2,h5ind),dsname,myid,nodes(1:3),idimfd(:,FRH),h5tab, &
         &            real(rhoav(nspec*2+1,0:nxsd,0:nysd,0:nzsd)/renrho,PRECI), &
         &            istat1,istat2,int(odoms,HSIZE_T))
        else if(func.eq.2) then
          call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
          if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
            write(*,*) h5ind, istats(:,h5ind), ' : rho close h5 file'
          end if
        end if
      end if
!
      do is=1,nspec
        if(irhsp(is).ge.1) then
          h5ind = h5ind + 1
          if(func.eq.0) then
            write(grname,'(a,i1.1,a)') 'nd',is,'p'
            write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
            call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
            if(id_sd(1,h5ind).le.0) then
              write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : ndp', is, 'create h5 file'
            end if
          else if(func.eq.1.and.ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) then
            call wrt3mp(id_sd(2,h5ind),dsname,myid,nodes(1:3),idimfd(:,FRD),h5tab, &
           &            real(abs(rhoav(nspec*0+is,0:nxsd,0:nysd,0:nzsd)/renrho/rho0),PRECI), &
           &            istat1,istat2,int(odoms,HSIZE_T))
          else if(func.eq.2) then
            call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
            if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
              write(*,*) h5ind, istats(:,h5ind), ' : ndp', is, 'close h5 file'
            end if
          end if
        end if
      end do
!
      if(ifxyz(7).ge.1) then
        h5ind = h5ind + 1
        if(func.eq.0) then
          write(grname,'(a)') 'rhobk'
          write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
          call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
          if(id_sd(1,h5ind).le.0) then
            write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : rhobk create h5 file'
          end if
        else if(func.eq.1.and.ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) then
          call wrt3mp(id_sd(2,h5ind),dsname,myid,nodes(1:3),idimfd(:,FRH),h5tab, &
         &            real(rhoav(nspec*2+2,0:nxsd,0:nysd,0:nzsd)/renrho,PRECI), &
         &            istat1,istat2,int(odoms,HSIZE_T))
        else if(func.eq.2) then
          call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
          if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
            write(*,*) h5ind, istats(:,h5ind), ' : rhobk close h5 file'
          end if
        end if
      end if
!
      if(ifxyz(7).ge.1) then
        h5ind = h5ind + 1
        if(func.eq.0) then
          write(grname,'(a)') 'phisp'
          write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
          call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
          if(id_sd(1,h5ind).le.0) then
            write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : phisp create h5 file'
          end if
        else if(func.eq.1.and.ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) then
          call wrt3mp(id_sd(2,h5ind),dsname,myid,nodes(1:3),idimfd(:,FPH),h5tab, &
         &            real(phiav(1,0:nxsd,0:nysd,0:nzsd)/renphi/coe_ph,PRECI), &
         &            istat1,istat2,int(odoms,HSIZE_T))
        else if(func.eq.2) then
          call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
          if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
            write(*,*) h5ind, istats(:,h5ind), ' : phisp close h5 file'
          end if
        end if
      end if


!     --------------- current density:
      if(func.eq.1) write(dsname,'(i4.4)') ljdiag
      if(ijxyz(JX).ge.1) then
        h5ind = h5ind + 1
        if(func.eq.0) then
          write(grname,'(a)') 'jx'
          write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
          call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
          if(id_sd(1,h5ind).le.0) then
            write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : jx create h5 file'
          end if
        else if(func.eq.1.and.ijdiag.ne.0.and.mod(idiag,ijdiag).eq.0) then
          call wrt3mp(id_sd(2,h5ind),dsname,myid,nodes(1:3),idimfd(:,FAJ),h5tab, &
         &            real(ajav(nspec*3+JX,0:nxsd,0:nysd,0:nzsd)/renj/coe_j,PRECI), &
         &            istat1,istat2,int(odoms,HSIZE_T))
        else if(func.eq.2) then
          call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
          if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
            write(*,*) h5ind, istats(:,h5ind), ' : jx close h5 file'
          end if
        end if
!
        do is=1,nspec
          h5ind = h5ind + 1
          if(func.eq.0) then
            write(grname,'(a,i1.1,a)') 'j',is,'x'
            write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
            call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
            if(id_sd(1,h5ind).le.0) then
              write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : jx', is, 'create h5 file'
            end if
          else if(func.eq.1.and.ijdiag.ne.0.and.mod(idiag,ijdiag).eq.0) then
            call wrt3mp(id_sd(2,h5ind),dsname,myid,nodes(1:3),idimfd(:,FAJ),h5tab, &
           &            real(ajav((is-1)*3+JX,0:nxsd,0:nysd,0:nzsd)/renj/coe_j,PRECI), &
           &            istat1,istat2,int(odoms,HSIZE_T))
          else if(func.eq.2) then
            call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
            if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
              write(*,*) h5ind, istats(:,h5ind), ' : jx', is, 'close h5 file'
            end if
          end if
        end do
      end if
!
      if(ijxyz(JY).ge.1) then
        h5ind = h5ind + 1
        if(func.eq.0) then
          write(grname,'(a)') 'jy'
          write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
          call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
          if(id_sd(1,h5ind).le.0) then
            write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : jy create h5 file'
          end if
        else if(func.eq.1.and.ijdiag.ne.0.and.mod(idiag,ijdiag).eq.0) then
          call wrt3mp(id_sd(2,h5ind),dsname,myid,nodes(1:3),idimfd(:,FAJ),h5tab, &
         &            real(ajav(nspec*3+JY,0:nxsd,0:nysd,0:nzsd)/renj/coe_j,PRECI), &
         &            istat1,istat2,int(odoms,HSIZE_T))
        else if(func.eq.2) then
          call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
          if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
            write(*,*) h5ind, istats(:,h5ind), ' : jy close h5 file'
          end if
        end if
!
        do is=1,nspec
          h5ind = h5ind + 1
          if(func.eq.0) then
            write(grname,'(a,i1.1,a)') 'j',is,'y'
            write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
            call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
            if(id_sd(1,h5ind).le.0) then
              write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : jy', is, 'create h5 file'
            end if
          else if(func.eq.1.and.ijdiag.ne.0.and.mod(idiag,ijdiag).eq.0) then
            call wrt3mp(id_sd(2,h5ind),dsname,myid,nodes(1:3),idimfd(:,FAJ),h5tab, &
           &            real(ajav((is-1)*3+JY,0:nxsd,0:nysd,0:nzsd)/renj/coe_j,PRECI), &
           &            istat1,istat2,int(odoms,HSIZE_T))
          else if(func.eq.2) then
            call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
            if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
              write(*,*) h5ind, istats(:,h5ind), ' : jy', is, 'close h5 file'
            end if
          end if
        end do
      end if
!
      if(ijxyz(JZ).ge.1) then
        h5ind = h5ind + 1
        if(func.eq.0) then
          write(grname,'(a)') 'jz'
          write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
          call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
          if(id_sd(1,h5ind).le.0) then
            write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : jz create h5 file'
          end if
        else if(func.eq.1.and.ijdiag.ne.0.and.mod(idiag,ijdiag).eq.0) then
          call wrt3mp(id_sd(2,h5ind),dsname,myid,nodes(1:3),idimfd(:,FAJ),h5tab, &
         &            real(ajav(nspec*3+JZ,0:nxsd,0:nysd,0:nzsd)/renj/coe_j,PRECI), &
         &            istat1,istat2,int(odoms,HSIZE_T))
        else if(func.eq.2) then
          call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
          if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
            write(*,*) h5ind, istats(:,h5ind), ' : jz close h5 file'
          end if
        end if
!
        do is=1,nspec
          h5ind = h5ind + 1
          if(func.eq.0) then
            write(grname,'(a,i1.1,a)') 'j',is,'z'
            write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
            call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
            if(id_sd(1,h5ind).le.0) then
              write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : jz', is, 'create h5 file'
            end if
          else if(func.eq.1.and.ijdiag.ne.0.and.mod(idiag,ijdiag).eq.0) then
            call wrt3mp(id_sd(2,h5ind),dsname,myid,nodes(1:3),idimfd(:,FAJ),h5tab, &
           &            real(ajav((is-1)*3+JZ,0:nxsd,0:nysd,0:nzsd)/renj/coe_j,PRECI), &
           &            istat1,istat2,int(odoms,HSIZE_T))
          else if(func.eq.2) then
            call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
            if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
              write(*,*) h5ind, istats(:,h5ind), ' : jz', is, 'close h5 file'
            end if
          end if
        end do
      end if
!
      if(ijxyz(JS).ge.1) then
        h5ind = h5ind + 1
        if(func.eq.0) then
          write(grname,'(a)') 'js'
          write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
          call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
          if(id_sd(1,h5ind).le.0) then
            write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : js create h5 file'
          end if
        else if(func.eq.1.and.ijdiag.ne.0.and.mod(idiag,ijdiag).eq.0) then
          call wrt3mp(id_sd(2,h5ind),dsname,myid,nodes(1:3),idimfd(:,FRD),h5tab, &
         &            real(sum(rhoav(nspec*1+1:nspec*1+minsp,0:nxsd,0:nysd,0:nzsd),dim=1) &
         &                     *0.5d0/renj/coe_j,PRECI), &
         &            istat1,istat2,int(odoms,HSIZE_T))
        else if(func.eq.2) then
          call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
          if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
            write(*,*) h5ind, istats(:,h5ind), ' : js close h5 file'
          end if
        end if
!
        do is=1,nspec
          if(irhsp(is).ge.1) then
          h5ind = h5ind + 1
          if(func.eq.0) then
            write(grname,'(a,i1.1,a)') 'j',is,'s'
            write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
            call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
            if(id_sd(1,h5ind).le.0) then
              write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : js', is, 'create h5 file'
            end if
          else if(func.eq.1.and.ijdiag.ne.0.and.mod(idiag,ijdiag).eq.0) then
            call wrt3mp(id_sd(2,h5ind),dsname,myid,nodes(1:3),idimfd(:,FRD),h5tab, &
           &            real(abs(rhoav(nspec*1+is,0:nxsd,0:nysd,0:nzsd)*0.5d0/renj/coe_j), &
           &                 PRECI), &
           &            istat1,istat2,int(odoms,HSIZE_T))
          else if(func.eq.2) then
            call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
            if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
              write(*,*) h5ind, istats(:,h5ind), ' : js', is, 'close h5 file'
            end if
          end if
          end if
        end do
      end if
!
      if(ijxyz(JS).ge.2) then
        h5ind = h5ind + 1
        if(func.eq.0) then
          write(grname,'(a)') 'jcx'
          write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
          call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
          if(id_sd(1,h5ind).le.0) then
            write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : jcx create h5 file'
          end if
        else if(func.eq.1.and.ijdiag.ne.0.and.mod(idiag,ijdiag).eq.0) then
          call wrt3mp(id_sd(2,h5ind),dsname,myid,nodes(1:3),idimfd(:,FAJ),h5tab, &
         &            real(ajav((nspec+1)*3+JX,0:nxsd,0:nysd,0:nzsd)/renj/coe_j,PRECI), &
         &            istat1,istat2,int(odoms,HSIZE_T))
        else if(func.eq.2) then
          call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
          if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
            write(*,*) h5ind, istats(:,h5ind), ' : jcx close h5 file'
          end if
        end if
!
        h5ind = h5ind + 1
        if(func.eq.0) then
          write(grname,'(a)') 'jcy'
          write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
          call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
          if(id_sd(1,h5ind).le.0) then
            write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : jcy create h5 file'
          end if
        else if(func.eq.1.and.ijdiag.ne.0.and.mod(idiag,ijdiag).eq.0) then
          call wrt3mp(id_sd(2,h5ind),dsname,myid,nodes(1:3),idimfd(:,FAJ),h5tab, &
         &            real(ajav((nspec+1)*3+JY,0:nxsd,0:nysd,0:nzsd)/renj/coe_j,PRECI), &
         &            istat1,istat2,int(odoms,HSIZE_T))
        else if(func.eq.2) then
          call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
          if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
            write(*,*) h5ind, istats(:,h5ind), ' : jcy close h5 file'
          end if
        end if
!
        h5ind = h5ind + 1
        if(func.eq.0) then
          write(grname,'(a)') 'jcz'
          write(plhdfn,'(a,i2.2,a,i4.4,a)') trim(grname),istep/i1hdf,'_',nfsnap,'.h5'
          call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
          if(id_sd(1,h5ind).le.0) then
            write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : jcz create h5 file'
          end if
        else if(func.eq.1.and.ijdiag.ne.0.and.mod(idiag,ijdiag).eq.0) then
          call wrt3mp(id_sd(2,h5ind),dsname,myid,nodes(1:3),idimfd(:,FAJ),h5tab, &
         &            real(ajav((nspec+1)*3+JZ,0:nxsd,0:nysd,0:nzsd)/renj/coe_j,PRECI), &
         &            istat1,istat2,int(odoms,HSIZE_T))
        else if(func.eq.2) then
          call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(2,h5ind))
          if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
            write(*,*) h5ind, istats(:,h5ind), ' : jcz close h5 file'
          end if
        end if
      end if


!-------------------- 
!     --------------- clear rho/phiave
      if(func.eq.1.and.ifdiag.ne.0.and.mod(idiag,ifdiag).eq.0) then
        lfdiag = lfdiag + 1
        if(daverg.eq.1) then
          if(myid.eq.0) print*,"h5fcount =",h5fcount
          h5fcount = 0
          ebav(:,:,:,:,:) = 0.0d0
          rhoav(:,:,:,:) = 0.0d0
          phiav(:,:,:,:) = 0.0d0
        end if
      end if
!
!     --------------- clear ajav
      if(func.eq.1.and.ijdiag.ne.0.and.mod(idiag,ijdiag).eq.0) then
        ljdiag = ljdiag + 1
        if(daverg.eq.1) then
          if(myid.eq.0) print*,"h5jcount =",h5jcount
          h5jcount = 0
          ajav(:,:,:,:) = 0.0d0
        end if
      end if


!-------------------- other components
!     --------------- antenna surface current/e-field:
      if(iadiag.ge.1) then
!        write(plhdfn,'(a,i2.2,a)') 'surfj',istep/i2hdf,'.h5'
!        call hdfopen_pg(plhdfn,grname,id_sd(24),DFACC_CREATE,MCW)
!        if(id_sd(24).le.0) then
!          write(*,*) id_sd(24), plhdfn, ' : surfj create h5 file'
!        end if
!
!        write(plhdfn,'(a,i2.2,a)') 'surfe',istep/i2hdf,'.h5'
!        call hdfopen_pg(plhdfn,grname,id_sd(25),DFACC_CREATE,MCW)
!        if(id_sd(25).le.0) then
!          write(*,*) id_sd(25), plhdfn, ' : surfe create h5 file'
!        end if
!
!        write(plhdfn,'(a,i2.2,a)') 'phi1d',istep/i2hdf,'.h5'
!        call hdfopen_pg(plhdfn,grname,id_sd(26),DFACC_CREATE,MCW)
!        if(id_sd(26).le.0) then
!          write(*,*) id_sd(26), plhdfn, ' : phi1d create h5 file'
!        end if
!
!        write(plhdfn,'(a,i2.2,a)') 'surfeU',istep/i2hdf,'.h5'
!        call hdfopen_pg(plhdfn,grname,id_sd(27),DFACC_CREATE,MCW)
!        if(id_sd(27).le.0) then
!          write(*,*) id_sd(27), plhdfn, ' : surfeU create h5 file'
!        end if
!
!        write(plhdfn,'(a,i2.2,a)') 'surfeT',istep/i2hdf,'.h5'
!        call hdfopen_pg(plhdfn,grname,id_sd(28),DFACC_CREATE,MCW)
!        if(id_sd(28).le.0) then
!          write(*,*) id_sd(28), plhdfn, ' : surfeT create h5 file'
!        end if
      end if


!-------------------- particles
!     --------------- particle positions + velocities:
      if(func.eq.1.and. &
!     &   any(ipahdf(1:nspec).ne.0.and.mod(idiag,ipahdf(1:nspec)).eq.0)) then
     &   product(mod(idiag,ipahdf(1:nspec)), &
     &           MASK=(ipahdf(1:nspec).ne.0)).eq.0) then
        lnoutp(:) = 0
        noutp(:) = 0
        poffset(:) = 0
        neeP = pbase(1)
        neeS = pbase(2)
        do is=1,nspec
          nsP = neeP + 1
          neP = neeP + totalp(is,1)
          neeP = neeP + totalp(is,1)
          nsS = neeS + 1
          neS = neeS + totalp(is,2)
          neeS = neeS + totalp(is,2)
!
          lnoutp(is) = int((neP - nsP + ipadig(is))/ipadig(is),8) &
         &           + int((neS - nsS + ipadig(is))/ipadig(is),8)
        end do
        call MPI_Allreduce(lnoutp(1:nspec),noutp(1:nspec),nspec,MPI_INTEGER8,MPI_SUM,MCW,ierr)
        call MPI_Exscan(lnoutp(1:nspec),poffset(1:nspec),nspec,MPI_INTEGER8,MPI_SUM,MCW,ierr)
        nne = 0
        do is=1,nspec
          nns = nne + 1
          nne = nne + npr(is)
!
          if(noutp(is).le.0.and.nflag_emit(is).ge.1.and.idiag.eq.0) then
            lnoutp(is) = int((nne - nns + ipadig(is))/ipadig(is),8)
          else if(noutp(is).le.0) then
            if(myid.eq.0) then
              lnoutp(is) = 1
            else
              lnoutp(is) = 0
            end if
          end if
        end do
        call MPI_Allreduce(lnoutp(1:nspec),noutp(1:nspec),nspec,MPI_INTEGER8,MPI_SUM,MCW,ierr)
        call MPI_Exscan(lnoutp(1:nspec),poffset(1:nspec),nspec,MPI_INTEGER8,MPI_SUM,MCW,ierr)
!        do iid=0,nnode
!          if(iid.eq.myid) then
!            print*, "myid:", myid
!            print*, "  lnoutp", lnoutp(1:nspec)
!            print*, "  noutp", noutp(1:nspec)
!            print*, "  poffset", poffset(1:nspec)
!            flush(6)
!          end if
!          call MPI_Barrier(MCW,ierr)
!        end do
      end if


      neeP = pbase(1)
      neeS = pbase(2)
      nne = 0
      do is=1,nspec
        if(func.eq.1.and.ipahdf(is).ne.0.and.mod(idiag,ipahdf(is)).eq.0) then
          nsP = neeP + 1
          neP = neeP + totalp(is,1)
          neeP = neeP + totalp(is,1)
          nsS = neeS + 1
          neS = neeS + totalp(is,2)
          neeS = neeS + totalp(is,2)
          nns = nne + 1
          nne = nne + npr(is)
          write(dsname,'(i4.4)') lpdiag(is)
!          if(myid.eq.0) print*, "ipadig:noutp", noutp(1:nspec)
        end if
!
        if(ipaxyz(1,is).ge.1) then
          h5ind = h5ind + 1
          if(func.eq.0) then
            write(grname,'(a)') 'xe'
            write(plhdfn,'(a,i1,a,i2.2,a,i4.4,a)') 'p',is,trim(grname),istep/i3hdf,'_',nfsnap,'.h5'
            call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
            if(id_sd(1,h5ind).le.0) then
              write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : px', is, 'create h5 file'
            end if
          else if(func.eq.1.and.ipahdf(is).ne.0.and.mod(idiag,ipahdf(is)).eq.0) then
            if(noutp(is).gt.0) then
              call wrt1vmp &
             &     (id_sd(2,h5ind),dsname,myid, &
             &      lnoutp(is:is),noutp(is:is),poffset(is:is), &
             &      real((/pbuf(nsP:neP:ipadig(is))%x/renr, &
             &             pbuf(nsS:neS:ipadig(is))%x/renr/),PRECI), &
             &      istat1,istat2)
            else if(nflag_emit(is).ge.1.and.idiag.eq.0) then
              call wrt1vmp &
             &     (id_sd(2,h5ind),dsname,myid, &
             &      lnoutp(is:is),noutp(is:is),poffset(is:is), &
             &      real(vxr(nns:nne:ipadig(is))/renv,PRECI), &
             &      istat1,istat2)
            else
              if(myid.eq.0) pout(1) = -10.0d0
              call wrt1vmp &
             &     (id_sd(2,h5ind),dsname,myid, &
             &      lnoutp(is:is),noutp(is:is),poffset(is:is), &
             &      real(pout(1:lnoutp(is)),PRECI), &
             &      istat1,istat2)
            end if
          else if(func.eq.2) then
            call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(1,h5ind))
            if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
              write(*,*) h5ind, istats(:,h5ind), ' : px', is, 'close h5 file'
            end if
          end if
        end if
!
        if(ipaxyz(2,is).ge.1) then
          h5ind = h5ind + 1
          if(func.eq.0) then
            write(grname,'(a)') 'ye'
            write(plhdfn,'(a,i1,a,i2.2,a,i4.4,a)') 'p',is,trim(grname),istep/i3hdf,'_',nfsnap,'.h5'
            call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
            if(id_sd(1,h5ind).le.0) then
              write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : py', is, 'create h5 file'
            end if
          else if(func.eq.1.and.ipahdf(is).ne.0.and.mod(idiag,ipahdf(is)).eq.0) then
            if(noutp(is).gt.0) then
              call wrt1vmp &
             &     (id_sd(2,h5ind),dsname,myid, &
             &      lnoutp(is:is),noutp(is:is),poffset(is:is), &
             &      real((/pbuf(nsP:neP:ipadig(is))%y/renr, &
             &             pbuf(nsS:neS:ipadig(is))%y/renr/),PRECI), &
             &      istat1,istat2)
            else if(nflag_emit(is).ge.1.and.idiag.eq.0) then
              call wrt1vmp &
             &     (id_sd(2,h5ind),dsname,myid, &
             &      lnoutp(is:is),noutp(is:is),poffset(is:is), &
             &      real(vyr(nns:nne:ipadig(is))/renv,PRECI), &
             &      istat1,istat2)
            else
              if(myid.eq.0) pout(1) = -10.0d0
              call wrt1vmp &
             &     (id_sd(2,h5ind),dsname,myid, &
             &      lnoutp(is:is),noutp(is:is),poffset(is:is), &
             &      real(pout(1:lnoutp(is)),PRECI), &
             &      istat1,istat2)
            end if
          else if(func.eq.2) then
            call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(1,h5ind))
            if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
              write(*,*) h5ind, istats(:,h5ind), ' : py', is, 'close h5 file'
            end if
          end if
        end if
!
        if(ipaxyz(3,is).ge.1) then
          h5ind = h5ind + 1
          if(func.eq.0) then
            write(grname,'(a)') 'ze'
            write(plhdfn,'(a,i1,a,i2.2,a,i4.4,a)') 'p',is,trim(grname),istep/i3hdf,'_',nfsnap,'.h5'
            call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
            if(id_sd(1,h5ind).le.0) then
              write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : pz', is, 'create h5 file'
            end if
          else if(func.eq.1.and.ipahdf(is).ne.0.and.mod(idiag,ipahdf(is)).eq.0) then
            if(noutp(is).gt.0) then
              call wrt1vmp &
             &     (id_sd(2,h5ind),dsname,myid, &
             &      lnoutp(is:is),noutp(is:is),poffset(is:is), &
             &      real((/pbuf(nsP:neP:ipadig(is))%z/renr, &
             &             pbuf(nsS:neS:ipadig(is))%z/renr/),PRECI), &
             &      istat1,istat2)
            else if(nflag_emit(is).ge.1.and.idiag.eq.0) then
              call wrt1vmp &
             &     (id_sd(2,h5ind),dsname,myid, &
             &      lnoutp(is:is),noutp(is:is),poffset(is:is), &
             &      real(vzr(nns:nne:ipadig(is))/renv,PRECI), &
             &      istat1,istat2)
            else
              if(myid.eq.0) pout(1) = -10.0d0
              call wrt1vmp &
             &     (id_sd(2,h5ind),dsname,myid, &
             &      lnoutp(is:is),noutp(is:is),poffset(is:is), &
             &      real(pout(1:lnoutp(is)),PRECI), &
             &      istat1,istat2)
            end if
          else if(func.eq.2) then
            call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(1,h5ind))
            if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
              write(*,*) h5ind, istats(:,h5ind), ' : pz', is, 'close h5 file'
            end if
          end if
        end if
!
        if(ipaxyz(4,is).ge.1) then
          h5ind = h5ind + 1
          if(func.eq.0) then
            write(grname,'(a)') 'vxe'
            write(plhdfn,'(a,i1,a,i2.2,a,i4.4,a)') 'p',is,trim(grname),istep/i3hdf,'_',nfsnap,'.h5'
            call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
            if(id_sd(1,h5ind).le.0) then
              write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : pvx', is, 'create h5 file'
            end if
          else if(func.eq.1.and.ipahdf(is).ne.0.and.mod(idiag,ipahdf(is)).eq.0) then
            if(noutp(is).gt.0) then
              call wrt1vmp &
             &     (id_sd(2,h5ind),dsname,myid, &
             &      lnoutp(is:is),noutp(is:is),poffset(is:is), &
             &      real((/pbuf(nsP:neP:ipadig(is))%vx/renv, &
             &             pbuf(nsS:neS:ipadig(is))%vx/renv/),PRECI), &
             &      istat1,istat2)
            else if(nflag_emit(is).ge.1.and.idiag.eq.0) then
              call wrt1vmp &
             &     (id_sd(2,h5ind),dsname,myid, &
             &      lnoutp(is:is),noutp(is:is),poffset(is:is), &
             &      real(vxf(nns:nne:ipadig(is))/renv,PRECI), &
             &      istat1,istat2)
            else
              if(myid.eq.0) pout(1) = -10.0d0
              call wrt1vmp &
             &     (id_sd(2,h5ind),dsname,myid, &
             &      lnoutp(is:is),noutp(is:is),poffset(is:is), &
             &      real(pout(1:lnoutp(is)),PRECI), &
             &      istat1,istat2)
            end if
          else if(func.eq.2) then
            call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(1,h5ind))
            if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
              write(*,*) h5ind, istats(:,h5ind), ' : pvx', is, 'close h5 file'
            end if
          end if
        end if
!
        if(ipaxyz(5,is).ge.1) then
          h5ind = h5ind + 1
          if(func.eq.0) then
            write(grname,'(a)') 'vye'
            write(plhdfn,'(a,i1,a,i2.2,a,i4.4,a)') 'p',is,trim(grname),istep/i3hdf,'_',nfsnap,'.h5'
            call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
            if(id_sd(1,h5ind).le.0) then
              write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : pvy', is, 'create h5 file'
            end if
          else if(func.eq.1.and.ipahdf(is).ne.0.and.mod(idiag,ipahdf(is)).eq.0) then
            if(noutp(is).gt.0) then
              call wrt1vmp &
             &     (id_sd(2,h5ind),dsname,myid, &
             &      lnoutp(is:is),noutp(is:is),poffset(is:is), &
             &      real((/pbuf(nsP:neP:ipadig(is))%vy/renv, &
             &             pbuf(nsS:neS:ipadig(is))%vy/renv/),PRECI), &
             &      istat1,istat2)
            else if(nflag_emit(is).ge.1.and.idiag.eq.0) then
              call wrt1vmp &
             &     (id_sd(2,h5ind),dsname,myid, &
             &      lnoutp(is:is),noutp(is:is),poffset(is:is), &
             &      real(vyf(nns:nne:ipadig(is))/renv,PRECI), &
             &      istat1,istat2)
            else
              if(myid.eq.0) pout(1) = -10.0d0
              call wrt1vmp &
             &     (id_sd(2,h5ind),dsname,myid, &
             &      lnoutp(is:is),noutp(is:is),poffset(is:is), &
             &      real(pout(1:lnoutp(is)),PRECI), &
             &      istat1,istat2)
            end if
          else if(func.eq.2) then
            call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(1,h5ind))
            if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
              write(*,*) h5ind, istats(:,h5ind), ' : pvy', is, 'close h5 file'
            end if
          end if
        end if
!
        if(ipaxyz(6,is).ge.1) then
          h5ind = h5ind + 1
          if(func.eq.0) then
            write(grname,'(a)') 'vze'
            write(plhdfn,'(a,i1,a,i2.2,a,i4.4,a)') 'p',is,trim(grname),istep/i3hdf,'_',nfsnap,'.h5'
            call hdfopen_pg(plhdfn,grname,id_sd(1,h5ind),id_sd(2,h5ind),DFACC_CREATE,MCW)
            if(id_sd(1,h5ind).le.0) then
              write(*,*) h5ind, id_sd(1,h5ind), plhdfn, ' : pvz', is, 'create h5 file'
            end if
          else if(func.eq.1.and.ipahdf(is).ne.0.and.mod(idiag,ipahdf(is)).eq.0) then
            if(noutp(is).gt.0) then
              call wrt1vmp &
             &     (id_sd(2,h5ind),dsname,myid, &
             &      lnoutp(is:is),noutp(is:is),poffset(is:is), &
             &      real((/pbuf(nsP:neP:ipadig(is))%vz/renv, &
             &             pbuf(nsS:neS:ipadig(is))%vz/renv/),PRECI), &
             &      istat1,istat2)
            else if(nflag_emit(is).ge.1.and.idiag.eq.0) then
              call wrt1vmp &
             &     (id_sd(2,h5ind),dsname,myid, &
             &      lnoutp(is:is),noutp(is:is),poffset(is:is), &
             &      real(vzf(nns:nne:ipadig(is))/renv,PRECI), &
             &      istat1,istat2)
            else
              if(myid.eq.0) pout(1) = -10.0d0
              call wrt1vmp &
             &     (id_sd(2,h5ind),dsname,myid, &
             &      lnoutp(is:is),noutp(is:is),poffset(is:is), &
             &      real(pout(1:lnoutp(is)),PRECI), &
             &      istat1,istat2)
            end if
          else if(func.eq.2) then
            call hdfclose_g(id_sd(1,h5ind),id_sd(2,h5ind),istats(1,h5ind),istats(1,h5ind))
            if(istats(1,h5ind).lt.0.or.istats(2,h5ind).lt.0) then
              write(*,*) h5ind, istats(:,h5ind), ' : pvz', is, 'close h5 file'
            end if
          end if
        end if
!
        if(func.eq.1.and.ipahdf(is).ne.0.and.mod(idiag,ipahdf(is)).eq.0) then
          lpdiag(is) = lpdiag(is) + 1
        end if
      end do


  return
  end subroutine hdfdig
