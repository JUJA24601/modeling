#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine ohinit
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   O H I N I T
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
  use hdf
#define MSS MPI_STATUS_SIZE
  implicit none
!
  integer(kind=8) :: npmax, nranmax
  integer(kind=4) :: n
  integer(kind=4) :: is
!
  integer(kind=4) :: inttmp(2,3)
  integer(kind=HID_T) :: fileid
  integer(kind=4) :: stats0,stats1
  integer(kind=8) :: dims(2)
  character(len=30) :: filename,dsname


!------------------------------ 
      pcoord(:) = nodes(1:3)
      n = pcoord(1)*pcoord(2)*pcoord(3)
      minsp = max(nspec,1)
      allocate(nphgram(n,minsp,4))
      allocate(sdoms(2,OH_DIMENSION,n))
      allocate(bounds(2,OH_DIMENSION,n))
      allocate(bared(2,OH_DIMENSION+1,n))
      allocate(medges(2,OH_DIMENSION,n))
      allocate(famind(n+1),fammbr(2*n))
      allocate(slcs(n,3))
      allocate(scnts(n),rcnts(n))
      allocate(ncpmxs(n),scpmxs(n))
!
      allocate(totalp(minsp,4))
!
      npmax = 0
      do is=1,nspec
        npmax = npmax + np(is)
        if(myid.eq.0) print*, "is,np =", is,np(is)
      end do
      npmax = max(npmax,1000)
      if(myid.eq.0) print*, "npmax =", npmax
      maxlocalp = oh_max_local_particles(npmax,MAXFRAC,OH_IPBUF_SIZE)
      if(myid.eq.0) print*, "maxlocalp =", maxlocalp
      allocate(pbuf(maxlocalp))
!
      nranmax = 0
      do is=1,nspec
        if(nranmax.lt.np(is)/nnode+1) nranmax = np(is)/nnode + 1
      end do
      nranmax = min(max(nranmax*(100+MAXFRAC)/100+1,100000),maxlocalp)
      allocate(dranu(nranmax))
      if(myid.eq.0) print*, "nranmax =", nranmax


!------------------------------ 
      nbor(1,1,1) = -1
      sdoms(1,1,1) = 0;  sdoms(2,1,1) = -1
      scoord(1,1) = 0; scoord(2,1) = nx
      scoord(1,2) = 0; scoord(2,2) = ny
      scoord(1,3) = 0; scoord(2,3) = nz
      idxsd = mod(myid,nodes(1))
      idysd = mod(int(myid/nodes(1)),nodes(2))
      idzsd = int(myid/(nodes(1)*nodes(2)))
      if(myid.eq.0) print*, "scoord",scoord
      nbound = BNR+1
      if(nfbnd(1).ne.1.and.nx.gt.1) then
        bcond(:,1) = 1
      else
        bcond(:,1) = 2
      end if
      if(nfbnd(2).ne.1.and.ny.gt.1) then
        bcond(:,2) = 1
      else
        bcond(:,2) = 2
      end if
      if(nfbnd(3).ne.1.and.nz.gt.1) then
        bcond(:,3) = 1
      else
        bcond(:,3) = 2
      end if
!
      ftypes(:,FEB) = (/6,        -1,2, -1,2,  0,0/)             ! for eb()
      ftypes(:,FDB) = (/3,        -1,2, -1,2,  0,0/)             ! for db()
      ftypes(:,FAJ) = (/3,        -1,2,  0,0, -1,2/)             ! for aj()
      ftypes(:,FRH) = (/1,         0,1,  0,0,  0,1/)             ! for rho()
      ftypes(:,FRB) = (/2,         0,1,  0,0,  0,1/)             ! for rhobk()
      ftypes(:,FPH) = (/1,        -1,2, -1,2,  0,0/)             ! for phi()
      ftypes(:,FJD) = (/minsp*3,  -1,2,  0,0, -1,2/)             ! for ajdg()
      ftypes(:,FRD) = (/minsp*2+1, 0,1,  0,0,  0,1/)             ! for rhodg()
      ftypes(:,FSC) = (/12,       -1,2, -1,2,  0,0/)             ! for ebsc()
      ftypes(1,FNR+1) = 0                                         ! terminator
!
      cfields(:) = (/FEB,FAJ,FRH,FRB,FJD,FRD,0/)
!     ctypes(:,:,*,***) = reshape((/downward,  upward/), (/3,2))) ! sample
!     ctypes(:,:,*,***) = reshape((/   s,r,n,  s, r,n/), (/3,2))) ! sample
      ctypes(:,:,1,CEB) = reshape((/   0,0,2, -1,-1,1/), (/3,2/)) ! for eb()
      ctypes(:,:,1,CAJ) = reshape((/  -1,2,3, -1,-4,3/), (/3,2/)) ! for aj()
      ctypes(:,:,1,CRH) = reshape((/   0,1,1,  0,-1,1/), (/3,2/)) ! for rho()
      ctypes(:,:,1,CRB) = reshape((/   0,1,1,  0,-1,1/), (/3,2/)) ! for rhobk()
      ctypes(:,:,1,CJD) = reshape((/  -1,2,3, -1,-4,3/), (/3,2/)) ! for ajdg()
      ctypes(:,:,1,CRD) = reshape((/   0,1,1,  0,-1,1/), (/3,2/)) ! for rhodg()
!
!      ctypes(:,:,2,CEB) = reshape((/   0,0,0,  0, 0,0/), (/3,2/)) ! for eb()
!      ctypes(:,:,2,CAJ) = reshape((/   0,0,0,  0, 0,0/), (/3,2/)) ! for aj()
!      ctypes(:,:,2,CRH) = reshape((/   0,0,0,  0, 0,0/), (/3,2/)) ! for rho()
!      ctypes(:,:,2,CRB) = reshape((/   0,0,0,  0, 0,0/), (/3,2/)) ! for rhobk()
!      ctypes(:,:,2,CJD) = reshape((/   0,0,0,  0, 0,0/), (/3,2/)) ! for ajdg()
!      ctypes(:,:,2,CRD) = reshape((/   0,0,0,  0, 0,0/), (/3,2/)) ! for rhodg()
      ctypes(:,:,2,CEB) = reshape((/   0,0,2, -1,-1,1/), (/3,2/)) ! for eb()
      ctypes(:,:,2,CAJ) = reshape((/  -1,2,3, -1,-4,3/), (/3,2/)) ! for aj()
      ctypes(:,:,2,CRH) = reshape((/   0,1,1,  0,-1,1/), (/3,2/)) ! for rho()
      ctypes(:,:,2,CRB) = reshape((/   0,1,1,  0,-1,1/), (/3,2/)) ! for rhobk()
      ctypes(:,:,2,CJD) = reshape((/  -1,2,3, -1,-4,3/), (/3,2/)) ! for ajdg()
      ctypes(:,:,2,CRD) = reshape((/   0,1,1,  0,-1,1/), (/3,2/)) ! for rhodg()


!------------------------------ 
      call oh3_init(sdid(:), minsp, MAXFRAC, nphgram(:,:,:), totalp(:,:), &
     &              pbuf(:), pbase(:), maxlocalp, mycomm, nbor(:,:,:), &
     &              pcoord(:), sdoms(:,:,:), scoord(:,:), nbound, bcond(:,:), &
     &              bounds(:,:,:), ftypes(:,:), cfields(:), ctypes(:,:,:,:), &
     &              fsizes(:,:,:), OHHELP_stats, OHHELP_repiter, &
     &              OHHELP_verbose)


!------------------------------ 
!      allocate(medi(1,          fsizes(1,1,FEB):fsizes(2,1,FEB), &
!     &                          fsizes(1,2,FEB):fsizes(2,2,FEB), &
!     &                          fsizes(1,3,FEB):fsizes(2,3,FEB)))
!      allocate(medi(1,          -1:nx+1, &
!     &                          -1:ny+1, &
!     &                          -1:nz+1))
      allocate( eb(6,           fsizes(1,1,FEB):fsizes(2,1,FEB), &
     &                          fsizes(1,2,FEB):fsizes(2,2,FEB), &
     &                          fsizes(1,3,FEB):fsizes(2,3,FEB), (nebfld+1)*2))
      allocate( ebsc(12,        fsizes(1,1,FSC):fsizes(2,1,FSC), &
     &                          fsizes(1,2,FSC):fsizes(2,2,FSC), &
     &                          fsizes(1,3,FSC):fsizes(2,3,FSC), (nebfld)*2))
      allocate( ebav(6,         fsizes(1,1,FEB):fsizes(2,1,FEB), &
     &                          fsizes(1,2,FEB):fsizes(2,2,FEB), &
     &                          fsizes(1,3,FEB):fsizes(2,3,FEB), nebfld+1))
      allocate( mp(6,           fsizes(1,1,FEB):fsizes(2,1,FEB), &
     &                          fsizes(1,2,FEB):fsizes(2,2,FEB), &
     &                          fsizes(1,3,FEB):fsizes(2,3,FEB), 4))
      allocate( db(3,           fsizes(1,1,FDB):fsizes(2,1,FDB), &
     &                          fsizes(1,2,FDB):fsizes(2,2,FDB), &
     &                          fsizes(1,3,FDB):fsizes(2,3,FDB), 2))
      allocate( aj(3,           fsizes(1,1,FAJ):fsizes(2,1,FAJ), &
     &                          fsizes(1,2,FAJ):fsizes(2,2,FAJ), &
     &                          fsizes(1,3,FAJ):fsizes(2,3,FAJ), 4))
      allocate(ajdg(minsp*3,    fsizes(1,1,FJD):fsizes(2,1,FJD), &
     &                          fsizes(1,2,FJD):fsizes(2,2,FJD), &
     &                          fsizes(1,3,FJD):fsizes(2,3,FJD), 2))
      allocate(ajav(minsp*3+6,  fsizes(1,1,FJD):fsizes(2,1,FJD), &
     &                          fsizes(1,2,FJD):fsizes(2,2,FJD), &
     &                          fsizes(1,3,FJD):fsizes(2,3,FJD)))
      allocate(rho(1,           fsizes(1,1,FRH):fsizes(2,1,FRH), &
     &                          fsizes(1,2,FRH):fsizes(2,2,FRH), &
     &                          fsizes(1,3,FRH):fsizes(2,3,FRH), 4))
      allocate(rhobk(2,         fsizes(1,1,FRB):fsizes(2,1,FRB), &
     &                          fsizes(1,2,FRB):fsizes(2,2,FRB), &
     &                          fsizes(1,3,FRB):fsizes(2,3,FRB), 3))
      allocate(rhodg(minsp*2+1, fsizes(1,1,FRD):fsizes(2,1,FRD), &
     &                          fsizes(1,2,FRD):fsizes(2,2,FRD), &
     &                          fsizes(1,3,FRD):fsizes(2,3,FRD), 2))
      allocate(rhoav(minsp*2+3, fsizes(1,1,FRD):fsizes(2,1,FRD), &
     &                          fsizes(1,2,FRD):fsizes(2,2,FRD), &
     &                          fsizes(1,3,FRD):fsizes(2,3,FRD)))
      allocate(phi(1,           fsizes(1,1,FPH):fsizes(2,1,FPH), &
     &                          fsizes(1,2,FPH):fsizes(2,2,FPH), &
     &                          fsizes(1,3,FPH):fsizes(2,3,FPH), 2))
      allocate(phiav(1,         fsizes(1,1,FPH):fsizes(2,1,FPH), &
     &                          fsizes(1,2,FPH):fsizes(2,2,FPH), &
     &                          fsizes(1,3,FPH):fsizes(2,3,FPH)))
      allocate(wrk(9,           fsizes(1,1,FAJ):fsizes(2,1,FAJ), &
     &                          fsizes(1,2,FAJ):fsizes(2,2,FAJ), &
     &                          fsizes(1,3,FAJ):fsizes(2,3,FAJ)))
      allocate(colf(1,          fsizes(1,1,FPH):fsizes(2,1,FPH), &
     &                          fsizes(1,2,FPH):fsizes(2,2,FPH), &
     &                          fsizes(1,3,FPH):fsizes(2,3,FPH), 2))
      allocate(pmlc(6,          -1:max(nx,ny,nz)+1,2))


!------------------------------ sub-domain/array dimensions
      nxsd = sdoms(2,1,sdid(1)+1) - sdoms(1,1,sdid(1)+1)
      nysd = sdoms(2,2,sdid(1)+1) - sdoms(1,2,sdid(1)+1)
      nzsd = sdoms(2,3,sdid(1)+1) - sdoms(1,3,sdid(1)+1)


!------------------------------ 
      do n=1,FNR
        lxsd(n) = fsizes(2,1,n) - fsizes(1,1,n) + 1
        lysd(n) = fsizes(2,2,n) - fsizes(1,2,n) + 1
        lzsd(n) = fsizes(2,3,n) - fsizes(1,3,n) + 1
      end do


!------------------------------ 
      if(jobnum(1).gt.0) then
        if(jobnum(1).eq.1) then
          write(filename,'(a,i4.4,a)') './SNAPSHOT0/esdat', myid, '.h5'
        else
          write(filename,'(a,i4.4,a)') './SNAPSHOT0/emdat', myid, '.h5'
        end if
        call hdfopen(filename,fileid,DFACC_READ)
!
        dsname = 'sdoms'
        dims(1) = 2; dims(2) = 3
        call read2i(fileid,dsname,dims(1:2),inttmp(1:2,1:3),stats0,stats1)
        if(inttmp(1,1).ne.sdoms(1,1,1).or.inttmp(2,1).ne.sdoms(2,1,1).or. &
       &   inttmp(1,2).ne.sdoms(1,2,1).or.inttmp(2,2).ne.sdoms(2,2,1).or. &
       &   inttmp(1,3).ne.sdoms(1,3,1).or.inttmp(2,3).ne.sdoms(2,3,1)) then
          if(myid.eq.0) print*, "sdoms is inconsistent with continued job data: STOP"
          stop
        end if
!
        call hdfclose(fileid,stats0)
      end if

!------------------------------ 
      if(myid.eq.0) print*,'maxlocalp=',maxlocalp
!
      if(myid.eq.0) print*,'EB:', (size(eb) + size(ebav) + size(mp))*8*1e-6
      if(myid.eq.0) print*,'  ',fsizes(1,1,FEB),fsizes(2,1,FEB),lxsd(FEB)
      if(myid.eq.0) print*,'  ',fsizes(1,2,FEB),fsizes(2,2,FEB),lysd(FEB)
      if(myid.eq.0) print*,'  ',fsizes(1,3,FEB),fsizes(2,3,FEB),lzsd(FEB)
      if(myid.eq.0) print*,'DB:', (size(db))*8*1e-6
      if(myid.eq.0) print*,'  ',fsizes(1,1,FDB),fsizes(2,1,FDB),lxsd(FDB)
      if(myid.eq.0) print*,'  ',fsizes(1,2,FDB),fsizes(2,2,FDB),lysd(FDB)
      if(myid.eq.0) print*,'  ',fsizes(1,3,FDB),fsizes(2,3,FDB),lzsd(FDB)
      if(myid.eq.0) print*,'AJ:', (size(aj) + size(wrk))*8*1e-6
      if(myid.eq.0) print*,'  ',fsizes(1,1,FAJ),fsizes(2,1,FAJ),lxsd(FAJ)
      if(myid.eq.0) print*,'  ',fsizes(1,2,FAJ),fsizes(2,2,FAJ),lysd(FAJ)
      if(myid.eq.0) print*,'  ',fsizes(1,3,FAJ),fsizes(2,3,FAJ),lzsd(FAJ)
      if(myid.eq.0) print*,'JD:', (size(ajdg) + size(ajav))*8*1e-6
      if(myid.eq.0) print*,'  ',fsizes(1,1,FJD),fsizes(2,1,FJD),lxsd(FJD)
      if(myid.eq.0) print*,'  ',fsizes(1,2,FJD),fsizes(2,2,FJD),lysd(FJD)
      if(myid.eq.0) print*,'  ',fsizes(1,3,FJD),fsizes(2,3,FJD),lzsd(FJD)
      if(myid.eq.0) print*,'RH:', (size(rho))*8*1e-6
      if(myid.eq.0) print*,'  ',fsizes(1,1,FRH),fsizes(2,1,FRH),lxsd(FRH)
      if(myid.eq.0) print*,'  ',fsizes(1,2,FRH),fsizes(2,2,FRH),lysd(FRH)
      if(myid.eq.0) print*,'  ',fsizes(1,3,FRH),fsizes(2,3,FRH),lzsd(FRH)
      if(myid.eq.0) print*,'RB:', (size(rhobk))*8*1e-6
      if(myid.eq.0) print*,'  ',fsizes(1,1,FRB),fsizes(2,1,FRB),lxsd(FRB)
      if(myid.eq.0) print*,'  ',fsizes(1,2,FRB),fsizes(2,2,FRB),lysd(FRB)
      if(myid.eq.0) print*,'  ',fsizes(1,3,FRB),fsizes(2,3,FRB),lzsd(FRB)
      if(myid.eq.0) print*,'RD:', (size(rhodg) + size(rhoav))*8*1e-6
      if(myid.eq.0) print*,'  ',fsizes(1,1,FRD),fsizes(2,1,FRD),lxsd(FRD)
      if(myid.eq.0) print*,'  ',fsizes(1,2,FRD),fsizes(2,2,FRD),lysd(FRD)
      if(myid.eq.0) print*,'  ',fsizes(1,3,FRD),fsizes(2,3,FRD),lzsd(FRD)
      if(myid.eq.0) print*,'PH:', (size(phi) + size(phiav) + size(colf))*8*1e-6
      if(myid.eq.0) print*,'  ',fsizes(1,1,FPH),fsizes(2,1,FPH),lxsd(FPH)
      if(myid.eq.0) print*,'  ',fsizes(1,2,FPH),fsizes(2,2,FPH),lysd(FPH)
      if(myid.eq.0) print*,'  ',fsizes(1,3,FPH),fsizes(2,3,FPH),lzsd(FPH)
      if(myid.eq.0) print*,'Total allocated field-array size (MB/proc):', &
     &  (size(eb) + size(ebsc) + size(ebav) + size(mp) + size(pmlc) + size(db) &
     & + size(aj) + size(wrk) + size(ajdg) + size(ajav) + size(rho) + size(rhobk) &
     & + size(rhodg) + size(rhoav) + size(phi) + size(phiav) + size(colf))*8*1e-6
      if(myid.eq.0) print*,'Total allocated particle-buffer size (MB/proc):', &
     &  maxlocalp*8*1e-6

!      if(myid.eq.0) print*,'bounds:'
!      if(myid.eq.0) print*,'  ',bounds


  return
  end subroutine ohinit
