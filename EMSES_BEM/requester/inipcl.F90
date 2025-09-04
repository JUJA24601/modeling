#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine inipcl
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   I N I P C L
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .  this subroutine gives an initial setting of particles.  .
!   ............................................................
!
!-------------------- parameter and common blocks
  use oh_type
  use paramt
  use allcom
  use hdf
! #define MCW MPI_COMM_WORLD
#define MCW CTCA_subcomm
  implicit none
!
  integer(kind=8) :: m,mm, ns,ne
  integer(kind=4) :: is
  integer(kind=4) :: iran, nran, nranmax
  integer(kind=4) :: ipc, jpc
  integer(kind=4) :: icon
  integer(kind=4) :: ps, rid
  integer(kind=4) :: ierr
  integer(kind=4) :: xl,xu, yl,yu, zl,zu
  integer(kind=4) :: ngx,ngy,ngz
!  integer(kind=4) :: oh3_map_region_to_node
  real(kind=8) :: dxl,dxu, dyl,dyu, dzl,dzu
  real(kind=8) :: gustep
  real(kind=8) :: vpesq
  real(kind=8) :: vxtemp, vytemp, vztemp
  real(kind=8) :: disp1,disp2,disp3
  real(kind=8) :: radsq(inpc),radsqi(inpc)
  real(kind=8) :: rbwlsq(3),rbwlsqi(3), rdomsq
  real(kind=8) :: xmove,ymove,zmove, xmovei,ymovei,zmovei
  real(kind=8) :: xsepa,ysepa,zsepa
  real(kind=8) :: xpast0,ypast0,zpast0, xpast1,ypast1,zpast1
  real(kind=8) :: vxyz, vxyztmp(ispec)
  real(kind=8) :: tfrac, tfrac1,tfrac2,tfrac3,tfrac4,tfrac5,tfrac6
  real(kind=8) :: termA,termB
  type(particle3) :: ptmp
!
  logical :: pcond(inpc)
!
  integer(kind=4) :: inttmp(ispec)
  integer(kind=HID_T) :: fileid
  integer(kind=4) :: stats0,stats1
  integer(kind=8) :: dims(1)
  real(kind=8) :: rens(4)
  character(len=30) :: filename,dsname


!-------------------- 
      ps = 3


!-------------------- 
      xl = sdoms(1,1,sdid(1)+1); xu = sdoms(2,1,sdid(1)+1)
      yl = sdoms(1,2,sdid(1)+1); yu = sdoms(2,2,sdid(1)+1)
      zl = sdoms(1,3,sdid(1)+1); zu = sdoms(2,3,sdid(1)+1)
      ngx = xu - xl; ngy = yu - yl; ngz = zu - zl
      do ipc=1,npc
        pcond(ipc) = (ncond(ipc).gt.0)
      end do


!-------------------- 
      do ipc=1,npc
        if(geotype(ipc).eq.2) then
          if(cylinder(ipc)%radius.gt.0.0d0) then
            radsq(ipc) = cylinder(ipc)%radius*cylinder(ipc)%radius
            radsqi(ipc) = 1.0d0/radsq(ipc)
          else
            radsq(ipc) = 0.0d0
            radsqi(ipc) = 1.0d100
          end if
        else if(geotype(ipc).eq.3) then
          if(sphere(ipc)%radius.gt.0.0d0) then
            radsq(ipc) = sphere(ipc)%radius*sphere(ipc)%radius
            radsqi(ipc) = 1.0d0/radsq(ipc)
          else
            radsq(ipc) = 0.0d0
            radsqi(ipc) = 1.0d100
          end if
        end if
      end do
      if(rbowl(1).gt.0.0d0) then
        rbwlsq(1:3) = rbowl(1:3)*rbowl(1:3)
        rbwlsqi(1:3) = 1.0d0/rbwlsq(1:3)
      else
        rbwlsq(1:3) = 0.0d0
        rbwlsqi(1:3) = 1.0d100
      end if
      if(rdome.gt.0.0d0) then
        rdomsq = rdome*rdome
      else
        rdomsq = 0.0d0
      end if
      if(rhole.gt.0.0d0) then
        rholesq = rhole*rhole
        xlwall = xholec - rhole
        xuwall = xholec + rhole
      else
        rholesq = 0.0d0
        xlwall = -abs(dmiss)
        xuwall = +abs(dmiss)
      end if


    if(jobnum(1).eq.0) then
!****************************** top of species loop 1
      nphgram(:,:,:) = 0
      ne = 0
ISL1: do is=1,nspec
        totalp(is,1) = npin(is)/nnode
        if(myid.lt.mod(npin(is),nnode)) then
          totalp(is,1) = totalp(is,1) + 1
        end if
        if(is.eq.4.and.myid.ge.0.and.myid.le.3) then
          totalp(is,1) = 1
        else if(is.eq.4.and.myid.gt.3) then
          totalp(is,1) = 0
        end if
        ns = ne + 1
        ne = ne + totalp(is,1)
!
        nranmax = int(size(dranu))
        if(type_rdist(is).eq.0) then
!          call RANU0(dranu,totalp(is,1),icon)
!          if(icon.ne.0) print*, "Warning(RANU0): myid,icon=",myid,icon
!          do m=ns,ne
!            pbuf(m)%vx = dranu(m-ns+1)
!          end do
!          call RANU0(dranu,totalp(is,1),icon)
!          if(icon.ne.0) print*, "Warning(RANU0): myid,icon=",myid,icon
!          do m=ns,ne
!            pbuf(m)%vy = dranu(m-ns+1)
!          end do
!          call RANU0(dranu,totalp(is,1),icon)
!          if(icon.ne.0) print*, "Warning(RANU0): myid,icon=",myid,icon
!          do m=ns,ne
!            pbuf(m)%vz = dranu(m-ns+1)
!          end do
!
          m = ns
          do mm=0,totalp(is,1)-1
!           if(is.ne.4) then
            iran = mod(mm,nranmax) + 1
            if(iran.eq.1) then
              nran = min(totalp(is,1),nranmax)
              call RANU0(dranu,nran,icon)
              pbuf(ns:ns+nran-1)%vx = dranu(1:nran)
              call RANU0(dranu,nran,icon)
              pbuf(ns:ns+nran-1)%vy = dranu(1:nran)
              call RANU0(dranu,nran,icon)
              pbuf(ns:ns+nran-1)%vz = dranu(1:nran)
            end if
            pbuf(m)%x = slx*pbuf(ns+iran-1)%vx
            pbuf(m)%y = sly*pbuf(ns+iran-1)%vy
            pbuf(m)%z = slz*pbuf(ns+iran-1)%vz
!           else
!            if(myid.eq.0) then!if(mm.eq.0) then
!              pbuf(m)%vx = +0.0d0
!              pbuf(m)%vy = +0.0d0
!              pbuf(m)%vz = -1.5d0
!              pbuf(m)%x = 16.0d0
!              pbuf(m)%y = 32.0d0
!              pbuf(m)%z = 101.d0
!            else if(myid.eq.1) then!else if(mm.eq.1) then
!              pbuf(m)%vx = +0.0d0
!              pbuf(m)%vy = +0.0d0
!              pbuf(m)%vz = -1.5d0
!              pbuf(m)%x = 48.0d0
!              pbuf(m)%y = 32.0d0
!              pbuf(m)%z = 101.d0
!            else if(myid.eq.2) then!else if(mm.eq.2) then
!              pbuf(m)%vx = +1.5d0
!              pbuf(m)%vy = +0.0d0
!              pbuf(m)%vz = +0.0d0
!              pbuf(m)%x = 37.0d0
!              pbuf(m)%y = 32.0d0
!              pbuf(m)%z = 90.0d0
!            else if(myid.eq.3) then!else if(mm.eq.3) then
!              pbuf(m)%vx = +0.0d0
!              pbuf(m)%vy = +0.0d0
!              pbuf(m)%vz = -1.5d0
!              pbuf(m)%x = 32.0d0
!              pbuf(m)%y = 32.0d0
!              pbuf(m)%z = 77.0d0
!            end if
!            print*, "Special v:", pbuf(m)%x,pbuf(m)%y,pbuf(m)%z,pbuf(m)%vx,pbuf(m)%vy,pbuf(m)%vz
!           end if
            ptmp%x = pbuf(m)%x/dr
            ptmp%y = pbuf(m)%y/dr
            ptmp%z = pbuf(m)%z/dr
            rid = oh3_map_particle_to_subdomain &
           &        (ptmp%x,ptmp%y,ptmp%z)
            pbuf(m)%nid = rid
!
            gustep = 2.0d0
            xmove = 0.1d0; ymove = 0.1d0; zmove = 0.1d0
#include "defsurf2.fnc"
#include "defbody.fnc"
!
            if(radring(2).gt.0.0d0) then
              if((pbuf(m)%x-rring(1))**2+(pbuf(m)%y-rring(2))**2.gt. &
             &    radring(2)*radring(2)) then
                pbuf(m)%x = 0.0d0
                pbuf(m)%y = 0.0d0
                pbuf(m)%z = 0.0d0
                pbuf(m)%nid = -1
              end if
            end if
!
            if(pbuf(m)%preside.lt.0) then
              pbuf(m)%x = 0.0d0
              pbuf(m)%y = 0.0d0
              pbuf(m)%z = 0.0d0
              pbuf(m)%nid = -1
              pbuf(m)%preside = 0
            end if
!
            if(pbuf(m)%nid.ge.0) then
              nphgram(pbuf(m)%nid+1,is,1) = nphgram(pbuf(m)%nid+1,is,1) + 1
              m = m + 1
            end if
          end do
!
        else if(type_rdist(is).eq.1) then
!          call RANU0(dranu,totalp(is,1),icon)
!          if(icon.ne.0) print*, "Warning(RANU0): myid,icon=",myid,icon
!          do m=ns,ne
!            pbuf(m)%vx = dranu(m-ns+1)
!          end do
!          call RANU0(dranu,totalp(is,1),icon)
!          if(icon.ne.0) print*, "Warning(RANU0): myid,icon=",myid,icon
!          do m=ns,ne
!            pbuf(m)%vy = dranu(m-ns+1)
!          end do
!          call RANU0(dranu,totalp(is,1),icon)
!          if(icon.ne.0) print*, "Warning(RANU0): myid,icon=",myid,icon
!          do m=ns,ne
!            pbuf(m)%vz = dranu(m-ns+1)
!          end do
!
          m = ns
          do mm=0,totalp(is,1)-1
            iran = mod(mm,nranmax) + 1
            if(iran.eq.1) then
              nran = min(totalp(is,1),nranmax)
              call RANU0(dranu,nran,icon)
              pbuf(ns:ns+nran-1)%vx = dranu(1:nran)
              call RANU0(dranu,nran,icon)
              pbuf(ns:ns+nran-1)%vy = dranu(1:nran)
              call RANU0(dranu,nran,icon)
              pbuf(ns:ns+nran-1)%vz = dranu(1:nran)
            end if
            pbuf(m)%x = ngx*pbuf(ns+iran-1)%vx + xl
            pbuf(m)%y = ngy*pbuf(ns+iran-1)%vy + yl
            pbuf(m)%z = ngz*pbuf(ns+iran-1)%vz + zl
            ptmp%x = pbuf(m)%x
            ptmp%y = pbuf(m)%y
            ptmp%z = pbuf(m)%z
            rid = oh3_map_particle_to_subdomain &
           &        (ptmp%x,ptmp%y,ptmp%z)
            pbuf(m)%nid = rid
!
            xmove = 0.1d0; ymove = 0.1d0; zmove = 0.1d0
#include "defsurf2.fnc"
#include "defbody.fnc"
!
            if(radring(2).gt.0.0d0) then
              if((pbuf(m)%x-rring(1))**2+(pbuf(m)%y-rring(2))**2.gt. &
             &    radring(2)*radring(2)) then
                pbuf(m)%x = 0.0d0
                pbuf(m)%y = 0.0d0
                pbuf(m)%z = 0.0d0
                pbuf(m)%nid = -1
              end if
            end if
!
            if(pbuf(m)%nid.eq.-1.or.pbuf(m)%preside.lt.0) then
              pbuf(m)%x = 0.0d0
              pbuf(m)%y = 0.0d0
              pbuf(m)%z = 0.0d0
              pbuf(m)%nid = -1
              pbuf(m)%preside = 0
            end if
!
            if(pbuf(m)%nid.ge.0) then
              nphgram(pbuf(m)%nid+1,is,1) = nphgram(pbuf(m)%nid+1,is,1) + 1
              m = m + 1
            end if
          end do
!
        else if(type_rdist(is).eq.2) then
!          if(sprd_rdist(1,is).ge.1.0d0) then
!            call RANN0(0.0d0,1.0d0,dranu,totalp(is,1),icon)
!            if(icon.ne.0) print*, "Warning(RANN0): myid,icon=",myid,icon 
!          else
!            call RANU0(dranu,totalp(is,1),icon)
!            if(icon.ne.0) print*, "Warning(RANU0): myid,icon=",myid,icon
!          end if
!          do m=ns,ne
!            pbuf(m)%vx = dranu(m-ns+1)
!          end do
!          if(sprd_rdist(2,is).ge.1.0d0) then
!            call RANN0(0.0d0,1.0d0,dranu,totalp(is,1),icon)
!            if(icon.ne.0) print*, "Warning(RANN0): myid,icon=",myid,icon 
!          else
!            call RANU0(dranu,totalp(is,1),icon)
!            if(icon.ne.0) print*, "Warning(RANU0): myid,icon=",myid,icon
!          end if
!          do m=ns,ne
!            pbuf(m)%vy = dranu(m-ns+1)
!          end do
!          if(sprd_rdist(3,is).ge.1.0d0) then
!            call RANN0(0.0d0,1.0d0,dranu,totalp(is,1),icon)
!            if(icon.ne.0) print*, "Warning(RANN0): myid,icon=",myid,icon 
!          else
!            call RANU0(dranu,totalp(is,1),icon)
!            if(icon.ne.0) print*, "Warning(RANU0): myid,icon=",myid,icon
!          end if
!          do m=ns,ne
!            pbuf(m)%vz = dranu(m-ns+1)
!          end do
!
          m = ns
          do mm=0,totalp(is,1)-1
            iran = mod(mm,nranmax) + 1
            if(iran.eq.1) then
              nran = min(totalp(is,1),nranmax)
              if(sprd_rdist(1,is).ge.1.0d0) then
                call RANN0(0.0d0,1.0d0,dranu,nran,icon)
              else
                call RANU0(dranu,nran,icon)
              end if
              pbuf(ns:ns+nran-1)%vx = dranu(1:nran)
              if(sprd_rdist(2,is).ge.1.0d0) then
                call RANN0(0.0d0,1.0d0,dranu,nran,icon)
              else
                call RANU0(dranu,nran,icon)
              end if
              pbuf(ns:ns+nran-1)%vy = dranu(1:nran)
              if(sprd_rdist(3,is).ge.1.0d0) then
                call RANN0(0.0d0,1.0d0,dranu,nran,icon)
              else
                call RANU0(dranu,nran,icon)
              end if
              pbuf(ns:ns+nran-1)%vz = dranu(1:nran)
            end if

            if(sprd_rdist(1,is).ge.1.0d0) then
              pbuf(m)%x = sprd_rdist(1,is)*pbuf(ns+mm)%vx &
             &          + cntr_rdist(1,is)
            else
              pbuf(m)%x = slx*pbuf(ns+mm)%vx
            end if
            if(sprd_rdist(2,is).ge.1.0d0) then
              pbuf(m)%y = sprd_rdist(2,is)*pbuf(ns+mm)%vy &
             &          + cntr_rdist(2,is)
            else
              pbuf(m)%y = sly*pbuf(ns+mm)%vy
            end if
            if(sprd_rdist(3,is).ge.1.0d0) then
              pbuf(m)%z = sprd_rdist(3,is)*pbuf(ns+mm)%vz &
             &          + cntr_rdist(3,is)
            else
              pbuf(m)%z = slz*pbuf(ns+mm)%vz
            end if
            ptmp%x = pbuf(m)%x/dr
            ptmp%y = pbuf(m)%y/dr
            ptmp%z = pbuf(m)%z/dr
            rid = oh3_map_particle_to_subdomain &
           &        (ptmp%x,ptmp%y,ptmp%z)
            pbuf(m)%nid = rid
!
            xmove = 0.1d0; ymove = 0.1d0; zmove = 0.1d0
#include "defsurf2.fnc"
#include "defbody.fnc"
!
            if(radring(2).gt.0.0d0) then
              if((pbuf(m)%x-rring(1))**2+(pbuf(m)%y-rring(2))**2.gt. &
             &    radring(2)*radring(2)) then
                pbuf(m)%x = 0.0d0
                pbuf(m)%y = 0.0d0
                pbuf(m)%z = 0.0d0
                pbuf(m)%nid = -1
              end if
            end if
!
            if(pbuf(m)%nid.eq.-1.or.pbuf(m)%preside.lt.0) then
              pbuf(m)%x = 0.0d0
              pbuf(m)%y = 0.0d0
              pbuf(m)%z = 0.0d0
              pbuf(m)%nid = -1
              pbuf(m)%preside = 0
            end if
!
            if(pbuf(m)%nid.ge.0) then
              nphgram(pbuf(m)%nid+1,is,1) = nphgram(pbuf(m)%nid+1,is,1) + 1
              m = m + 1
            end if
          end do
        end if
        totalp(is,1) = m - ns
        ne = m - 1
!
        if(lcgamma(is).ge.1.0d0.or.lcgamma(is).lt.0.0d0.or. &
       &   lcbeta(is).le.0.0d0) then
          call RANN0(spe(is)*dcos(speth(is)/180.0d0*pi),peth(is), &
         &            dranu,totalp(is,1),icon)
          if(icon.ne.0) print*, "Warning(RANN0): myid,icon=",myid,icon 
          do m=ns,ne
            pbuf(m)%vx = dranu(m-ns+1)
          end do
          call RANN0(spe(is)*dsin(speth(is)/180.0d0*pi),peth(is), &
         &            dranu,totalp(is,1),icon)
          if(icon.ne.0) print*, "Warning(RANN0): myid,icon=",myid,icon 
          do m=ns,ne
            pbuf(m)%vy = dranu(m-ns+1)
          end do
        else
          print*, "WARNING!!"
          call RANN0(spe(is)*dcos(speth(is)/180.0d0*pi),peth(is), &
         &            dranu,int(totalp(is,1)*lcgamma(is)),icon)
          if(icon.ne.0) print*, "Warning(RANN0): myid,icon=",myid,icon
          do m=ns,ns-1+int(totalp(is,1)*lcgamma(is))
            pbuf(m)%vx = dranu(m-ns+1)
          end do
          call RANN0(spe(is)*dsin(speth(is)/180.0d0*pi),peth(is), &
         &            dranu,int(totalp(is,1)*lcgamma(is)),icon)
          if(icon.ne.0) print*, "Warning(RANN0): myid,icon=",myid,icon
          do m=ns,ns-1+int(totalp(is,1)*lcgamma(is))
            pbuf(m)%vy = dranu(m-ns+1)
          end do
!
          m = ns + int(totalp(is,1)*lcgamma(is))
          do while(m.le.ne)
            call RANN0(spe(is),peth(is),dranu(1:1),1,icon)
            call RANN0(spe(is),peth(is),dranu(2:2),1,icon)
            vpesq = dranu(1)*dranu(1) + dranu(2)*dranu(2)
            call RANU0(dranu(3:3),1,icon)
            if(exp(-vpesq/lcbeta(is)/peth(is)/peth(is)).le.dranu(3)) then
              pbuf(m)%vx = dranu(1) + spe(is)*dcos(speth(is)/180.0d0*pi)
              pbuf(m)%vy = dranu(2) + spe(is)*dsin(speth(is)/180.0d0*pi)
              m = m + 1
            end if
          end do
        end if

!       if(is.ne.4) then
        call RANN0(spa(is),path(is),dranu,totalp(is,1),icon)
        if(icon.ne.0) print*, "Warning(RANN0): myid,icon=",myid,icon 
        do m=ns,ne
          pbuf(m)%vz = dranu(m-ns+1)
        end do
!       else
!        if(myid.eq.0) then
!          pbuf(ns+0)%vx = +0.0d0
!          pbuf(ns+0)%vy = +0.0d0
!          pbuf(ns+0)%vz = -1.5d0
!        else if(myid.eq.1) then
!          pbuf(ns+0)%vx = +0.0d0
!          pbuf(ns+0)%vy = +0.0d0
!          pbuf(ns+0)%vz = -1.5d0
!        else if(myid.eq.2) then
!          pbuf(ns+0)%vx = +1.5d0
!          pbuf(ns+0)%vy = +0.0d0
!          pbuf(ns+0)%vz = +0.0d0
!        else if(myid.eq.3) then
!          pbuf(ns+0)%vx = +0.0d0
!          pbuf(ns+0)%vy = +0.0d0
!          pbuf(ns+0)%vz = -1.5d0
!        end if
!        print*, "Special v:", pbuf(m)%x,pbuf(m)%y,pbuf(m)%z,pbuf(m)%vx,pbuf(m)%vy,pbuf(m)%vz
!       end if
!
        do m=ns,ne
          vxtemp = pbuf(m)%vx*t11 + pbuf(m)%vy*t12 + pbuf(m)%vz*t13
          vytemp = pbuf(m)%vx*t21 + pbuf(m)%vy*t22 + pbuf(m)%vz*t23
          vztemp = pbuf(m)%vx*t31 + pbuf(m)%vy*t32 + pbuf(m)%vz*t33
          vxtemp = vxtemp &
         &       + vdri(is)*sin(vdthz(is)/180.0d0*pi)*cos(vdthxy(is)/180.0d0*pi)
          vytemp = vytemp &
         &       + vdri(is)*sin(vdthz(is)/180.0d0*pi)*sin(vdthxy(is)/180.0d0*pi)
          vztemp = vztemp &
         &       + vdri(is)*cos(vdthz(is)/180.0d0*pi)
          pbuf(m)%vx = vxtemp + vdx(is)
          pbuf(m)%vy = vytemp + vdy(is)
          pbuf(m)%vz = vztemp + vdz(is)
        end do
!
!        if((ewmodel.eq.1.or.ewmodel.eq.2).and.is.eq.1.and.phiz.eq.0.0d0) then
!          do m=ns,ne
!            pbuf(m)%vx = pbuf(m)%vx + Ew(1,is,1)*qm(is)*omegaw(1)*id2omega*sin(0.0d0) &
!           &                        + Ew(2,is,1)*qm(is)*wc*id2omega*sin(0.0d0)
!            pbuf(m)%vy = pbuf(m)%vy - Ew(1,is,1)*qm(is)*wc*id2omega*cos(0.0d0) &
!           &                        - Ew(2,is,1)*qm(is)*omegaw(1)*id2omega*cos(0.0d0)
!            pbuf(m)%vz = pbuf(m)%vz
!          end do
!        end if
!
        do m=ns,ne
          pbuf(m)%spec = is
          pbuf(m)%preside = 0
          pbuf(m)%pid = 0
        end do
!
        if(ne.ge.ns) then
          do m=ns,ne
            vxyz = sqrt((pbuf(m)%vx-vdtx(is))*(pbuf(m)%vx-vdtx(is)) &
           &          + (pbuf(m)%vy-vdty(is))*(pbuf(m)%vy-vdty(is)) &
           &          + (pbuf(m)%vz-vdtz(is))*(pbuf(m)%vz-vdtz(is)))
            if(vxyztmp(is).lt.vxyz) vxyztmp(is) = vxyz
          end do
        else
          vxyztmp(is) = sqrt(path(is)*path(is)*25.0d0 &
           &               + peth(is)*peth(is)*25.0d0 &
           &               + peth(is)*peth(is)*25.0d0 )
        end if
      end do ISL1
!
      call MPI_Allreduce(vxyztmp,vxyzmax,minsp,MPI_REAL8,MPI_MAX,MCW,ierr)
      if(myid.eq.0) print*, "vxyzmax =",vxyzmax(1:minsp)
!
      do is=1,nspec
        vxyzmax(is) = vxyzmax(is)*1.25d0
        ivxyzmax(is) = 1.0d0/vxyzmax(is)
      end do
!
    else
!****************************** top of species loop 2
      nphgram(:,:,:) = 0
      if(jobnum(1).eq.1) then
        write(filename,'(a,i4.4,a)') './SNAPSHOT0/esdat', myid, '.h5'
      else
        write(filename,'(a,i4.4,a)') './SNAPSHOT0/emdat', myid, '.h5'
      end if
      call hdfopen(filename,fileid,DFACC_READ)
!
      dsname = 'np'
      dims(1) = nspec
      call read1i(fileid,dsname,dims(1:1),inttmp(1:nspec),stats0,stats1)
!
      dsname = 'rens'
      dims(1) = 4
      call read1d(fileid,dsname,dims(1:1),rens(1:4),stats0,stats1)
!
      dsname = 'px'
      dims(1) = sum(inttmp(1:nspec))
      call read1d(fileid,dsname,dims(1:1),pbuf(1:dims(1))%x,stats0,stats1)
!
      dsname = 'py'
      dims(1) = sum(inttmp(1:nspec))
      call read1d(fileid,dsname,dims(1:1),pbuf(1:dims(1))%y,stats0,stats1)
!
      dsname = 'pz'
      dims(1) = sum(inttmp(1:nspec))
      call read1d(fileid,dsname,dims(1:1),pbuf(1:dims(1))%z,stats0,stats1)
!
      dsname = 'pvx'
      dims(1) = sum(inttmp(1:nspec))
      call read1d(fileid,dsname,dims(1:1),pbuf(1:dims(1))%vx,stats0,stats1)
!
      dsname = 'pvy'
      dims(1) = sum(inttmp(1:nspec))
      call read1d(fileid,dsname,dims(1:1),pbuf(1:dims(1))%vy,stats0,stats1)
!
      dsname = 'pvz'
      dims(1) = sum(inttmp(1:nspec))
      call read1d(fileid,dsname,dims(1:1),pbuf(1:dims(1))%vz,stats0,stats1)
!
      dsname = 'pres'
      dims(1) = sum(inttmp(1:nspec))
      call read1i(fileid,dsname,dims(1:1),pbuf(1:dims(1))%preside,stats0,stats1)
!
      dsname = 'pis'
      dims(1) = sum(inttmp(1:nspec))
      call read1i(fileid,dsname,dims(1:1),pbuf(1:dims(1))%spec,stats0,stats1)
!
      totalp(1:nspec,1) = inttmp(1:nspec)
!
      call specsort
!
      do m=1,sum(inttmp(1:nspec))
        ptmp%x = pbuf(m)%x/rens(1)/dr
        ptmp%y = pbuf(m)%y/rens(1)/dr
        ptmp%z = pbuf(m)%z/rens(1)/dr
        rid = oh3_map_particle_to_subdomain &
       &        (ptmp%x,ptmp%y,ptmp%z)
        pbuf(m)%nid = rid
        if(rid.ge.0) &
       &  nphgram(rid+1,pbuf(m)%spec,1) = nphgram(rid+1,pbuf(m)%spec,1) + 1
        pbuf(m)%pid = 0
      end do
!
      call hdfclose(fileid,stats0)
    end if


  return
  end subroutine
