#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine psuper3(ps,ustep,sstep)
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   P S U P E R
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .  this subroutine supervises particle injection from      .
!   .  the outer/inner boundaries of the simulation box.       .
!   ............................................................
!
!-------------------- parameter and common block
  use oh_type
  use paramt
  use allcom
! #define MCW MPI_COMM_WORLD
#define MCW CTCA_subcomm
#define MSS MPI_STATUS_SIZE
  implicit none
! 
  integer(kind=4) :: m, l
  integer(kind=8) :: nprs,npre
  integer(kind=8) :: nofluxt(inpc)
  integer(kind=8) :: nofluxsum
  integer(kind=4) :: i,j,k, i1,j1,k1
  integer(kind=4) :: ib,jb,kb, i1b,j1b,k1b
  integer(kind=4) :: ia,ja,ka, i1a,j1a,k1a
  integer(kind=4) :: nns,nne, isw
  integer(kind=4) :: ie, nemitw
  integer(kind=4) :: nranmax
  integer(kind=4) :: iran,jran,kran, ipc, iepl, is, iss
  integer(kind=4) :: ustep, sstep
  integer(kind=4) :: addr0,addr1,addr2
  integer(kind=4) :: icon
  integer(kind=4) :: omni, algn
  integer(kind=4) :: xl,xu, yl,yu, zl,zu
  integer(kind=4) :: func
  integer(kind=4) :: idim, ps, rid
  integer(kind=4) :: ierr
!  integer(kind=4) :: oh3_map_region_to_node
  real(kind=8) :: fnemit
  real(kind=8) :: ground, radius, axis(2), cntr(3)
  real(kind=8) :: lxlej,lxuej, lylej,lyuej, lzlej,lzuej
  real(kind=8) :: dxl,dxu, dyl,dyu, dzl,dzu
  real(kind=8) :: xlocalb,ylocalb,zlocalb
  real(kind=8) :: xb,yb,zb
  real(kind=8) :: xlocala,ylocala,zlocala
  real(kind=8) :: xa,ya,za
  real(kind=8) :: xr,yr,zr
  real(kind=8) :: xd1,yd1,zd1
  real(kind=8) :: xd2,yd2,zd2
  real(kind=8) :: tslx,tsly, dustep, rustep, qs
  real(kind=8) :: vvx1,vvy1,vvz1
  real(kind=8) :: vvx2,vvy2,vvz2
  real(kind=8) :: vx1w1,vx1w2,vx1w3,vx1w4
  real(kind=8) :: vy1w1,vy1w2,vy1w3,vy1w4
  real(kind=8) :: vz1w1,vz1w2,vz1w3,vz1w4
  real(kind=8) :: vx2w1,vx2w2,vx2w3,vx2w4
  real(kind=8) :: vy2w1,vy2w2,vy2w3,vy2w4
  real(kind=8) :: vz2w1,vz2w2,vz2w3,vz2w4
  real(kind=8) :: csz,snz, csxy,snxy, xew,yew,zew
  real(kind=8) :: te11,te12,te13,te21,te22,te23,te31,te32,te33
  real(kind=8) :: arearatio(12)
  real(kind=8) :: betav, zdepth, rtmp
  real(kind=8) :: xsepa,ysepa,zsepa, xysepa
  real(kind=8) :: tan1,cos1
  real(kind=8) :: xlocal,ylocal,zlocal
  real(kind=8) :: x1,y1,z1, z2, xy1,xz1,yz1, xz2,yz2
  real(kind=8) :: v1,v2,v3,v4,v5,v6,v7,v8
  real(kind=8) :: v11,v12,v13,v21,v22,v23,v31,v32,v33,iv31
  real(kind=8) :: lbdome,lcdome,lddome,lbbowl,lcbowl,ldbowl
  real(kind=8) :: rdomsq,tdome,rbwlsq(3),rbwlsqi(3),tbowl
  real(kind=8) :: vxtmp,vytmp,vztmp
!  type(oh_particle) :: pinj
!  integer(kind=4) :: ninjct

  logical :: npbnd_x_injct,npbnd_y_injct,npbnd_z_injct
  logical :: pcond


!-------------------- 
      xl = sdoms(1,1,sdid(1)+1); xu = sdoms(2,1,sdid(1)+1)
      yl = sdoms(1,2,sdid(1)+1); yu = sdoms(2,2,sdid(1)+1)
      zl = sdoms(1,3,sdid(1)+1); zu = sdoms(2,3,sdid(1)+1)
      dxl = xl; dxu = xu; dyl = yl; dyu = yu; dzl = zl; dzu = zu


!-------------------- 
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
      else
        rholesq = 0.0d0
      end if


!-------------------- supervision 0
      if(sstep.eq.0) then
        nne = 0
  ISL2: do is=1,nspec
          nns = nne + 1
          nne = nne + nepl(is)
   IEPLL: do iepl=nns,nne
            ipc = ipc_ej(iepl)
            if(ncond(ipc).gt.0) then
              gcount(1)%chgacm(1:2,ipc) = &
             &  gcount(1)%chgacm(1:2,ipc) &
             &  - q(is)*nemit(iepl)*sqdscaled(ipc)
            end if
            gcount(1)%outflux(ipc,is) = &
           &  gcount(1)%outflux(ipc,is) + nemit(iepl)
            nemit(iepl) = 0
          end do IEPLL
        end do ISL2
        gcount(1)%chgacm(1:2,1) = &
             &  gcount(1)%chgacm(1:2,1) &
             &  - q(3)*nsecemit(1)*sqdscaled(1)
        gcount(1)%outflux(1,3) = &
           &  gcount(1)%outflux(1,3) + nsecemit(1)
        gcount(1)%chgacm(1:2,2) = &
             &  gcount(1)%chgacm(1:2,2) &
             &  - q(3)*nsecemit(2)*sqdscaled(1)
        gcount(1)%outflux(2,3) = &
           &  gcount(1)%outflux(2,3) + nsecemit(2)
        nsecemit(:) = 0
        sey(:,:) = 0.0d0
        return
      else if(sstep.eq.1) then
  ISL3: do is=1,nspec
          do ipc=1,npc
            if(ncond(ipc).gt.0) then
              gcount(1)%chgacm(1:2,ipc) = gcount(1)%chgacm(1:2,ipc) &
             &  + q(is)*gcount(1)%influx(ipc,is)*sqdscaled(ipc)
            end if
          end do
        end do ISL3
      end if


!-------------------- 
      if(istep.ge.nstep) return


!-------------------- supervision 1
!
!     injection of particles from inner/outer boundaries
!
!      do is=1,nspec
!        injctisw(is) = injct(is)*ustep
!      end do
!
!****************************** top of species loop 1
      npre = 0
      nne = 0
      dustep = ustep
      rustep = 1.0d0/dustep
      tslx = 2.0d0*slx
      tsly = 2.0d0*sly
ISL1: do is=1,nspec
        nprs = npre + 1
        npre = npre + npr(is)
        nns = nne + 1
        nne = nne + nepl(is)
!       ------------- 
        npbnd_x_injct = &
       &  (emmode.and.npbnd(1,is).eq.2.and.mtd_vbnd(1).eq.2)
        npbnd_y_injct = &
       &  (emmode.and.npbnd(2,is).eq.2.and.mtd_vbnd(2).eq.2)
        npbnd_z_injct = &
       &  (emmode.and.npbnd(3,is).eq.2.and.mtd_vbnd(3).eq.2)
!
!-------------------- 
        if(sstep.eq.1.and.nflag_emit(is).eq.0.and. &
       &   (npbnd(1,is).eq.2.or.npbnd(2,is).eq.2.or.npbnd(3,is).eq.2)) then
!
!============================== from outer boundary
!         ----------- x-bottom
          if(npbnd(1,is).eq.2.and.bared(1,1,sdid(1)+1).eq.1) then
            narrxf(is) = (int((istep + lstep)*arrxf(is),8) &
           &            - int((istep + lstep - injct(is))*arrxf(is),8))*ustep
            if(narrxf(is).gt.0) then
              call RANU0(dranu,narrxf(is)*5,icon)
              iran = 0
              do m=1,narrxf(is)
                iran = iran + 5
!               ----- pick up from velocity distribution
                addr1 = nprs + npr(is)*dranu(iran-4)
                addr2 = nprs + npr(is)*dranu(iran-3)
!               ----- assign position and velocity
                pinj(1)%x = xl
                pinj(1)%y = yl + (yu - yl)*dranu(iran-2)
                pinj(1)%z = zl + (zu - zl)*dranu(iran-1)
                if(pinj(1)%z.le.zssurf.or.pinj(1)%z.le.ubhole) cycle
                pinj(1)%vx = vxf(addr1)
                pinj(1)%vy = vyr(addr2)
                pinj(1)%vz = vzr(addr2)
!                pinj(1)%x = pinj(1)%x &
!               &          - pinj(1)%vx*dustep*dranu(iran)
!                pinj(1)%y = pinj(1)%y &
!               &          - pinj(1)%vy*dustep*dranu(iran)
!                pinj(1)%z = pinj(1)%z &
!               &          - pinj(1)%vz*dustep*dranu(iran)
                pinj(1)%nid = sdid(1)
                pinj(1)%spec = is
                pinj(1)%pid = 0
                pinj(1)%preside = -9
                call oh2_inject_particle(pinj(1))
                totalp(is,3) = totalp(is,3) + 1
!
                if(npbnd_x_injct) then
                  pinj(1)%nid = sdid(1)
                  pinj(1)%preside = REBD
                  call oh2_inject_particle(pinj(1))
                  totalp(is,3) = totalp(is,3) + 1
                end if
              end do
            end if
          end if
!         ----------- x-top
          if(npbnd(1,is).eq.2.and.bared(2,1,sdid(1)+1).eq.1) then
            narrxb(is) = (int((istep + lstep)*arrxb(is),8) &
           &            - int((istep + lstep - injct(is))*arrxb(is),8))*ustep
            if(narrxb(is).gt.0) then
              call RANU0(dranu,narrxb(is)*5,icon)
              iran = 0
              do m=1,narrxb(is)
                iran = iran + 5
!               ----- pick up from velocity distribution
                addr1 = nprs + npr(is)*dranu(iran-4)
                addr2 = nprs + npr(is)*dranu(iran-3)
!               ----- assign position and velocity
                pinj(1)%x = xu
                pinj(1)%y = yl + (yu - yl)*dranu(iran-2)
                pinj(1)%z = zl + (zu - zl)*dranu(iran-1)
                if(pinj(1)%z.le.zssurf.or.pinj(1)%z.le.ubhole) cycle
                pinj(1)%vx = vxb(addr1)
                pinj(1)%vy = vyr(addr2)
                pinj(1)%vz = vzr(addr2)
!                pinj(1)%x = pinj(1)%x &
!               &          - pinj(1)%vx*dustep*dranu(iran)
!                pinj(1)%y = pinj(1)%y &
!               &          - pinj(1)%vy*dustep*dranu(iran)
!                pinj(1)%z = pinj(1)%z &
!               &          - pinj(1)%vz*dustep*dranu(iran)
                pinj(1)%nid = sdid(1)
                pinj(1)%spec = is
                pinj(1)%preside = -9
                pinj(1)%pid = 0
                call oh2_inject_particle(pinj(1))
                totalp(is,3) = totalp(is,3) + 1
!
                if(npbnd_x_injct) then
                  pinj(1)%nid = sdid(1)
                  pinj(1)%preside = REBD
                  call oh2_inject_particle(pinj(1))
                  totalp(is,3) = totalp(is,3) + 1
                end if
              end do
            end if
          end if
!         ----------- y-bottom
          if(npbnd(2,is).eq.2.and.bared(1,2,sdid(1)+1).eq.1) then
            narryf(is) = (int((istep + lstep)*arryf(is),8) &
           &            - int((istep + lstep - injct(is))*arryf(is),8))*ustep
            if(narryf(is).gt.0) then
              call RANU0(dranu,narryf(is)*5,icon)
              iran = 0
              do m=1,narryf(is)
                iran = iran + 5
!               ----- pick up from velocity distribution
                addr1 = nprs + npr(is)*dranu(iran-4)
                addr2 = nprs + npr(is)*dranu(iran-3)
!               ----- assign position and velocity
                pinj(1)%x = xl + (xu - xl)*dranu(iran-2)
                pinj(1)%y = yl
                pinj(1)%z = zl + (zu - zl)*dranu(iran-1)
                if(pinj(1)%z.le.zssurf.or.pinj(1)%z.le.ubhole) cycle
                pinj(1)%vx = vxr(addr2)
                pinj(1)%vy = vyf(addr1)
                pinj(1)%vz = vzr(addr2)
!                pinj(1)%x = pinj(1)%x &
!               &          - pinj(1)%vx*dustep*dranu(iran)
!                pinj(1)%y = pinj(1)%y &
!               &          - pinj(1)%vy*dustep*dranu(iran)
!                pinj(1)%z = pinj(1)%z &
!               &          - pinj(1)%vz*dustep*dranu(iran)
                pinj(1)%nid = sdid(1)
                pinj(1)%spec = is
                pinj(1)%preside = -9
                pinj(1)%pid = 0
                call oh2_inject_particle(pinj(1))
                totalp(is,3) = totalp(is,3) + 1
!
                if(npbnd_y_injct) then
                  pinj(1)%nid = sdid(1)
                  pinj(1)%preside = REBD
                  call oh2_inject_particle(pinj(1))
                  totalp(is,3) = totalp(is,3) + 1
                end if
              end do
            end if
          end if
!         ----------- y-top
          if(npbnd(2,is).eq.2.and.bared(2,2,sdid(1)+1).eq.1) then
            narryb(is) = (int((istep + lstep)*arryb(is),8) &
           &            - int((istep + lstep - injct(is))*arryb(is),8))*ustep
            if(narryb(is).gt.0) then
              call RANU0(dranu,narryb(is)*5,icon)
              iran = 0
              do m=1,narryb(is)
                iran = iran + 5
!               ----- pick up from velocity distribution
                addr1 = nprs + npr(is)*dranu(iran-4)
                addr2 = nprs + npr(is)*dranu(iran-3)
!               ----- assign position and velocity
                pinj(1)%x = xl + (xu - xl)*dranu(iran-2)
                pinj(1)%y = yu
                pinj(1)%z = zl + (zu - zl)*dranu(iran-1)
                if(pinj(1)%z.le.zssurf.or.pinj(1)%z.le.ubhole) cycle
                pinj(1)%vx = vxr(addr2)
                pinj(1)%vy = vyb(addr1)
                pinj(1)%vz = vzr(addr2)
!                pinj(1)%x = pinj(1)%x &
!               &          - pinj(1)%vx*dustep*dranu(iran)
!                pinj(1)%y = pinj(1)%y &
!               &          - pinj(1)%vy*dustep*dranu(iran)
!                pinj(1)%z = pinj(1)%z &
!               &          - pinj(1)%vz*dustep*dranu(iran)
                pinj(1)%nid = sdid(1)
                pinj(1)%spec = is
                pinj(1)%preside = -9
                pinj(1)%pid = 0
                call oh2_inject_particle(pinj(1))
                totalp(is,3) = totalp(is,3) + 1
!
                if(npbnd_y_injct) then
                  pinj(1)%nid = sdid(1)
                  pinj(1)%preside = REBD
                  call oh2_inject_particle(pinj(1))
                  totalp(is,3) = totalp(is,3) + 1
                end if
              end do
            end if
          end if
!         ----------- z-bottom
          if(npbnd(3,is).eq.2.and.bared(1,3,sdid(1)+1).eq.1) then
            narrzf(is) = (int((istep + lstep)*arrzf(is),8) &
           &            - int((istep + lstep - injct(is))*arrzf(is),8))*ustep
            if(narrzf(is).gt.0) then
              call RANU0(dranu,narrzf(is)*5,icon)
              iran = 0
              do m=1,narrzf(is)
                iran = iran + 5
!               ----- pick up from velocity distribution
                addr1 = nprs + npr(is)*dranu(iran-4)
                addr2 = nprs + npr(is)*dranu(iran-3)
!               ----- assign position and velocity
                pinj(1)%x = xl + (xu - xl)*dranu(iran-2)
                pinj(1)%y = yl + (yu - yl)*dranu(iran-1)
                pinj(1)%z = zl
                if(pinj(1)%z.le.zssurf.or.pinj(1)%z.le.ubhole) cycle
                pinj(1)%vx = vxr(addr2)
                pinj(1)%vy = vyr(addr2)
                pinj(1)%vz = vzf(addr1)
!                pinj(1)%x = pinj(1)%x &
!               &          - pinj(1)%vx*dustep*dranu(iran)
!                pinj(1)%y = pinj(1)%y &
!               &          - pinj(1)%vy*dustep*dranu(iran)
!                pinj(1)%z = pinj(1)%z &
!               &          - pinj(1)%vz*dustep*dranu(iran)
                pinj(1)%nid = sdid(1)
                pinj(1)%spec = is
                pinj(1)%preside = -9
                pinj(1)%pid = 0
                call oh2_inject_particle(pinj(1))
                totalp(is,3) = totalp(is,3) + 1
!
                if(npbnd_z_injct) then
                  pinj(1)%nid = sdid(1)
                  pinj(1)%preside = REBD
                  call oh2_inject_particle(pinj(1))
                  totalp(is,3) = totalp(is,3) + 1
                end if
              end do
            end if
          end if
!         ----------- z-top
          if(npbnd(3,is).eq.2.and.bared(2,3,sdid(1)+1).eq.1) then
            narrzb(is) = (int((istep + lstep)*arrzb(is),8) &
           &            - int((istep + lstep - injct(is))*arrzb(is),8))*ustep
            if(narrzb(is).gt.0) then
              call RANU0(dranu,narrzb(is)*5,icon)
              iran = 0
              do m=1,narrzb(is)
                iran = iran + 5
!               ------- pick up from velocity distribution
                addr1 = nprs + npr(is)*dranu(iran-4)
                addr2 = nprs + npr(is)*dranu(iran-3)
!               ------- assign position and velocity
                pinj(1)%x = xl + (xu - xl)*dranu(iran-2)
                pinj(1)%y = yl + (yu - yl)*dranu(iran-1)
                pinj(1)%z = zu
                if(pinj(1)%z.le.zssurf.or.pinj(1)%z.le.ubhole) cycle
                pinj(1)%vx = vxr(addr2)
                pinj(1)%vy = vyr(addr2)
                pinj(1)%vz = vzb(addr1)
!                pinj(1)%x = pinj(1)%x &
!               &          - pinj(1)%vx*dustep*dranu(iran)
!                pinj(1)%y = pinj(1)%y &
!               &          - pinj(1)%vy*dustep*dranu(iran)
!                pinj(1)%z = pinj(1)%z &
!               &          - pinj(1)%vz*dustep*dranu(iran)
                pinj(1)%nid = sdid(1)
                pinj(1)%spec = is
                pinj(1)%preside = -9
                pinj(1)%pid = 0
                call oh2_inject_particle(pinj(1))
                totalp(is,3) = totalp(is,3) + 1
!
                if(npbnd_z_injct) then
                  pinj(1)%nid = sdid(1)
                  pinj(1)%preside = REBD
                  call oh2_inject_particle(pinj(1))
                  totalp(is,3) = totalp(is,3) + 1
                end if
              end do
            end if
          end if

!-------------------- 
        else if(nflag_emit(is).ge.1.and.nflag_emit(is).le.2) then
          if(injct(is).eq.0) cycle
          if(istep.lt.staendinj(1,is).or.istep.gt.staendinj(2,is)) cycle
!          ie = mod(istep,injct(is))
!          if(ie.ne.0.or.inpf(is).eq.0) cycle
!
!============================== from prescribed surfaces
          do iepl=nns,nne
            if(sstep.eq.1) nemit(iepl) = 0
            nemitw = &
           &  (int((istep + lstep + sstep*dtinj)*fluxf(iepl),8) &
           & - int((istep + lstep + (sstep - 1)*dtinj)*fluxf(iepl),8))*ustep
            ipc = ipc_ej(iepl)
            omni = omni_ej(iepl)
            algn = algn_ej(iepl)
            ground = peject(iepl)%grd
            lxlej = peject(iepl)%xl; lxuej = peject(iepl)%xu
            lylej = peject(iepl)%yl; lyuej = peject(iepl)%yu
            lzlej = peject(iepl)%zl; lzuej = peject(iepl)%zu
            radius = radi_ej(iepl)
            cntr(1:3) = cntr_ej(1:3,iepl)
            pcond = (ncond(ipc).gt.0)
            if(istep.lt.staendinjs(1,iepl).or.istep.gt.staendinjs(2,iepl)) cycle
            nranmax = int(size(dranu)/6)
!
!           ------------------- pos. & vel. assignment
            if(geom_ej(iepl).eq.1) then
        IPL1: do m=1,nemitw
                iran = mod(m-1,nranmax) + 1
                jran = iran*6 - 5
                if(iran.eq.1) then
                  call RANU0(dranu,min(nemitw*6,nranmax*6),icon)
                end if
!
                if(nemd(iepl).eq.+1) then
                  call emit_from_rectangular_surf &
                 &  (pinj(1)%vx,pinj(1)%vy,pinj(1)%vz, &
                 &   pinj(1)%x,pinj(1)%y,pinj(1)%z,+1, &
                 &   ground,lylej,lyuej,lzlej,lzuej,jran)
                else if(nemd(iepl).eq.-1) then
                  call emit_from_rectangular_surf &
                 &  (pinj(1)%vx,pinj(1)%vy,pinj(1)%vz, &
                 &   pinj(1)%x,pinj(1)%y,pinj(1)%z,-1, &
                 &   ground,lylej,lyuej,lzlej,lzuej,jran)
                else if(nemd(iepl).eq.+2) then
                  call emit_from_rectangular_surf &
                 &  (pinj(1)%vy,pinj(1)%vz,pinj(1)%vx, &
                 &   pinj(1)%y,pinj(1)%z,pinj(1)%x,+1, &
                 &   ground,lzlej,lzuej,lxlej,lxuej,jran)
                else if(nemd(iepl).eq.-2) then
                  call emit_from_rectangular_surf &
                 &  (pinj(1)%vy,pinj(1)%vz,pinj(1)%vx, &
                 &   pinj(1)%y,pinj(1)%z,pinj(1)%x,-1, &
                 &   ground,lzlej,lzuej,lxlej,lxuej,jran)
                else if(nemd(iepl).eq.+3) then
                  call emit_from_rectangular_surf &
                 &  (pinj(1)%vz,pinj(1)%vx,pinj(1)%vy, &
                 &   pinj(1)%z,pinj(1)%x,pinj(1)%y,+1, &
                 &   ground,lxlej,lxuej,lylej,lyuej,jran)
                else if(nemd(iepl).eq.-3) then
                  call emit_from_rectangular_surf &
                 &  (pinj(1)%vz,pinj(1)%vx,pinj(1)%vy, &
                 &   pinj(1)%z,pinj(1)%x,pinj(1)%y,-1, &
                 &   ground,lxlej,lxuej,lylej,lyuej,jran)
                end if
!
                pinj(1)%nid = &
               &  oh3_map_particle_to_subdomain(pinj(1)%x,pinj(1)%y,pinj(1)%z)
                if(pinj(1)%nid.ne.sdid(1)) cycle IPL1
                pinj(1)%spec = is
                pinj(1)%pid = 0
                pinj(1)%preside = -9
                call oh2_inject_particle(pinj(1))
                nemit(iepl) = nemit(iepl) + 1
                totalp(is,3) = totalp(is,3) + 1
!
                if(pcond) then
                  pinj(1)%preside = RPCB
                  call oh2_inject_particle(pinj(1))
                  totalp(is,3) = totalp(is,3) + 1
                else
                  pinj(1)%preside = RNCB
                  call oh2_inject_particle(pinj(1))
                  totalp(is,3) = totalp(is,3) + 1
                end if
              end do IPL1
            else if(geom_ej(iepl).eq.2) then

              if(abs(nemd(iepl)).eq.algn) then
          IPL2: do m=1,nemitw
                  iran = mod(m-1,nranmax) + 1
                  jran = iran*6 - 5
                  if(iran.eq.1) then
                    call RANU0(dranu,min(nemitw*6,nranmax*6),icon)
                  end if
!
                  if(nemd(iepl).eq.+1) then
                    call emit_from_circular_surf &
                   &  (pinj(1)%vx,pinj(1)%vy,pinj(1)%vz, &
                   &   pinj(1)%x,pinj(1)%y,pinj(1)%z,+1, &
                   &   ground,cntr(1),cntr(2),jran)
                  else if(nemd(iepl).eq.-1) then
                    call emit_from_circular_surf &
                   &  (pinj(1)%vx,pinj(1)%vy,pinj(1)%vz, &
                   &   pinj(1)%x,pinj(1)%y,pinj(1)%z,-1, &
                   &   ground,cntr(1),cntr(2),jran)
                  else if(nemd(iepl).eq.+2) then
                    call emit_from_circular_surf &
                   &  (pinj(1)%vy,pinj(1)%vz,pinj(1)%vx, &
                   &   pinj(1)%y,pinj(1)%z,pinj(1)%x,+1, &
                   &   ground,cntr(1),cntr(2),jran)
                  else if(nemd(iepl).eq.-2) then
                    call emit_from_circular_surf &
                   &  (pinj(1)%vy,pinj(1)%vz,pinj(1)%vx, &
                   &   pinj(1)%y,pinj(1)%z,pinj(1)%x,-1, &
                   &   ground,cntr(1),cntr(2),jran)
                  else if(nemd(iepl).eq.+3) then
                    call emit_from_circular_surf &
                   &  (pinj(1)%vz,pinj(1)%vx,pinj(1)%vy, &
                   &   pinj(1)%z,pinj(1)%x,pinj(1)%y,+1, &
                   &   ground,cntr(1),cntr(2),jran)
                  else if(nemd(iepl).eq.-3) then
                    call emit_from_circular_surf &
                   &  (pinj(1)%vz,pinj(1)%vx,pinj(1)%vy, &
                   &   pinj(1)%z,pinj(1)%x,pinj(1)%y,-1, &
                   &   ground,cntr(1),cntr(2),jran)
                  end if
!
                  pinj(1)%nid = &
                 &  oh3_map_particle_to_subdomain(pinj(1)%x,pinj(1)%y,pinj(1)%z)
                  if(pinj(1)%nid.ne.sdid(1)) cycle IPL2
                  pinj(1)%spec = is
                  pinj(1)%pid = 0
                  pinj(1)%preside = -9
                  call oh2_inject_particle(pinj(1))
                  nemit(iepl) = nemit(iepl) + 1
                  totalp(is,3) = totalp(is,3) + 1
!
                  if(pcond) then
                    pinj(1)%preside = RPCB
                    call oh2_inject_particle(pinj(1))
                    totalp(is,3) = totalp(is,3) + 1
                  else
                    pinj(1)%preside = RNCB
                    call oh2_inject_particle(pinj(1))
                    totalp(is,3) = totalp(is,3) + 1
                  end if
                end do IPL2
              else !if(abs(nemd(iepl)).ne.algn) then
          IPL3: do m=1,nemitw
                  iran = mod(m-1,nranmax) + 1
                  jran = iran*6 - 5
                  kran = iran*5 - 4
                  if(iran.eq.1) then
                    call RANU0(dranu,min(nemitw*6,nranmax*6),icon)
                  end if
!
                  if(algn.eq.1.and.omni.eq.0.and.radius.ge.0.1d0) then
                    call emit_from_cylindrical_surf1 &
                   &  (pinj(1)%vz,pinj(1)%vy,pinj(1)%vx, &
                   &   pinj(1)%z,pinj(1)%y,pinj(1)%x,psizy, &
                   &   cntr(2),cntr(1),lxlej,lxuej,jran)
                  else if(algn.eq.1.and.omni.eq.0.and.radius.lt.0.1d0) then
                    call emit_from_cylindrical_surf2 &
                   &  (pinj(1)%vz,pinj(1)%vy,pinj(1)%vx, &
                   &   pinj(1)%z,pinj(1)%y,pinj(1)%x,psizy, &
                   &   cntr(2),cntr(1),lxlej,lxuej,kran)
                  else if(algn.eq.1.and.omni.eq.1) then
                    call emit_from_cylindrical_surf3 &
                   &  (pinj(1)%vz,pinj(1)%vy,pinj(1)%vx, &
                   &   pinj(1)%z,pinj(1)%y,pinj(1)%x,psizy, &
                   &   cntr(2),cntr(1),lxlej,lxuej,jran)
                  else if(algn.eq.2.and.omni.eq.0.and.radius.ge.0.1d0) then
                    call emit_from_cylindrical_surf1 &
                   &  (pinj(1)%vz,pinj(1)%vx,pinj(1)%vy, &
                   &   pinj(1)%z,pinj(1)%x,pinj(1)%y,psizx, &
                   &   cntr(1),cntr(2),lylej,lyuej,jran)
                  else if(algn.eq.2.and.omni.eq.0.and.radius.lt.0.1d0) then
                    call emit_from_cylindrical_surf2 &
                   &  (pinj(1)%vz,pinj(1)%vx,pinj(1)%vy, &
                   &   pinj(1)%z,pinj(1)%x,pinj(1)%y,psizx, &
                   &   cntr(1),cntr(2),lylej,lyuej,kran)
                  else if(algn.eq.2.and.omni.eq.1) then
                    call emit_from_cylindrical_surf3 &
                   &  (pinj(1)%vz,pinj(1)%vx,pinj(1)%vy, &
                   &   pinj(1)%z,pinj(1)%x,pinj(1)%y,psizx, &
                   &   cntr(1),cntr(2),lylej,lyuej,jran)
                  else if(algn.eq.3.and.omni.eq.0.and.radius.ge.0.1d0) then
                    call emit_from_cylindrical_surf1 &
                   &  (pinj(1)%vx,pinj(1)%vy,pinj(1)%vz, &
                   &   pinj(1)%x,pinj(1)%y,pinj(1)%z,psixy, &
                   &   cntr(1),cntr(2),lzlej,lzuej,jran)
                  else if(algn.eq.3.and.omni.eq.0.and.radius.lt.0.1d0) then
                    call emit_from_cylindrical_surf2 &
                   &  (pinj(1)%vx,pinj(1)%vy,pinj(1)%vz, &
                   &   pinj(1)%x,pinj(1)%y,pinj(1)%z,psixy, &
                   &   cntr(1),cntr(2),lzlej,lzuej,kran)
                  else if(algn.eq.3.and.omni.eq.1) then
                    call emit_from_cylindrical_surf3 &
                   &  (pinj(1)%vx,pinj(1)%vy,pinj(1)%vz, &
                   &   pinj(1)%x,pinj(1)%y,pinj(1)%z,psixy, &
                   &   cntr(1),cntr(2),lzlej,lzuej,jran)
                  end if
!
                  pinj(1)%nid = &
                 &  oh3_map_particle_to_subdomain(pinj(1)%x,pinj(1)%y,pinj(1)%z)
                  if(pinj(1)%nid.ne.sdid(1)) cycle IPL3
                  pinj(1)%spec = is
                  pinj(1)%pid = 0
                  pinj(1)%preside = -9
                  call oh2_inject_particle(pinj(1))
                  nemit(iepl) = nemit(iepl) + 1
                  totalp(is,3) = totalp(is,3) + 1
!
                  if(pcond) then
                    pinj(1)%preside = RPCB
                    call oh2_inject_particle(pinj(1))
                    totalp(is,3) = totalp(is,3) + 1
                  else
                    pinj(1)%preside = RNCB
                    call oh2_inject_particle(pinj(1))
                    totalp(is,3) = totalp(is,3) + 1
                  end if
                end do IPL3
              end if
            else if(geom_ej(iepl).eq.3) then
        IPL4: do m=1,nemitw
                iran = mod(m-1,nranmax) + 1
                kran = iran*5 - 4
                if(iran.eq.1) then
                  call RANU0(dranu,min(nemitw*5,nranmax*5),icon)
                end if
!
                if(omni.eq.0) then
                  call emit_from_spherical_surf1 &
                 &  (pinj(1)%vx,pinj(1)%vy,pinj(1)%vz, &
                 &   pinj(1)%x,pinj(1)%y,pinj(1)%z, &
                 &   radius,cntr(1),cntr(2),cntr(3),kran)
                else !if(omni.ne.0) then
                  call emit_from_spherical_surf2 &
                 &  (pinj(1)%vx,pinj(1)%vy,pinj(1)%vz, &
                 &   pinj(1)%x,pinj(1)%y,pinj(1)%z, &
                 &   radius,cntr(1),cntr(2),cntr(3),kran)
                end if
!
                pinj(1)%nid = &
               &  oh3_map_particle_to_subdomain(pinj(1)%x,pinj(1)%y,pinj(1)%z)
                if(pinj(1)%nid.ne.sdid(1)) cycle IPL4
                pinj(1)%spec = is
                pinj(1)%pid = 0
                pinj(1)%preside = -9
                call oh2_inject_particle(pinj(1))
                nemit(iepl) = nemit(iepl) + 1
                totalp(is,3) = totalp(is,3) + 1
!
                if(pcond) then
                  pinj(1)%preside = RPCB
                  call oh2_inject_particle(pinj(1))
                  totalp(is,3) = totalp(is,3) + 1
                else
                  pinj(1)%preside = RNCB
                  call oh2_inject_particle(pinj(1))
                  totalp(is,3) = totalp(is,3) + 1
                end if
              end do IPL4
            end if
!
!           --------- calculate accumulated charge
!            if(ncond(ipc).gt.0) then
!              gcount(1)%chgacm(:,ipc) = gcount(1)%chgacm(:,ipc) &
!             &                             - q(is)*nemitw*sqdscaled(ipc)
!            end if
!            gcount(1)%outflux(ipc,is) = gcount(1)%outflux(ipc,is) &
!           &                               + nemitw
!            if(nemit(iepl).gt.0.and.sstep.eq.nscycinj) &
!           &  print*, "PS3c: myid,nemit",myid,nemit(iepl),sstep
          end do


!============================== from solid surface
          if(zl.le.zssurf.and.zu.gt.zssurf) then
            nemitw = (int((istep + lstep)*fluxub,8) &
           &       -  int((istep + lstep - injct(is))*fluxub,8))*ustep
            nranmax = int(size(dranu)/5)
            if(nemitw.gt.0) then
              do m=1,nemitw
                iran = mod(m-1,nranmax) + 1
                if(iran.eq.1) then
                  call RANU0(dranu,min(nemitw*5,nranmax*5),icon)
                end if
!               ----- pick up from velocity distribution
                addr0 = nprs + npr(is)*dranu(iran*5-4)
                addr1 = nprs + npr(is)*dranu(iran*5-3)
                addr2 = nprs + npr(is)*dranu(iran*5-2)
!               ----- assign position and velocity
                pinj(1)%vx = vtangf(addr0)
                pinj(1)%vy = vtangf(addr1)
                pinj(1)%vz = vnormf(addr2)
                pinj(1)%x = xl + (xu - xl)*dranu(iran*5-1)
                pinj(1)%y = yl + (yu - yl)*dranu(iran*5  )
                pinj(1)%z = zssurf
!
#if ssurf==2
                xsepa = pinj(1)%x - xdomec
                ysepa = pinj(1)%y - ydomec
                zsepa = pinj(1)%z - zdomec
!
                lbdome = 2.0d0*(pinj(1)%x*dray(1) + pinj(1)%y*dray(2) + &
               &                pinj(1)%z*dray(3) - dcdome)
                lcdome = xsepa*xsepa+ysepa*ysepa+zsepa*zsepa - rdomsq
                lddome = lbdome*lbdome - 4.0d0*laray*lcdome
                if(lddome.ge.0.0d0) then
                  tdome = 0.5d0*(-lbdome-sqrt(lddome))/laray
                  if(tdome.lt.0.0d0) then
                    pinj(1)%x = pinj(1)%x + tdome*dray(1)
                    pinj(1)%y = pinj(1)%y + tdome*dray(2)
                    pinj(1)%z = pinj(1)%z + tdome*dray(3)
                    v13 = (pinj(1)%x - xdomec)/rdome
                    v23 = (pinj(1)%y - ydomec)/rdome
                    v33 = (pinj(1)%z - zdomec)/rdome
                    v31 = -sqrt(1.0d0 - v33*v33)
                    iv31 = 1.0d0/v31
                    v12 = +v23*iv31
                    v22 = -v13*iv31
                    v11 = +v33*v22
                    v21 = -v33*v12
                    vxtmp = pinj(1)%vx; vytmp = pinj(1)%vy; vztmp = pinj(1)%vz
                    pinj(1)%vx = v11*vxtmp + v12*vytmp + v13*vztmp
                    pinj(1)%vy = v21*vxtmp + v22*vytmp + v23*vztmp
                    pinj(1)%vz = v31*vxtmp             + v33*vztmp
                  end if
                end if
!
                xsepa = (pinj(1)%x - xbowlc)*maxval(rbowl(1:3))/rbowl(1)
                ysepa = (pinj(1)%y - ybowlc)*maxval(rbowl(1:3))/rbowl(2)
                zsepa = (pinj(1)%z - zbowlc)*maxval(rbowl(1:3))/rbowl(3)
!
                if(xsepa*xsepa+ysepa*ysepa+zsepa*zsepa.lt.maxval(rbwlsq(1:3))) then
                  lcbowl = xsepa*xsepa+ysepa*ysepa+zsepa*zsepa - maxval(rbwlsq(1:3))
                  lbbowl = 2.0d0*(pinj(1)%x*dray(1) + pinj(1)%y*dray(2) + &
                 &                pinj(1)%z*dray(3) - dcbowl)
                  ldbowl = lbbowl*lbbowl - 4.0d0*laray*lcbowl
                  if(ldbowl.ge.0.0d0) then
                    tbowl = 0.5d0*(-lbbowl+sqrt(ldbowl))/laray
                    xsepa = (xsepa + tbowl*dray(1))*rbowl(1)/maxval(rbowl(1:3))
                    ysepa = (ysepa + tbowl*dray(2))*rbowl(2)/maxval(rbowl(1:3))
                    zsepa = (zsepa + tbowl*dray(3))*rbowl(3)/maxval(rbowl(1:3))
                    pinj(1)%x = xbowlc + xsepa
                    pinj(1)%y = ybowlc + ysepa
                    pinj(1)%z = zbowlc + zsepa
                    xysepa = sqrt(xsepa*xsepa + ysepa*ysepa)
                    tan1 = -rbwlsq(1)*rbwlsqi(3)*zsepa*xysepa
                    cos1 = sqrt(1.0d0/(tan1*tan1+1.0d0))
                    v13 = -cos1*xsepa/xysepa
                    v23 = -cos1*ysepa/xysepa
                    v33 = sqrt(1-cos1*cos1)
                    v31 = -sqrt(1.0d0 - v33*v33)
                    iv31 = 1.0d0/v31
                    v12 = +v23*iv31
                    v22 = -v13*iv31
                    v11 = +v33*v22
                    v21 = -v33*v12
                    vxtmp = pinj(1)%vx; vytmp = pinj(1)%vy; vztmp = pinj(1)%vz
                    pinj(1)%vx = v11*vxtmp + v12*vytmp + v13*vztmp
                    pinj(1)%vy = v21*vxtmp + v22*vytmp + v23*vztmp
                    pinj(1)%vz = v31*vxtmp             + v33*vztmp
                  end if
                end if

!                lcbowl = xsepa2*xsepa2+ysepa2*ysepa2+zsepa2*zsepa2 - rbwlsq

!                if(lcbowl.lt.0.0d0) then
!                  lbbowl = 2.0d0*(pinj(1)%x*dray(1) + pinj(1)%y*dray(2) + &
!                 &                pinj(1)%z*dray(3) - dcbowl)
!                  ldbowl = lbbowl*lbbowl - 4.0d0*laray*lcbowl
!                  if(ldbowl.ge.0.0d0) then
!                    tbowl = 0.5d0*(-lbbowl+sqrt(ldbowl))/laray
!                    pinj(1)%x = pinj(1)%x + tbowl*dray(1)
!                    pinj(1)%y = pinj(1)%y + tbowl*dray(2)
!                    pinj(1)%z = pinj(1)%z + tbowl*dray(3)
!                    v13 = (xbowlc - pinj(1)%x)/rbowl
!                    v23 = (ybowlc - pinj(1)%y)/rbowl
!                    v33 = (zbowlc - pinj(1)%z)/rbowl
!                    v31 = -sqrt(1.0d0 - v33*v33)
!                    iv31 = 1.0d0/v31
!                    v12 = +v23*iv31
!                    v22 = -v13*iv31
!                    v11 = +v33*v22
!                    v21 = -v33*v12
!                    vxtmp = pinj(1)%vx; vytmp = pinj(1)%vy; vztmp = pinj(1)%vz
!                    pinj(1)%vx = v11*vxtmp + v12*vytmp + v13*vztmp
!                    pinj(1)%vy = v21*vxtmp + v22*vytmp + v23*vztmp
!                    pinj(1)%vz = v31*vxtmp             + v33*vztmp
!                  end if
!                end if
#elif ssurf==3
                if(zlrechole(1).le.zssurf.and.zssurf.le.zurechole(1).and. &
               &   xlrechole(1).le.pinj(1)%x.and.pinj(1)%x.le.xurechole(1).and. &
               &   ylrechole(1).le.pinj(1)%y.and.pinj(1)%y.le.yurechole(1)) then
                  cycle
                end if
#endif
!
                pinj(1)%nid = &
               &  oh3_map_particle_to_subdomain(pinj(1)%x,pinj(1)%y,pinj(1)%z)
                pinj(1)%spec = is
                pinj(1)%pid = 0
                pinj(1)%preside = -9
                call oh2_inject_particle(pinj(1))
                totalp(is,3) = totalp(is,3) + 1
!
                pinj(1)%preside = RNCB
                call oh2_inject_particle(pinj(1))
                totalp(is,3) = totalp(is,3) + 1
              end do
            end if
          end if

!-------------------- from hole basement


!-------------------- from hole flank


        end if
!============================== end of plane loop

      end do ISL1
!****************************** end of species loop 1

      pbase(4) = pbase(3) + sum(totalp(1:nspec,3))
      pbase(5) = pbase(4)


    return


    contains


    subroutine emit_from_rectangular_surf &
   &  (vnorm,vtang1,vtang2,rnorm,rtang1,rtang2,ud,surfp,hl1,hu1,hl2,hu2,ir)
!-------------------- args & vars
      integer(kind=4),intent(in) :: ud,ir
      real(kind=8),intent(in)    :: surfp,hl1,hu1,hl2,hu2
      real(kind=8),intent(out)   :: vnorm,vtang1,vtang2,rnorm,rtang1,rtang2
      integer(kind=4) :: addr0,addr1,addr2

!-------------------- pick up from velocity distribution
      addr0 = nprs + npr(is)*dranu(ir  )
      addr1 = nprs + npr(is)*dranu(ir+1)
      addr2 = nprs + npr(is)*dranu(ir+2)

!-------------------- assign position and velocity
      vnorm  = ud*vnormf(addr0)
      vtang1 = vtangf(addr1)
      vtang2 = vtangf(addr2)
      rnorm  = surfp
      rtang1 = hl1 + (hu1 - hl1)*dranu(ir+4)
      rtang2 = hl2 + (hu2 - hl2)*dranu(ir+5)

      return
    end subroutine emit_from_rectangular_surf


    subroutine emit_from_rectangular_surf_with_hole &
   &  (vnorm,vtang1,vtang2,rnorm,rtang1,rtang2,ud,surfp,hl1,hu1,hl2,hu2, &
   &   holel1,holeu1,holel2,holeu2,ir)
!-------------------- args & vars
      integer(kind=4),intent(in) :: ud,ir
      real(kind=8),intent(in)    :: surfp,hl1,hu1,hl2,hu2
      real(kind=8),intent(in)    :: holel1,holeu1,holel2,holeu2
      real(kind=8),intent(out)   :: vnorm,vtang1,vtang2,rnorm,rtang1,rtang2
      integer(kind=4) :: addr0,addr1,addr2, icon
      real(kind=8) :: dran(2)

!-------------------- pick up from velocity distribution
      addr0 = nprs + npr(is)*dranu(ir  )
      addr1 = nprs + npr(is)*dranu(ir+1)
      addr2 = nprs + npr(is)*dranu(ir+2)

!-------------------- assign position and velocity
      vnorm  = ud*vnormf(addr0)
      vtang1 = vtangf(addr1)
      vtang2 = vtangf(addr2)
      rnorm  = surfp
      do!while(.true.)
        call RANU0(dran,2,icon)
        rtang1 = hl1 + (hu1 - hl1)*dran(1)
        rtang2 = hl2 + (hu2 - hl2)*dran(2)
        if(rtang1.lt.holel1.or.rtang1.gt.holeu1.or. &
       &   rtang2.lt.holel2.or.rtang2.gt.holeu2) then
          rtang1 = rtang1
          rtang2 = rtang2
          exit
        end if
      end do

      return
    end subroutine emit_from_rectangular_surf_with_hole


    subroutine emit_from_circular_surf &
   &  (vnorm,vtang1,vtang2,rnorm,rtang1,rtang2,ud,surfp,axis1,axis2,ir)
!-------------------- args & vars
      integer(kind=4),intent(in) :: ud,ir
      real(kind=8),intent(in)    :: surfp,axis1,axis2
      real(kind=8),intent(out)   :: vnorm,vtang1,vtang2,rnorm,rtang1,rtang2
      integer(kind=4) :: addr0,addr1,addr2,j3
      real(kind=8) :: radial, angular

!-------------------- pick up from velocity distribution
      addr0 = nprs + npr(is)*dranu(ir  )
      addr1 = nprs + npr(is)*dranu(ir+1)
      addr2 = nprs + npr(is)*dranu(ir+2)
      j3 = nprs + npr(is)*dranu(ir+3)

!-------------------- assign position and velocity
      vnorm  = ud*vnormf(addr0)
      vtang1 = vtangf(addr1)
      vtang2 = vtangf(addr2)
      radial = radius*runiform(j3)
      angular = pi2*dranu(ir+4)
      rnorm  = surfp
      rtang1 = axis1 + radial*dcos(angular)
      rtang2 = axis2 + radial*dsin(angular)

      return
    end subroutine emit_from_circular_surf


    subroutine emit_from_cylindrical_surf1 &
   &  (vnorm,vtang1,vtang2,rnorm,rtang1,rtang2,psi,axis1,axis2,hl,hu,ir)
!-------------------- args & vars
      integer(kind=4),intent(in) :: ir
      real(kind=8),intent(in)    :: psi,axis1,axis2,hl,hu
      real(kind=8),intent(out)   :: vnorm,vtang1,vtang2,rnorm,rtang1,rtang2
      integer(kind=4) :: addr0,addr1,addr2,j3

!-------------------- pick up from velocity distribution
      addr0 = nprs + npr(is)*dranu(ir  )
      addr1 = nprs + npr(is)*dranu(ir+1)
      addr2 = nprs + npr(is)*dranu(ir+2)
      j3 = nprs + npr(is)*dranu(ir+3)

!-------------------- assign position and velocity
      vnorm  = vnormc(addr0)*dcos(psicos(addr2)+psi) &
     &       - vtangc(addr1)*dsin(psicos(addr2)+psi)
      vtang1 = vnormc(addr0)*dsin(psicos(addr2)+psi) &
     &       + vtangc(addr1)*dcos(psicos(addr2)+psi)
      vtang2 = vtangf(j3)
      rnorm  = axis1 &
     &       + radius*dcos(psicos(addr2)+psi)
      rtang1 = axis2 &
     &       + radius*dsin(psicos(addr2)+psi)
      rtang2 = hl + (hu - hl)*dranu(ir+5)
     
      return
    end subroutine emit_from_cylindrical_surf1


    subroutine emit_from_cylindrical_surf2 &
   &  (vnorm,vtang1,vtang2,rnorm,rtang1,rtang2,psi,axis1,axis2,hl,hu,ir)
!-------------------- args & vars
      integer(kind=4),intent(in) :: ir
      real(kind=8),intent(in)    :: psi,axis1,axis2,hl,hu
      real(kind=8),intent(out)   :: vnorm,vtang1,vtang2,rnorm,rtang1,rtang2
      integer(kind=4) :: addr0,addr1,addr2
      real(kind=8) :: v2normi

!-------------------- pick up from velocity distribution
      addr0 = nprs + npr(is)*dranu(ir  )
      addr1 = nprs + npr(is)*dranu(ir+1)
      addr2 = nprs + npr(is)*dranu(ir+2)

!-------------------- assign position and velocity
      vnorm  = vnormc(addr0)*dcos(psi) &
     &       - vtangc(addr1)*dsin(psi)
      vtang1 = vnormc(addr0)*dsin(psi) &
     &       + vtangc(addr1)*dcos(psi)
      vtang2 = vtangf(addr2)
      v2normi = radius/dsqrt(vnorm*vnorm + vtang1*vtang1)
      rnorm  = axis1 &
     &       + vnorm*v2normi
      rtang1 = axis2 &
     &       + vtang1*v2normi
      rtang2 = hl + (hu - hl)*dranu(ir+4)
     
      return
    end subroutine emit_from_cylindrical_surf2


    subroutine emit_from_cylindrical_surf3 &
   &  (vnorm,vtang1,vtang2,rnorm,rtang1,rtang2,psi,axis1,axis2,hl,hu,ir)
!-------------------- args & vars
      integer(kind=4),intent(in) :: ir
      real(kind=8),intent(in)    :: psi,axis1,axis2,hl,hu
      real(kind=8),intent(out)   :: vnorm,vtang1,vtang2,rnorm,rtang1,rtang2
      integer(kind=4) :: addr0,addr1,addr2,j3
      real(kind=8) :: eazimuth

!-------------------- pick up from velocity distribution
      addr0 = nprs + npr(is)*dranu(ir  )
      addr1 = nprs + npr(is)*dranu(ir+1)
      eazimuth = pi2*dranu(ir+2)
      j3 = nprs + npr(is)*dranu(ir+3)

!-------------------- assign position and velocity
      vnorm  = vnormc(addr0)*dcos(eazimuth) &
     &       - vtangc(addr1)*dsin(eazimuth)
      vtang1 = vnormc(addr0)*dsin(eazimuth) &
     &       + vtangc(addr1)*dcos(eazimuth)
      vtang2 = vtangf(j3)
      rnorm  = axis1 &
     &       + radius*dcos(eazimuth)
      rtang1 = axis2 &
     &       + radius*dsin(eazimuth)
      rtang2 = hl + (hu - hl)*dranu(ir+5)

      return
    end subroutine emit_from_cylindrical_surf3


    subroutine emit_from_spherical_surf1 &
   &  (vx,vy,vz,rx,ry,rz,rad,x0,y0,z0,ir)
!-------------------- args & vars
      integer(kind=4),intent(in) :: ir
      real(kind=8),intent(in)    :: rad, x0,y0,z0
      real(kind=8),intent(out)   :: vx,vy,vz,rx,ry,rz
      integer(kind=4) :: addr0,addr1,addr2
      real(kind=8) :: rxy,csrz,snrz,csrxy,snrxy
      real(kind=8) :: ttr11,ttr12,ttr13,ttr21,ttr22,ttr23,ttr31,ttr32,ttr33
      real(kind=8) :: xx,yy,zz,vrx,vry,vrz

!-------------------- pick up from velocity distribution
      addr0 = nprs + npr(is)*dranu(ir  )
      addr1 = nprs + npr(is)*dranu(ir+1)
      addr2 = nprs + npr(is)*dranu(ir+2)

!-------------------- assign position and velocity
      snrz = dranu(ir+3)
      csrz = sqrt(1.0d0 - snrz)
      snrz = sqrt(snrz)
      rxy = 2.0d0*pi*dranu(ir+4)
      snrxy = sin(rxy)
      csrxy = cos(rxy)

      ttr11 = csrz*csrxy
      ttr12 = -snrxy
      ttr13 = snrz*csrxy
      ttr21 = csrz*snrxy
      ttr22 = csrxy
      ttr23 = snrz*snrxy
      ttr31 = -snrz
      ttr32 = 0.0d0
      ttr33 = csrz

      xx = snrz*csrxy
      yy = snrz*snrxy
      zz = csrz

      rx = (xx*tt11 + yy*tt12 + zz*tt13)*rad + x0
      ry = (xx*tt21 + yy*tt22 + zz*tt23)*rad + y0
      rz = (xx*tt31 + yy*tt32 + zz*tt33)*rad + z0

      vx = vtangf(addr0)
      vy = vtangf(addr1)
      vz = vnormf(addr2)

      vrx = vx*ttr11 + vy*ttr12 + vz*ttr13
      vry = vx*ttr21 + vy*ttr22 + vz*ttr23
      vrz = vx*ttr31 + vy*ttr32 + vz*ttr33

      vx = vrx*tt11 + vry*tt12 + vrz*tt13
      vy = vrx*tt21 + vry*tt22 + vrz*tt23
      vz = vrx*tt31 + vry*tt32 + vrz*tt33
     
      return
    end subroutine emit_from_spherical_surf1


    subroutine emit_from_spherical_surf2 &
   &  (vx,vy,vz,rx,ry,rz,rad,x0,y0,z0,ir)
!-------------------- args & vars
      integer(kind=4),intent(in) :: ir
      real(kind=8),intent(in)    :: rad, x0,y0,z0
      real(kind=8),intent(out)   :: vx,vy,vz,rx,ry,rz
      integer(kind=4) :: addr0,addr1,addr2
      real(kind=8) :: rxy,csrz,snrz,csrxy,snrxy
      real(kind=8) :: ttr11,ttr12,ttr13,ttr21,ttr22,ttr23,ttr31,ttr32,ttr33
      real(kind=8) :: xx,yy,zz,vrx,vry,vrz

!-------------------- pick up from velocity distribution
      addr0 = nprs + npr(is)*dranu(ir  )
      addr1 = nprs + npr(is)*dranu(ir+1)
      addr2 = nprs + npr(is)*dranu(ir+2)

!-------------------- assign position and velocity
      csrz = 1.0d0-2.0d0*dranu(ir+3)
      snrz = sqrt(1.0d0-csrz*csrz)
      rxy = 2.0d0*pi*dranu(ir+4)
      snrxy = sin(rxy)
      csrxy = cos(rxy)

      ttr11 = csrz*csrxy
      ttr12 = -snrxy
      ttr13 = snrz*csrxy
      ttr21 = csrz*snrxy
      ttr22 = csrxy
      ttr23 = snrz*snrxy
      ttr31 = -snrz
      ttr32 = 0.0d0
      ttr33 = csrz

      rx = snrz*csrxy*rad + x0
      ry = snrz*snrxy*rad + y0
      rz = csrz*rad + z0

      vrx = vtangf(addr0)
      vry = vtangf(addr1)
      vrz = vnormf(addr2)

      vx = vrx*ttr11 + vry*ttr12 + vrz*ttr13
      vy = vrx*ttr21 + vry*ttr22 + vrz*ttr23
      vz = vrx*ttr31 + vry*ttr32 + vrz*ttr33
     
      return
    end subroutine emit_from_spherical_surf2


  end subroutine psuper3
