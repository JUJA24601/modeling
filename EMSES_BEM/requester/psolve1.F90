#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine psolve1(pstq,ustep,func)
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   P S O L V E
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .      this subroutine gives a current distribution to     .
!   .      grid points from particle locations by the          .
!   .      charge conservation method based on the "zigzag"    .
!   .      scheme. [Umeda, T., 2003]                           .
!   ............................................................
!
!
!-------------------- parameter and common block
  use oh_type
  use paramt
  use allcom
  implicit none
!
  integer(kind=4),parameter :: n_emission_pattern=2
  integer(kind=8) :: m,mm, ns,ne, nee
  integer(kind=8) :: nprs,npre
  integer(kind=8) :: isfluxt(inpc,ispec)
  integer(kind=8) :: influxsum
  integer(kind=4) :: i,j,k
  integer(kind=4) :: ib,jb,kb, i1b,j1b,k1b
  integer(kind=4) :: ia,ja,ka, i1a,j1a,k1a
  integer(kind=4) :: nns,nne
  integer(kind=4) :: is,iis, iss, ieb, iepl, ustep, ibdy, ipc,jpc, itch
  integer(kind=4) :: xl,xu, yl,yu, zl,zu
  integer(kind=4) :: ngx,ngy,ngz
  integer(kind=4) :: nb(3)
  integer(kind=4) :: func
  integer(kind=4) :: pstq, ps, psd, rid, nud
  integer(kind=4) :: nyield(n_emission_pattern)
  integer(kind=4) :: icon,iran
  integer(kind=4) :: addr0,addr1,addr2
!  integer(kind=4) :: oh3_map_region_to_node
  real(kind=8) :: dxl,dxu, dyl,dyu, dzl,dzu
  real(kind=8) :: dngx,dngy,dngz
  real(kind=8) :: xlc(inpc),ylc(inpc),zlc(inpc)
  real(kind=8) :: xuc(inpc),yuc(inpc),zuc(inpc)
  real(kind=8) :: xlocalb,ylocalb,zlocalb
  real(kind=8) :: xb,yb,zb
  real(kind=8) :: xlocala,ylocala,zlocala
  real(kind=8) :: xa,ya,za
  real(kind=8) :: xr,yr,zr
  real(kind=8) :: xd1,yd1,zd1, xd2,yd2,zd2
  real(kind=8) :: x1,y1,z1, z2, xy1,xz1,yz1, xz2,yz2
  real(kind=8) :: tslx,tsly,tslz, gustep,rustep, qs,qmp
  real(kind=8) :: v1,v2,v3,v4,v5,v6,v7,v8
  real(kind=8) :: eex,eey,eez, bbx,bby,bbz, boris, vxt,vyt,vzt, vxyz
  real(kind=8) :: vvx1,vvy1,vvz1, vvx2,vvy2,vvz2
  real(kind=8) :: vx1w1,vx1w2,vx1w3,vx1w4, vx2w1,vx2w2,vx2w3,vx2w4
  real(kind=8) :: vy1w1,vy1w2,vy1w3,vy1w4, vy2w1,vy2w2,vy2w3,vy2w4
  real(kind=8) :: vz1w1,vz1w2,vz1w3,vz1w4, vz2w1,vz2w2,vz2w3,vz2w4
  real(kind=8) :: perg, pergfrac, pemaxdble(ispec),pemaxinvh(ispec), costhi
  real(kind=8) :: xpast,ypast,zpast
  real(kind=8) :: yield,weightr,dltamax1114
!
  real(kind=8) :: displx,disply,displz,displ3i
!
  real(kind=8) :: disp1,disp2,disp3
  real(kind=8) :: radsq(inpc),radsqi(inpc)
  real(kind=8) :: rbwlsq(3),rbwlsqi(3), rdomsq
  real(kind=8) :: xpast0,ypast0,zpast0, xpast1,ypast1,zpast1
  real(kind=8) :: xmove,ymove,zmove, xmovei,ymovei,zmovei
  real(kind=8) :: xsepa,ysepa,zsepa, termA,termB
  real(kind=8) :: tfrac, tfrac1,tfrac2,tfrac3,tfrac4,tfrac5,tfrac6
  real(kind=8) :: xpast1a,ypast1a,xpast1b,ypast1b
  real(kind=8) :: xpast2a,ypast2a,xpast2b,ypast2b
  real(kind=8) :: nsigmadt, colfp, uppb
  real(kind=8) :: ewx,ewy,tew

  real(kind=8) :: keprt, ikeprt, keprt_xy, dotprod, pcosth, ipcosth
  real(kind=8) :: rkeprt, irkeprt, irkeprt035, rkeprt135
  real(kind=8) :: vxtemp,vytemp,vztemp

  logical :: npbnd_x_periodic,npbnd_y_periodic,npbnd_z_periodic
  logical :: npbnd_x_reflect,npbnd_y_reflect,npbnd_z_reflect
  logical :: npbnd_x_erase,npbnd_y_erase,npbnd_z_erase
  logical :: singridx,singridy,singridz
  logical :: pout_xl,pout_xu, pout_yl,pout_yu, pout_zl,pout_zu
  logical :: flag_collr(2)
  logical :: flag_tch
  logical :: flag_tq
  logical :: pcond(inpc)
  logical :: sec_emission


!-------------------- 
      ps = modulo(pstq-1,2) + 1


!-------------------- 
      xl = sdoms(1,1,sdid(ps)+1); xu = sdoms(2,1,sdid(ps)+1)
      yl = sdoms(1,2,sdid(ps)+1); yu = sdoms(2,2,sdid(ps)+1)
      zl = sdoms(1,3,sdid(ps)+1); zu = sdoms(2,3,sdid(ps)+1)
      ngx = xu - xl; ngy = yu - yl; ngz = zu - zl
      dxl = xl; dxu = xu; dyl = yl; dyu = yu; dzl = zl; dzu = zu
      dngx = ngx; dngy = ngy; dngz = ngz
      pcond(1:npc) = (ncond(1:npc).gt.0)
      xlc(1:npc) = xlpc(1:npc) - dxl; xuc(1:npc) = xupc(1:npc) - dxl
      ylc(1:npc) = ylpc(1:npc) - dyl; yuc(1:npc) = yupc(1:npc) - dyl
      zlc(1:npc) = zlpc(1:npc) - dzl; zuc(1:npc) = zupc(1:npc) - dzl


!--------------------
      pemaxdble(1:nspec) = 2.0d0*pemax(1:nspec)
      where(pemaxdble(1:nspec).ne.0.0d0)
        pemaxinvh(1:nspec) = 1.0d0/pemaxdble(1:nspec)
      elsewhere
        pemaxinvh(1:nspec) = 1.0d10
      endwhere


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
        rholeinv = 1.0d0/rhole
        xlwall = xholec - rhole
        xuwall = xholec + rhole
      else
        rholesq = 0.0d0
        rholeinv = 1.0d0/epsilon(1.0d0)
        xlwall = -abs(dmiss)
        xuwall = +abs(dmiss)
      end if


!-------------------- for three dimensional system
      do ieb=0,nebfld-1
        psd = ieb*2 + ps
!       ------------- zero clear aj
        aj(:,:,:,:,psd) = 0.0d0
      end do
!     --------------- zero clear ajdg
      if(func.eq.1) &
     &   ajdg(:,:,:,:,ps) = 0.0d0


!============================== species loop
      nee = pbase(pstq)
      rustep = 1.0d0/ustep
      tslx = 2.0d0*slx
      tsly = 2.0d0*sly
      tslz = 2.0d0*slz
      singridx = (nx.eq.1)
      singridy = (ny.eq.1)
      singridz = (nz.eq.1)
      flag_tch = (ntch.gt.0)
      flag_tq = (pstq.ge.3)
ISL1: do is=1,nspec
        ns = nee + 1
        ne = nee + totalp(is,pstq)
        nee = ne
        npre = 0
        do iis=1,isse(is)
          nprs = npre + 1
          npre = npre + npr(iis)
        end do
!       ------------- 
        npbnd_x_periodic = (npbnd(1,is).eq.0)
        npbnd_y_periodic = (npbnd(2,is).eq.0)
        npbnd_z_periodic = (npbnd(3,is).eq.0)
        npbnd_x_reflect = (npbnd(1,is).eq.1)
        npbnd_y_reflect = (npbnd(2,is).eq.1)
        npbnd_z_reflect = (npbnd(3,is).eq.1)
        npbnd_x_erase = &
       &  (emmode.and.npbnd(1,is).eq.2.and.mtd_vbnd(1).eq.2)
        npbnd_y_erase = &
       &  (emmode.and.npbnd(2,is).eq.2.and.mtd_vbnd(2).eq.2)
        npbnd_z_erase = &
       &  (emmode.and.npbnd(3,is).eq.2.and.mtd_vbnd(3).eq.2)
        sec_emission = &
       &  (isse(is).ge.1.and.isse(is).le.nspec)
        if(sec_emission) sec_emission = &
       &  (abs(q(isse(is))).gt.0.0d0.and. &
           deltaemax(is).gt.0.0d0.and.pemax(is).gt.0.0d0)
!       -------------
        if(sec_emission) then
          weightr = abs(q(is)/q(isse(is)))
          dltamax1114 = 1.114d0*weightr*deltaemax(is)
        end if
!       ------------- 
        if(ewmodel.eq.1.or.ewmodel.eq.2) then
          tew = t - dt*nretard
          ewx = Ew(1,is,1)*cos(omegaw(1)*tew)
          ewy = Ew(2,is,1)*sin(omegaw(1)*tew)
        else
          ewx = 0.0d0
          ewy = 0.0d0
        end if
!       -------------
        flag_collr(1) = (mfpath(is).gt.1.0d0)
        flag_collr(2) = (mfpath(is).lt.-1.0d0)
        if(flag_collr(1)) then
          nsigmadt = ustep/mfpath(is)
        else if(flag_collr(2)) then
          nsigmadt = ustep/mfpath(is)
        else
          nsigmadt = 0.0d0
        end if
        uppb = 1.0d0 - exp(-vxyzmax(is)*abs(nsigmadt))
!       ------------- zero clear work array
        wrk(TJX:TJZ,:,:,:) = 0.0d0
!       ------------- inner loop
   ML1: do m=ns,ne
          if(pbuf(m)%preside.le.-9) &
         &  pbuf(m)%preside = pbuf(m)%preside + 10
!
          if(pbuf(m)%nid.eq.-1.or.pbuf(m)%preside.lt.0) then
            cycle ML1
          else if(pbuf(m)%preside.eq.1) then
            call RANU0(dranu,1,icon)
            gustep = ustep*dranu(1)
!            qmp = qm(is)*ustep*max(0.0d0,dranu(1)-0.5d0)*0.5d0
            qmp = 0.0d0
          else
            gustep = ustep
            qmp = qm(is)*gustep*0.5d0
          end if
!
!          if(is.eq.4) then
!            print*, "A:pbuf(m)@is=4", pbuf(m)
!          end if
!
          xlocalb = pbuf(m)%x - dxl
          ylocalb = pbuf(m)%y - dyl
          zlocalb = pbuf(m)%z - dzl
!
          ib = floor(xlocalb)
          jb = floor(ylocalb)
          kb = floor(zlocalb)
!
          xb = ib
          yb = jb
          zb = kb
!
          i1b = ib + 1
          j1b = jb + 1
          k1b = kb + 1
!
          x1 = xlocalb - xb
          y1 = ylocalb - yb
          z1 = zlocalb - zb
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
          eex = wrk(EX,ib ,jb ,k1b)*v1 + wrk(EX,i1b,jb ,k1b)*v2 &
         &    + wrk(EX,i1b,j1b,k1b)*v3 + wrk(EX,ib ,j1b,k1b)*v4 &
         &    + wrk(EX,ib ,jb ,kb )*v5 + wrk(EX,i1b,jb ,kb )*v6 &
         &    + wrk(EX,i1b,j1b,kb )*v7 + wrk(EX,ib ,j1b,kb )*v8
          eey = wrk(EY,ib ,jb ,k1b)*v1 + wrk(EY,i1b,jb ,k1b)*v2 &
         &    + wrk(EY,i1b,j1b,k1b)*v3 + wrk(EY,ib ,j1b,k1b)*v4 &
         &    + wrk(EY,ib ,jb ,kb )*v5 + wrk(EY,i1b,jb ,kb )*v6 &
         &    + wrk(EY,i1b,j1b,kb )*v7 + wrk(EY,ib ,j1b,kb )*v8
          eez = wrk(EZ,ib ,jb ,k1b)*v1 + wrk(EZ,i1b,jb ,k1b)*v2 &
         &    + wrk(EZ,i1b,j1b,k1b)*v3 + wrk(EZ,ib ,j1b,k1b)*v4 &
         &    + wrk(EZ,ib ,jb ,kb )*v5 + wrk(EZ,i1b,jb ,kb )*v6 &
         &    + wrk(EZ,i1b,j1b,kb )*v7 + wrk(EZ,ib ,j1b,kb )*v8
          bbx = wrk(BX,ib ,jb ,k1b)*v1 + wrk(BX,i1b,jb ,k1b)*v2 &
         &    + wrk(BX,i1b,j1b,k1b)*v3 + wrk(BX,ib ,j1b,k1b)*v4 &
         &    + wrk(BX,ib ,jb ,kb )*v5 + wrk(BX,i1b,jb ,kb )*v6 &
         &    + wrk(BX,i1b,j1b,kb )*v7 + wrk(BX,ib ,j1b,kb )*v8
          bby = wrk(BY,ib ,jb ,k1b)*v1 + wrk(BY,i1b,jb ,k1b)*v2 &
         &    + wrk(BY,i1b,j1b,k1b)*v3 + wrk(BY,ib ,j1b,k1b)*v4 &
         &    + wrk(BY,ib ,jb ,kb )*v5 + wrk(BY,i1b,jb ,kb )*v6 &
         &    + wrk(BY,i1b,j1b,kb )*v7 + wrk(BY,ib ,j1b,kb )*v8
          bbz = wrk(BZ,ib ,jb ,k1b)*v1 + wrk(BZ,i1b,jb ,k1b)*v2 &
         &    + wrk(BZ,i1b,j1b,k1b)*v3 + wrk(BZ,ib ,j1b,k1b)*v4 &
         &    + wrk(BZ,ib ,jb ,kb )*v5 + wrk(BZ,i1b,jb ,kb )*v6 &
         &    + wrk(BZ,i1b,j1b,kb )*v7 + wrk(BZ,ib ,j1b,kb )*v8
!
!         ----------- 
          if(flag_tch) then
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
          end if
!
!         ----------- charge-to-mass ratio is taken into consideration
          eex = (eex + ewx)*qmp
          eey = (eey + ewy)*qmp
          eez = eez*qmp
          bbx = bbx*qmp
          bby = bby*qmp
          bbz = bbz*qmp
!
!         ----------- update particle velocities (Buneman-Boris method)
          boris = 2.0d0/(1.0d0 + bbx*bbx + bby*bby + bbz*bbz)
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

!         ----------- velocity modification associated with collisions
          if((flag_collr(1)).and.pbuf(m)%preside.ne.1.and. &
         &   (pbuf(m)%x.le.xlcol(1).or.pbuf(m)%x.ge.xucol(1).or. &
         &    pbuf(m)%y.le.ylcol(1).or.pbuf(m)%y.ge.yucol(1).or. &
         &    pbuf(m)%z.le.zlcol(1).or.pbuf(m)%z.ge.zucol(1))) then
            vxt = pbuf(m)%vx - vdtx(is)
            vyt = pbuf(m)%vy - vdty(is)
            vzt = pbuf(m)%vz - vdtz(is)
            vxyz = sqrt(vxt*vxt + vyt*vyt + vzt*vzt)
            call RANU0(dranu,1,icon)
            colfp = colf(1,ib ,jb ,k1b,ps)*v1 + colf(1,i1b,jb ,k1b,ps)*v2 &
           &      + colf(1,i1b,j1b,k1b,ps)*v3 + colf(1,ib ,j1b,k1b,ps)*v4 &
           &      + colf(1,ib ,jb ,kb ,ps)*v5 + colf(1,i1b,jb ,kb ,ps)*v6 &
           &      + colf(1,i1b,j1b,kb ,ps)*v7 + colf(1,ib ,j1b,kb ,ps)*v8
            if(exp(-vxyz*nsigmadt).lt.dranu(1)*colfp) then
              vxyz = 0.0d0
              call RANU0(dranu,1,icon)
              do while(exp(-vxyz*nsigmadt).ge.1.0d0-dranu(1)*uppb)
                call RANN0(0.0d0,peth(is),dranu,1,icon)
                vxt = dranu(1)
                call RANN0(0.0d0,peth(is),dranu,1,icon)
                vyt = dranu(1)
                call RANN0(0.0d0,path(is),dranu,1,icon)
                vzt = dranu(1)
                vxyz = sqrt(vxt*vxt + vyt*vyt + vzt*vzt)
                call RANU0(dranu,1,icon)
              end do
              pbuf(m)%vx = vxt + vdtx(is)
              pbuf(m)%vy = vyt + vdty(is)
              pbuf(m)%vz = vzt + vdtz(is)
            end if
          else if((flag_collr(2)).and.pbuf(m)%preside.ne.1.and. &
         &   (pbuf(m)%x.le.xlcol(1).or.pbuf(m)%x.ge.xucol(1).or. &
         &    pbuf(m)%y.le.ylcol(1).or.pbuf(m)%y.ge.yucol(1).or. &
         &    pbuf(m)%z.le.zlcol(1).or.pbuf(m)%z.ge.zucol(1))) then
            vxt = peth(is)
            vyt = peth(is)
            vzt = path(is)
            vxyz = sqrt(vxt*vxt + vyt*vyt + vzt*vzt)
            call RANU0(dranu,1,icon)
            colfp = colf(1,ib ,jb ,k1b,ps)*v1 + colf(1,i1b,jb ,k1b,ps)*v2 &
           &      + colf(1,i1b,j1b,k1b,ps)*v3 + colf(1,ib ,j1b,k1b,ps)*v4 &
           &      + colf(1,ib ,jb ,kb ,ps)*v5 + colf(1,i1b,jb ,kb ,ps)*v6 &
           &      + colf(1,i1b,j1b,kb ,ps)*v7 + colf(1,ib ,j1b,kb ,ps)*v8
            if(exp(vxyz*nsigmadt).lt.dranu(1)*colfp) then
              call RANN0(0.0d0,peth(is),dranu,1,icon)
              vxt = dranu(1)
              call RANN0(0.0d0,peth(is),dranu,1,icon)
              vyt = dranu(1)
              call RANN0(0.0d0,path(is),dranu,1,icon)
              vzt = dranu(1)
              pbuf(m)%vx = vxt + vdtx(is)
              pbuf(m)%vy = vyt + vdty(is)
              pbuf(m)%vz = vzt + vdtz(is)
            end if
          end if

!         ----------- update particle positions
          pbuf(m)%preside = 0
          xmove = pbuf(m)%vx*gustep
          ymove = pbuf(m)%vy*gustep
          zmove = pbuf(m)%vz*gustep
          pbuf(m)%x = pbuf(m)%x + xmove
          pbuf(m)%y = pbuf(m)%y + ymove
          pbuf(m)%z = pbuf(m)%z + zmove

!         ----------- internal boundary treatment
          nyield(1:2) = 0
#include "defsurf.fnc"
if(pbuf(m)%preside.ge.0) then
#include "defbody.fnc"
end if

!          if(is.eq.4) then
!            print*, "B:pbuf(m)@is=4", pbuf(m)
!          end if

!         ----------- external boundary treatment 1
          pout_xl = (pbuf(m)%x.lt.0.0d0)
          pout_xu = (pbuf(m)%x.gt.slx)
          pout_yl = (pbuf(m)%y.lt.0.0d0)
          pout_yu = (pbuf(m)%y.gt.sly)
          pout_zl = (pbuf(m)%z.lt.0.0d0)
          pout_zu = (pbuf(m)%z.gt.slz)
!
          if(pout_xl) then
            if(npbnd_x_reflect) then
              pbuf(m)%x = -pbuf(m)%x
              pbuf(m)%vx = -pbuf(m)%vx
            end if
          else if(pout_xu) then
            if(npbnd_x_reflect) then
              pbuf(m)%x = tslx - pbuf(m)%x
              pbuf(m)%vx = -pbuf(m)%vx
            end if
          end if
          if(pout_yl) then
            if(npbnd_y_reflect) then
              pbuf(m)%y = -pbuf(m)%y
              pbuf(m)%vy = -pbuf(m)%vy
            end if
          else if(pout_yu) then
            if(npbnd_y_reflect) then
              pbuf(m)%y = tsly - pbuf(m)%y
              pbuf(m)%vy = -pbuf(m)%vy
            end if
          end if
          if(pout_zl) then
            if(npbnd_z_reflect) then
              pbuf(m)%z = -pbuf(m)%z
              pbuf(m)%vz = -pbuf(m)%vz
            end if
          else if(pout_zu) then
            if(npbnd_z_reflect) then
              pbuf(m)%z = tslz - pbuf(m)%z
              pbuf(m)%vz = -pbuf(m)%vz
            end if
          end if

!         ----------- 
          xlocala = pbuf(m)%x - dxl
          ylocala = pbuf(m)%y - dyl
          zlocala = pbuf(m)%z - dzl
!
          ia = floor(xlocala)
          ja = floor(ylocala)
          ka = floor(zlocala)
!
          xa = ia
          ya = ja
          za = ka
!
          i1a = ia + 1
          j1a = ja + 1
          k1a = ka + 1
!
          if(ib.eq.ia) then
            xr = (xlocalb + xlocala)*0.5d0
          else
            xr = max(xb,xa)
          end if
!
          if(jb.eq.ja) then
            yr = (ylocalb + ylocala)*0.5d0
          else
            yr = max(yb,ya)
          end if
!
          if(kb.eq.ka) then
            zr = (zlocalb + zlocala)*0.5d0
          else
            zr = max(zb,za)
          end if

!         ----------- external boundary treatment 2
          nud = 14
          if(pbuf(m)%x.lt.dxl) then
            nud = nud - 1
          else if(pbuf(m)%x.gt.dxu) then
            nud = nud + 1
          end if
          if(pbuf(m)%y.lt.dyl) then
            nud = nud - 3
          else if(pbuf(m)%y.gt.dyu) then
            nud = nud + 3
          end if
          if(pbuf(m)%z.lt.dzl) then
            nud = nud - 9
          else if(pbuf(m)%z.gt.dzu) then
            nud = nud + 9
          end if
!
          if(pbuf(m)%nid.ne.-1) then
            if(nud.ne.14) then
              if(flag_tq) then
                call oh2_remove_injected_particle(pbuf(m))
                pbuf(m)%nid = &
               &  oh3_map_particle_to_subdomain &
               &  (pbuf(m)%x,pbuf(m)%y,pbuf(m)%z)
              else
                nphgram(pbuf(m)%nid+1,is,pstq) = &
               &  nphgram(pbuf(m)%nid+1,is,pstq) - 1
                pbuf(m)%nid = nborps(nud,is,ps)
              end if
!
              if(pout_xl) then
                if(singridx) then
                  pbuf(m)%x = modulo(pbuf(m)%x,1.0d0)
                else if(npbnd_x_periodic) then
                  pbuf(m)%x = pbuf(m)%x + slx
                else if(npbnd_x_erase) then
                  pbuf(m)%nid = sdid(ps)
                  pbuf(m)%preside = IEBD
                end if
              else if(pout_xu) then
                if(singridx) then
                  pbuf(m)%x = modulo(pbuf(m)%x,1.0d0)
                else if(npbnd_x_periodic) then
                  pbuf(m)%x = pbuf(m)%x - slx
                else if(npbnd_x_erase) then
                  pbuf(m)%nid = sdid(ps)
                  pbuf(m)%preside = IEBD
                end if
              end if
              if(pout_yl) then
                if(singridy) then
                  pbuf(m)%y = modulo(pbuf(m)%y,1.0d0)
                else if(npbnd_y_periodic) then
                  pbuf(m)%y = pbuf(m)%y + sly
                else if(npbnd_y_erase) then
                  pbuf(m)%nid = sdid(ps)
                  pbuf(m)%preside = IEBD
                end if
              else if(pout_yu) then
                if(singridy) then
                  pbuf(m)%y = modulo(pbuf(m)%y,1.0d0)
                else if(npbnd_y_periodic) then
                  pbuf(m)%y = pbuf(m)%y - sly
                else if(npbnd_y_erase) then
                  pbuf(m)%nid = sdid(ps)
                  pbuf(m)%preside = IEBD
                end if
              end if
              if(pout_zl) then
                if(singridz) then
                  pbuf(m)%z = modulo(pbuf(m)%z,1.0d0)
                else if(npbnd_z_periodic) then
                  pbuf(m)%z = pbuf(m)%z + slz
                else if(npbnd_z_erase) then
                  pbuf(m)%nid = sdid(ps)
                  pbuf(m)%preside = IEBD
                end if
              else if(pout_zu) then
                if(singridz) then
                  pbuf(m)%z = modulo(pbuf(m)%z,1.0d0)
                else if(npbnd_z_periodic) then
                  pbuf(m)%z = pbuf(m)%z - slz
                else if(npbnd_z_erase) then
                  pbuf(m)%nid = sdid(ps)
                  pbuf(m)%preside = IEBD
                end if
              end if
!
              if(pbuf(m)%nid.eq.-1) then
                gcount(1)%nesc(is) = gcount(1)%nesc(is) + 1
              else
                if(flag_tq) then
                  call oh2_remap_injected_particle(pbuf(m))
                else
                  nphgram(pbuf(m)%nid+1,is,pstq) = &
                 &  nphgram(pbuf(m)%nid+1,is,pstq) + 1
                end if
              end if
            end if
!
#include "secemit.fnc"
          end if


!-------------------- 
          xd1 = (xr + xlocalb)*0.5d0 - xb
          yd1 = (yr + ylocalb)*0.5d0 - yb
          zd1 = (zr + zlocalb)*0.5d0 - zb
!
          vvx1 = (xr - xlocalb)*rustep
          vvy1 = (yr - ylocalb)*rustep
          vvz1 = (zr - zlocalb)*rustep
!
          vx1w4 = yd1*zd1
          vx1w2 = yd1 - vx1w4
          vx1w3 = zd1 - vx1w4
          vx1w1 = 1.0d0 - yd1 - vx1w3
          vy1w4 = xd1*zd1
          vy1w2 = xd1 - vy1w4
          vy1w3 = zd1 - vy1w4
          vy1w1 = 1.0d0 - xd1 - vy1w3
          vz1w4 = xd1*yd1
          vz1w2 = xd1 - vz1w4
          vz1w3 = yd1 - vz1w4
          vz1w1 = 1.0d0 - xd1 - vz1w3
!         ----------- 
          xd2 = (xlocala + xr)*0.5d0 - xa
          yd2 = (ylocala + yr)*0.5d0 - ya
          zd2 = (zlocala + zr)*0.5d0 - za
!
          vvx2 = (xlocala - xr)*rustep
          vvy2 = (ylocala - yr)*rustep
          vvz2 = (zlocala - zr)*rustep
!
          vx2w4 = yd2*zd2
          vx2w2 = yd2 - vx2w4
          vx2w3 = zd2 - vx2w4
          vx2w1 = 1.0d0 - yd2 - vx2w3
          vy2w4 = xd2*zd2
          vy2w2 = xd2 - vy2w4
          vy2w3 = zd2 - vy2w4
          vy2w1 = 1.0d0 - xd2 - vy2w3
          vz2w4 = xd2*yd2
          vz2w2 = xd2 - vz2w4
          vz2w3 = yd2 - vz2w4
          vz2w1 = 1.0d0 - xd2 - vz2w3


!-------------------- 
          wrk(TJX,ib ,jb ,kb ) = wrk(TJX,ib ,jb ,kb ) + vvx1*vx1w1
          wrk(TJY,ib ,jb ,kb ) = wrk(TJY,ib ,jb ,kb ) + vvy1*vy1w1
          wrk(TJZ,ib ,jb ,kb ) = wrk(TJZ,ib ,jb ,kb ) + vvz1*vz1w1

          wrk(TJX,ib ,j1b,kb ) = wrk(TJX,ib ,j1b,kb ) + vvx1*vx1w2
          wrk(TJY,i1b,jb ,kb ) = wrk(TJY,i1b,jb ,kb ) + vvy1*vy1w2
          wrk(TJZ,i1b,jb ,kb ) = wrk(TJZ,i1b,jb ,kb ) + vvz1*vz1w2

          wrk(TJX,ib ,jb ,k1b) = wrk(TJX,ib ,jb ,k1b) + vvx1*vx1w3
          wrk(TJY,ib ,jb ,k1b) = wrk(TJY,ib ,jb ,k1b) + vvy1*vy1w3
          wrk(TJZ,ib ,j1b,kb ) = wrk(TJZ,ib ,j1b,kb ) + vvz1*vz1w3

          wrk(TJX,ib ,j1b,k1b) = wrk(TJX,ib ,j1b,k1b) + vvx1*vx1w4
          wrk(TJY,i1b,jb ,k1b) = wrk(TJY,i1b,jb ,k1b) + vvy1*vy1w4
          wrk(TJZ,i1b,j1b,kb ) = wrk(TJZ,i1b,j1b,kb ) + vvz1*vz1w4

          wrk(TJX,ia ,ja ,ka ) = wrk(TJX,ia ,ja ,ka ) + vvx2*vx2w1
          wrk(TJY,ia ,ja ,ka ) = wrk(TJY,ia ,ja ,ka ) + vvy2*vy2w1
          wrk(TJZ,ia ,ja ,ka ) = wrk(TJZ,ia ,ja ,ka ) + vvz2*vz2w1

          wrk(TJX,ia ,j1a,ka ) = wrk(TJX,ia ,j1a,ka ) + vvx2*vx2w2
          wrk(TJY,i1a,ja ,ka ) = wrk(TJY,i1a,ja ,ka ) + vvy2*vy2w2
          wrk(TJZ,i1a,ja ,ka ) = wrk(TJZ,i1a,ja ,ka ) + vvz2*vz2w2

          wrk(TJX,ia ,ja ,k1a) = wrk(TJX,ia ,ja ,k1a) + vvx2*vx2w3
          wrk(TJY,ia ,ja ,k1a) = wrk(TJY,ia ,ja ,k1a) + vvy2*vy2w3
          wrk(TJZ,ia ,j1a,ka ) = wrk(TJZ,ia ,j1a,ka ) + vvz2*vz2w3

          wrk(TJX,ia ,j1a,k1a) = wrk(TJX,ia ,j1a,k1a) + vvx2*vx2w4
          wrk(TJY,i1a,ja ,k1a) = wrk(TJY,i1a,ja ,k1a) + vvy2*vy2w4
          wrk(TJZ,i1a,j1a,ka ) = wrk(TJZ,i1a,j1a,ka ) + vvz2*vz2w4
        end do ML1


!---------------- store current value from work to ajx
        qs = q(is)
        do k=-1,zu-zl+1
        do j=-1,yu-yl+1
        do i=-1,xu-xl+1
          aj(JX,i,j,k,ps) = aj(JX,i,j,k,ps) + wrk(TJX,i,j,k)*qs
          aj(JY,i,j,k,ps) = aj(JY,i,j,k,ps) + wrk(TJY,i,j,k)*qs
          aj(JZ,i,j,k,ps) = aj(JZ,i,j,k,ps) + wrk(TJZ,i,j,k)*qs
        end do
        end do
        end do
        if(func.eq.1) then
          iss = (is - 1)*3
          do k=-1,zu-zl+1
          do j=-1,yu-yl+1
          do i=-1,xu-xl+1
            ajdg(iss+JX,i,j,k,ps) = ajdg(iss+JX,i,j,k,ps) &
           &                      + wrk(TJX,i,j,k)*qs
            ajdg(iss+JY,i,j,k,ps) = ajdg(iss+JY,i,j,k,ps) &
           &                      + wrk(TJY,i,j,k)*qs
            ajdg(iss+JZ,i,j,k,ps) = ajdg(iss+JZ,i,j,k,ps) &
           &                      + wrk(TJZ,i,j,k)*qs
          end do
          end do
          end do
        end if
!
      end do ISL1
!============================== species loop end


!-------------------- boundary treatment for jx,jy,jz
!      call crntbd
!      call fbound(5)


!-------------------- filtering jx, jy, jz in x, y and z
!      call filtr3(ajx,ix,iy,iz,nx1p2,ny1p2,nz1p2,jxfltr)
!      call filtr3(ajy,ix,iy,iz,nx1p2,ny1p2,nz1p2,jyfltr)
!      call filtr3(ajz,ix,iy,iz,nx1p2,ny1p2,nz1p2,jzfltr)


!-------------------- uniform component cancellation
!      if(juncan.ge.1) call curcrt


!-------------------- masking current component
!      call fsmask(7)


!-------------------- increments of current
!      do k=1,nzm
!      do j=1,nym
!      do i=1,nxm
!        bjx(i,j,k) = rimlt*(ajx(i,j,k)-cjx(i,j,k))
!        bjy(i,j,k) = rimlt*(ajy(i,j,k)-cjy(i,j,k))
!        bjz(i,j,k) = rimlt*(ajz(i,j,k)-cjz(i,j,k))
!      end do
!      end do
!      end do


!-------------------- diagnostics
!    if(myid.eq.0) then
!      if(intfoc.ne.0) then
!        influxsum = 0
!        do ipc=1,npc
!          influx(ipc) = 0
!          do is=1,nspec
!            influx(ipc) = influx(ipc) + isfluxt(ipc,is)
!          end do
!          influxsum = influxsum + influx(ipc)
!        end do
!        if(istep.eq.0.or.mod(istep-1,intfoc).eq.0) then
!          open(90,file='influx',position='append')
!          open(91,file='isflux',position='append')
!          open(92,file='nesc',position='append')
!        end if
!        write(90,*) t, (influx(ipc),ipc=1,npc), influxsum
!        write(91,*) t, ((isflux(ipc,is),ipc=1,npc),is=1,nspec)
!        write(92,*) t, (nesc(is),is=1,nspec)
!        if(mod(istep,intfoc).eq.0.or.istep.eq.nstep) then
!          close(90)
!          close(91)
!          close(92)
!        end if
!      end if
!    end if


  return


  contains


    function round_prob(num,rand_number) result(rounded)
      real(kind=8), intent(in) :: num
      real(kind=8), intent(in) :: rand_number
      integer :: rounded
      real(kind=8) :: fractional_part

      fractional_part = num - nint(num)
      ! stochastic rounding up or down depending on the size of the decimal point
      if(fractional_part > 0.0d0) then
        if(rand_number < fractional_part) then
          rounded = nint(num) + 1.0d0 ! rounding up
        else
          rounded = nint(num)         ! rounding down
        end if
      else
        rounded = nint(num)           ! if only integer part, return as is
      end if
    end function round_prob


    function sey_model_1(delta,efrac,costh)
      real(kind=8) :: sey_model_1
      real(kind=8),intent(in) :: delta, efrac, costh

      sey_model_1 = delta*efrac &
     &                   *exp(2.0d0 - 2.0d0*sqrt(efrac)) &
     &                   *exp(2.0d0*(1.0d0 - costh))

    return
    end function sey_model_1


    function sey_model_2(delta,efrac,costh)
      real(kind=8) :: sey_model_2
      real(kind=8),intent(in) :: delta, efrac, costh

      sey_model_2 = 1.114d0*delta/costh*efrac**(-0.35d0) &
     &                     *(1.0d0 - exp(-2.28d0*costh*efrac**1.35d0))

    return
    end function sey_model_2


  end subroutine psolve1



!
  subroutine add_boundary_current1(ps)
!
!   ____________________________________________________________
!
!                       S U B R O U T I N E
!                 B O U N D A R Y _ C U R R E N T
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
  integer(kind=4) :: xu,yu,zu
  integer(kind=4) :: sl,su, dl,du, nl,nu
  integer(kind=4) :: ps


!-------------------- 
      xu = sdoms(2,1,sdid(ps)+1) - sdoms(1,1,sdid(ps)+1)
      yu = sdoms(2,2,sdid(ps)+1) - sdoms(1,2,sdid(ps)+1)
      zu = sdoms(2,3,sdid(ps)+1) - sdoms(1,3,sdid(ps)+1)


!-------------------- 
      sl = ctypes(2,2,1,CAJ)	!=-4
      su = ctypes(2,1,1,CAJ)	!=+2
!
      nl = ctypes(3,2,1,CAJ)	!= 3
      nu = ctypes(3,1,1,CAJ)	!= 3
!
      dl = sl + nl		!=-1
      du = su - nu		!=-1

!-------------------- 
      if(bounds(1,3,sdid(ps)+1).eq.1) then
        do k=0,nl-1		!(i.e., do k=0,2)
        do j=0,yu+(su+nu-sl)-1	!(i.e., do j=0,yu+8)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+8)
          aj(JX,sl+i,sl+j,dl+k,ps) = aj(JX,sl+i,sl+j,dl+k,ps) &
         &                         + aj(JX,sl+i,sl+j,sl+k,ps)
          aj(JY,sl+i,sl+j,dl+k,ps) = aj(JY,sl+i,sl+j,dl+k,ps) &
         &                         + aj(JY,sl+i,sl+j,sl+k,ps)
          aj(JZ,sl+i,sl+j,dl+k,ps) = aj(JZ,sl+i,sl+j,dl+k,ps) &
         &                         + aj(JZ,sl+i,sl+j,sl+k,ps)
        end do
        end do
        end do
      else if(bounds(1,3,sdid(ps)+1).eq.2.and.nfbnd(3).eq.0.and.nz.eq.1) then
        do k=0,nl-1		!(i.e., do k=0,2)
        do j=0,yu+(su+nu-sl)-1	!(i.e., do j=0,yu+8)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+8)
          aj(JX,sl+i,sl+j,dl+k,ps) = aj(JX,sl+i,sl+j,dl+k,ps) &
         &                         + aj(JX,sl+i,sl+j,sl+k,ps)
          aj(JY,sl+i,sl+j,dl+k,ps) = aj(JY,sl+i,sl+j,dl+k,ps) &
         &                         + aj(JY,sl+i,sl+j,sl+k,ps)
          aj(JZ,sl+i,sl+j,dl+k,ps) = aj(JZ,sl+i,sl+j,dl+k,ps) &
         &                         + aj(JZ,sl+i,sl+j,sl+k,ps)
        end do
        end do
        end do
        do k=0,nl-2		!(i.e., do k=0,1)
        do j=0,yu+(su+nu-sl)-1	!(i.e., do j=0,yu+8)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+8)
          aj(JX,sl+i,sl+j,dl+k,ps) = aj(JX,sl+i,sl+j,dl  +k,ps) &
         &                         + aj(JX,sl+i,sl+j,sl+1+k,ps)
          aj(JY,sl+i,sl+j,dl+k,ps) = aj(JY,sl+i,sl+j,dl  +k,ps) &
         &                         + aj(JY,sl+i,sl+j,sl+1+k,ps)
          aj(JZ,sl+i,sl+j,dl+k,ps) = aj(JZ,sl+i,sl+j,dl  +k,ps) &
         &                         + aj(JZ,sl+i,sl+j,sl+1+k,ps)
        end do
        end do
        end do
      end if

!-------------------- 
      if(bounds(2,3,sdid(ps)+1).eq.1) then
        do k=0,nu-1		!(i.e., do k=0,2)
        do j=0,yu+(su+nu-sl)-1	!(i.e., do j=0,yu+8)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+8)
          aj(JX,sl+i,sl+j,zu+du+k,ps) = aj(JX,sl+i,sl+j,zu+du+k,ps) &
         &                            + aj(JX,sl+i,sl+j,zu+su+k,ps)
          aj(JY,sl+i,sl+j,zu+du+k,ps) = aj(JY,sl+i,sl+j,zu+du+k,ps) &
         &                            + aj(JY,sl+i,sl+j,zu+su+k,ps)
          aj(JZ,sl+i,sl+j,zu+du+k,ps) = aj(JZ,sl+i,sl+j,zu+du+k,ps) &
         &                            + aj(JZ,sl+i,sl+j,zu+su+k,ps)
        end do
        end do
        end do
      else if(bounds(2,3,sdid(ps)+1).eq.2.and.nfbnd(3).eq.0.and.nz.eq.1) then
        do k=0,nu-1		!(i.e., do k=0,2)
        do j=0,yu+(su+nu-sl)-1	!(i.e., do j=0,yu+8)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+8)
          aj(JX,sl+i,sl+j,zu+du+k,ps) = aj(JX,sl+i,sl+j,zu+du+k,ps) &
         &                            + aj(JX,sl+i,sl+j,zu+su+k,ps)
          aj(JY,sl+i,sl+j,zu+du+k,ps) = aj(JY,sl+i,sl+j,zu+du+k,ps) &
         &                            + aj(JY,sl+i,sl+j,zu+su+k,ps)
          aj(JZ,sl+i,sl+j,zu+du+k,ps) = aj(JZ,sl+i,sl+j,zu+du+k,ps) &
         &                            + aj(JZ,sl+i,sl+j,zu+su+k,ps)
        end do
        end do
        end do
        do k=0,nu-2		!(i.e., do k=0,1)
        do j=0,yu+(su+nu-sl)-1	!(i.e., do j=0,yu+8)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+8)
          aj(JX,sl+i,sl+j,zu   +k,ps) = aj(JX,sl+i,sl+j,zu   +k,ps) &
         &                            + aj(JX,sl+i,sl+j,zu+su+k,ps)
          aj(JY,sl+i,sl+j,zu   +k,ps) = aj(JY,sl+i,sl+j,zu   +k,ps) &
         &                            + aj(JY,sl+i,sl+j,zu+su+k,ps)
          aj(JZ,sl+i,sl+j,zu   +k,ps) = aj(JZ,sl+i,sl+j,zu   +k,ps) &
         &                            + aj(JZ,sl+i,sl+j,zu+su+k,ps)
        end do
        end do
        end do
      end if

!-------------------- 
      if(bounds(1,2,sdid(ps)+1).eq.1) then
        do k=0,zu+(du+nu-dl)-1	!(i.e., do k=0,zu+2)
        do j=0,nl-1		!(i.e., do j=0,2)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+8)
          aj(JX,sl+i,dl+j,dl+k,ps) = aj(JX,sl+i,dl+j,dl+k,ps) &
         &                         + aj(JX,sl+i,sl+j,dl+k,ps)
          aj(JY,sl+i,dl+j,dl+k,ps) = aj(JY,sl+i,dl+j,dl+k,ps) &
         &                         + aj(JY,sl+i,sl+j,dl+k,ps)
          aj(JZ,sl+i,dl+j,dl+k,ps) = aj(JZ,sl+i,dl+j,dl+k,ps) &
         &                         + aj(JZ,sl+i,sl+j,dl+k,ps)
        end do
        end do
        end do
      else if(bounds(1,2,sdid(ps)+1).eq.2.and.nfbnd(2).eq.0.and.ny.eq.1) then
        do k=0,zu+(du+nu-dl)-1	!(i.e., do k=0,zu+2)
        do j=0,nl-1		!(i.e., do j=0,2)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+8)
          aj(JX,sl+i,dl+j,dl+k,ps) = aj(JX,sl+i,dl+j,dl+k,ps) &
         &                         + aj(JX,sl+i,sl+j,dl+k,ps)
          aj(JY,sl+i,dl+j,dl+k,ps) = aj(JY,sl+i,dl+j,dl+k,ps) &
         &                         + aj(JY,sl+i,sl+j,dl+k,ps)
          aj(JZ,sl+i,dl+j,dl+k,ps) = aj(JZ,sl+i,dl+j,dl+k,ps) &
         &                         + aj(JZ,sl+i,sl+j,dl+k,ps)
        end do
        end do
        end do
        do k=0,zu+(du+nu-dl)-1	!(i.e., do k=0,zu+2)
        do j=0,nl-2		!(i.e., do j=0,1)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+8)
          aj(JX,sl+i,dl+j,dl+k,ps) = aj(JX,sl+i,dl+j,dl  +k,ps) &
         &                         + aj(JX,sl+i,sl+j,dl+1+k,ps)
          aj(JY,sl+i,dl+j,dl+k,ps) = aj(JY,sl+i,dl+j,dl  +k,ps) &
         &                         + aj(JY,sl+i,sl+j,dl+1+k,ps)
          aj(JZ,sl+i,dl+j,dl+k,ps) = aj(JZ,sl+i,dl+j,dl  +k,ps) &
         &                         + aj(JZ,sl+i,sl+j,dl+1+k,ps)
        end do
        end do
        end do
      end if

!-------------------- 
      if(bounds(2,2,sdid(ps)+1).eq.1) then
        do k=0,zu+(du+nu-dl)-1	!(i.e., do k=0,zu+2)
        do j=0,nu-1		!(i.e., do j=0,2)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+8)
          aj(JX,sl+i,yu+du+j,dl+k,ps) = aj(JX,sl+i,yu+du+j,dl+k,ps) &
         &                            + aj(JX,sl+i,yu+su+j,dl+k,ps)
          aj(JY,sl+i,yu+du+j,dl+k,ps) = aj(JY,sl+i,yu+du+j,dl+k,ps) &
         &                            + aj(JY,sl+i,yu+su+j,dl+k,ps)
          aj(JZ,sl+i,yu+du+j,dl+k,ps) = aj(JZ,sl+i,yu+du+j,dl+k,ps) &
         &                            + aj(JZ,sl+i,yu+su+j,dl+k,ps)
        end do
        end do
        end do
      else if(bounds(2,2,sdid(ps)+1).eq.2.and.nfbnd(2).eq.0.and.ny.eq.1) then
        do k=0,zu+(du+nu-dl)-1	!(i.e., do k=0,zu+2)
        do j=0,nu-1		!(i.e., do j=0,2)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+8)
          aj(JX,sl+i,yu+du+j,dl+k,ps) = aj(JX,sl+i,yu+du+j,dl+k,ps) &
         &                            + aj(JX,sl+i,yu+su+j,dl+k,ps)
          aj(JY,sl+i,yu+du+j,dl+k,ps) = aj(JY,sl+i,yu+du+j,dl+k,ps) &
         &                            + aj(JY,sl+i,yu+su+j,dl+k,ps)
          aj(JZ,sl+i,yu+du+j,dl+k,ps) = aj(JZ,sl+i,yu+du+j,dl+k,ps) &
         &                            + aj(JZ,sl+i,yu+su+j,dl+k,ps)
        end do
        end do
        end do
        do k=0,zu+(du+nu-dl)-1	!(i.e., do k=0,zu+2)
        do j=0,nu-2		!(i.e., do j=0,1)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+8)
          aj(JX,sl+i,yu   +j,dl+k,ps) = aj(JX,sl+i,yu   +j,dl+k,ps) &
         &                            + aj(JX,sl+i,yu+su+j,dl+k,ps)
          aj(JY,sl+i,yu   +j,dl+k,ps) = aj(JY,sl+i,yu   +j,dl+k,ps) &
         &                            + aj(JY,sl+i,yu+su+j,dl+k,ps)
          aj(JZ,sl+i,yu   +j,dl+k,ps) = aj(JZ,sl+i,yu   +j,dl+k,ps) &
         &                            + aj(JZ,sl+i,yu+su+j,dl+k,ps)
        end do
        end do
        end do
      end if

!-------------------- 
      if(bounds(1,1,sdid(ps)+1).eq.1) then
        do k=0,zu+(du+nu-dl)-1	!(i.e., do k=0,zu+2)
        do j=0,yu+(du+nu-dl)-1	!(i.e., do j=0,yu+2)
        do i=0,nl-1		!(i.e., do i=0,2)
          aj(JX,dl+i,dl+j,dl+k,ps) = aj(JX,dl+i,dl+j,dl+k,ps) &
         &                         + aj(JX,sl+i,dl+j,dl+k,ps)
          aj(JY,dl+i,dl+j,dl+k,ps) = aj(JY,dl+i,dl+j,dl+k,ps) &
         &                         + aj(JY,sl+i,dl+j,dl+k,ps)
          aj(JZ,dl+i,dl+j,dl+k,ps) = aj(JZ,dl+i,dl+j,dl+k,ps) &
         &                         + aj(JZ,sl+i,dl+j,dl+k,ps)
        end do
        end do
        end do
      else if(bounds(1,1,sdid(ps)+1).eq.2.and.nfbnd(1).eq.0.and.nx.eq.1) then
        do k=0,zu+(du+nu-dl)-1	!(i.e., do k=0,zu+2)
        do j=0,yu+(du+nu-dl)-1	!(i.e., do j=0,yu+2)
        do i=0,nl-1		!(i.e., do i=0,2)
          aj(JX,dl+i,dl+j,dl+k,ps) = aj(JX,dl+i,dl+j,dl+k,ps) &
         &                         + aj(JX,sl+i,dl+j,dl+k,ps)
          aj(JY,dl+i,dl+j,dl+k,ps) = aj(JY,dl+i,dl+j,dl+k,ps) &
         &                         + aj(JY,sl+i,dl+j,dl+k,ps)
          aj(JZ,dl+i,dl+j,dl+k,ps) = aj(JZ,dl+i,dl+j,dl+k,ps) &
         &                         + aj(JZ,sl+i,dl+j,dl+k,ps)
        end do
        end do
        end do
        do k=0,zu+(du+nu-dl)-1	!(i.e., do k=0,zu+2)
        do j=0,yu+(du+nu-dl)-1	!(i.e., do j=0,yu+2)
        do i=0,nl-2		!(i.e., do i=0,1)
          aj(JX,dl+i,dl+j,dl+k,ps) = aj(JX,dl  +i,dl+j,dl+k,ps) &
         &                         + aj(JX,sl+1+i,dl+j,dl+k,ps)
          aj(JY,dl+i,dl+j,dl+k,ps) = aj(JY,dl  +i,dl+j,dl+k,ps) &
         &                         + aj(JY,sl+1+i,dl+j,dl+k,ps)
          aj(JZ,dl+i,dl+j,dl+k,ps) = aj(JZ,dl  +i,dl+j,dl+k,ps) &
         &                         + aj(JZ,sl+1+i,dl+j,dl+k,ps)
        end do
        end do
        end do
      end if

!-------------------- 
      if(bounds(2,1,sdid(ps)+1).eq.1) then
        do k=0,zu+(du+nu-dl)-1	!(i.e., do k=0,zu+2)
        do j=0,yu+(du+nu-dl)-1	!(i.e., do j=0,yu+2)
        do i=0,nu-1		!(i.e., do i=0,2)
          aj(JX,xu+du+i,dl+j,dl+k,ps) = aj(JX,xu+du+i,dl+j,dl+k,ps) &
         &                            + aj(JX,xu+su+i,dl+j,dl+k,ps)
          aj(JY,xu+du+i,dl+j,dl+k,ps) = aj(JY,xu+du+i,dl+j,dl+k,ps) &
         &                            + aj(JY,xu+su+i,dl+j,dl+k,ps)
          aj(JZ,xu+du+i,dl+j,dl+k,ps) = aj(JZ,xu+du+i,dl+j,dl+k,ps) &
         &                            + aj(JZ,xu+su+i,dl+j,dl+k,ps)
        end do
        end do
        end do
      else if(bounds(2,1,sdid(ps)+1).eq.1.and.nfbnd(1).eq.0.and.nx.eq.1) then
        do k=0,zu+(du+nu-dl)-1	!(i.e., do k=0,zu+2)
        do j=0,yu+(du+nu-dl)-1	!(i.e., do j=0,yu+2)
        do i=0,nu-1		!(i.e., do i=0,2)
          aj(JX,xu+du+i,dl+j,dl+k,ps) = aj(JX,xu+du+i,dl+j,dl+k,ps) &
         &                            + aj(JX,xu+su+i,dl+j,dl+k,ps)
          aj(JY,xu+du+i,dl+j,dl+k,ps) = aj(JY,xu+du+i,dl+j,dl+k,ps) &
         &                            + aj(JY,xu+su+i,dl+j,dl+k,ps)
          aj(JZ,xu+du+i,dl+j,dl+k,ps) = aj(JZ,xu+du+i,dl+j,dl+k,ps) &
         &                            + aj(JZ,xu+su+i,dl+j,dl+k,ps)
        end do
        end do
        end do
        do k=0,zu+(du+nu-dl)-1	!(i.e., do k=0,zu+2)
        do j=0,yu+(du+nu-dl)-1	!(i.e., do j=0,yu+2)
        do i=0,nu-2		!(i.e., do i=0,1)
          aj(JX,xu   +i,dl+j,dl+k,ps) = aj(JX,xu   +i,dl+j,dl+k,ps) &
         &                            + aj(JX,xu+su+i,dl+j,dl+k,ps)
          aj(JY,xu   +i,dl+j,dl+k,ps) = aj(JY,xu   +i,dl+j,dl+k,ps) &
         &                            + aj(JY,xu+su+i,dl+j,dl+k,ps)
          aj(JZ,xu   +i,dl+j,dl+k,ps) = aj(JZ,xu   +i,dl+j,dl+k,ps) &
         &                            + aj(JZ,xu+su+i,dl+j,dl+k,ps)
        end do
        end do
        end do
      end if


  return
  end subroutine add_boundary_current1


!
  subroutine add_boundary_current2(ps)
!
!   ____________________________________________________________
!
!                       S U B R O U T I N E
!                 B O U N D A R Y _ C U R R E N T
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
  integer(kind=4) :: xu,yu,zu
  integer(kind=4) :: sl,su, dl,du, nl,nu
  integer(kind=4) :: is
  integer(kind=4) :: ps


!-------------------- 
      xu = sdoms(2,1,sdid(ps)+1) - sdoms(1,1,sdid(ps)+1)
      yu = sdoms(2,2,sdid(ps)+1) - sdoms(1,2,sdid(ps)+1)
      zu = sdoms(2,3,sdid(ps)+1) - sdoms(1,3,sdid(ps)+1)


!-------------------- 
      sl = ctypes(2,2,1,CJD)	!=-4
      su = ctypes(2,1,1,CJD)	!=+2
!
      nl = ctypes(3,2,1,CJD)	!= 3
      nu = ctypes(3,1,1,CJD)	!= 3
!
      dl = sl + nl		!=-1
      du = su - nu		!=-1

!-------------------- 
      if(bounds(1,3,sdid(ps)+1).eq.1) then
        do k=0,nl-1		!(i.e., do k=0,2)
        do j=0,yu+(su+nu-sl)-1	!(i.e., do j=0,yu+8)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+8)
          ajdg(:,sl+i,sl+j,dl+k,ps) = ajdg(:,sl+i,sl+j,dl+k,ps) &
         &                          + ajdg(:,sl+i,sl+j,sl+k,ps)
        end do
        end do
        end do
      else if(bounds(1,3,sdid(ps)+1).eq.2.and.nfbnd(3).eq.0.and.nz.eq.1) then
        do k=0,nl-1		!(i.e., do k=0,2)
        do j=0,yu+(su+nu-sl)-1	!(i.e., do j=0,yu+8)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+8)
          ajdg(:,sl+i,sl+j,dl+k,ps) = ajdg(:,sl+i,sl+j,dl+k,ps) &
         &                          + ajdg(:,sl+i,sl+j,sl+k,ps)
        end do
        end do
        end do
        do k=0,nl-2		!(i.e., do k=0,1)
        do j=0,yu+(su+nu-sl)-1	!(i.e., do j=0,yu+8)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+8)
          ajdg(:,sl+i,sl+j,dl+k,ps) = ajdg(:,sl+i,sl+j,dl  +k,ps) &
         &                          + ajdg(:,sl+i,sl+j,sl+1+k,ps)
        end do
        end do
        end do
      end if

!-------------------- 
      if(bounds(2,3,sdid(ps)+1).eq.1) then
        do k=0,nu-1		!(i.e., do k=0,2)
        do j=0,yu+(su+nu-sl)-1	!(i.e., do j=0,yu+8)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+8)
          ajdg(:,sl+i,sl+j,zu+du+k,ps) = ajdg(:,sl+i,sl+j,zu+du+k,ps) &
         &                             + ajdg(:,sl+i,sl+j,zu+su+k,ps)
        end do
        end do
        end do
      else if(bounds(2,3,sdid(ps)+1).eq.2.and.nfbnd(3).eq.0.and.nz.eq.1) then
        do k=0,nu-1		!(i.e., do k=0,2)
        do j=0,yu+(su+nu-sl)-1	!(i.e., do j=0,yu+8)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+8)
          ajdg(:,sl+i,sl+j,zu+du+k,ps) = ajdg(:,sl+i,sl+j,zu+du+k,ps) &
         &                             + ajdg(:,sl+i,sl+j,zu+su+k,ps)
        end do
        end do
        end do
        do k=0,nu-2		!(i.e., do k=0,1)
        do j=0,yu+(su+nu-sl)-1	!(i.e., do j=0,yu+8)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+8)
          ajdg(:,sl+i,sl+j,zu   +k,ps) = ajdg(:,sl+i,sl+j,zu   +k,ps) &
         &                             + ajdg(:,sl+i,sl+j,zu+su+k,ps)
        end do
        end do
        end do
      end if

!-------------------- 
      if(bounds(1,2,sdid(ps)+1).eq.1) then
        do k=0,zu+(du+nu-dl)-1	!(i.e., do k=0,zu+2)
        do j=0,nl-1		!(i.e., do j=0,2)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+8)
          ajdg(:,sl+i,dl+j,dl+k,ps) = ajdg(:,sl+i,dl+j,dl+k,ps) &
         &                          + ajdg(:,sl+i,sl+j,dl+k,ps)
        end do
        end do
        end do
      else if(bounds(1,2,sdid(ps)+1).eq.2.and.nfbnd(2).eq.0.and.ny.eq.1) then
        do k=0,zu+(du+nu-dl)-1	!(i.e., do k=0,zu+2)
        do j=0,nl-1		!(i.e., do j=0,2)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+8)
          ajdg(:,sl+i,dl+j,dl+k,ps) = ajdg(:,sl+i,dl+j,dl+k,ps) &
         &                          + ajdg(:,sl+i,sl+j,dl+k,ps)
        end do
        end do
        end do
        do k=0,zu+(du+nu-dl)-1	!(i.e., do k=0,zu+2)
        do j=0,nl-2		!(i.e., do j=0,1)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+8)
          ajdg(:,sl+i,dl+j,dl+k,ps) = ajdg(:,sl+i,dl  +j,dl+k,ps) &
         &                          + ajdg(:,sl+i,sl+1+j,dl+k,ps)
        end do
        end do
        end do
      end if

!-------------------- 
      if(bounds(2,2,sdid(ps)+1).eq.1) then
        do k=0,zu+(du+nu-dl)-1	!(i.e., do k=0,zu+2)
        do j=0,nu-1		!(i.e., do j=0,2)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+8)
          ajdg(:,sl+i,yu+du+j,dl+k,ps) = ajdg(:,sl+i,yu+du+j,dl+k,ps) &
         &                             + ajdg(:,sl+i,yu+su+j,dl+k,ps)
        end do
        end do
        end do
      else if(bounds(2,2,sdid(ps)+1).eq.2.and.nfbnd(2).eq.0.and.ny.eq.1) then
        do k=0,zu+(du+nu-dl)-1	!(i.e., do k=0,zu+2)
        do j=0,nu-1		!(i.e., do j=0,2)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+8)
          ajdg(:,sl+i,yu+du+j,dl+k,ps) = ajdg(:,sl+i,yu+du+j,dl+k,ps) &
         &                             + ajdg(:,sl+i,yu+su+j,dl+k,ps)
        end do
        end do
        end do
        do k=0,zu+(du+nu-dl)-1	!(i.e., do k=0,zu+2)
        do j=0,nu-2		!(i.e., do j=0,1)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+8)
          ajdg(:,sl+i,yu   +j,dl+k,ps) = ajdg(:,sl+i,yu   +j,dl+k,ps) &
         &                             + ajdg(:,sl+i,yu+su+j,dl+k,ps)
        end do
        end do
        end do
      end if

!-------------------- 
      if(bounds(1,1,sdid(ps)+1).eq.1) then
        do k=0,zu+(du+nu-dl)-1	!(i.e., do k=0,zu+2)
        do j=0,yu+(du+nu-dl)-1	!(i.e., do j=0,yu+2)
        do i=0,nl-1		!(i.e., do i=0,2)
          ajdg(:,dl+i,dl+j,dl+k,ps) = ajdg(:,dl+i,dl+j,dl+k,ps) &
         &                          + ajdg(:,sl+i,dl+j,dl+k,ps)
        end do
        end do
        end do
      else if(bounds(1,1,sdid(ps)+1).eq.2.and.nfbnd(1).eq.0.and.nx.eq.1) then
        do k=0,zu+(du+nu-dl)-1	!(i.e., do k=0,zu+2)
        do j=0,yu+(du+nu-dl)-1	!(i.e., do j=0,yu+2)
        do i=0,nl-1		!(i.e., do i=0,2)
          ajdg(:,dl+i,dl+j,dl+k,ps) = ajdg(:,dl+i,dl+j,dl+k,ps) &
         &                          + ajdg(:,sl+i,dl+j,dl+k,ps)
        end do
        end do
        end do
        do k=0,zu+(du+nu-dl)-1	!(i.e., do k=0,zu+2)
        do j=0,yu+(du+nu-dl)-1	!(i.e., do j=0,yu+2)
        do i=0,nl-2		!(i.e., do i=0,1)
          ajdg(:,dl+i,dl+j,dl+k,ps) = ajdg(:,dl  +i,dl+j,dl+k,ps) &
         &                          + ajdg(:,sl+1+i,dl+j,dl+k,ps)
        end do
        end do
        end do
      end if

!-------------------- 
      if(bounds(2,1,sdid(ps)+1).eq.1) then
        do k=0,zu+(du+nu-dl)-1	!(i.e., do k=0,zu+2)
        do j=0,yu+(du+nu-dl)-1	!(i.e., do j=0,yu+2)
        do i=0,nu-1		!(i.e., do i=0,2)
          ajdg(:,xu+du+i,dl+j,dl+k,ps) = ajdg(:,xu+du+i,dl+j,dl+k,ps) &
         &                             + ajdg(:,xu+su+i,dl+j,dl+k,ps)
        end do
        end do
        end do
      else if(bounds(2,1,sdid(ps)+1).eq.2.and.nfbnd(1).eq.0.and.nx.eq.1) then
        do k=0,zu+(du+nu-dl)-1	!(i.e., do k=0,zu+2)
        do j=0,yu+(du+nu-dl)-1	!(i.e., do j=0,yu+2)
        do i=0,nu-1		!(i.e., do i=0,2)
          ajdg(:,xu+du+i,dl+j,dl+k,ps) = ajdg(:,xu+du+i,dl+j,dl+k,ps) &
         &                             + ajdg(:,xu+su+i,dl+j,dl+k,ps)
        end do
        end do
        end do
        do k=0,zu+(du+nu-dl)-1	!(i.e., do k=0,zu+2)
        do j=0,yu+(du+nu-dl)-1	!(i.e., do j=0,yu+2)
        do i=0,nu-2		!(i.e., do i=0,1)
          ajdg(:,xu   +i,dl+j,dl+k,ps) = ajdg(:,xu   +i,dl+j,dl+k,ps) &
         &                             + ajdg(:,xu+su+i,dl+j,dl+k,ps)
        end do
        end do
        end do
      end if


  return
  end subroutine add_boundary_current2



!
  subroutine exchange_lcur_borders1(ps)
!
!   ____________________________________________________________
!
!                       S U B R O U T I N E
!            E X C H A N G E _ L C U R _ B O R D E R S
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
  integer(kind=4),intent(in) :: ps
  integer(kind=4) :: i,j,k
  integer(kind=4) :: xu,yu,zu
  integer(kind=4) :: sl,su, dl,du, nl,nu



!-------------------- 
      xu = sdoms(2,1,sdid(ps)+1) - sdoms(1,1,sdid(ps)+1)
      yu = sdoms(2,2,sdid(ps)+1) - sdoms(1,2,sdid(ps)+1)
      zu = sdoms(2,3,sdid(ps)+1) - sdoms(1,3,sdid(ps)+1)


!-------------------- 
      sl = ctypes(2,2,1,CAJ)	!=-4
      su = ctypes(2,1,1,CAJ)	!=+2
!
      nl = ctypes(3,2,1,CAJ)	!= 3
      nu = ctypes(3,1,1,CAJ)	!= 3
!
      dl = sl + nl		!=-1
      du = su - nu		!=-1


!-------------------- 
      if(bounds(1,3,sdid(ps)+1).eq.2.and.nfbnd(3).eq.0.and.nz.eq.1) then
        do k=0,nl-1		!(i.e., do k=0,2)
        do j=0,yu+(su+nu-sl)-1	!(i.e., do j=0,yu+8)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+8)
          aj(JX:JZ,sl+i,sl+j,zu+su+k,ps) = aj(JX:JZ,sl+i,sl+j,dl+k,ps)
        end do
        end do
        end do
      end if

!-------------------- 
      if(bounds(2,3,sdid(ps)+1).eq.2.and.nfbnd(3).eq.0.and.nz.eq.1) then
        do k=0,nu-1		!(i.e., do k=0,2)
        do j=0,yu+(su+nu-sl)-1	!(i.e., do j=0,yu+8)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+8)
          aj(JX:JZ,sl+i,sl+j,sl+k,ps) = aj(JX:JZ,sl+i,sl+j,zu+du+k,ps)
        end do
        end do
        end do
      end if

!-------------------- 
      if(bounds(1,2,sdid(ps)+1).eq.2.and.nfbnd(2).eq.0.and.ny.eq.1) then
        do k=0,zu+(su+nu-sl)-1	!(i.e., do k=0,zu+8)
        do j=0,nl-1		!(i.e., do j=0,2)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+8)
          aj(JX:JZ,sl+i,yu+su+j,sl+k,ps) = aj(JX:JZ,sl+i,dl+j,sl+k,ps)
        end do
        end do
        end do
      end if

!-------------------- 
      if(bounds(2,2,sdid(ps)+1).eq.2.and.nfbnd(2).eq.0.and.ny.eq.1) then
        do k=0,zu+(su+nu-sl)-1	!(i.e., do k=0,zu+8)
        do j=0,nu-1		!(i.e., do j=0,2)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+8)
          aj(JX:JZ,sl+i,sl+j,sl+k,ps) = aj(JX:JZ,sl+i,yu+du+j,sl+k,ps)
        end do
        end do
        end do
      end if

!-------------------- 
      if(bounds(1,1,sdid(ps)+1).eq.2.and.nfbnd(1).eq.0.and.nx.eq.1) then
        do k=0,zu+(su+nu-sl)-1	!(i.e., do k=0,zu+8)
        do j=0,yu+(su+nu-sl)-1	!(i.e., do j=0,yu+8)
        do i=0,nl-1		!(i.e., do i=0,2)
          aj(JX:JZ,xu+su+i,sl+j,sl+k,ps) = aj(JX:JZ,dl+i,sl+j,sl+k,ps)
        end do
        end do
        end do
      end if

!-------------------- 
      if(bounds(2,1,sdid(ps)+1).eq.2.and.nfbnd(1).eq.0.and.nx.eq.1) then
        do k=0,zu+(su+nu-sl)-1	!(i.e., do k=0,zu+8)
        do j=0,yu+(su+nu-sl)-1	!(i.e., do j=0,yu+8)
        do i=0,nu-1		!(i.e., do i=0,2)
          aj(JX:JZ,sl+i,sl+j,sl+k,ps) = aj(JX:JZ,xu+du+i,sl+j,sl+k,ps)
        end do
        end do
        end do
      end if


  return
  end subroutine exchange_lcur_borders1



!
  subroutine exchange_lcur_borders2(ps)
!
!   ____________________________________________________________
!
!                       S U B R O U T I N E
!            E X C H A N G E _ L C U R _ B O R D E R S
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
  integer(kind=4),intent(in) :: ps
  integer(kind=4) :: i,j,k
  integer(kind=4) :: xu,yu,zu
  integer(kind=4) :: sl,su, dl,du, nl,nu


!-------------------- 
      xu = sdoms(2,1,sdid(ps)+1) - sdoms(1,1,sdid(ps)+1)
      yu = sdoms(2,2,sdid(ps)+1) - sdoms(1,2,sdid(ps)+1)
      zu = sdoms(2,3,sdid(ps)+1) - sdoms(1,3,sdid(ps)+1)


!-------------------- 
      sl = ctypes(2,2,1,CJD)	!=-4
      su = ctypes(2,1,1,CJD)	!=+2
!
      nl = ctypes(3,2,1,CJD)	!= 3
      nu = ctypes(3,1,1,CJD)	!= 3
!
      dl = sl + nl		!=-1
      du = su - nu		!=-1


!-------------------- 
      if(bounds(1,3,sdid(ps)+1).eq.2.and.nfbnd(3).eq.0.and.nz.eq.1) then
        do k=0,nl-1		!(i.e., do k=0,2)
        do j=0,yu+(su+nu-sl)-1	!(i.e., do j=0,yu+8)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+8)
          ajdg(:,sl+i,sl+j,zu+su+k,ps) = ajdg(:,sl+i,sl+j,dl+k,ps)
        end do
        end do
        end do
      end if

!-------------------- 
      if(bounds(2,3,sdid(ps)+1).eq.2.and.nfbnd(3).eq.0.and.nz.eq.1) then
        do k=0,nu-1		!(i.e., do k=0,2)
        do j=0,yu+(su+nu-sl)-1	!(i.e., do j=0,yu+8)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+8)
          ajdg(:,sl+i,sl+j,sl+k,ps) = ajdg(:,sl+i,sl+j,zu+du+k,ps)
        end do
        end do
        end do
      end if

!-------------------- 
      if(bounds(1,2,sdid(ps)+1).eq.2.and.nfbnd(2).eq.0.and.ny.eq.1) then
        do k=0,zu+(su+nu-sl)-1	!(i.e., do k=0,zu+8)
        do j=0,nl-1		!(i.e., do j=0,2)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+8)
          ajdg(:,sl+i,yu+su+j,sl+k,ps) = ajdg(:,sl+i,dl+j,sl+k,ps)
        end do
        end do
        end do
      end if

!-------------------- 
      if(bounds(2,2,sdid(ps)+1).eq.2.and.nfbnd(2).eq.0.and.ny.eq.1) then
        do k=0,zu+(su+nu-sl)-1	!(i.e., do k=0,zu+8)
        do j=0,nu-1		!(i.e., do j=0,2)
        do i=0,xu+(su+nu-sl)-1	!(i.e., do i=0,xu+8)
          ajdg(:,sl+i,sl+j,sl+k,ps) = ajdg(:,sl+i,yu+du+j,sl+k,ps)
        end do
        end do
        end do
      end if

!-------------------- 
      if(bounds(1,1,sdid(ps)+1).eq.2.and.nfbnd(1).eq.0.and.nx.eq.1) then
        do k=0,zu+(su+nu-sl)-1	!(i.e., do k=0,zu+8)
        do j=0,yu+(su+nu-sl)-1	!(i.e., do j=0,yu+8)
        do i=0,nl-1		!(i.e., do i=0,2)
          ajdg(:,xu+su+i,sl+j,sl+k,ps) = ajdg(:,dl+i,sl+j,sl+k,ps)
        end do
        end do
        end do
      end if

!-------------------- 
      if(bounds(2,1,sdid(ps)+1).eq.2.and.nfbnd(1).eq.0.and.nx.eq.1) then
        do k=0,zu+(su+nu-sl)-1	!(i.e., do k=0,zu+8)
        do j=0,yu+(su+nu-sl)-1	!(i.e., do j=0,yu+8)
        do i=0,nu-1		!(i.e., do i=0,2)
          ajdg(:,sl+i,sl+j,sl+k,ps) = ajdg(:,xu+du+i,sl+j,sl+k,ps)
        end do
        end do
        end do
      end if


  return
  end subroutine exchange_lcur_borders2



!
  subroutine add_source_current(ps)
!
!   ____________________________________________________________
!
!                       S U B R O U T I N E
!              A D D _ S O U R C E _ C U R R E N T
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
  integer(kind=4) :: i,j,k, ijs, jcomp
  integer(kind=4) :: xl,yl,zl, xu,yu,zu
  integer(kind=4) :: ps


!-------------------- 
      xl = sdoms(1,1,sdid(ps)+1); xu = sdoms(2,1,sdid(ps)+1)
      yl = sdoms(1,2,sdid(ps)+1); yu = sdoms(2,2,sdid(ps)+1)
      zl = sdoms(1,3,sdid(ps)+1); zu = sdoms(2,3,sdid(ps)+1)


!-------------------- 
      do ijs=1,njs
        i = rjs(1,ijs); j = rjs(2,ijs); k = rjs(3,ijs)
        if(i.ge.xl.and.i.lt.xu.and. &
       &   j.ge.yl.and.j.lt.yu.and. &
       &   k.ge.zl.and.k.lt.zu) then
          aj(JX:JZ,i-xl,j-yl,k-zl,ps) &
         &  = aj(JX:JZ,i-xl,j-yl,k-zl,ps) &
         &  + ajs(JX:JZ,ijs) &
         &  *sin(wjs(ijs)*(2.0d0*istep + 1.0d0) + th0js(ijs))
        end if
      end do


  return
  end subroutine add_source_current
