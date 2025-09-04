#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine frmttd(ustep)
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   F R M T T D
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
  integer(kind=4) :: nintgrtdstep = 4000
  integer(kind=4) :: ustep, ndnstep
  integer(kind=4) :: ipc,is
  integer(kind=4) :: minpc
  real(kind=8) :: icur(3,inpc,ispec), ocur(3,inpc,ispec)
  real(kind=8) :: expema
  real(kind=8) :: tew
  character(len=40) :: sformat

  integer(kind=4) :: i,j,k
  integer(kind=4) :: xl,xu, yl,yu, zl,zu
  real(kind=8) :: eai0,eas0,eas1,eas2
  real(kind=8) :: bai0,bas0,bas1,bas2
  real(kind=8) :: bzi0,bzs0,bzs1,bzs2


!-------------------- 
      xl = sdoms(1,1,sdid(1)+1); xu = sdoms(2,1,sdid(1)+1)
      yl = sdoms(1,2,sdid(1)+1); yu = sdoms(2,2,sdid(1)+1)
      zl = sdoms(1,3,sdid(1)+1); zu = sdoms(2,3,sdid(1)+1)


!-------------------- 
      if(nstep.gt.0) then
        ndnstep = int(log10(real(nstep*2))) + 1
      else
        ndnstep = 6
      end if


!-------------------- 
      if(istep.gt.0.and.emarlxt.gt.0) then
        expema = exp(-1.0d0/emarlxt)
      else
        expema = 0.0d0
      end if


!-------------------- 
      minpc = max(npc,1)


!-------------------- diagnostics
      if(myid.eq.0) then
        if(istep.eq.0.and.nintgrtdstep.gt.nstep) then
          if(nstep.gt.0) then
            nintgrtdstep = nstep
          else
            nintgrtdstep = 1
          end if
        end if
        if(istep.eq.0) then
          icur(:,:,:) = 0.0d0
          ocur(:,:,:) = 0.0d0
        end if
        if(istep.eq.nstep-nintgrtdstep) then
!          selfpintgrtd(1:npc+1) = 0.0d0
!          outfluxintgrtd(1:npc,1:nspec) = 0
!          infhistintgrtd(0:npc,1:npc,1:nspec) = 0
!          rhoindintgrtd(1:npc) = 0.0d0
!          vidata(1:10) = 0.0d0
        end if
        if(istep.gt.nstep-nintgrtdstep) then
!          selfpintgrtd(1:npc) = selfpintgrtd(1:npc) + selfp(1:npc)
!          selfpintgrtd(npc+1) = selfpintgrtd(npc+1) + phiref(2)
!          outfluxintgrtd(1:npc,1:nspec) = &
!         &  outfluxintgrtd(1:npc,1:nspec) &
!         &  + gcount(2)%outflux(1:npc,1:nspec)
!          infhistintgrtd(0:npc,1:npc,1:nspec) = &
!         &  infhistintgrtd(0:npc,1:npc,1:nspec) &
!         &  + gcount(2)%infhist(0:npc,1:npc,1:nspec)
!          rhoindintgrtd(1:npc) = &
!         &  rhoindintgrtd(1:npc) &
!         &  + rhoind(1:npc)/wrelax - rhoindh(1:npc)/wrelax - gcount(2)%chgacm(2,1:npc)
!          vidata( 1) = vidata( 1) + (selfp(1)+selfp(npc))*0.5d0
!          vidata( 2) = vidata( 2) + (biasc(1)%val+biasc(2)%val)*0.5d0
!          vidata( 3) = vidata( 3) + (sum(gcount(2)%outflux(1:2,nspec)) &
!         &                          +sum(gcount(2)%outflux(npc-1:npc,nspec)) &
!         &                          -sum(gcount(2)%infhist(1,1:npc,nspec)) &
!         &                          -sum(gcount(2)%infhist(2,1:npc,nspec)) &
!         &                          -sum(gcount(2)%infhist(npc-1,1:npc,nspec)) &
!         &                          -sum(gcount(2)%infhist(npc,1:npc,nspec)))*q(3)*0.5d0
!          vidata( 4) = vidata( 4) - (sum(gcount(2)%infhist(0,1:2,1)) &
!         &                          +sum(gcount(2)%infhist(0,npc-1:npc,1)))*q(1)*0.5d0
!          vidata( 5) = vidata( 5) - (sum(gcount(2)%infhist(0,1:2,2)) &
!         &                          +sum(gcount(2)%infhist(0,npc-1:npc,2)))*q(2)*0.5d0
!          vidata( 6) = vidata( 6) + (sum(gcount(2)%infhist(1,3:npc,nspec)) &
!         &                          +sum(gcount(2)%infhist(2,3:npc,nspec)) &
!         &                          +sum(gcount(2)%infhist(npc-1,1:npc-2,nspec)) &
!         &                          +sum(gcount(2)%infhist(npc,1:npc-2,nspec)))*q(3)*0.5d0
!          vidata( 7) = vidata( 7) - (sum(gcount(2)%infhist(3:npc,1,nspec)) &
!         &                          +sum(gcount(2)%infhist(3:npc,2,nspec)) &
!         &                          +sum(gcount(2)%infhist(1:npc-2,npc-1,nspec)) &
!         &                          +sum(gcount(2)%infhist(1:npc-2,npc,nspec)))*q(3)*0.5d0
!          vidata( 8) = vidata( 8) - (   (gcount(2)%infhist(3,1,nspec)) &
!         &                          +   (gcount(2)%infhist(3,2,nspec)) &
!         &                          +   (gcount(2)%infhist(npc-2,npc-1,nspec)) &
!         &                          +   (gcount(2)%infhist(npc-2,npc,nspec)))*q(3)*0.5d0
!          vidata( 9) = vidata( 9) - (   (gcount(2)%infhist(4,1,nspec)) &
!         &                          +   (gcount(2)%infhist(4,2,nspec)) &
!         &                          +   (gcount(2)%infhist(npc-3,npc-1,nspec)) &
!         &                          +   (gcount(2)%infhist(npc-3,npc,nspec)))*q(3)*0.5d0
!          vidata(10) = vidata(10) - (sum(gcount(2)%infhist(5:npc-4,1,nspec)) &
!         &                          +sum(gcount(2)%infhist(5:npc-4,2,nspec)) &
!         &                          +sum(gcount(2)%infhist(5:npc-4,npc-1,nspec)) &
!         &                          +sum(gcount(2)%infhist(5:npc-4,npc,nspec)))*q(3)*0.5d0
        end if
        if(istep.eq.nstep) then
!          selfpintgrtd(1:npc+1) = selfpintgrtd(1:npc+1)/nintgrtdstep/renphi
!          outfluxintgrtd(1:npc,1:nspec) = &
!         &  outfluxintgrtd(1:npc,1:nspec)/nintgrtdstep
!          infhistintgrtd(0:npc,1:npc,1:nspec) = &
!         &  infhistintgrtd(0:npc,1:npc,1:nspec)/nintgrtdstep
!          rhoindintgrtd(1:npc) = &
!         &  rhoindintgrtd(1:npc)/nintgrtdstep/renq*rent
!          vidata(1) = vidata(1)/nintgrtdstep/renphi
!          vidata(2:10) = vidata(2:10)/nintgrtdstep/renq
        end if
!
        if(intfoc.ne.0.and.mod(istep,intfoc).eq.0) then
          print*,"  globalp =",gcount(2)%globalp(1:nspec,:)
        end if
        if(intfoc.gt.0) then
          if(istep.eq.0.or.mod(istep-1,intfoc).eq.0) then
            open(FILEEN,file='energy',position='append')
            open(FILEPP,file='pprobe',position='append')
            open(FILEPB,file='pbody',position='append')
            open(FILEIC,file='icur',position='append')
            open(FILEOC,file='ocur',position='append')
!            open(79,file='energy',position='append')
!            open(86,file='pbody',position='append')
!            open(89,file='influx',position='append')
!            open(90,file='noflux',position='append')
!            open(91,file='isflux',position='append')
!            open(92,file='nesc',position='append')
            open(93,file='chgacm1',position='append')
            open(94,file='chgacm2',position='append')
!            open(95,file='chgmov',position='append')
            open(96,file='seyield',position='append')
!            open(97,file='icur',position='append')
!            open(98,file='ocur',position='append')
!            open(99,file='ewave',position='append')
          end if
!         ----------- 
          sformat = '(Ix.x,1X,xxxx(ES14.7e2,1X))'
          write(sformat(3:5),'(I1,A,I1)') ndnstep, ".", ndnstep
          write(sformat(10:13),'(I4.4)') 2+3*nspec
          write(FILEEN,sformat) istep, engebg(1:2)/rened, engkg(1:3,1:nspec)/rened
!         ----------- 
          if(npprb.gt.0) then
            sformat = '(Ix.x,1X,xxxx(ES14.7e2,1X))'
            write(sformat(10:13),'(I4.4)') npprb
          else
            sformat = '(Ix.x)'
          end if
          write(sformat(3:5),'(I1,A,I1)') ndnstep, ".", ndnstep
          write(FILEPP,sformat) istep, pprb(1:npprb,2)/renphi
!         ----------- 
          sformat = '(Ix.x,1X,xxxx(ES14.7e2,1X))'
          write(sformat(3:5),'(I1,A,I1)') ndnstep, ".", ndnstep
          write(sformat(10:13),'(I4.4)') npc+1
          write(FILEPB,sformat) istep, selfp(1:npc)/renphi, phiref(2)/renphi
!         ----------- 
!          sformat = '(Ix.x,1X,xxxx(I,1X))'
!          write(sformat(3:5),'(I1,A,I1)') ndnstep, ".", ndnstep
!          write(sformat(10:13),'(I4.4)') minpc*minsp
!          write(89,sformat) istep, int(gcount(2)%influx(1:minpc,1:minsp)*2/ustep)
!         -----------
          if(nspec.gt.0.and.npc.gt.0) then
            sformat = '(xxxx(Ix.x,xxxx(1X,ES9.4e1),2X))'
            write(sformat(2:5),'(I4.4)') nspec
            write(sformat(8:10),'(I1,A,I1)') ndnstep, ".", ndnstep
            write(sformat(12:15),'(I4.4)') npc*2
            do is=1,nspec
            do ipc=1,npc
              icur(1,ipc,is) = abs(gcount(2)%influx(ipc,is)*2/ustep*q(is)/renq/dt)
              icur(1,ipc,is) = icur(1,ipc,is)*sqdscaled(ipc)
              icur(3,ipc,is) = icur(3,ipc,is)*expema &
             &               + icur(1,ipc,is)*(1.0d0-expema)
              if(abs(icur(3,ipc,is)).ge.1.0d-9) then
                icur(2,ipc,is) = icur(3,ipc,is)
              else
                icur(2,ipc,is) = 0.0d0
              end if
            end do
            end do
            write(FILEIC,sformat) (istep, icur(1:2,1:npc,is), is=1,nspec)
          else
            sformat = '(Ix.x)'
            write(sformat(3:5),'(I1,A,I1)') ndnstep, ".", ndnstep
            write(FILEIC,sformat) istep
          end if
!         ----------- 
!          sformat = '(Ix.x,1X,xxxx(I,1X))'
!          write(sformat(3:5),'(I1,A,I1)') ndnstep, ".", ndnstep
!          write(sformat(10:13),'(I4.4)') minpc*minsp
!          write(90,sformat) istep, int(gcount(2)%outflux(1:minpc,1:minsp)*2/ustep)
!         ----------- 
          if(nspec.gt.0.and.npc.gt.0) then
            sformat = '(xxxx(Ix.x,xxxx(1X,ES9.4e1),2X))'
            write(sformat(2:5),'(I4.4)') nspec
            write(sformat(8:10),'(I1,A,I1)') ndnstep, ".", ndnstep
            write(sformat(12:15),'(I4.4)') npc*2
            do is=1,nspec
            do ipc=1,npc
              ocur(1,ipc,is) = abs(gcount(2)%outflux(ipc,is)*2/ustep*q(is)/renq/dt)
              ocur(1,ipc,is) = ocur(1,ipc,is)*sqdscaled(ipc)
              ocur(3,ipc,is) = ocur(3,ipc,is)*expema &
             &               + ocur(1,ipc,is)*(1.0d0-expema)
              if(abs(ocur(3,ipc,is)).ge.1.0d-9) then
                ocur(2,ipc,is) = ocur(3,ipc,is)
              else
                ocur(2,ipc,is) = 0.0d0
              end if
            end do
            end do
            write(FILEOC,sformat) (istep, ocur(1:2,1:minpc,is), is=1,minsp)
          else
            sformat = '(Ix.x)'
            write(sformat(3:5),'(I1,A,I1)') ndnstep, ".", ndnstep
            write(FILEOC,sformat) istep
          end if
!         -----------
!          sformat = '(Ix.x,1X,xxxx(I,1X))'
!          write(sformat(3:5),'(I1,A,I1)') ndnstep, ".", ndnstep
!          write(sformat(10:13),'(I4.4)') (1+npc)*minpc*minsp
!          write(91,sformat) istep, int(gcount(2)%infhist(0:npc,1:minpc,1:minsp))
!         ----------- 
!          sformat = '(Ix.x,1X,xxxx(I,1X))'
!          write(sformat(3:5),'(I1,A,I1)') ndnstep, ".", ndnstep
!          write(sformat(10:13),'(I4.4)') nspec
!          write(92,sformat) istep, int(gcount(2)%nesc(1:nspec)*2/ustep)
!         ----------- 
          sformat = '(Ix.x,1X,xxxx(ES14.7e2,1X))'
          write(sformat(3:5),'(I1,A,I1)') ndnstep, ".", ndnstep
          write(sformat(10:13),'(I4.4)') minpc
          write(93,sformat) istep, gcount(2)%chgacm(1,1:minpc)/renq
!         ----------- 
          sformat = '(Ix.x,1X,xxxx(ES14.7e2,1X))'
          write(sformat(3:5),'(I1,A,I1)') ndnstep, ".", ndnstep
          write(sformat(10:13),'(I4.4)') minpc
          write(94,sformat) istep, gcount(2)%chgacm(2,1:minpc)/renq
!         ----------- 



!         ----------- 
!          sformat = '(Ix.x,1X,xxxx(ES14.7e2,1X))'
!          write(sformat(3:5),'(I1,A,I1)') ndnstep, ".", ndnstep
!          write(sformat(10:13),'(I4.4)') minpc
!          write(95,sformat) istep, (rhoind(1:npc)/wrelax-rhoindh(1:npc)/wrelax-gcount(2)%chgacm(2,1:npc))/renq*rent*2.0d0/ustep
!         ----------- 
          sformat = '(Ix.x,1X,xxxx(ES14.7e2,1X))'
          write(sformat(3:5),'(I1,A,I1)') ndnstep, ".", ndnstep
          write(sformat(10:13),'(I4.4)') 3
          write(96,sformat) istep, abs(seygl(1:3,3)*2/ustep/renq/dt)
!          write(96,'(f7.2,1X,36(f9.2,1X))') &
!         &  t, seygl(1:11,1),seygl(12,1)*4.0d0,seygl(1:11,3), &
!         &     seygl(12,3)*4.0d0,seygl(1:11,4),seygl(12,4)*4.0d0
!         ----------- 
!          sformat = '(Ix.x,xxxx(1X,xxxx(1X,ES9.4e1)))'
!          write(sformat(3:5),'(I1,A,I1)') ndnstep, ".", ndnstep
!          write(sformat(7:10),'(I4.4)') nspec + 1
!          write(sformat(15:18),'(I4.4)') max(2,1)
!          tew = t - dt*nretard
!          write(99,sformat) istep,(Ew(1,is,1)*cos(omegaw(1)*tew)/rene, &
!         &                         Ew(2,is,1)*sin(omegaw(1)*tew)/rene,is=0,nspec)
!         ----------- 
!          sformat = '(Ix.x,1X,xxxx(ES14.7e2,1X))'
!          write(sformat(3:5),'(I1,A,I1)') ndnstep, ".", ndnstep
!          write(sformat(10:13),'(I4.4)') 10
!          write(99,sformat) istep, (selfp(1)+selfp(npc))*0.5d0/renphi, &

!         &               (biasc(1)%val+biasc(2)%val)*0.5d0/renq, &

!         &              +(sum(gcount(2)%outflux(1:2,nspec)) &
!         &               +sum(gcount(2)%outflux(npc-1:npc,nspec)) &
!         &               -sum(gcount(2)%infhist(1,1:npc,nspec)) &
!         &               -sum(gcount(2)%infhist(2,1:npc,nspec)) &
!         &               -sum(gcount(2)%infhist(npc-1,1:npc,nspec)) &
!         &               -sum(gcount(2)%infhist(npc,1:npc,nspec)))*q(3)*0.5d0/renq, &

!         &              -(sum(gcount(2)%infhist(0,1:2,1)) &
!         &               +sum(gcount(2)%infhist(0,npc-1:npc,1)))*q(1)*0.5d0/renq, &

!         &              -(sum(gcount(2)%infhist(0,1:2,2)) &
!         &               +sum(gcount(2)%infhist(0,npc-1:npc,2)))*q(2)*0.5d0/renq, &

!         &              +(sum(gcount(2)%infhist(1,3:npc,nspec)) &
!         &               +sum(gcount(2)%infhist(2,3:npc,nspec)) &
!         &               +sum(gcount(2)%infhist(npc-1,1:npc-2,nspec)) &
!         &               +sum(gcount(2)%infhist(npc,1:npc-2,nspec)))*q(3)*0.5d0/renq, &

!         &              -(sum(gcount(2)%infhist(3:npc,1,nspec)) &
!         &               +sum(gcount(2)%infhist(3:npc,2,nspec)) &
!         &               +sum(gcount(2)%infhist(1:npc-2,npc-1,nspec)) &
!         &               +sum(gcount(2)%infhist(1:npc-2,npc,nspec)))*q(3)*0.5d0/renq, &

!         &              -(   (gcount(2)%infhist(3,1,nspec)) &
!         &               +   (gcount(2)%infhist(3,2,nspec)) &
!         &               +   (gcount(2)%infhist(npc-2,npc-1,nspec)) &
!         &               +   (gcount(2)%infhist(npc-2,npc,nspec)))*q(3)*0.5d0/renq, &

!         &              -(   (gcount(2)%infhist(4,1,nspec)) &
!         &               +   (gcount(2)%infhist(4,2,nspec)) &
!         &               +   (gcount(2)%infhist(npc-3,npc-1,nspec)) &
!         &               +   (gcount(2)%infhist(npc-3,npc,nspec)))*q(3)*0.5d0/renq, &

!         &              -(sum(gcount(2)%infhist(5:npc-4,1,nspec)) &
!         &               +sum(gcount(2)%infhist(5:npc-4,2,nspec)) &
!         &               +sum(gcount(2)%infhist(5:npc-4,npc-1,nspec)) &
!         &               +sum(gcount(2)%infhist(5:npc-4,npc,nspec)))*q(3)*0.5d0/renq
!          if(istep.eq.nstep) &
!         &  write(99,sformat) istep*2, vidata(1:10)
          if(mod(istep,intfoc).eq.0.or.istep.eq.nstep) then
            close(FILEEN)
            close(FILEPP)
            close(FILEPB)
            close(FILEIC)
            close(FILEOC)
!            close(79)
!            close(86)
!            close(89)
!            close(90)
!            close(91)
!            close(92)
            close(93)
            close(94)
!            close(95)
            close(96)
!            close(97)
!            close(98)
!            close(99)

            close(501)
            close(502)
            close(503)
            close(504)
            close(505)
            close(506)
            close(507)
            close(508)
          end if
        end if
      end if
!
!-------------------- clear particle counter
      gcount(1)%influx = 0
      gcount(1)%outflux = 0
      gcount(1)%infhist = 0
      gcount(1)%nesc = 0
      gcount(1)%chgacm(2,:) = 0.0d0


!-------------------- 
      i = 2560 - 250 - xl
      j = 2560 - yl
      k = 0 - zl
      if(i.ge.0.and.i.lt.xu-xl.and. &
         j.ge.0.and.j.lt.yu-yl.and. &
         k.ge.0.and.k.lt.zu-zl) then
        if(istep.eq.0) then
          open(501,file='eabs_i0_0250',status='replace')
          open(502,file='eabs_s0_0250',status='replace')
          open(503,file='eabs_s1_0250',status='replace')
          open(504,file='eabs_s2_0250',status='replace')
          open(505,file='babs_i0_0250',status='replace')
          open(506,file='babs_s0_0250',status='replace')
          open(507,file='babs_s1_0250',status='replace')
          open(508,file='babs_s2_0250',status='replace')
          open(509,file='bz_i0_0250',status='replace')
          open(510,file='bz_s0_0250',status='replace')
          open(511,file='bz_s1_0250',status='replace')
          open(512,file='bz_s2_0250',status='replace')
        end if
!
        eai0 = sqrt(0.25d0* &
       &            (eb(EX,i,j,k,5) + eb(EX,i-1,j,k,5))* &
       &            (eb(EX,i,j,k,5) + eb(EX,i-1,j,k,5)) &
       &          + 0.25d0* &
       &            (eb(EY,i,j,k,5) + eb(EY,i,j-1,k,5))* &
       &            (eb(EY,i,j,k,5) + eb(EY,i,j-1,k,5)) &
       &          + 0.25d0* &
       &            (eb(EZ,i,j,k,5) + eb(EZ,i,j,k-1,5))* &
       &            (eb(EZ,i,j,k,5) + eb(EZ,i,j,k-1,5)) )
        eas0 = sqrt(0.25d0* &
       &            (eb(EX,i,j,k,1) + eb(EX,i-1,j,k,1) + &
       &             eb(EX,i,j,k,3) + eb(EX,i-1,j,k,3))* &
       &            (eb(EX,i,j,k,1) + eb(EX,i-1,j,k,1) + &
       &             eb(EX,i,j,k,3) + eb(EX,i-1,j,k,3)) &
       &          + 0.25d0* &
       &            (eb(EY,i,j,k,1) + eb(EY,i,j-1,k,1) + &
       &             eb(EY,i,j,k,3) + eb(EY,i,j-1,k,3))* &
       &            (eb(EY,i,j,k,1) + eb(EY,i,j-1,k,1) + &
       &             eb(EY,i,j,k,3) + eb(EY,i,j-1,k,3)) &
       &          + 0.25d0* &
       &            (eb(EZ,i,j,k,1) + eb(EZ,i,j,k-1,1) + &
       &             eb(EZ,i,j,k,3) + eb(EZ,i,j,k-1,3))* &
       &            (eb(EZ,i,j,k,1) + eb(EZ,i,j,k-1,1) + &
       &             eb(EZ,i,j,k,3) + eb(EZ,i,j,k-1,3)) )
        eas1 = sqrt(0.25d0* &
       &            (eb(EX,i,j,k,1) + eb(EX,i-1,j,k,1))* &
       &            (eb(EX,i,j,k,1) + eb(EX,i-1,j,k,1)) &
       &          + 0.25d0* &
       &            (eb(EY,i,j,k,1) + eb(EY,i,j-1,k,1))* &
       &            (eb(EY,i,j,k,1) + eb(EY,i,j-1,k,1)) &
       &          + 0.25d0* &
       &            (eb(EZ,i,j,k,1) + eb(EZ,i,j,k-1,1))* &
       &            (eb(EZ,i,j,k,1) + eb(EZ,i,j,k-1,1)) )
        eas2 = sqrt(0.25d0* &
       &            (eb(EX,i,j,k,3) + eb(EX,i-1,j,k,3))* &
       &            (eb(EX,i,j,k,3) + eb(EX,i-1,j,k,3)) &
       &          + 0.25d0* &
       &            (eb(EY,i,j,k,3) + eb(EY,i,j-1,k,3))* &
       &            (eb(EY,i,j,k,3) + eb(EY,i,j-1,k,3)) &
       &          + 0.25d0* &
       &            (eb(EZ,i,j,k,3) + eb(EZ,i,j,k-1,3))* &
       &            (eb(EZ,i,j,k,3) + eb(EZ,i,j,k-1,3)) )
!
        bai0 = sqrt(0.0625d0* &
       &            (eb(BX,i,j,k  ,5) + eb(BX,i,j-1,k  ,5) + &
       &             eb(BX,i,j,k-1,5) + eb(BX,i,j-1,k-1,5))* &
       &            (eb(BX,i,j,k  ,5) + eb(BX,i,j-1,k  ,5) + &
       &             eb(BX,i,j,k-1,5) + eb(BX,i,j-1,k-1,5)) &
       &          + 0.0625d0* &
       &            (eb(BY,i,j,k  ,5) + eb(BY,i-1,j,k  ,5) + &
       &             eb(BY,i,j,k-1,5) + eb(BY,i-1,j,k-1,5))* &
       &            (eb(BY,i,j,k  ,5) + eb(BY,i-1,j,k  ,5) + &
       &             eb(BY,i,j,k-1,5) + eb(BY,i-1,j,k-1,5)) &
       &          + 0.0625d0* &
       &            (eb(BZ,i,j  ,k,5) + eb(BZ,i-1,j  ,k,5) + &
       &             eb(BZ,i,j-1,k,5) + eb(BZ,i-1,j-1,k,5))* &
       &            (eb(BZ,i,j  ,k,5) + eb(BZ,i-1,j  ,k,5) + &
       &             eb(BZ,i,j-1,k,5) + eb(BZ,i-1,j-1,k,5)) )
        bas0 = sqrt(0.0625d0* &
       &            (eb(BX,i,j,k  ,1) + eb(BX,i,j-1,k  ,1) + &
       &             eb(BX,i,j,k-1,1) + eb(BX,i,j-1,k-1,1) + &
       &             eb(BX,i,j,k  ,3) + eb(BX,i,j-1,k  ,3) + &
       &             eb(BX,i,j,k-1,3) + eb(BX,i,j-1,k-1,3))* &
       &            (eb(BX,i,j,k  ,1) + eb(BX,i,j-1,k  ,1) + &
       &             eb(BX,i,j,k-1,1) + eb(BX,i,j-1,k-1,1) + &
       &             eb(BX,i,j,k  ,3) + eb(BX,i,j-1,k  ,3) + &
       &             eb(BX,i,j,k-1,3) + eb(BX,i,j-1,k-1,3)) &
       &          + 0.0625d0* &
       &            (eb(BY,i,j,k  ,1) + eb(BY,i-1,j,k  ,1) + &
       &             eb(BY,i,j,k-1,1) + eb(BY,i-1,j,k-1,1) + &
       &             eb(BY,i,j,k  ,3) + eb(BY,i-1,j,k  ,3) + &
       &             eb(BY,i,j,k-1,3) + eb(BY,i-1,j,k-1,3))* &
       &            (eb(BY,i,j,k  ,1) + eb(BY,i-1,j,k  ,1) + &
       &             eb(BY,i,j,k-1,1) + eb(BY,i-1,j,k-1,1) + &
       &             eb(BY,i,j,k  ,3) + eb(BY,i-1,j,k  ,3) + &
       &             eb(BY,i,j,k-1,3) + eb(BY,i-1,j,k-1,3)) &
       &          + 0.0625d0* &
       &            (eb(BZ,i,j  ,k,1) + eb(BZ,i-1,j  ,k,1) + &
       &             eb(BZ,i,j-1,k,1) + eb(BZ,i-1,j-1,k,1) + &
       &             eb(BZ,i,j  ,k,3) + eb(BZ,i-1,j  ,k,3) + &
       &             eb(BZ,i,j-1,k,3) + eb(BZ,i-1,j-1,k,3))* &
       &            (eb(BZ,i,j  ,k,1) + eb(BZ,i-1,j  ,k,1) + &
       &             eb(BZ,i,j-1,k,1) + eb(BZ,i-1,j-1,k,1) + &
       &             eb(BZ,i,j  ,k,3) + eb(BZ,i-1,j  ,k,3) + &
       &             eb(BZ,i,j-1,k,3) + eb(BZ,i-1,j-1,k,3)) )
        bas1 = sqrt(0.0625d0* &
       &            (eb(BX,i,j,k  ,1) + eb(BX,i,j-1,k  ,1) + &
       &             eb(BX,i,j,k-1,1) + eb(BX,i,j-1,k-1,1))* &
       &            (eb(BX,i,j,k  ,1) + eb(BX,i,j-1,k  ,1) + &
       &             eb(BX,i,j,k-1,1) + eb(BX,i,j-1,k-1,1)) &
       &          + 0.0625d0* &
       &            (eb(BY,i,j,k  ,1) + eb(BY,i-1,j,k  ,1) + &
       &             eb(BY,i,j,k-1,1) + eb(BY,i-1,j,k-1,1))* &
       &            (eb(BY,i,j,k  ,1) + eb(BY,i-1,j,k  ,1) + &
       &             eb(BY,i,j,k-1,1) + eb(BY,i-1,j,k-1,1)) &
       &          + 0.0625d0* &
       &            (eb(BZ,i,j  ,k,1) + eb(BZ,i-1,j  ,k,1) + &
       &             eb(BZ,i,j-1,k,1) + eb(BZ,i-1,j-1,k,1))* &
       &            (eb(BZ,i,j  ,k,1) + eb(BZ,i-1,j  ,k,1) + &
       &             eb(BZ,i,j-1,k,1) + eb(BZ,i-1,j-1,k,1)) )
        bas2 = sqrt(0.0625d0* &
       &            (eb(BX,i,j,k  ,3) + eb(BX,i,j-1,k  ,3) + &
       &             eb(BX,i,j,k-1,3) + eb(BX,i,j-1,k-1,3))* &
       &            (eb(BX,i,j,k  ,3) + eb(BX,i,j-1,k  ,3) + &
       &             eb(BX,i,j,k-1,3) + eb(BX,i,j-1,k-1,3)) &
       &          + 0.0625d0* &
       &            (eb(BY,i,j,k  ,3) + eb(BY,i-1,j,k  ,3) + &
       &             eb(BY,i,j,k-1,3) + eb(BY,i-1,j,k-1,3))* &
       &            (eb(BY,i,j,k  ,3) + eb(BY,i-1,j,k  ,3) + &
       &             eb(BY,i,j,k-1,3) + eb(BY,i-1,j,k-1,3)) &
       &          + 0.0625d0* &
       &            (eb(BZ,i,j  ,k,3) + eb(BZ,i-1,j  ,k,3) + &
       &             eb(BZ,i,j-1,k,3) + eb(BZ,i-1,j-1,k,3))* &
       &            (eb(BZ,i,j  ,k,3) + eb(BZ,i-1,j  ,k,3) + &
       &             eb(BZ,i,j-1,k,3) + eb(BZ,i-1,j-1,k,3)) )
        bzi0 = 0.25d0* &
       &            (eb(BZ,i,j  ,k,5) + eb(BZ,i-1,j  ,k,5) + &
       &             eb(BZ,i,j-1,k,5) + eb(BZ,i-1,j-1,k,5))
        bzs1 = 0.25d0* &
       &            (eb(BZ,i,j  ,k,1) + eb(BZ,i-1,j  ,k,1) + &
       &             eb(BZ,i,j-1,k,1) + eb(BZ,i-1,j-1,k,1))
        bzs2 = 0.25d0* &
       &            (eb(BZ,i,j  ,k,3) + eb(BZ,i-1,j  ,k,3) + &
       &             eb(BZ,i,j-1,k,3) + eb(BZ,i-1,j-1,k,3))
        bzs0 = bzs1 + bzs2
!
        write(501,*) eai0/abs(rene)
        write(502,*) eas0/abs(rene)
        write(503,*) eas1/abs(rene)
        write(504,*) eas2/abs(rene)
        write(505,*) bai0/abs(renb)
        write(506,*) bas0/abs(renb)
        write(507,*) bas1/abs(renb)
        write(508,*) bas2/abs(renb)
        write(509,*) bzi0/renb
        write(510,*) bzs0/renb
        write(511,*) bzs1/renb
        write(512,*) bzs2/renb
!
        if(istep.eq.nstep) then
          close(501)
          close(502)
          close(503)
          close(504)
          close(505)
          close(506)
          close(507)
          close(508)
          close(509)
          close(510)
          close(511)
          close(512)
        end if
      end if


!-------------------- 
      i = 2560 - 500 - xl
      j = 2560 - yl
      k = 0 - zl
      if(i.ge.0.and.i.lt.xu-xl.and. &
         j.ge.0.and.j.lt.yu-yl.and. &
         k.ge.0.and.k.lt.zu-zl) then
        if(istep.eq.0) then
          open(601,file='eabs_i0_0500',status='replace')
          open(602,file='eabs_s0_0500',status='replace')
          open(603,file='eabs_s1_0500',status='replace')
          open(604,file='eabs_s2_0500',status='replace')
          open(605,file='babs_i0_0500',status='replace')
          open(606,file='babs_s0_0500',status='replace')
          open(607,file='babs_s1_0500',status='replace')
          open(608,file='babs_s2_0500',status='replace')
          open(609,file='bz_i0_0500',status='replace')
          open(610,file='bz_s0_0500',status='replace')
          open(611,file='bz_s1_0500',status='replace')
          open(612,file='bz_s2_0500',status='replace')
        end if
!
        eai0 = sqrt(0.25d0* &
       &            (eb(EX,i,j,k,5) + eb(EX,i-1,j,k,5))* &
       &            (eb(EX,i,j,k,5) + eb(EX,i-1,j,k,5)) &
       &          + 0.25d0* &
       &            (eb(EY,i,j,k,5) + eb(EY,i,j-1,k,5))* &
       &            (eb(EY,i,j,k,5) + eb(EY,i,j-1,k,5)) &
       &          + 0.25d0* &
       &            (eb(EZ,i,j,k,5) + eb(EZ,i,j,k-1,5))* &
       &            (eb(EZ,i,j,k,5) + eb(EZ,i,j,k-1,5)) )
        eas0 = sqrt(0.25d0* &
       &            (eb(EX,i,j,k,1) + eb(EX,i-1,j,k,1) + &
       &             eb(EX,i,j,k,3) + eb(EX,i-1,j,k,3))* &
       &            (eb(EX,i,j,k,1) + eb(EX,i-1,j,k,1) + &
       &             eb(EX,i,j,k,3) + eb(EX,i-1,j,k,3)) &
       &          + 0.25d0* &
       &            (eb(EY,i,j,k,1) + eb(EY,i,j-1,k,1) + &
       &             eb(EY,i,j,k,3) + eb(EY,i,j-1,k,3))* &
       &            (eb(EY,i,j,k,1) + eb(EY,i,j-1,k,1) + &
       &             eb(EY,i,j,k,3) + eb(EY,i,j-1,k,3)) &
       &          + 0.25d0* &
       &            (eb(EZ,i,j,k,1) + eb(EZ,i,j,k-1,1) + &
       &             eb(EZ,i,j,k,3) + eb(EZ,i,j,k-1,3))* &
       &            (eb(EZ,i,j,k,1) + eb(EZ,i,j,k-1,1) + &
       &             eb(EZ,i,j,k,3) + eb(EZ,i,j,k-1,3)) )
        eas1 = sqrt(0.25d0* &
       &            (eb(EX,i,j,k,1) + eb(EX,i-1,j,k,1))* &
       &            (eb(EX,i,j,k,1) + eb(EX,i-1,j,k,1)) &
       &          + 0.25d0* &
       &            (eb(EY,i,j,k,1) + eb(EY,i,j-1,k,1))* &
       &            (eb(EY,i,j,k,1) + eb(EY,i,j-1,k,1)) &
       &          + 0.25d0* &
       &            (eb(EZ,i,j,k,1) + eb(EZ,i,j,k-1,1))* &
       &            (eb(EZ,i,j,k,1) + eb(EZ,i,j,k-1,1)) )
        eas2 = sqrt(0.25d0* &
       &            (eb(EX,i,j,k,3) + eb(EX,i-1,j,k,3))* &
       &            (eb(EX,i,j,k,3) + eb(EX,i-1,j,k,3)) &
       &          + 0.25d0* &
       &            (eb(EY,i,j,k,3) + eb(EY,i,j-1,k,3))* &
       &            (eb(EY,i,j,k,3) + eb(EY,i,j-1,k,3)) &
       &          + 0.25d0* &
       &            (eb(EZ,i,j,k,3) + eb(EZ,i,j,k-1,3))* &
       &            (eb(EZ,i,j,k,3) + eb(EZ,i,j,k-1,3)) )
!
        bai0 = sqrt(0.0625d0* &
       &            (eb(BX,i,j,k  ,5) + eb(BX,i,j-1,k  ,5) + &
       &             eb(BX,i,j,k-1,5) + eb(BX,i,j-1,k-1,5))* &
       &            (eb(BX,i,j,k  ,5) + eb(BX,i,j-1,k  ,5) + &
       &             eb(BX,i,j,k-1,5) + eb(BX,i,j-1,k-1,5)) &
       &          + 0.0625d0* &
       &            (eb(BY,i,j,k  ,5) + eb(BY,i-1,j,k  ,5) + &
       &             eb(BY,i,j,k-1,5) + eb(BY,i-1,j,k-1,5))* &
       &            (eb(BY,i,j,k  ,5) + eb(BY,i-1,j,k  ,5) + &
       &             eb(BY,i,j,k-1,5) + eb(BY,i-1,j,k-1,5)) &
       &          + 0.0625d0* &
       &            (eb(BZ,i,j  ,k,5) + eb(BZ,i-1,j  ,k,5) + &
       &             eb(BZ,i,j-1,k,5) + eb(BZ,i-1,j-1,k,5))* &
       &            (eb(BZ,i,j  ,k,5) + eb(BZ,i-1,j  ,k,5) + &
       &             eb(BZ,i,j-1,k,5) + eb(BZ,i-1,j-1,k,5)) )
        bas0 = sqrt(0.0625d0* &
       &            (eb(BX,i,j,k  ,1) + eb(BX,i,j-1,k  ,1) + &
       &             eb(BX,i,j,k-1,1) + eb(BX,i,j-1,k-1,1) + &
       &             eb(BX,i,j,k  ,3) + eb(BX,i,j-1,k  ,3) + &
       &             eb(BX,i,j,k-1,3) + eb(BX,i,j-1,k-1,3))* &
       &            (eb(BX,i,j,k  ,1) + eb(BX,i,j-1,k  ,1) + &
       &             eb(BX,i,j,k-1,1) + eb(BX,i,j-1,k-1,1) + &
       &             eb(BX,i,j,k  ,3) + eb(BX,i,j-1,k  ,3) + &
       &             eb(BX,i,j,k-1,3) + eb(BX,i,j-1,k-1,3)) &
       &          + 0.0625d0* &
       &            (eb(BY,i,j,k  ,1) + eb(BY,i-1,j,k  ,1) + &
       &             eb(BY,i,j,k-1,1) + eb(BY,i-1,j,k-1,1) + &
       &             eb(BY,i,j,k  ,3) + eb(BY,i-1,j,k  ,3) + &
       &             eb(BY,i,j,k-1,3) + eb(BY,i-1,j,k-1,3))* &
       &            (eb(BY,i,j,k  ,1) + eb(BY,i-1,j,k  ,1) + &
       &             eb(BY,i,j,k-1,1) + eb(BY,i-1,j,k-1,1) + &
       &             eb(BY,i,j,k  ,3) + eb(BY,i-1,j,k  ,3) + &
       &             eb(BY,i,j,k-1,3) + eb(BY,i-1,j,k-1,3)) &
       &          + 0.0625d0* &
       &            (eb(BZ,i,j  ,k,1) + eb(BZ,i-1,j  ,k,1) + &
       &             eb(BZ,i,j-1,k,1) + eb(BZ,i-1,j-1,k,1) + &
       &             eb(BZ,i,j  ,k,3) + eb(BZ,i-1,j  ,k,3) + &
       &             eb(BZ,i,j-1,k,3) + eb(BZ,i-1,j-1,k,3))* &
       &            (eb(BZ,i,j  ,k,1) + eb(BZ,i-1,j  ,k,1) + &
       &             eb(BZ,i,j-1,k,1) + eb(BZ,i-1,j-1,k,1) + &
       &             eb(BZ,i,j  ,k,3) + eb(BZ,i-1,j  ,k,3) + &
       &             eb(BZ,i,j-1,k,3) + eb(BZ,i-1,j-1,k,3)) )
        bas1 = sqrt(0.0625d0* &
       &            (eb(BX,i,j,k  ,1) + eb(BX,i,j-1,k  ,1) + &
       &             eb(BX,i,j,k-1,1) + eb(BX,i,j-1,k-1,1))* &
       &            (eb(BX,i,j,k  ,1) + eb(BX,i,j-1,k  ,1) + &
       &             eb(BX,i,j,k-1,1) + eb(BX,i,j-1,k-1,1)) &
       &          + 0.0625d0* &
       &            (eb(BY,i,j,k  ,1) + eb(BY,i-1,j,k  ,1) + &
       &             eb(BY,i,j,k-1,1) + eb(BY,i-1,j,k-1,1))* &
       &            (eb(BY,i,j,k  ,1) + eb(BY,i-1,j,k  ,1) + &
       &             eb(BY,i,j,k-1,1) + eb(BY,i-1,j,k-1,1)) &
       &          + 0.0625d0* &
       &            (eb(BZ,i,j  ,k,1) + eb(BZ,i-1,j  ,k,1) + &
       &             eb(BZ,i,j-1,k,1) + eb(BZ,i-1,j-1,k,1))* &
       &            (eb(BZ,i,j  ,k,1) + eb(BZ,i-1,j  ,k,1) + &
       &             eb(BZ,i,j-1,k,1) + eb(BZ,i-1,j-1,k,1)) )
        bas2 = sqrt(0.0625d0* &
       &            (eb(BX,i,j,k  ,3) + eb(BX,i,j-1,k  ,3) + &
       &             eb(BX,i,j,k-1,3) + eb(BX,i,j-1,k-1,3))* &
       &            (eb(BX,i,j,k  ,3) + eb(BX,i,j-1,k  ,3) + &
       &             eb(BX,i,j,k-1,3) + eb(BX,i,j-1,k-1,3)) &
       &          + 0.0625d0* &
       &            (eb(BY,i,j,k  ,3) + eb(BY,i-1,j,k  ,3) + &
       &             eb(BY,i,j,k-1,3) + eb(BY,i-1,j,k-1,3))* &
       &            (eb(BY,i,j,k  ,3) + eb(BY,i-1,j,k  ,3) + &
       &             eb(BY,i,j,k-1,3) + eb(BY,i-1,j,k-1,3)) &
       &          + 0.0625d0* &
       &            (eb(BZ,i,j  ,k,3) + eb(BZ,i-1,j  ,k,3) + &
       &             eb(BZ,i,j-1,k,3) + eb(BZ,i-1,j-1,k,3))* &
       &            (eb(BZ,i,j  ,k,3) + eb(BZ,i-1,j  ,k,3) + &
       &             eb(BZ,i,j-1,k,3) + eb(BZ,i-1,j-1,k,3)) )
        bzi0 = 0.25d0* &
       &            (eb(BZ,i,j  ,k,5) + eb(BZ,i-1,j  ,k,5) + &
       &             eb(BZ,i,j-1,k,5) + eb(BZ,i-1,j-1,k,5))
        bzs1 = 0.25d0* &
       &            (eb(BZ,i,j  ,k,1) + eb(BZ,i-1,j  ,k,1) + &
       &             eb(BZ,i,j-1,k,1) + eb(BZ,i-1,j-1,k,1))
        bzs2 = 0.25d0* &
       &            (eb(BZ,i,j  ,k,3) + eb(BZ,i-1,j  ,k,3) + &
       &             eb(BZ,i,j-1,k,3) + eb(BZ,i-1,j-1,k,3))
        bzs0 = bzs1 + bzs2
!
        write(601,*) eai0/abs(rene)
        write(602,*) eas0/abs(rene)
        write(603,*) eas1/abs(rene)
        write(604,*) eas2/abs(rene)
        write(605,*) bai0/abs(renb)
        write(606,*) bas0/abs(renb)
        write(607,*) bas1/abs(renb)
        write(608,*) bas2/abs(renb)
        write(609,*) bzi0/renb
        write(610,*) bzs0/renb
        write(611,*) bzs1/renb
        write(612,*) bzs2/renb
!
        if(istep.eq.nstep) then
          close(601)
          close(602)
          close(603)
          close(604)
          close(605)
          close(606)
          close(607)
          close(608)
          close(609)
          close(610)
          close(611)
          close(612)
        end if
      end if


!-------------------- 
      i = 2560 - 750 - xl
      j = 2560 - yl
      k = 0 - zl
      if(i.ge.0.and.i.lt.xu-xl.and. &
         j.ge.0.and.j.lt.yu-yl.and. &
         k.ge.0.and.k.lt.zu-zl) then
        if(istep.eq.0) then
          open(701,file='eabs_i0_0750',status='replace')
          open(702,file='eabs_s0_0750',status='replace')
          open(703,file='eabs_s1_0750',status='replace')
          open(704,file='eabs_s2_0750',status='replace')
          open(705,file='babs_i0_0750',status='replace')
          open(706,file='babs_s0_0750',status='replace')
          open(707,file='babs_s1_0750',status='replace')
          open(708,file='babs_s2_0750',status='replace')
          open(709,file='bz_i0_0750',status='replace')
          open(710,file='bz_s0_0750',status='replace')
          open(711,file='bz_s1_0750',status='replace')
          open(712,file='bz_s2_0750',status='replace')
        end if
!
        eai0 = sqrt(0.25d0* &
       &            (eb(EX,i,j,k,5) + eb(EX,i-1,j,k,5))* &
       &            (eb(EX,i,j,k,5) + eb(EX,i-1,j,k,5)) &
       &          + 0.25d0* &
       &            (eb(EY,i,j,k,5) + eb(EY,i,j-1,k,5))* &
       &            (eb(EY,i,j,k,5) + eb(EY,i,j-1,k,5)) &
       &          + 0.25d0* &
       &            (eb(EZ,i,j,k,5) + eb(EZ,i,j,k-1,5))* &
       &            (eb(EZ,i,j,k,5) + eb(EZ,i,j,k-1,5)) )
        eas0 = sqrt(0.25d0* &
       &            (eb(EX,i,j,k,1) + eb(EX,i-1,j,k,1) + &
       &             eb(EX,i,j,k,3) + eb(EX,i-1,j,k,3))* &
       &            (eb(EX,i,j,k,1) + eb(EX,i-1,j,k,1) + &
       &             eb(EX,i,j,k,3) + eb(EX,i-1,j,k,3)) &
       &          + 0.25d0* &
       &            (eb(EY,i,j,k,1) + eb(EY,i,j-1,k,1) + &
       &             eb(EY,i,j,k,3) + eb(EY,i,j-1,k,3))* &
       &            (eb(EY,i,j,k,1) + eb(EY,i,j-1,k,1) + &
       &             eb(EY,i,j,k,3) + eb(EY,i,j-1,k,3)) &
       &          + 0.25d0* &
       &            (eb(EZ,i,j,k,1) + eb(EZ,i,j,k-1,1) + &
       &             eb(EZ,i,j,k,3) + eb(EZ,i,j,k-1,3))* &
       &            (eb(EZ,i,j,k,1) + eb(EZ,i,j,k-1,1) + &
       &             eb(EZ,i,j,k,3) + eb(EZ,i,j,k-1,3)) )
        eas1 = sqrt(0.25d0* &
       &            (eb(EX,i,j,k,1) + eb(EX,i-1,j,k,1))* &
       &            (eb(EX,i,j,k,1) + eb(EX,i-1,j,k,1)) &
       &          + 0.25d0* &
       &            (eb(EY,i,j,k,1) + eb(EY,i,j-1,k,1))* &
       &            (eb(EY,i,j,k,1) + eb(EY,i,j-1,k,1)) &
       &          + 0.25d0* &
       &            (eb(EZ,i,j,k,1) + eb(EZ,i,j,k-1,1))* &
       &            (eb(EZ,i,j,k,1) + eb(EZ,i,j,k-1,1)) )
        eas2 = sqrt(0.25d0* &
       &            (eb(EX,i,j,k,3) + eb(EX,i-1,j,k,3))* &
       &            (eb(EX,i,j,k,3) + eb(EX,i-1,j,k,3)) &
       &          + 0.25d0* &
       &            (eb(EY,i,j,k,3) + eb(EY,i,j-1,k,3))* &
       &            (eb(EY,i,j,k,3) + eb(EY,i,j-1,k,3)) &
       &          + 0.25d0* &
       &            (eb(EZ,i,j,k,3) + eb(EZ,i,j,k-1,3))* &
       &            (eb(EZ,i,j,k,3) + eb(EZ,i,j,k-1,3)) )
!
        bai0 = sqrt(0.0625d0* &
       &            (eb(BX,i,j,k  ,5) + eb(BX,i,j-1,k  ,5) + &
       &             eb(BX,i,j,k-1,5) + eb(BX,i,j-1,k-1,5))* &
       &            (eb(BX,i,j,k  ,5) + eb(BX,i,j-1,k  ,5) + &
       &             eb(BX,i,j,k-1,5) + eb(BX,i,j-1,k-1,5)) &
       &          + 0.0625d0* &
       &            (eb(BY,i,j,k  ,5) + eb(BY,i-1,j,k  ,5) + &
       &             eb(BY,i,j,k-1,5) + eb(BY,i-1,j,k-1,5))* &
       &            (eb(BY,i,j,k  ,5) + eb(BY,i-1,j,k  ,5) + &
       &             eb(BY,i,j,k-1,5) + eb(BY,i-1,j,k-1,5)) &
       &          + 0.0625d0* &
       &            (eb(BZ,i,j  ,k,5) + eb(BZ,i-1,j  ,k,5) + &
       &             eb(BZ,i,j-1,k,5) + eb(BZ,i-1,j-1,k,5))* &
       &            (eb(BZ,i,j  ,k,5) + eb(BZ,i-1,j  ,k,5) + &
       &             eb(BZ,i,j-1,k,5) + eb(BZ,i-1,j-1,k,5)) )
        bas0 = sqrt(0.0625d0* &
       &            (eb(BX,i,j,k  ,1) + eb(BX,i,j-1,k  ,1) + &
       &             eb(BX,i,j,k-1,1) + eb(BX,i,j-1,k-1,1) + &
       &             eb(BX,i,j,k  ,3) + eb(BX,i,j-1,k  ,3) + &
       &             eb(BX,i,j,k-1,3) + eb(BX,i,j-1,k-1,3))* &
       &            (eb(BX,i,j,k  ,1) + eb(BX,i,j-1,k  ,1) + &
       &             eb(BX,i,j,k-1,1) + eb(BX,i,j-1,k-1,1) + &
       &             eb(BX,i,j,k  ,3) + eb(BX,i,j-1,k  ,3) + &
       &             eb(BX,i,j,k-1,3) + eb(BX,i,j-1,k-1,3)) &
       &          + 0.0625d0* &
       &            (eb(BY,i,j,k  ,1) + eb(BY,i-1,j,k  ,1) + &
       &             eb(BY,i,j,k-1,1) + eb(BY,i-1,j,k-1,1) + &
       &             eb(BY,i,j,k  ,3) + eb(BY,i-1,j,k  ,3) + &
       &             eb(BY,i,j,k-1,3) + eb(BY,i-1,j,k-1,3))* &
       &            (eb(BY,i,j,k  ,1) + eb(BY,i-1,j,k  ,1) + &
       &             eb(BY,i,j,k-1,1) + eb(BY,i-1,j,k-1,1) + &
       &             eb(BY,i,j,k  ,3) + eb(BY,i-1,j,k  ,3) + &
       &             eb(BY,i,j,k-1,3) + eb(BY,i-1,j,k-1,3)) &
       &          + 0.0625d0* &
       &            (eb(BZ,i,j  ,k,1) + eb(BZ,i-1,j  ,k,1) + &
       &             eb(BZ,i,j-1,k,1) + eb(BZ,i-1,j-1,k,1) + &
       &             eb(BZ,i,j  ,k,3) + eb(BZ,i-1,j  ,k,3) + &
       &             eb(BZ,i,j-1,k,3) + eb(BZ,i-1,j-1,k,3))* &
       &            (eb(BZ,i,j  ,k,1) + eb(BZ,i-1,j  ,k,1) + &
       &             eb(BZ,i,j-1,k,1) + eb(BZ,i-1,j-1,k,1) + &
       &             eb(BZ,i,j  ,k,3) + eb(BZ,i-1,j  ,k,3) + &
       &             eb(BZ,i,j-1,k,3) + eb(BZ,i-1,j-1,k,3)) )
        bas1 = sqrt(0.0625d0* &
       &            (eb(BX,i,j,k  ,1) + eb(BX,i,j-1,k  ,1) + &
       &             eb(BX,i,j,k-1,1) + eb(BX,i,j-1,k-1,1))* &
       &            (eb(BX,i,j,k  ,1) + eb(BX,i,j-1,k  ,1) + &
       &             eb(BX,i,j,k-1,1) + eb(BX,i,j-1,k-1,1)) &
       &          + 0.0625d0* &
       &            (eb(BY,i,j,k  ,1) + eb(BY,i-1,j,k  ,1) + &
       &             eb(BY,i,j,k-1,1) + eb(BY,i-1,j,k-1,1))* &
       &            (eb(BY,i,j,k  ,1) + eb(BY,i-1,j,k  ,1) + &
       &             eb(BY,i,j,k-1,1) + eb(BY,i-1,j,k-1,1)) &
       &          + 0.0625d0* &
       &            (eb(BZ,i,j  ,k,1) + eb(BZ,i-1,j  ,k,1) + &
       &             eb(BZ,i,j-1,k,1) + eb(BZ,i-1,j-1,k,1))* &
       &            (eb(BZ,i,j  ,k,1) + eb(BZ,i-1,j  ,k,1) + &
       &             eb(BZ,i,j-1,k,1) + eb(BZ,i-1,j-1,k,1)) )
        bas2 = sqrt(0.0625d0* &
       &            (eb(BX,i,j,k  ,3) + eb(BX,i,j-1,k  ,3) + &
       &             eb(BX,i,j,k-1,3) + eb(BX,i,j-1,k-1,3))* &
       &            (eb(BX,i,j,k  ,3) + eb(BX,i,j-1,k  ,3) + &
       &             eb(BX,i,j,k-1,3) + eb(BX,i,j-1,k-1,3)) &
       &          + 0.0625d0* &
       &            (eb(BY,i,j,k  ,3) + eb(BY,i-1,j,k  ,3) + &
       &             eb(BY,i,j,k-1,3) + eb(BY,i-1,j,k-1,3))* &
       &            (eb(BY,i,j,k  ,3) + eb(BY,i-1,j,k  ,3) + &
       &             eb(BY,i,j,k-1,3) + eb(BY,i-1,j,k-1,3)) &
       &          + 0.0625d0* &
       &            (eb(BZ,i,j  ,k,3) + eb(BZ,i-1,j  ,k,3) + &
       &             eb(BZ,i,j-1,k,3) + eb(BZ,i-1,j-1,k,3))* &
       &            (eb(BZ,i,j  ,k,3) + eb(BZ,i-1,j  ,k,3) + &
       &             eb(BZ,i,j-1,k,3) + eb(BZ,i-1,j-1,k,3)) )
        bzi0 = 0.25d0* &
       &            (eb(BZ,i,j  ,k,5) + eb(BZ,i-1,j  ,k,5) + &
       &             eb(BZ,i,j-1,k,5) + eb(BZ,i-1,j-1,k,5))
        bzs1 = 0.25d0* &
       &            (eb(BZ,i,j  ,k,1) + eb(BZ,i-1,j  ,k,1) + &
       &             eb(BZ,i,j-1,k,1) + eb(BZ,i-1,j-1,k,1))
        bzs2 = 0.25d0* &
       &            (eb(BZ,i,j  ,k,3) + eb(BZ,i-1,j  ,k,3) + &
       &             eb(BZ,i,j-1,k,3) + eb(BZ,i-1,j-1,k,3))
        bzs0 = bzs1 + bzs2
!
        write(701,*) eai0/abs(rene)
        write(702,*) eas0/abs(rene)
        write(703,*) eas1/abs(rene)
        write(704,*) eas2/abs(rene)
        write(705,*) bai0/abs(renb)
        write(706,*) bas0/abs(renb)
        write(707,*) bas1/abs(renb)
        write(708,*) bas2/abs(renb)
        write(709,*) bzi0/renb
        write(710,*) bzs0/renb
        write(711,*) bzs1/renb
        write(712,*) bzs2/renb
!
        if(istep.eq.nstep) then
          close(701)
          close(702)
          close(703)
          close(704)
          close(705)
          close(706)
          close(707)
          close(708)
          close(709)
          close(710)
          close(711)
          close(712)
        end if
      end if


!-------------------- 
      i = 2560 - 1000 - xl
      j = 2560 - yl
      k = 0 - zl
      if(i.ge.0.and.i.lt.xu-xl.and. &
         j.ge.0.and.j.lt.yu-yl.and. &
         k.ge.0.and.k.lt.zu-zl) then
        if(istep.eq.0) then
          open(801,file='eabs_i0_1000',status='replace')
          open(802,file='eabs_s0_1000',status='replace')
          open(803,file='eabs_s1_1000',status='replace')
          open(804,file='eabs_s2_1000',status='replace')
          open(805,file='babs_i0_1000',status='replace')
          open(806,file='babs_s0_1000',status='replace')
          open(807,file='babs_s1_1000',status='replace')
          open(808,file='babs_s2_1000',status='replace')
          open(809,file='bz_i0_1000',status='replace')
          open(810,file='bz_s0_1000',status='replace')
          open(811,file='bz_s1_1000',status='replace')
          open(812,file='bz_s2_1000',status='replace')
        end if
!
        eai0 = sqrt(0.25d0* &
       &            (eb(EX,i,j,k,5) + eb(EX,i-1,j,k,5))* &
       &            (eb(EX,i,j,k,5) + eb(EX,i-1,j,k,5)) &
       &          + 0.25d0* &
       &            (eb(EY,i,j,k,5) + eb(EY,i,j-1,k,5))* &
       &            (eb(EY,i,j,k,5) + eb(EY,i,j-1,k,5)) &
       &          + 0.25d0* &
       &            (eb(EZ,i,j,k,5) + eb(EZ,i,j,k-1,5))* &
       &            (eb(EZ,i,j,k,5) + eb(EZ,i,j,k-1,5)) )
        eas0 = sqrt(0.25d0* &
       &            (eb(EX,i,j,k,1) + eb(EX,i-1,j,k,1) + &
       &             eb(EX,i,j,k,3) + eb(EX,i-1,j,k,3))* &
       &            (eb(EX,i,j,k,1) + eb(EX,i-1,j,k,1) + &
       &             eb(EX,i,j,k,3) + eb(EX,i-1,j,k,3)) &
       &          + 0.25d0* &
       &            (eb(EY,i,j,k,1) + eb(EY,i,j-1,k,1) + &
       &             eb(EY,i,j,k,3) + eb(EY,i,j-1,k,3))* &
       &            (eb(EY,i,j,k,1) + eb(EY,i,j-1,k,1) + &
       &             eb(EY,i,j,k,3) + eb(EY,i,j-1,k,3)) &
       &          + 0.25d0* &
       &            (eb(EZ,i,j,k,1) + eb(EZ,i,j,k-1,1) + &
       &             eb(EZ,i,j,k,3) + eb(EZ,i,j,k-1,3))* &
       &            (eb(EZ,i,j,k,1) + eb(EZ,i,j,k-1,1) + &
       &             eb(EZ,i,j,k,3) + eb(EZ,i,j,k-1,3)) )
        eas1 = sqrt(0.25d0* &
       &            (eb(EX,i,j,k,1) + eb(EX,i-1,j,k,1))* &
       &            (eb(EX,i,j,k,1) + eb(EX,i-1,j,k,1)) &
       &          + 0.25d0* &
       &            (eb(EY,i,j,k,1) + eb(EY,i,j-1,k,1))* &
       &            (eb(EY,i,j,k,1) + eb(EY,i,j-1,k,1)) &
       &          + 0.25d0* &
       &            (eb(EZ,i,j,k,1) + eb(EZ,i,j,k-1,1))* &
       &            (eb(EZ,i,j,k,1) + eb(EZ,i,j,k-1,1)) )
        eas2 = sqrt(0.25d0* &
       &            (eb(EX,i,j,k,3) + eb(EX,i-1,j,k,3))* &
       &            (eb(EX,i,j,k,3) + eb(EX,i-1,j,k,3)) &
       &          + 0.25d0* &
       &            (eb(EY,i,j,k,3) + eb(EY,i,j-1,k,3))* &
       &            (eb(EY,i,j,k,3) + eb(EY,i,j-1,k,3)) &
       &          + 0.25d0* &
       &            (eb(EZ,i,j,k,3) + eb(EZ,i,j,k-1,3))* &
       &            (eb(EZ,i,j,k,3) + eb(EZ,i,j,k-1,3)) )
!
        bai0 = sqrt(0.0625d0* &
       &            (eb(BX,i,j,k  ,5) + eb(BX,i,j-1,k  ,5) + &
       &             eb(BX,i,j,k-1,5) + eb(BX,i,j-1,k-1,5))* &
       &            (eb(BX,i,j,k  ,5) + eb(BX,i,j-1,k  ,5) + &
       &             eb(BX,i,j,k-1,5) + eb(BX,i,j-1,k-1,5)) &
       &          + 0.0625d0* &
       &            (eb(BY,i,j,k  ,5) + eb(BY,i-1,j,k  ,5) + &
       &             eb(BY,i,j,k-1,5) + eb(BY,i-1,j,k-1,5))* &
       &            (eb(BY,i,j,k  ,5) + eb(BY,i-1,j,k  ,5) + &
       &             eb(BY,i,j,k-1,5) + eb(BY,i-1,j,k-1,5)) &
       &          + 0.0625d0* &
       &            (eb(BZ,i,j  ,k,5) + eb(BZ,i-1,j  ,k,5) + &
       &             eb(BZ,i,j-1,k,5) + eb(BZ,i-1,j-1,k,5))* &
       &            (eb(BZ,i,j  ,k,5) + eb(BZ,i-1,j  ,k,5) + &
       &             eb(BZ,i,j-1,k,5) + eb(BZ,i-1,j-1,k,5)) )
        bas0 = sqrt(0.0625d0* &
       &            (eb(BX,i,j,k  ,1) + eb(BX,i,j-1,k  ,1) + &
       &             eb(BX,i,j,k-1,1) + eb(BX,i,j-1,k-1,1) + &
       &             eb(BX,i,j,k  ,3) + eb(BX,i,j-1,k  ,3) + &
       &             eb(BX,i,j,k-1,3) + eb(BX,i,j-1,k-1,3))* &
       &            (eb(BX,i,j,k  ,1) + eb(BX,i,j-1,k  ,1) + &
       &             eb(BX,i,j,k-1,1) + eb(BX,i,j-1,k-1,1) + &
       &             eb(BX,i,j,k  ,3) + eb(BX,i,j-1,k  ,3) + &
       &             eb(BX,i,j,k-1,3) + eb(BX,i,j-1,k-1,3)) &
       &          + 0.0625d0* &
       &            (eb(BY,i,j,k  ,1) + eb(BY,i-1,j,k  ,1) + &
       &             eb(BY,i,j,k-1,1) + eb(BY,i-1,j,k-1,1) + &
       &             eb(BY,i,j,k  ,3) + eb(BY,i-1,j,k  ,3) + &
       &             eb(BY,i,j,k-1,3) + eb(BY,i-1,j,k-1,3))* &
       &            (eb(BY,i,j,k  ,1) + eb(BY,i-1,j,k  ,1) + &
       &             eb(BY,i,j,k-1,1) + eb(BY,i-1,j,k-1,1) + &
       &             eb(BY,i,j,k  ,3) + eb(BY,i-1,j,k  ,3) + &
       &             eb(BY,i,j,k-1,3) + eb(BY,i-1,j,k-1,3)) &
       &          + 0.0625d0* &
       &            (eb(BZ,i,j  ,k,1) + eb(BZ,i-1,j  ,k,1) + &
       &             eb(BZ,i,j-1,k,1) + eb(BZ,i-1,j-1,k,1) + &
       &             eb(BZ,i,j  ,k,3) + eb(BZ,i-1,j  ,k,3) + &
       &             eb(BZ,i,j-1,k,3) + eb(BZ,i-1,j-1,k,3))* &
       &            (eb(BZ,i,j  ,k,1) + eb(BZ,i-1,j  ,k,1) + &
       &             eb(BZ,i,j-1,k,1) + eb(BZ,i-1,j-1,k,1) + &
       &             eb(BZ,i,j  ,k,3) + eb(BZ,i-1,j  ,k,3) + &
       &             eb(BZ,i,j-1,k,3) + eb(BZ,i-1,j-1,k,3)) )
        bas1 = sqrt(0.0625d0* &
       &            (eb(BX,i,j,k  ,1) + eb(BX,i,j-1,k  ,1) + &
       &             eb(BX,i,j,k-1,1) + eb(BX,i,j-1,k-1,1))* &
       &            (eb(BX,i,j,k  ,1) + eb(BX,i,j-1,k  ,1) + &
       &             eb(BX,i,j,k-1,1) + eb(BX,i,j-1,k-1,1)) &
       &          + 0.0625d0* &
       &            (eb(BY,i,j,k  ,1) + eb(BY,i-1,j,k  ,1) + &
       &             eb(BY,i,j,k-1,1) + eb(BY,i-1,j,k-1,1))* &
       &            (eb(BY,i,j,k  ,1) + eb(BY,i-1,j,k  ,1) + &
       &             eb(BY,i,j,k-1,1) + eb(BY,i-1,j,k-1,1)) &
       &          + 0.0625d0* &
       &            (eb(BZ,i,j  ,k,1) + eb(BZ,i-1,j  ,k,1) + &
       &             eb(BZ,i,j-1,k,1) + eb(BZ,i-1,j-1,k,1))* &
       &            (eb(BZ,i,j  ,k,1) + eb(BZ,i-1,j  ,k,1) + &
       &             eb(BZ,i,j-1,k,1) + eb(BZ,i-1,j-1,k,1)) )
        bas2 = sqrt(0.0625d0* &
       &            (eb(BX,i,j,k  ,3) + eb(BX,i,j-1,k  ,3) + &
       &             eb(BX,i,j,k-1,3) + eb(BX,i,j-1,k-1,3))* &
       &            (eb(BX,i,j,k  ,3) + eb(BX,i,j-1,k  ,3) + &
       &             eb(BX,i,j,k-1,3) + eb(BX,i,j-1,k-1,3)) &
       &          + 0.0625d0* &
       &            (eb(BY,i,j,k  ,3) + eb(BY,i-1,j,k  ,3) + &
       &             eb(BY,i,j,k-1,3) + eb(BY,i-1,j,k-1,3))* &
       &            (eb(BY,i,j,k  ,3) + eb(BY,i-1,j,k  ,3) + &
       &             eb(BY,i,j,k-1,3) + eb(BY,i-1,j,k-1,3)) &
       &          + 0.0625d0* &
       &            (eb(BZ,i,j  ,k,3) + eb(BZ,i-1,j  ,k,3) + &
       &             eb(BZ,i,j-1,k,3) + eb(BZ,i-1,j-1,k,3))* &
       &            (eb(BZ,i,j  ,k,3) + eb(BZ,i-1,j  ,k,3) + &
       &             eb(BZ,i,j-1,k,3) + eb(BZ,i-1,j-1,k,3)) )
        bzi0 = 0.25d0* &
       &            (eb(BZ,i,j  ,k,5) + eb(BZ,i-1,j  ,k,5) + &
       &             eb(BZ,i,j-1,k,5) + eb(BZ,i-1,j-1,k,5))
        bzs1 = 0.25d0* &
       &            (eb(BZ,i,j  ,k,1) + eb(BZ,i-1,j  ,k,1) + &
       &             eb(BZ,i,j-1,k,1) + eb(BZ,i-1,j-1,k,1))
        bzs2 = 0.25d0* &
       &            (eb(BZ,i,j  ,k,3) + eb(BZ,i-1,j  ,k,3) + &
       &             eb(BZ,i,j-1,k,3) + eb(BZ,i-1,j-1,k,3))
        bzs0 = bzs1 + bzs2
!
        write(801,*) eai0/abs(rene)
        write(802,*) eas0/abs(rene)
        write(803,*) eas1/abs(rene)
        write(804,*) eas2/abs(rene)
        write(805,*) bai0/abs(renb)
        write(806,*) bas0/abs(renb)
        write(807,*) bas1/abs(renb)
        write(808,*) bas2/abs(renb)
        write(809,*) bzi0/renb
        write(810,*) bzs0/renb
        write(811,*) bzs1/renb
        write(812,*) bzs2/renb
!
        if(istep.eq.nstep) then
          close(801)
          close(802)
          close(803)
          close(804)
          close(805)
          close(806)
          close(807)
          close(808)
          close(809)
          close(810)
          close(811)
          close(812)
        end if
      end if


!-------------------- 
      i = 2560 - 2000 - xl
      j = 2560 - yl
      k = 0 - zl
      if(i.ge.0.and.i.lt.xu-xl.and. &
         j.ge.0.and.j.lt.yu-yl.and. &
         k.ge.0.and.k.lt.zu-zl) then
        if(istep.eq.0) then
          open(901,file='eabs_i0_2000',status='replace')
          open(902,file='eabs_s0_2000',status='replace')
          open(903,file='eabs_s1_2000',status='replace')
          open(904,file='eabs_s2_2000',status='replace')
          open(905,file='babs_i0_2000',status='replace')
          open(906,file='babs_s0_2000',status='replace')
          open(907,file='babs_s1_2000',status='replace')
          open(908,file='babs_s2_2000',status='replace')
          open(909,file='bz_i0_2000',status='replace')
          open(910,file='bz_s0_2000',status='replace')
          open(911,file='bz_s1_2000',status='replace')
          open(912,file='bz_s2_2000',status='replace')
        end if
!
        eai0 = sqrt(0.25d0* &
       &            (eb(EX,i,j,k,5) + eb(EX,i-1,j,k,5))* &
       &            (eb(EX,i,j,k,5) + eb(EX,i-1,j,k,5)) &
       &          + 0.25d0* &
       &            (eb(EY,i,j,k,5) + eb(EY,i,j-1,k,5))* &
       &            (eb(EY,i,j,k,5) + eb(EY,i,j-1,k,5)) &
       &          + 0.25d0* &
       &            (eb(EZ,i,j,k,5) + eb(EZ,i,j,k-1,5))* &
       &            (eb(EZ,i,j,k,5) + eb(EZ,i,j,k-1,5)) )
        eas0 = sqrt(0.25d0* &
       &            (eb(EX,i,j,k,1) + eb(EX,i-1,j,k,1) + &
       &             eb(EX,i,j,k,3) + eb(EX,i-1,j,k,3))* &
       &            (eb(EX,i,j,k,1) + eb(EX,i-1,j,k,1) + &
       &             eb(EX,i,j,k,3) + eb(EX,i-1,j,k,3)) &
       &          + 0.25d0* &
       &            (eb(EY,i,j,k,1) + eb(EY,i,j-1,k,1) + &
       &             eb(EY,i,j,k,3) + eb(EY,i,j-1,k,3))* &
       &            (eb(EY,i,j,k,1) + eb(EY,i,j-1,k,1) + &
       &             eb(EY,i,j,k,3) + eb(EY,i,j-1,k,3)) &
       &          + 0.25d0* &
       &            (eb(EZ,i,j,k,1) + eb(EZ,i,j,k-1,1) + &
       &             eb(EZ,i,j,k,3) + eb(EZ,i,j,k-1,3))* &
       &            (eb(EZ,i,j,k,1) + eb(EZ,i,j,k-1,1) + &
       &             eb(EZ,i,j,k,3) + eb(EZ,i,j,k-1,3)) )
        eas1 = sqrt(0.25d0* &
       &            (eb(EX,i,j,k,1) + eb(EX,i-1,j,k,1))* &
       &            (eb(EX,i,j,k,1) + eb(EX,i-1,j,k,1)) &
       &          + 0.25d0* &
       &            (eb(EY,i,j,k,1) + eb(EY,i,j-1,k,1))* &
       &            (eb(EY,i,j,k,1) + eb(EY,i,j-1,k,1)) &
       &          + 0.25d0* &
       &            (eb(EZ,i,j,k,1) + eb(EZ,i,j,k-1,1))* &
       &            (eb(EZ,i,j,k,1) + eb(EZ,i,j,k-1,1)) )
        eas2 = sqrt(0.25d0* &
       &            (eb(EX,i,j,k,3) + eb(EX,i-1,j,k,3))* &
       &            (eb(EX,i,j,k,3) + eb(EX,i-1,j,k,3)) &
       &          + 0.25d0* &
       &            (eb(EY,i,j,k,3) + eb(EY,i,j-1,k,3))* &
       &            (eb(EY,i,j,k,3) + eb(EY,i,j-1,k,3)) &
       &          + 0.25d0* &
       &            (eb(EZ,i,j,k,3) + eb(EZ,i,j,k-1,3))* &
       &            (eb(EZ,i,j,k,3) + eb(EZ,i,j,k-1,3)) )
!
        bai0 = sqrt(0.0625d0* &
       &            (eb(BX,i,j,k  ,5) + eb(BX,i,j-1,k  ,5) + &
       &             eb(BX,i,j,k-1,5) + eb(BX,i,j-1,k-1,5))* &
       &            (eb(BX,i,j,k  ,5) + eb(BX,i,j-1,k  ,5) + &
       &             eb(BX,i,j,k-1,5) + eb(BX,i,j-1,k-1,5)) &
       &          + 0.0625d0* &
       &            (eb(BY,i,j,k  ,5) + eb(BY,i-1,j,k  ,5) + &
       &             eb(BY,i,j,k-1,5) + eb(BY,i-1,j,k-1,5))* &
       &            (eb(BY,i,j,k  ,5) + eb(BY,i-1,j,k  ,5) + &
       &             eb(BY,i,j,k-1,5) + eb(BY,i-1,j,k-1,5)) &
       &          + 0.0625d0* &
       &            (eb(BZ,i,j  ,k,5) + eb(BZ,i-1,j  ,k,5) + &
       &             eb(BZ,i,j-1,k,5) + eb(BZ,i-1,j-1,k,5))* &
       &            (eb(BZ,i,j  ,k,5) + eb(BZ,i-1,j  ,k,5) + &
       &             eb(BZ,i,j-1,k,5) + eb(BZ,i-1,j-1,k,5)) )
        bas0 = sqrt(0.0625d0* &
       &            (eb(BX,i,j,k  ,1) + eb(BX,i,j-1,k  ,1) + &
       &             eb(BX,i,j,k-1,1) + eb(BX,i,j-1,k-1,1) + &
       &             eb(BX,i,j,k  ,3) + eb(BX,i,j-1,k  ,3) + &
       &             eb(BX,i,j,k-1,3) + eb(BX,i,j-1,k-1,3))* &
       &            (eb(BX,i,j,k  ,1) + eb(BX,i,j-1,k  ,1) + &
       &             eb(BX,i,j,k-1,1) + eb(BX,i,j-1,k-1,1) + &
       &             eb(BX,i,j,k  ,3) + eb(BX,i,j-1,k  ,3) + &
       &             eb(BX,i,j,k-1,3) + eb(BX,i,j-1,k-1,3)) &
       &          + 0.0625d0* &
       &            (eb(BY,i,j,k  ,1) + eb(BY,i-1,j,k  ,1) + &
       &             eb(BY,i,j,k-1,1) + eb(BY,i-1,j,k-1,1) + &
       &             eb(BY,i,j,k  ,3) + eb(BY,i-1,j,k  ,3) + &
       &             eb(BY,i,j,k-1,3) + eb(BY,i-1,j,k-1,3))* &
       &            (eb(BY,i,j,k  ,1) + eb(BY,i-1,j,k  ,1) + &
       &             eb(BY,i,j,k-1,1) + eb(BY,i-1,j,k-1,1) + &
       &             eb(BY,i,j,k  ,3) + eb(BY,i-1,j,k  ,3) + &
       &             eb(BY,i,j,k-1,3) + eb(BY,i-1,j,k-1,3)) &
       &          + 0.0625d0* &
       &            (eb(BZ,i,j  ,k,1) + eb(BZ,i-1,j  ,k,1) + &
       &             eb(BZ,i,j-1,k,1) + eb(BZ,i-1,j-1,k,1) + &
       &             eb(BZ,i,j  ,k,3) + eb(BZ,i-1,j  ,k,3) + &
       &             eb(BZ,i,j-1,k,3) + eb(BZ,i-1,j-1,k,3))* &
       &            (eb(BZ,i,j  ,k,1) + eb(BZ,i-1,j  ,k,1) + &
       &             eb(BZ,i,j-1,k,1) + eb(BZ,i-1,j-1,k,1) + &
       &             eb(BZ,i,j  ,k,3) + eb(BZ,i-1,j  ,k,3) + &
       &             eb(BZ,i,j-1,k,3) + eb(BZ,i-1,j-1,k,3)) )
        bas1 = sqrt(0.0625d0* &
       &            (eb(BX,i,j,k  ,1) + eb(BX,i,j-1,k  ,1) + &
       &             eb(BX,i,j,k-1,1) + eb(BX,i,j-1,k-1,1))* &
       &            (eb(BX,i,j,k  ,1) + eb(BX,i,j-1,k  ,1) + &
       &             eb(BX,i,j,k-1,1) + eb(BX,i,j-1,k-1,1)) &
       &          + 0.0625d0* &
       &            (eb(BY,i,j,k  ,1) + eb(BY,i-1,j,k  ,1) + &
       &             eb(BY,i,j,k-1,1) + eb(BY,i-1,j,k-1,1))* &
       &            (eb(BY,i,j,k  ,1) + eb(BY,i-1,j,k  ,1) + &
       &             eb(BY,i,j,k-1,1) + eb(BY,i-1,j,k-1,1)) &
       &          + 0.0625d0* &
       &            (eb(BZ,i,j  ,k,1) + eb(BZ,i-1,j  ,k,1) + &
       &             eb(BZ,i,j-1,k,1) + eb(BZ,i-1,j-1,k,1))* &
       &            (eb(BZ,i,j  ,k,1) + eb(BZ,i-1,j  ,k,1) + &
       &             eb(BZ,i,j-1,k,1) + eb(BZ,i-1,j-1,k,1)) )
        bas2 = sqrt(0.0625d0* &
       &            (eb(BX,i,j,k  ,3) + eb(BX,i,j-1,k  ,3) + &
       &             eb(BX,i,j,k-1,3) + eb(BX,i,j-1,k-1,3))* &
       &            (eb(BX,i,j,k  ,3) + eb(BX,i,j-1,k  ,3) + &
       &             eb(BX,i,j,k-1,3) + eb(BX,i,j-1,k-1,3)) &
       &          + 0.0625d0* &
       &            (eb(BY,i,j,k  ,3) + eb(BY,i-1,j,k  ,3) + &
       &             eb(BY,i,j,k-1,3) + eb(BY,i-1,j,k-1,3))* &
       &            (eb(BY,i,j,k  ,3) + eb(BY,i-1,j,k  ,3) + &
       &             eb(BY,i,j,k-1,3) + eb(BY,i-1,j,k-1,3)) &
       &          + 0.0625d0* &
       &            (eb(BZ,i,j  ,k,3) + eb(BZ,i-1,j  ,k,3) + &
       &             eb(BZ,i,j-1,k,3) + eb(BZ,i-1,j-1,k,3))* &
       &            (eb(BZ,i,j  ,k,3) + eb(BZ,i-1,j  ,k,3) + &
       &             eb(BZ,i,j-1,k,3) + eb(BZ,i-1,j-1,k,3)) )
        bzi0 = 0.25d0* &
       &            (eb(BZ,i,j  ,k,5) + eb(BZ,i-1,j  ,k,5) + &
       &             eb(BZ,i,j-1,k,5) + eb(BZ,i-1,j-1,k,5))
        bzs1 = 0.25d0* &
       &            (eb(BZ,i,j  ,k,1) + eb(BZ,i-1,j  ,k,1) + &
       &             eb(BZ,i,j-1,k,1) + eb(BZ,i-1,j-1,k,1))
        bzs2 = 0.25d0* &
       &            (eb(BZ,i,j  ,k,3) + eb(BZ,i-1,j  ,k,3) + &
       &             eb(BZ,i,j-1,k,3) + eb(BZ,i-1,j-1,k,3))
        bzs0 = bzs1 + bzs2
!
        write(901,*) eai0/abs(rene)
        write(902,*) eas0/abs(rene)
        write(903,*) eas1/abs(rene)
        write(904,*) eas2/abs(rene)
        write(905,*) bai0/abs(renb)
        write(906,*) bas0/abs(renb)
        write(907,*) bas1/abs(renb)
        write(908,*) bas2/abs(renb)
        write(909,*) bzi0/renb
        write(910,*) bzs0/renb
        write(911,*) bzs1/renb
        write(912,*) bzs2/renb
!
        if(istep.eq.nstep) then
          close(901)
          close(902)
          close(903)
          close(904)
          close(905)
          close(906)
          close(907)
          close(908)
          close(909)
          close(910)
          close(911)
          close(912)
        end if
      end if


  return
  end subroutine
