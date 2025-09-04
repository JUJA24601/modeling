!***************************************************************
! ElectroMagnetic Spacecraft Environment Simulator named EMSES
! ~      ~        ~          ~           ~
!      load-balanced by One-handed Help (OhHelp) Algorithm
!
!           controle program for emses
!
!***************************************************************
!
                  subroutine emses
!
! **** emses standard version ****
!
!---      Three Dimensional ElectroMagnetic Simulation Code     --
!   This code is based on
!        Kyoto university ElectroMagnetic Particle cOde (KEMPO)
!
!========================= History of kempo =========================
!
!       Programed by H.Matsumoto          July, 1980
!                 by H.tahara             February, 1981
!                 by I.Ayukawa            February, 1982
!                 by M.Ohashi             April, 1982
!                 by Y.Omura              September, 1982
!                 by Y.Omura, K.Fukuchi
!                          and T.Yamada   September, 1983
!           Optimized for vector processor
!                 by Y.Omura              February, 1984
!                 by Y.Omura, T.Kimura    June, 1984
!           Revised as version 7 & 8
!                 by Y.Omura and N.Komori November, 1985
!           Revised as version 9
!                 by Y.Omura & T.Tanaka   April, 1986
!           Data output format is revised as version 10
!                 by Y.Omura & K.Inagaki  April,1986
!           Velocity distribution diagnostics is added as version 11
!                 by Y.Omura              July, 1986
!
!           Expanded to three dimensional system
!                 by H.Yokoyama           January, 1992
!           Free boundary is added as version 2 (3D)
!                 by M.Yamane             May, 1993
!
!========================= History of emses =========================
!
!       Programed by Y.Miyake             October, 2008
!           OhHelp Library (version 0.9) is developed
!                 by H.Nakashima          August, 2009
!           OhHelp is applied to EMSES
!                 by Y.Miyake             July, 2010
!
!=====================================================================
!
!            Radio Atomospheric Science Center (-2004)
!        Research Institute for Sustainable Humanosphere (2004-)
!          Academic Center for Computing and Media Studies
!              Kyoto University, Kyoto, Japan.
!
!
!-------------------- parameter and common blocks
  use oh_type
  use paramt
  use allcom
#include "oh_stats.h"
! #define MCW MPI_COMM_WORLD
#define MCW CTCA_subcomm
#define MSS MPI_STATUS_SIZE
  implicit none
!
  integer(kind=4) :: ustep
  integer(kind=4) :: fajdg
  integer(kind=8) :: lsumphgram,gsumphgram
  integer(kind=4) :: mpierr
  integer(kind=4) :: nprocs,myrank
  integer(kind=4) :: m,is,ps,sd,ieb
  integer(kind=4) :: sstep_injct
  real(kind=8) :: smltime,emltime
  logical :: rebalance


!--- MPI Initialize ---
!              -------------------------------------------------
                ! call MPI_Init(mpierr)
                call CTCAR_init()
                call MPI_Comm_size(MPI_COMM_WORLD,nprocs,mpierr)
                call MPI_Comm_rank(MPI_COMM_WORLD,myrank,mpierr)
                call MPI_Comm_size(MCW,nnode,mpierr)
                call MPI_Comm_rank(MCW,myid,mpierr)
                call oh1_fam_comm(CTCA_subcomm)
                call hdfinit()

                ! if(myid.eq.0) print*, "EMSES: MCW, np, rk, Scom, nn, id", MPI_COMM_WORLD, nprocs, myrank, CTCA_subcomm, nnode, myid
!              -------------------------------------------------
!--- start ---
!                              -------------
                                call  input
          if(emflag.ge.1) then
                                emmode = .true.
          else
                                emmode = .false.
          end if
          if(emflag.ge.2) then
                                implic = .true.
          else
                                implic = .false.
          end if
!                              -------------
!--- initialization ---              *
!                              -------------
                                call ohinit
                                call inital
!                              -------------
!--- particle initialization ---     *
!                              -------------
                                call inipcl
                                call inipin
!                              -------------
!--- field initialization ---        *
!                              -------------
                                call inifld
!                              -------------
!--- medium properties ---           *
!                              -------------
                                call medium
!                              -------------
!--- loadbalancing ---               *
!                          --------------------------------
          if(nspec.gt.0) then
                            currmode = oh3_transbound(0,0)
            do ps=1,2; do is=1,nspec
              if(ps.eq.1.or.sdid(ps).ge.0) &
             &              nphgram(sdid(ps)+1,is,ps) = totalp(is,ps)
            end do; end do
            do ps=3,4
              do is=1,nspec
                            totalp(is,ps) = 0
              end do
                            pbase(ps+1) = pbase(3)
            end do
                            famind = -1; fammbr = -1;
                            call oh1_families(famind,fammbr)
!        if(myid.eq.0) print*, "fam{ind,mbr}",famind,fammbr
                            call create_nborps(1)
            if(currmode.lt.0) then
                            call create_nborps(2)
            do ieb=0,(nebfld-1)*2,2
                            call oh3_bcast_field &
                           &(eb(1,0,0,0,ieb+1),eb(1,0,0,0,ieb+2),FEB)
                            call oh3_bcast_field &
                           &(ebsc(1,0,0,0,ieb+1),ebsc(1,0,0,0,ieb+2),FSC)
            end do
                            call oh3_bcast_field &
                           &(mp(1,0,0,0,1),mp(1,0,0,0,2),FEB)
                            call oh3_bcast_field &
                           &(mp(1,0,0,0,3),mp(1,0,0,0,4),FEB)
                            call oh3_bcast_field &
                           &(colf(1,0,0,0,1),colf(1,0,0,0,2),FPH)
                            currmode = 1
            end if
                            gcount(1)%globalp(1:nspec,:) = totalp(1:nspec,:)
                            call MPI_Allreduce(gcount(1),gcount(2),1, &
                           &mpi_type_count,mpi_sum_count,MCW,mpierr)
          end if
!                          --------------------------------
!--- variable transformation ---     *
!                              -------------
                                call rescal
!                              -------------
!--- capacity matrix ---             *
!                              -------------
                                call getccl
                                call capmtx
!                              -------------
!--- charge neutrality ---           *
!                              -------------
                                call chgntr
!                              -------------
!                                    *
          if(jobnum(1).eq.1.and.jobnum(3).eq.1) then
!--- v(0.0) back to v(-0.5) ---      *
!                              -------------
                                call grdcrt(1)
                                call vpushh(1,2,-1,0,1)
              if(sdid(2).ge.0)  call grdcrt(2)
              if(sdid(2).ge.0)  call vpushh(2,2,-1,0,1)
!                              -------------
          end if
!--- loadbalancing ---               *
!    (with injected particles preloaded)
!                          --------------------------------
          if(nspec.gt.0) then
            rebalance = .false.
            do sstep_injct=1,nscycinj
                            call psuper3(1,2,sstep_injct)
                            currmode = oh3_transbound(currmode,0)
              do ps=1,2; do is=1,nspec
                if(ps.eq.1.or.sdid(ps).ge.0) &
               &            nphgram(sdid(ps)+1,is,ps) = totalp(is,ps)
              end do; end do
              do ps=3,4
                do is=1,nspec
                            totalp(is,ps) = 0
                end do
                            pbase(ps+1) = pbase(3)
              end do
                            call oh1_families(famind,fammbr)
!                            call chkres(1,0,1," chkresA1:")
!                            call chkres(2,0,1," chkresA2:")
              if(currmode.lt.0) then
                            rebalance = .true.
                            currmode = 1
              end if
            end do
            if(rebalance) then
                            call create_nborps(2)
              do ieb=0,(nebfld-1)*2,2
                            call oh3_bcast_field &
                           &(eb(1,0,0,0,ieb+1),eb(1,0,0,0,ieb+2),FEB)
                            call oh3_bcast_field &
                           &(ebsc(1,0,0,0,ieb+1),ebsc(1,0,0,0,ieb+2),FSC)
              end do
                            call oh3_bcast_field &
                           &(mp(1,0,0,0,1),mp(1,0,0,0,2),FEB)
!                            call oh3_bcast_field &
!                           &(mp(1,0,0,0,3),mp(1,0,0,0,4),FEB)
                            call oh3_bcast_field &
                           &(colf(1,0,0,0,1),colf(1,0,0,0,2),FPH)
            end if
                            gcount(1)%globalp(1:nspec,:) = totalp(1:nspec,:)
                            call MPI_Allreduce(gcount(1),gcount(2),1, &
                           &mpi_type_count,mpi_sum_count,MCW,mpierr)
          end if
!                          --------------------------------
!--- rho(0.5+m/2) ---                *
!                              -------------
            if(ifdiag.ne.0) then
                                call charge(1,1)
              if(sdid(2).ge.0)  call charge(2,1)
              if(currmode.ne.0) call oh3_reduce_field &
                               &(rho(1,0,0,0,1),rho(1,0,0,0,2),FRH)
                                call oh3_exchange_borders &
                               &(rho(1,0,0,0,1),rho(1,0,0,0,2),CRH,0)
                                call exchange_lchg_borders(rho(:,:,:,:,1),1,FRH,CRH,1)
                                call add_boundary_charge(rho(:,:,:,:,1),1,FRH,CRH,1)
              if(currmode.ne.0) call oh3_reduce_field &
                               &(rhobk(1,0,0,0,1),rhobk(1,0,0,0,2),FRB)
                                call oh3_exchange_borders &
                               &(rhobk(1,0,0,0,1),rhobk(1,0,0,0,2),CRB,0)
                                call exchange_lchg_borders(rhobk(:,:,:,:,1),1,FRB,CRB,2)
                                call add_boundary_charge(rhobk(:,:,:,:,1),1,FRB,CRB,2)
              if(currmode.ne.0) call oh3_reduce_field &
                               &(rhodg(1,0,0,0,1),rhodg(1,0,0,0,2),FRD)
                                call oh3_exchange_borders &
                               &(rhodg(1,0,0,0,1),rhodg(1,0,0,0,2),CRD,0)
                                call exchange_lchg_borders(rhodg(:,:,:,:,1),1,FRD,CRD,nspec*2)
                                call add_boundary_charge(rhodg(:,:,:,:,1),1,FRD,CRD,nspec*2)
            else
                                call charge(1,0)
              if(sdid(2).ge.0)  call charge(2,0)
              if(currmode.ne.0) call oh3_reduce_field &
                               &(rho(1,0,0,0,1),rho(1,0,0,0,2),FRH)
                                call oh3_exchange_borders &
                               &(rho(1,0,0,0,1),rho(1,0,0,0,2),CRH,0)
                                call exchange_lchg_borders(rho(:,:,:,:,1),1,FRH,CRH,1)
                                call add_boundary_charge(rho(:,:,:,:,1),1,FRH,CRH,1)
              if(currmode.ne.0) call oh3_reduce_field &
                               &(rhobk(1,0,0,0,1),rhobk(1,0,0,0,2),FRB)
                                call oh3_exchange_borders &
                               &(rhobk(1,0,0,0,1),rhobk(1,0,0,0,2),CRB,0)
                                call exchange_lchg_borders(rhobk(:,:,:,:,1),1,FRB,CRB,2)
                                call add_boundary_charge(rhobk(:,:,:,:,1),1,FRB,CRB,2)
            end if
                                rhobk(:,:,:,:,3) = rhobk(:,:,:,:,3) &
                               &                 + rhobk(:,:,:,:,1)
                                call add_background_charge(1)
!                              -------------
!--- charge density masking ---      *
!                              -------------
              if(imask(5).eq.1) call fsmask(1,4)
!                              -------------
!--- e(0.5+m/2):ES ---               *
!                              -------------
              if(npc.ge.1)      call surchg
                                call esfld1
              if(currmode.ne.0) call oh3_bcast_field &
                               &(phi(1,0,0,0,1),phi(1,0,0,0,2),FPH)
                                call esfld3(1)
              if(sdid(2).ge.0)  call esfld3(2)
!                              -------------
!--- define external wave field ---  *
!                              -------------
                                call extfld(1)
              if(sdid(2).ge.0)  call extfld(2)
!                              -------------
!--- correction of em in body ---    *
!                              -------------
                                call bdyfld(1)
              if(sdid(2).ge.0)  call bdyfld(2)
!                              -------------
!--- exchanging boundary data ---    *
!                              -------------
            do ieb=0,(nebfld-1)*2,2
                                call oh3_exchange_borders &
                               &(eb(1,0,0,0,ieb+1),eb(1,0,0,0,ieb+2),CEB,currmode)
                                call boundary_emfld(ieb+1)
              if(sdid(2).ge.0)  call boundary_emfld(ieb+2)
            end do
!                              -------------
!--- check of parameters ---         *
!                              -------------
                                call chkprm
!                              -------------
!--- electrostatic potential ---     *
!                              -------------
!          if(nspec.ne.0)        call phispc(0)
!                              -------------
!--- energy diagnostics---           *
!                              -------------
!                                call energy
!                              -------------
!--- particle sort ---               *
!      ======================  -------------
            if(isort(1).ne.0) then
                                call spsort(1)
              if(sdid(2).ge.0)  call spsort(2)
            end if
!      ======================  -------------
!--- initial diagnostics ---         *
!                              -------------
                                call MPI_Barrier(MCW,mpierr)
                                call frmttd(1)
                                call digset
                                call hdfdig(0,0)
                                call hdfdig(1,0)
!                                call hdfdig(2,0)
!                              -------------
!                                    *
!--- main loop ---                   * ::::::::::<:::::::::
!                                    *                    :
                            mnmlocal = 0
                            mnmtotal = 0
                            call MPI_Barrier(MCW,mpierr)
!                            call gettod(smltime)
!                            smltime = smltime*1.0d-6
                            smltime = MPI_Wtime()
!                                    *                    :
!                        **************************       :
                  MAINL: do istep=1,nstep
!                        **************************       :
!--- time ---                        *                    :
!                              +++++++++++++              :
                                 t = t+dt
                               itime = istep
      if(myid.eq.0) write(6,*) '**** step --------- ',itime
!                         if(istep.ne.nstep) then
                                 ustep = 2
!                         else
!                                 ustep = 1
!                         end if
                 if(ijdiag.eq.0.or.istep.lt.hdfdigstart) then
                                 fajdg = 0
                 else
                                 fajdg = 1
                 end if
!                              +++++++++++++              :
!--- relocation of field pe and pb   *                    :
!                              -------------              :
                                call grdcrt(1)
!                              -------------              :
!--- v(n-0.5) to v(n+0.5) ---        *                    :
!--- r(n) to r(n+1) ---              *                    :
!--- j(n+0.5)_main ---               *                    :
!                              -------------              :
                                call psolve1(1,ustep,fajdg)
!                                call psolve1(3,ustep,fajdg)
!                              -------------              :
              if(sdid(2).ge.0) then
!--- relocation of field pe and pb   *                    :
!                              -------------              :
                                call grdcrt(2)
!                              -------------              :
!--- v(n-0.5) to v(n+0.5) ---        *                    :
!--- r(n) to r(n+1) ---              *                    :
!--- j(n+0.5)_main ---               *                    :
!                              -------------              :
                                call psolve1(2,ustep,fajdg)
!                              -------------              :
              end if
!--- j(n+0.5)_sub ---                *                    :
!      ======================  -------------              :
              if(currmode.ne.0) call oh3_allreduce_field &
                               &(aj(1,0,0,0,1),aj(1,0,0,0,2),FAJ)
                                call oh3_exchange_borders &
                               &(aj(1,0,0,0,1),aj(1,0,0,0,2),CAJ,currmode)
                                call exchange_lcur_borders1(1)
              if(sdid(2).ge.0)  call exchange_lcur_borders1(2)
                                call add_boundary_current1(1)
              if(sdid(2).ge.0)  call add_boundary_current1(2)
              if(currmode.ne.0) call oh3_allreduce_field &
                               &(ajdg(1,0,0,0,1),ajdg(1,0,0,0,2),FJD)
                                call oh3_exchange_borders &
                               &(ajdg(1,0,0,0,1),ajdg(1,0,0,0,2),CJD,currmode)
                                call exchange_lcur_borders2(1)
              if(sdid(2).ge.0)  call exchange_lcur_borders2(2)
                                call add_boundary_current2(1)
              if(sdid(2).ge.0)  call add_boundary_current2(2)
!      ======================  -------------              :

!--- current density masking ---     *                    :
!                              -------------              :
            if(imask(4).eq.1) then
                                call fsmask(1,3)
              if(sdid(2).ge.0)  call fsmask(2,3)
            end if
!                              -------------              :
!--- b(n) to b(n+0.5) ---            *                    :
!      ======================  -------------              :
            do ieb=0,(nebfld-1)*2,2
                                call bfield(ieb+1,1)
              if(sdid(2).ge.0)  call bfield(ieb+2,1)
            end do
!      ======================  -------------              :
!--- b-field masking          ---    *                    :
!                              -------------              :
            if(imask(1).eq.1) then
                                call fsmask(1,0)
              if(sdid(2).ge.0)  call fsmask(2,0)
            end if
!                              -------------              :
!--- e(n) to e(n+1) ---              *                    :
!      ======================  -------------              :
            do ieb=0,(nebfld-1)*2,2
                                call efield(ieb+1)
              if(sdid(2).ge.0)  call efield(ieb+2)
            end do
!      ======================  -------------              :
!--- exchanging boundary data ---                         :
!                              -------------              :
            do ieb=0,(nebfld-1)*2,2
                                call oh3_exchange_borders &
                               &(eb(1,0,0,0,ieb+1),eb(1,0,0,0,ieb+2),CEB,currmode)
                                call boundary_emfld(ieb+1)
              if(sdid(2).ge.0)  call boundary_emfld(ieb+2)
            end do
!                              -------------              :
!--- define external wave field ---  *                    :
!                              -------------              :
!                                call extfld(1)
!              if(sdid(2).ge.0)  call extfld(2)
!                              -------------              :
!--- correction of em in body ---    *                    :
!                              -------------              :
!                                call bdyfld(1)
!              if(sdid(2).ge.0)  call bdyfld(2)
!                              -------------              :
!--- loadbalancing ---               *                    :
!    (with injected particles preloaded)                  :
!                          --------------------------------
          if(nspec.gt.0) then
            rebalance = .false.
            do sstep_injct=1,nscycinj
                            call psuper3(1,2,sstep_injct)
                            currmode = oh3_transbound(currmode,1)
              do ps=1,2; do is=1,nspec
                if(ps.eq.1.or.sdid(ps).ge.0) &
               &            nphgram(sdid(ps)+1,is,ps) = totalp(is,ps)
              end do; end do
              do ps=3,4
                do is=1,nspec
                            totalp(is,ps) = 0
                end do
                            pbase(ps+1) = pbase(3)
              end do
                            call oh1_families(famind,fammbr)
!                            call chkres(1,0,1," chkresA3:")
!                            call chkres(2,0,1," chkresA4:")
              if(currmode.lt.0) then
                            rebalance = .true.
                            currmode = 1
              end if
            end do
            if(rebalance) then
                            call create_nborps(2)
              do ieb=0,(nebfld-1)*2,2
                            call oh3_bcast_field &
                           &(eb(1,0,0,0,ieb+1),eb(1,0,0,0,ieb+2),FEB)
                            call oh3_bcast_field &
                           &(ebsc(1,0,0,0,ieb+1),ebsc(1,0,0,0,ieb+2),FSC)
              end do
                            call oh3_bcast_field &
                           &(mp(1,0,0,0,1),mp(1,0,0,0,2),FEB)
!                            call oh3_bcast_field &
!                           &(mp(1,0,0,0,3),mp(1,0,0,0,4),FEB)
                            call oh3_bcast_field &
                           &(colf(1,0,0,0,1),colf(1,0,0,0,2),FPH)
            end if
                            gcount(1)%globalp(1:nspec,:) = totalp(1:nspec,:)
                            mnmlocal = mnmlocal + sum(totalp(:,:))
                            call MPI_Allreduce(gcount(1),gcount(2),1, &
                           &mpi_type_count,mpi_sum_count,MCW,mpierr)
          end if
!                          --------------------------------
!--- particle sort ---               *                    :
!      ======================  -------------              :
            if(isort(1).ne.0.and.mod(istep,isort(1)).eq.0.and. &
           &   istep.ne.nstep) then
                                call spsort(1)
              if(sdid(2).ge.0)  call spsort(2)
            end if
!      ======================  -------------              :
!--- rho(n+1) ---                    *                    :
!                              -------------              :
            if(ifdiag.ne.0.and.(mod(istep,ifdiag).eq.0.or.daverg.ge.1)) then
                                call charge(1,1)
              if(sdid(2).ge.0)  call charge(2,1)
              if(currmode.ne.0) call oh3_reduce_field &
                               &(rho(1,0,0,0,1),rho(1,0,0,0,2),FRH)
                                call oh3_exchange_borders &
                               &(rho(1,0,0,0,1),rho(1,0,0,0,2),CRH,0)
                                call exchange_lchg_borders(rho(:,:,:,:,1),1,FRH,CRH,1)
                                call add_boundary_charge(rho(:,:,:,:,1),1,FRH,CRH,1)
              if(currmode.ne.0) call oh3_reduce_field &
                               &(rhobk(1,0,0,0,1),rhobk(1,0,0,0,2),FRB)
                                call oh3_exchange_borders &
                               &(rhobk(1,0,0,0,1),rhobk(1,0,0,0,2),CRB,0)
                                call exchange_lchg_borders(rhobk(:,:,:,:,1),1,FRB,CRB,2)
                                call add_boundary_charge(rhobk(:,:,:,:,1),1,FRB,CRB,2)
              if(currmode.ne.0) call oh3_reduce_field &
                               &(rhodg(1,0,0,0,1),rhodg(1,0,0,0,2),FRD)
                                call oh3_exchange_borders &
                               &(rhodg(1,0,0,0,1),rhodg(1,0,0,0,2),CRD,0)
                                call exchange_lchg_borders(rhodg(:,:,:,:,1),1,FRD,CRD,nspec*2)
                                call add_boundary_charge(rhodg(:,:,:,:,1),1,FRD,CRD,nspec*2)
            else
                                call charge(1,0)
              if(sdid(2).ge.0)  call charge(2,0)
              if(currmode.ne.0) call oh3_reduce_field &
                               &(rho(1,0,0,0,1),rho(1,0,0,0,2),FRH)
                                call oh3_exchange_borders &
                               &(rho(1,0,0,0,1),rho(1,0,0,0,2),CRH,0)
                                call exchange_lchg_borders(rho(:,:,:,:,1),1,FRH,CRH,1)
                                call add_boundary_charge(rho(:,:,:,:,1),1,FRH,CRH,1)
              if(currmode.ne.0) call oh3_reduce_field &
                               &(rhobk(1,0,0,0,1),rhobk(1,0,0,0,2),FRB)
                                call oh3_exchange_borders &
                               &(rhobk(1,0,0,0,1),rhobk(1,0,0,0,2),CRB,0)
                                call exchange_lchg_borders(rhobk(:,:,:,:,1),1,FRB,CRB,2)
                                call add_boundary_charge(rhobk(:,:,:,:,1),1,FRB,CRB,2)
            end if
!                                call rhoupd(1)
                                rhobk(:,:,:,:,3) = rhobk(:,:,:,:,3) &
                               &                 + rhobk(:,:,:,:,1)
                                call add_background_charge(1)
!                              -------------              :
!--- charge density masking ---      *
!                              -------------
            if(imask(5).eq.1) then
                                call fsmask(1,4)
            end if
!                              -------------
!--- correction of e ---             *                    :
!                              -------------              :
              if(npc.ge.1) then
                                call surchg
              end if
                                call esfld1
              if(currmode.ne.0) call oh3_bcast_field &
                               &(phi(1,0,0,0,1),phi(1,0,0,0,2),FPH)
                                call esfld3(1)
              if(sdid(2).ge.0)  call esfld3(2)
!                              -------------              :
!--- electrostatic potential ---     *                    :
!                              -------------              :
!!            if(ifdiag.ne.0.and.mod(istep,ifdiag).eq.0) then
!!                                call esfld1(0)
!!            end if
!                              -------------              :
!--- define external wave field ---  *                    :
!                              -------------              :
                                call extfld(1)
              if(sdid(2).ge.0)  call extfld(2)
!                              -------------              :
!--- correction of em in body ---    *                    :
!                              -------------              :
                                call bdyfld(1)
              if(sdid(2).ge.0)  call bdyfld(2)
!                              -------------              :
!--- e-field masking          ---    *                    :
!                              -------------              :
            if(imask(2).eq.1) then
                                call fsmask(1,1)
              if(sdid(2).ge.0)  call fsmask(2,1)
            end if
            if(imask(3).eq.1) then
                                call fsmask(1,2)
              if(sdid(2).ge.0)  call fsmask(2,2)
            end if
!                              -------------              :
!--- exchanging boundary data ---                         :
!                              -------------              :
            do ieb=0,(nebfld-1)*2,2
                                call oh3_exchange_borders &
                               &(eb(1,0,0,0,ieb+1),eb(1,0,0,0,ieb+2),CEB,currmode)
                                call boundary_emfld(ieb+1)
              if(sdid(2).ge.0)  call boundary_emfld(ieb+2)
            end do
!                              -------------              :
!--- b(n-0.5) to b(n) ---            *                    :
!                              -------------              :
            do ieb=0,(nebfld-1)*2,2
                                call bfield(ieb+1,2)
              if(sdid(2).ge.0)  call bfield(ieb+2,2)
            end do
!                              -------------              :
!--- diagnostics ---                 *                    :
!                              -------------              :
                                call frmttd(ustep)
!                              -------------              :
!                                    *                    :
!                              -------------              :
!                                call hdfdig(0,istep)
!                              -------------              :
!                                    *                    :
!                              -------------              :
                                call hdfdig(1,istep)
!                              -------------              :
!                                    *                    :
!                              -------------              :
!                                call hdfdig(2,istep)
!                              -------------              :
!                                    *                    :
!                              -------------              :
!                                call energy
!                              -------------              :
!                                    *                    :
!      ================        -------------              :
!        if(neutr.eq.1)          call vcrrct
!      ================        -------------              :
!                                    *                    :
!      ================         >>>>>>>>>>>               :
!        if(icond.ne.0)           go to 999
!      ================         <<<<<<<<<<<               :
!                                    *                    :
!                              *************              :
                                end do MAINL
!                              *************              :
                            call MPI_Barrier(MCW,mpierr)
!                            call gettod(emltime)
!                            emltime=emltime*1.0d-6
                            emltime = MPI_Wtime()
!                                    *                    :
!                                    * ::::::::::>:::::::::
!                                    *
!                               +++++++++++
!                               lstep=istep-1
!                               +++++++++++
!                                    *
!                              -------------
          if(jobnum(2).eq.1.and.jobnum(3).eq.1) then
                                call grdcrt(1)
                                call vpushh(1,2,+1,1,0)
              if(sdid(2).ge.0)  call grdcrt(2)
              if(sdid(2).ge.0)  call vpushh(2,2,+1,1,0)
          end if
!                              -------------
!                                    *
!                          ---------------------------------------
                            currmode = oh3_transbound(currmode,1)
!                          --------------------------------------
!                                    *
!                              -------------
                                call hdfdig(2,istep)
!                              -------------
!                                    *
!--- system call ---                 *
        if(jobnum(2).eq.1) then
!                              ------------- 
          if(myid.eq.0) &
         &                      call system('mkdir SNAPSHOT1')
!                              -------------
!--- save data for next job ---      *
!                              -------------
                                call MPI_Barrier(MCW,mpierr)
                                call save_emsnap
!                              -------------
        end if
!                                    *
!                              -------------
!!!!                                call chkcst
                                call MPI_Barrier(MCW,mpierr)
 999                              continue
!                       print*,"time main loop", emltime-smltime
                        call MPI_Reduce(mnmlocal,mnmtotal,1, &
                       &                MPI_INTEGER8,MPI_SUM,0,MCW,mpierr)
                        if (myid.eq.0) then
                          print*,"time main loop", (emltime-smltime)
                          print*,"total N. of processed part.", mnmtotal
                          print*, &
                       &    "performance   ", &
                       &    dble(mnmtotal)/(emltime-smltime)*1.0d-6, &
                       &    "Mpart/s"
                        end if
!                              -------------
!                                    *
!                       ---------------------------
                         call hdffinalize()
                        !  call MPI_Finalize(mpierr)
                         call CTCAR_finalize()
!                       ---------------------------
!                                    *
!                                *********
                                   return
!                                *********
!                                    *
!                                 *******
                            end subroutine emses
