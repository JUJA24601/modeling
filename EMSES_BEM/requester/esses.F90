!***************************************************************
! ElectroMagnetic Spacecraft Environment Simulator named EMSES
! ~      ~        ~          ~           ~
!      load-balanced by One-handed Help (OhHelp) Algorithm
!
!           controle program for emses
!
!***************************************************************
!
                  subroutine esses
!
! **** "esses" is an electrostatic version of emses ****
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
!           OhHelp is applied to ESSES
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
  use bem
#include "oh_stats.h"
! #define MCW MPI_COMM_WORLD
#define MCW CTCA_subcomm
#define MSS MPI_STATUS_SIZE
  implicit none
!
  integer(kind=4) :: ustep
  integer(kind=4) :: fajdg
  integer(kind=4) :: mpierr
  integer(kind=4) :: nprocs,myrank
  integer(kind=4) :: is,ps,m
  integer(kind=4) :: sstep_injct
  real(kind=8) :: smltime,emltime
  real(kind=8) :: mltime=0.0d0,estime=0.0d0
  logical :: rebalance
  
  integer(kind=4) :: area_x, area_y, area_z, grid_id, x_grid, y_grid, z_grid, x_tick, y_tick, z_tick, i_bem, j_bem, k_bem, l_bem
  real(kind=8), allocatable :: phisnap_local(:,:,:), phisnap_root(:,:,:)
  character(len=100) :: filename

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
!--- CoToCoA initialization ---      *
!                              --------------
                                ! ファイル名
                                file_triangle = "./data_triangle.csv"
                                file_rectangle = "./data_rectangle.csv"
                                inquire(file=file_triangle, exist=triangle_exists)
                                inquire(file=file_rectangle, exist=rectangle_exists)
                                if (triangle_exists) call read_triangle(file_triangle)
                                if (rectangle_exists) call read_rectangle(file_rectangle)
                                num_of_faces = num_of_rectangles + num_of_triangles
                                ! print *, "num_of_faces: ", num_of_faces
                                allocate(coef(-1:nx/nodes(1)+1,-1:ny/nodes(2)+1,-1:nz/nodes(3)+1, num_of_faces))
                                allocate(omegas(num_of_faces))

                                ! 電位を計算するときの係数を計算しておく
                                do i_bem = -1, nx/nodes(1)+1
                                  do j_bem = -1, ny/nodes(2)+1
                                    do k_bem = -1, nz/nodes(3)+1
                                      do l_bem = 1, num_of_faces
                                        if (l_bem<=num_of_rectangles) then
                                          coef(i_bem,j_bem,k_bem,l_bem) = potential_elem_rect(                     &
                                            rectangles(l_bem),                                                     &
                                            [                                                                      &
                                              real(i_bem+mod(myid,nodes(1))*nx/nodes(1),8),                        &
                                              real(j_bem+mod(myid/nodes(1),nodes(2))*ny/nodes(2),8),               &
                                              real(k_bem+(myid/(nodes(1)*nodes(2)))*nz/nodes(3),8)                 & 
                                              ],                                                                   &
                                            1.0d0                                                                  &
                                            )
                                        else
                                          coef(i_bem,j_bem,k_bem,l_bem) = potential_elem_tri(                      &
                                            triangles(l_bem-num_of_rectangles),                                    &
                                            [                                                                      &
                                              real(i_bem+mod(myid,nodes(1))*nx/nodes(1),8),                        &
                                              real(j_bem+mod(myid/nodes(1),nodes(2))*ny/nodes(2),8),               &
                                              real(k_bem+(myid/(nodes(1)*nodes(2)))*nz/nodes(3),8)                 & 
                                              ],                                                                   &
                                            1.0d0                                                                  &
                                            )
                                        end if
                                      end do
                                    end do
                                  end do
                                end do
                                
                                ! print *, grid_point_rows
                                if (myid == 0) then
                                  call calc_gridpoint()
                                end if
                                call MPI_Bcast(grid_point_rows, 1, MPI_INTEGER4, 0, MCW, mpierr)
                                print *, "grid_point_rows: ", grid_point_rows

                                allocate(rtowdat(grid_point_rows+2), source=0.0d0)
                                allocate(phi_local(grid_point_rows), source=0.0d0)
                                allocate(phi_root(grid_point_rows), source=0.0d0)
                                allocate(phisnap_local(nx+1,ny+1,nz+1), source=0.0d0)
                                allocate(phisnap_root(nx+1,ny+1,nz+1), source=0.0d0)
                                allocate(wtordat(num_of_faces+1))
                                allocate(grid_point(grid_point_rows,3))
                                call CTCAR_regarea_real8(rtowdat, grid_point_rows+2, areaid(1))
                                call CTCAR_regarea_real8(wtordat, num_of_faces+1, areaid(2))
                                call CTCAR_regarea_int(grid_point, grid_point_rows*3, areaid(3))

                                dataint(1) = myid
                                dataint(2) = nstep

                                if (myid == 0) then
                                  print *, "grid_point_rows: ", grid_point_rows
                                  print *, "requester sendreq : ", myid
                                  print *, "nstep : ", nstep
                                  call CTCAR_sendreq(dataint, 2)    ! リクエスト送る
                                  print *, "requester sendreq done"
                                end if

                                call MPI_Barrier(MCW, mpierr)

                                call MPI_Bcast(grid_point, grid_point_rows*3, MPI_INTEGER4, 0, MCW, mpierr)

                                if (myid == 0) wtordat(1) = -1.0d0


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
!                            call oh1_families(famind,fammbr)
!            if(myid.eq.0) print*, "fam{ind,mbr}",famind,fammbr
                            call create_nborps(1)
            if(currmode.lt.0) then
                            call create_nborps(2)
                            call oh3_bcast_field &
                           &(eb(1,0,0,0,1),eb(1,0,0,0,2),FEB)
                            call oh3_bcast_field &
                           &(mp(1,0,0,0,1),mp(1,0,0,0,2),FEB)
!                            call oh3_bcast_field &
!                           &(mp(1,0,0,0,3),mp(1,0,0,0,4),FEB)
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
!                            call oh1_families(famind,fammbr)
!                            call chkres(1,0,1," chkresA1:")
!                            call chkres(2,0,1," chkresA2:")
              if(currmode.lt.0) then
                            rebalance = .true.
                            currmode = 1
              end if
            end do
            if(rebalance) then
                            call create_nborps(2)
                            call oh3_bcast_field &
                           &(eb(1,0,0,0,1),eb(1,0,0,0,2),FEB)
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
                                call esfld2(1)
              if(sdid(2).ge.0)  call esfld2(2)
!                              -------------
!--- check of parameters ---         *
!                              -------------
                                call chkprm
!                              -------------
!--- electrostatic potential ---     *
!                              -------------
!          if(nspec.ne.0)        call phispc(0)
!                              -------------
!                                    *
!---initial diagnostics---     -------------
                                call energy
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
                            mltime = mltime - MPI_Wtime()
!                                    *                    :
!                        **************************       :
                  MAINL: do istep=1,nstep
!                        **************************       :
!                                    *
!--- time ---                        *                    :
!                              +++++++++++++              :
                                 t = t+dt
                               itime = istep
      if(myid.eq.0) write(6,*) '**** step --------- ',itime
! check
                                ! ! 電位確認用
                                ! if (mod(istep,5000)==0) then
                                !   do x_tick = 0, nx
                                !     do y_tick = 0, ny
                                !       do z_tick = 0, nz
                                !         ! どの分割に含まれるかを計算
                                !         area_x = min(floor(real(x_tick)*nodes(1)/nx)+1, nodes(1))
                                !         area_y = min(floor(real(y_tick)*nodes(2)/ny)+1, nodes(2))
                                !         area_z = min(floor(real(z_tick)*nodes(3)/nz)+1, nodes(3))
                                !         ! 点がどのidのノードに含まれるかを計算
                                !         if (myid == nodes(1)*nodes(2)*(area_z-1) + nodes(1)*(area_y-1) + (area_x-1)) then
                                !           x_grid = x_tick - (nx/nodes(1))*(area_x-1)
                                !           y_grid = y_tick - (ny/nodes(2))*(area_y-1)
                                !           z_grid = z_tick - (nz/nodes(3))*(area_z-1)
                                !           phisnap_local(x_tick+1,y_tick+1,z_tick+1) = phi(1,x_grid,y_grid,z_grid,3)
                                !         end if
                                !       end do
                                !     end do
                                !   end do
                                !   call MPI_Reduce(phisnap_local, phisnap_root, (nx+1)*(ny+1)*(nz+1), MPI_REAL8, MPI_SUM, 0, MCW, mpierr)
                                !   if (myid == 0) then
                                !     write(filename, '(A,I5.5,A)') './phi_EMSES_init', istep, '.txt'
                                !     open (unit=24601,file=filename,action="write",status="replace")
                                !     ! write (out_unit,*) "istep", istep
                                !     do x_tick = 0, nx
                                !       do y_tick = 0, ny
                                !         do z_tick = 0, nz
                                !           write (24601,*) x_tick, y_tick, z_tick, phisnap_root(x_tick+1,y_tick+1,z_tick+1)
                                !         end do
                                !       end do
                                !     end do
                                !     close (24601)
                                !   end if
                                ! end if


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
!
!--- relocation of field pe and pb   *                    :
!                              -------------              :
                                call grdcrt(1)
!                              -------------              :
!---                                 *                    :
!                              -------------              :
                                call psuper3(1,2,0)
!                              -------------              :
!                                    *
!--- v(n+0.5) to v(n+0.5) ---        *                    :
!--- r(n) to r(n+1) ---              *                    :
!--- rho(n+1)_main --                *                    :
!                              -------------              :
                if(fajdg.eq.0) then
                                call psolve2(1,ustep)
!                                call psolve2(3,ustep)
                else
                                call psolve1(1,ustep,fajdg)
!                                call psolve1(3,ustep,fajdg)
                end if
!                              -------------              :
              if(sdid(2).ge.0) then
!--- relocation of field pe and pb   *                    :
!                              -------------              :
                                call grdcrt(2)
!                              -------------              :
!--- v(n+0.5) to v(n+0.5) ---        *                    :
!--- r(n) to r(n+1) ---              *                    :
!--- rho(n+1)_main ---               *                    :
!                              -------------              :
                if(fajdg.eq.0) then
                                call psolve2(2,ustep)
                else
                                call psolve1(2,ustep,fajdg)
                end if
!                              -------------              :
              end if
!--- j(n+0.5)_diag ---               *                    :
!      ======================  -------------              :
            if(fajdg.eq.1) then
              if(currmode.ne.0) call oh3_allreduce_field &
                               &(aj(1,0,0,0,1),aj(1,0,0,0,2),FAJ)
                                call oh3_exchange_borders &
                               &(aj(1,0,0,0,1),aj(1,0,0,0,2),CAJ,0)
                                call exchange_lcur_borders1(1)
                                call add_boundary_current1(1)
              if(currmode.ne.0) call oh3_allreduce_field &
                               &(ajdg(1,0,0,0,1),ajdg(1,0,0,0,2),FJD)
                                call oh3_exchange_borders &
                               &(ajdg(1,0,0,0,1),ajdg(1,0,0,0,2),CJD,0)
                                call exchange_lcur_borders2(1)
                                call add_boundary_current2(1)
            end if
!      ======================  -------------              :
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
!                            call oh1_families(famind,fammbr)
!                            call chkres(1,0,1," chkresA3:")
!                            call chkres(2,0,1," chkresA4:")
              if(currmode.lt.0) then
                            rebalance = .true.
                            currmode = 1
              end if
            end do
            if(rebalance) then
                            call create_nborps(2)
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

                            sey(1:3,3) = -nsecemit(1:3)*q(3)
                            call MPI_AllReduce(sey,seygl,12*ispec, &
                           &MPI_DOUBLE_PRECISION,MPI_SUM,MCW,mpierr)
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
!--- rho(n+0.5+m/2)_sub ---          *            :       :
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
                                rhobk(:,:,:,:,3) = rhobk(:,:,:,:,3) &
                               &                 + rhobk(:,:,:,:,1)
                                call add_background_charge(1)
!                              -------------
!--- charge density masking ---      *
!                              -------------
              if(imask(5).eq.1) call fsmask(1,4)
!                              -------------
!--- correction of e ---             *            v       :
!                              -------------      :       :
              if(npc.ge.1)      call surchg
                                estime = estime - MPI_Wtime()
                                call esfld1
                                estime = estime + MPI_Wtime()

                                ! ここで導体蓄積電荷を用いて計算
                                do grid_id = 1, grid_point_rows
                                  ! どの分割に含まれるかを計算
                                  area_x = min(floor(real(grid_point(grid_id,1))*nodes(1)/nx)+1, nodes(1))
                                  area_y = min(floor(real(grid_point(grid_id,2))*nodes(2)/ny)+1, nodes(2))
                                  area_z = min(floor(real(grid_point(grid_id,3))*nodes(3)/nz)+1, nodes(3))
                                  ! 点がどのidのノードに含まれるかを計算
                                  if (myid == nodes(1)*nodes(2)*(area_z-1) + nodes(1)*(area_y-1) + (area_x-1)) then
                                    x_grid = grid_point(grid_id,1) - (nx/nodes(1))*(area_x-1)
                                    y_grid = grid_point(grid_id,2) - (ny/nodes(2))*(area_y-1)
                                    z_grid = grid_point(grid_id,3) - (nz/nodes(3))*(area_z-1)
                                    phi_local(grid_id) = phi(1,x_grid,y_grid,z_grid,1)
                                  end if
                                end do
                                call MPI_Reduce(phi_local, phi_root, grid_point_rows, MPI_REAL8, MPI_SUM, 0, MCW, mpierr)

                                ! 通信
                                if (myid == 0) then
                                  ! print *, "req-step:", istep
                
                                  ! Req → Wrk 情報伝達
                                  ! R1  データ読込未完了@Wrkチェック
                                  do while (.true.)
                                    if (rtowdat(1) < 0.0d0) exit ! R1C rtowdat(1) 負がデータ読み込み完了
                                  end do
                
                                  ! R2  データ書込@Req実⾏
                                  rtowdat(2) = gcount(2)%chgacm(1,1)
                                  print *, "renq", renq
                                  print *, "remphi", renphi
                                  print *, "gcount(2)%chgacm(1,1)_esses", gcount(2)%chgacm(1,1)
                                  rtowdat(3:) = phi_root
                                  rtowdat(1) = istep
                
                                  ! Wrk → Req 情報伝達
                                  ! R3  データ書込@Wrk待ち
                                  do while (.true.)
                                    if (wtordat(1) >= 0.0d0) exit ! wtordat(1) 正がデータ書込完了
                                  end do
                                  omegas = wtordat(2:)
                                  ! R4  データ読込@Req
                                  ! if (istep == 4000) then
                                  !   open (unit=24601,file="./sigma_esses.txt",action="write",status="replace")
                                  !   do i = 1, num_of_faces
                                  !       write (24601,*), omegas(i)
                                  !   end do
                                  !   close (24601)
                                  ! end if
                                  
                                  ! R5  データ読込@Req完了フラグ
                                  wtordat(1) = -1.0d0
                
                                end if
                
                                call MPI_Bcast(omegas, num_of_faces, MPI_REAL8, 0, MCW, mpierr)
                
                                do i_bem = -1, nx/nodes(1)+1
                                  do j_bem = -1, ny/nodes(2)+1
                                    do k_bem = -1, nz/nodes(3)+1
                                      phi(1,i_bem,j_bem,k_bem,3) = 0.0d0
                                      do l_bem = 1, num_of_faces
                                        phi(1,i_bem,j_bem,k_bem,3) = phi(1,i_bem,j_bem,k_bem,3) + coef(i_bem,j_bem,k_bem,l_bem)*omegas(l_bem)
                                      end do
                                    end do
                                  end do
                                end do

                                ! 電位確認用
                                if (mod(istep,1000)==0) then
                                  do i_bem = 1,2
                                    do x_tick = 0, nx
                                      do y_tick = 0, ny
                                        do z_tick = 0, nz
                                          ! どの分割に含まれるかを計算
                                          area_x = min(floor(real(x_tick)*nodes(1)/nx)+1, nodes(1))
                                          area_y = min(floor(real(y_tick)*nodes(2)/ny)+1, nodes(2))
                                          area_z = min(floor(real(z_tick)*nodes(3)/nz)+1, nodes(3))
                                          ! 点がどのidのノードに含まれるかを計算
                                          if (myid == nodes(1)*nodes(2)*(area_z-1) + nodes(1)*(area_y-1) + (area_x-1)) then
                                            x_grid = x_tick - (nx/nodes(1))*(area_x-1)
                                            y_grid = y_tick - (ny/nodes(2))*(area_y-1)
                                            z_grid = z_tick - (nz/nodes(3))*(area_z-1)
                                            if (i_bem==1) then 
                                              phisnap_local(x_tick+1,y_tick+1,z_tick+1) = phi(1,x_grid,y_grid,z_grid,1)
                                            else
                                              phisnap_local(x_tick+1,y_tick+1,z_tick+1) = phi(1,x_grid,y_grid,z_grid,3) !+ phi(1,x_grid,y_grid,z_grid,1)
                                            end if
                                          end if
                                        end do
                                      end do
                                    end do
                                    call MPI_Reduce(phisnap_local, phisnap_root, (nx+1)*(ny+1)*(nz+1), MPI_REAL8, MPI_SUM, 0, MCW, mpierr)
                                    if (myid == 0) then
                                      if (i_bem==1) then
                                        write(filename, '(A,I6.6,A)') './phi_EMSES_esses', istep, '.txt'
                                        open (unit=24601,file=filename,action="write",status="replace")
                                          ! write (out_unit,*) "istep", istep
                                          do x_tick = 0, nx
                                            do y_tick = 0, ny
                                              do z_tick = 0, nz
                                                write (24601,*) x_tick, y_tick, z_tick, phisnap_root(x_tick+1,y_tick+1,z_tick+1)
                                              end do
                                            end do
                                          end do
                                        close (24601)
                                      else
                                        write(filename, '(A,I6.6,A)') './phi_BEM_esses', istep, '.txt'
                                        open (unit=24601,file=filename,action="write",status="replace")
                                          ! write (out_unit,*) "istep", istep
                                          do x_tick = 0, nx
                                            do y_tick = 0, ny
                                              do z_tick = 0, nz
                                                write (24601,*) x_tick, y_tick, z_tick, phisnap_root(x_tick+1,y_tick+1,z_tick+1)
                                              end do
                                            end do
                                          end do
                                        close (24601)
                                      end if
                                    end if
                                  end do
                                end if

                                phi(1,:,:,:,1) = phi(1,:,:,:,1) + phi(1,:,:,:,3)
                                
              if(currmode.ne.0) call oh3_bcast_field &
                               &(phi(1,0,0,0,1),phi(1,0,0,0,2),FPH)
                                call esfld2(1)
              if(sdid(2).ge.0)  call esfld2(2)
                                call poichk(1)
!                              -------------      :       :
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
!                                call oh3_exchange_borders &
!                               &(eb(1,0,0,0,1),eb(1,0,0,0,2),CEB,currmode)
!                                call boundary_emfld(1)
!              if(sdid(2).ge.0)  call boundary_emfld(2)
!                              -------------              :
!--- diagnostics ---                 *                    :
!                              -------------              :
                                call energy
!                              -------------              :
!                                    *                    :
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
                            mltime = mltime + MPI_Wtime()
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
                                call save_esdat
!                              -------------
        end if
!                                    *
!                              -------------
!!!!                                call chkcst
 999                              continue
!                       print*,"time main loop", emltime-smltime
                        call MPI_Reduce(mnmlocal,mnmtotal,1, &
                       &                MPI_INTEGER8,MPI_SUM,0,MCW,mpierr)
                        if (myid.eq.0) then
                          print*,"time main loop", (emltime-smltime)
                          print*,"time main loop2", mltime
                          print*,"time essolver", estime
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
                            end subroutine esses
