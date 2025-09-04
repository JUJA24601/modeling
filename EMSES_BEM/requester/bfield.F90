#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine bfield(psd,round)
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   B F I E L D
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .   subroutine for update of b-field for half time step    .
!   ............................................................
!
!-------------------- parameter and common blocks
  use oh_type
  use paramt
  use allcom
  implicit none
!
  integer(kind=4) :: i,j,k, x,y,z
  integer(kind=4) :: xl,xu, yl,yu, zl,zu
  integer(kind=4) :: ps, psd
  integer(kind=4) :: round


!-------------------- test particle simulation
      if(juncan.ge.1000) return


!-------------------- 
      ps = modulo(psd-1,2) + 1
      xl = sdoms(1,1,sdid(ps)+1); xu = sdoms(2,1,sdid(ps)+1)
      yl = sdoms(1,2,sdid(ps)+1); yu = sdoms(2,2,sdid(ps)+1)
      zl = sdoms(1,3,sdid(ps)+1); zu = sdoms(2,3,sdid(ps)+1)


!-------------------- loop for update of the b-field
      if(round.eq.3) then
        do j=-1,yu-yl+1
        do i=-1,xu-xl+1
          if(eb(EX,i,j,-1,psd).ne.eb(EX,i,j,0,psd).or. &
         &   eb(EX,i,j,+1,psd).ne.eb(EX,i,j,0,psd).or. &
         &   eb(EX,i,j,+2,psd).ne.eb(EX,i,j,0,psd)) then
            print*, "illEX:", myid,psd,i,j,eb(EX,i,j,-1,psd)/rene, &
           &        eb(EX,i,j, 0,psd)/rene,eb(EX,i,j,+1,psd)/rene, &
           &        eb(EX,i,j,+2,psd)/rene
          end if
          if(eb(EY,i,j,-1,psd).ne.eb(EY,i,j,0,psd).or. &
         &   eb(EY,i,j,+1,psd).ne.eb(EY,i,j,0,psd).or. &
         &   eb(EY,i,j,+2,psd).ne.eb(EY,i,j,0,psd)) then
            print*, "illEY:", myid,psd,i,j,eb(EY,i,j,-1,psd)/rene, &
           &        eb(EY,i,j, 0,psd)/rene,eb(EY,i,j,+1,psd)/rene, &
           &        eb(EY,i,j,+2,psd)/rene
          end if
          if(eb(EZ,i,j,-1,psd).ne.eb(EZ,i,j,0,psd).or. &
         &   eb(EZ,i,j,+1,psd).ne.eb(EZ,i,j,0,psd).or. &
         &   eb(EZ,i,j,+2,psd).ne.eb(EZ,i,j,0,psd)) then
            print*, "illEZ:", myid,psd,i,j,eb(EZ,i,j,-1,psd)/rene, &
           &        eb(EZ,i,j, 0,psd)/rene,eb(EZ,i,j,+1,psd)/rene, &
           &        eb(EZ,i,j,+2,psd)/rene
          end if
          if(eb(BX,i,j,-1,psd).ne.eb(BX,i,j,0,psd).or. &
         &   eb(BX,i,j,+1,psd).ne.eb(BX,i,j,0,psd).or. &
         &   eb(BX,i,j,+2,psd).ne.eb(BX,i,j,0,psd)) then
            print*, "illBX:", myid,psd,i,j,eb(BX,i,j,-1,psd)/rene, &
           &        eb(BX,i,j, 0,psd)/rene,eb(BX,i,j,+1,psd)/rene, &
           &        eb(BX,i,j,+2,psd)/rene
          end if
          if(eb(BY,i,j,-1,psd).ne.eb(BY,i,j,0,psd).or. &
         &   eb(BY,i,j,+1,psd).ne.eb(BY,i,j,0,psd).or. &
         &   eb(BY,i,j,+2,psd).ne.eb(BY,i,j,0,psd)) then
            print*, "illBY:", myid,psd,i,j,eb(BY,i,j,-1,psd)/rene, &
           &        eb(BY,i,j, 0,psd)/rene,eb(BY,i,j,+1,psd)/rene, &
           &        eb(BY,i,j,+2,psd)/rene
          end if
          if(eb(BZ,i,j,-1,psd).ne.eb(BZ,i,j,0,psd).or. &
         &   eb(BZ,i,j,+1,psd).ne.eb(BZ,i,j,0,psd).or. &
         &   eb(BZ,i,j,+2,psd).ne.eb(BZ,i,j,0,psd)) then
            print*, "illBZ:", myid,psd,i,j,eb(BZ,i,j,-1,psd)/rene, &
           &        eb(BZ,i,j, 0,psd)/rene,eb(BZ,i,j,+1,psd)/rene, &
           &        eb(BZ,i,j,+2,psd)/rene
          end if
        end do
        end do
      end if


!-------------------- loop for update of the b-field
      if(pmluse) then
        do k=-1,zu-zl
        do j=-1,yu-yl
        do i=-1,xu-xl
          x = i + xl; y = j + yl; z = k + zl
          if(x.ge.pmlxl.and.x.lt.pmlxu.and. &
         &   y.ge.pmlyl.and.y.lt.pmlyu.and. &
         &   z.ge.pmlzl.and.z.lt.pmlzu) then
            eb(BX,i,j,k,psd) = eb(BX,i  ,j  ,k  ,psd) &
           &                 + eb(EY,i  ,j  ,k+1,psd) - eb(EY,i,j,k,psd) &
           &                 - eb(EZ,i  ,j+1,k  ,psd) + eb(EZ,i,j,k,psd)
            eb(BY,i,j,k,psd) = eb(BY,i  ,j  ,k  ,psd) &
           &                 + eb(EZ,i+1,j  ,k  ,psd) - eb(EZ,i,j,k,psd) &
           &                 - eb(EX,i  ,j  ,k+1,psd) + eb(EX,i,j,k,psd)
            eb(BZ,i,j,k,psd) = eb(BZ,i  ,j  ,k  ,psd) &
           &                 + eb(EX,i  ,j+1,k  ,psd) - eb(EX,i,j,k,psd) &
           &                 - eb(EY,i+1,j  ,k  ,psd) + eb(EY,i,j,k,psd)
          else
            ebsc(BXY,i,j,k,psd) = pmlc(3,y,2)*ebsc(BXY,i,j,k,psd) &
           &                    - pmlc(4,y,2) &
           &                    *( eb(EZ,i  ,j+1,k  ,psd) - eb(EZ,i,j,k,psd) )
            ebsc(BXZ,i,j,k,psd) = pmlc(5,z,2)*ebsc(BXZ,i,j,k,psd) &
           &                    + pmlc(6,z,2) &
           &                    *( eb(EY,i  ,j  ,k+1,psd) - eb(EY,i,j,k,psd) )
            eb(BX,i,j,k,psd) = ebsc(BXY,i,j,k,psd) + ebsc(BXZ,i,j,k,psd)
!
            ebsc(BYZ,i,j,k,psd) = pmlc(5,z,2)*ebsc(BYZ,i,j,k,psd) &
           &                    - pmlc(6,z,2) &
           &                    *( eb(EX,i  ,j  ,k+1,psd) - eb(EX,i,j,k,psd) )
            ebsc(BYX,i,j,k,psd) = pmlc(1,x,2)*ebsc(BYX,i,j,k,psd) &
           &                    + pmlc(2,x,2) &
           &                    *( eb(EZ,i+1,j  ,k  ,psd) - eb(EZ,i,j,k,psd) )
            eb(BY,i,j,k,psd) = ebsc(BYZ,i,j,k,psd) + ebsc(BYX,i,j,k,psd)
!
            ebsc(BZX,i,j,k,psd) = pmlc(1,x,2)*ebsc(BZX,i,j,k,psd) &
           &                    - pmlc(2,x,2) &
           &                    *( eb(EY,i+1,j  ,k  ,psd) - eb(EY,i,j,k,psd) )
            ebsc(BZY,i,j,k,psd) = pmlc(3,y,2)*ebsc(BZY,i,j,k,psd) &
           &                    + pmlc(4,y,2) &
           &                    *( eb(EX,i  ,j+1,k  ,psd) - eb(EX,i,j,k,psd) )
            eb(BZ,i,j,k,psd) = ebsc(BZX,i,j,k,psd) + ebsc(BZY,i,j,k,psd)
          end if
        end do
        end do
        end do
      else
        do k=-1,zu-zl
        do j=-1,yu-yl
        do i=-1,xu-xl
          eb(BX,i,j,k,psd) = eb(BX,i  ,j  ,k  ,psd) &
         &                 + eb(EY,i  ,j  ,k+1,psd) - eb(EY,i,j,k,psd) &
         &                 - eb(EZ,i  ,j+1,k  ,psd) + eb(EZ,i,j,k,psd)
          eb(BY,i,j,k,psd) = eb(BY,i  ,j  ,k  ,psd) &
         &                 + eb(EZ,i+1,j  ,k  ,psd) - eb(EZ,i,j,k,psd) &
         &                 - eb(EX,i  ,j  ,k+1,psd) + eb(EX,i,j,k,psd)
          eb(BZ,i,j,k,psd) = eb(BZ,i  ,j  ,k  ,psd) &
         &                 + eb(EX,i  ,j+1,k  ,psd) - eb(EX,i,j,k,psd) &
         &                 - eb(EY,i+1,j  ,k  ,psd) + eb(EY,i,j,k,psd)
        end do
        end do
        end do
      end if


!-------------------- boundary treatment
!      call fsmask(1)
!      call fbound(1)


  return
  end subroutine bfield
