#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine grdcrt(ps)
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   G R D C R T
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .  this subroutine relocates electric and magnetic fields  .
!   .  defined at grids for interpolation to particle          .
!   .  positions.  (self-force of particle is eliminated by    .
!   .  the relocation)                                         .
!   ............................................................
!
!-------------------- parameter and common blocks
  use oh_type
  use paramt
  use allcom
  implicit none
!
  integer(kind=4) :: i,j,k, ieb, ipc, icap
  integer(kind=4) :: xl,yl,zl, xu,yu,zu
  integer(kind=4) :: ps,psd
  real(kind=8) :: eserg


!-------------------- 
      xl = sdoms(1,1,sdid(ps)+1); xu = sdoms(2,1,sdid(ps)+1)
      yl = sdoms(1,2,sdid(ps)+1); yu = sdoms(2,2,sdid(ps)+1)
      zl = sdoms(1,3,sdid(ps)+1); zu = sdoms(2,3,sdid(ps)+1)


!-------------------- set e0, b0
      do k=0,zu-zl
      do j=0,yu-yl
      do i=0,xu-xl
        wrk(EX:BZ,i,j,k) = (/e0x,e0y,e0z,b0x,b0y,b0z/)
      end do
      end do
      end do


!-------------------- particle pushing fields & relocation
      do ieb=0,nebfld
        psd = ieb*2 + ps
!       ------------- pe-fields
        do k=0,zu-zl
        do j=0,yu-yl
        do i=0,xu-xl
          wrk(EX,i,j,k) = wrk(EX,i,j,k) &
         &              + 0.5d0*(eb(EX,i,j,k,psd) + eb(EX,i-1,j,k,psd))
          wrk(EY,i,j,k) = wrk(EY,i,j,k) &
         &              + 0.5d0*(eb(EY,i,j,k,psd) + eb(EY,i,j-1,k,psd))
          wrk(EZ,i,j,k) = wrk(EZ,i,j,k) &
         &              + 0.5d0*(eb(EZ,i,j,k,psd) + eb(EZ,i,j,k-1,psd))
        end do
        end do
        end do
!
!       ------------- pb-fields
        do k=0,zu-zl
        do j=0,yu-yl
        do i=0,xu-xl
          wrk(BX,i,j,k) = wrk(BX,i,j,k) &
         &              + 0.25d0*(eb(BX,i,j,k  ,psd) + eb(BX,i,j-1,k  ,psd) &
         &                      + eb(BX,i,j,k-1,psd) + eb(BX,i,j-1,k-1,psd))
          wrk(BY,i,j,k) = wrk(BY,i,j,k) &
         &              + 0.25d0*(eb(BY,i,j,k  ,psd) + eb(BY,i-1,j,k  ,psd) &
         &                      + eb(BY,i,j,k-1,psd) + eb(BY,i-1,j,k-1,psd))
          wrk(BZ,i,j,k) = wrk(BZ,i,j,k) &
         &              + 0.25d0*(eb(BZ,i,j  ,k,psd) + eb(BZ,i-1,j  ,k,psd) &
         &                      + eb(BZ,i,j-1,k,psd) + eb(BZ,i-1,j-1,k,psd))
        end do
        end do
        end do
      end do


!-------------------- body surface electric field
      if(sfecrrct.ge.1) then
        eserg = path(1)*path(1)
        do ipc=1,npc
          do icap=nscpmx(ipc),nmxcpmx(ipc)-1
            i = bdygrid(1,icap) - xl
            j = bdygrid(2,icap) - yl
            k = bdygrid(3,icap) - zl
            if(i.ge.0.and.i.le.xu-xl.and. &
           &   j.ge.0.and.j.le.yu-yl.and. &
           &   k.ge.0.and.k.le.zu-zl) then
              if(abs(eb(EX,i,j,k,ps)).ge.1.0d-8*eserg.and. &
             &   abs(eb(EX,i-1,j,k,ps)).lt.1.0d-8*eserg) then
                wrk(EX,i,j,k) = eb(EX,i,j,k,ps) + e0x
              else if(abs(eb(EX,i-1,j,k,ps)).ge.1.0d-8*eserg.and. &
             &        abs(eb(EX,i,j,k,ps)).lt.1.0d-8*eserg) then
                wrk(EX,i,j,k) = eb(EX,i-1,j,k,ps) + e0x
              end if
              if(abs(eb(EY,i,j,k,ps)).ge.1.0d-8*eserg.and. & 
             &   abs(eb(EY,i,j-1,k,ps)).lt.1.0d-8*eserg) then
                wrk(EY,i,j,k) = eb(EY,i,j,k,ps) + e0y
              else if(abs(eb(EY,i,j-1,k,ps)).ge.1.0d-8*eserg.and. &
             &        abs(eb(EY,i,j,k,ps)).lt.1.0d-8*eserg) then
                wrk(EY,i,j,k) = eb(EY,i,j-1,k,ps) + e0y
              end if
              if(abs(eb(EZ,i,j,k,ps)).ge.1.0d-8*eserg.and. &
             &   abs(eb(EZ,i,j,k-1,ps)).lt.1.0d-8*eserg) then
                wrk(EZ,i,j,k) = eb(EZ,i,j,k,ps) + e0z
              else if(abs(eb(EZ,i,j,k-1,ps)).ge.1.0d-8*eserg.and. &
             &        abs(eb(EZ,i,j,k,ps)).lt.1.0d-8*eserg) then
                wrk(EZ,i,j,k) = eb(EZ,i,j,k-1,ps) + e0z
              end if
            end if
          end do
        end do
      end if


!--------------------- masking of longitudinal e-field
!      call fsmask(6)
!-------------------- boundary treatment
!      call fbound(4)
!      call fbound(3)


!-------------------- filtering pex, pey and pez in y and x respectively
!      call filtr3(pex,ix,iy,iz,nxm,nym,nzm,ipexfl)
!      call filtr3(pey,ix,iy,iz,nxm,nym,nzm,ipeyfl)
!      call filtr3(pez,ix,iy,iz,nxm,nym,nzm,ipezfl)
!      call fbound(3)


  return
  end subroutine
