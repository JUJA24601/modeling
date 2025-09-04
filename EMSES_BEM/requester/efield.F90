!
#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
!
  subroutine efield(psd)
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   E F I E L D
!   ____________________________________________________________
!
!   ............................................................
!   .                                                          .
!   .   subroutine to update e-field by current contribution   .
!   ............................................................
!
!-------------------- parameter and common blocks
  use oh_type
  use paramt
  use allcom
  implicit none
!
  integer(kind=4) :: i,j,k, x,y,z, ieb, ig
  integer(kind=4) :: xl,xu, yl,yu, zl,zu
  integer(kind=4) :: ps, psd
  real(kind=8) :: dlt
  real(kind=8) :: capc1


!-------------------- test particle simulation
      if(nflag_testp.ne.1) then
        dlt = 2.0d0
      else
        dlt = 0.0d0
      end if


!-------------------- 
      ps = modulo(psd-1,2) + 1
      xl = sdoms(1,1,sdid(ps)+1); xu = sdoms(2,1,sdid(ps)+1)
      yl = sdoms(1,2,sdid(ps)+1); yu = sdoms(2,2,sdid(ps)+1)
      zl = sdoms(1,3,sdid(ps)+1); zu = sdoms(2,3,sdid(ps)+1)


!-------------------- 
      if(psd.eq.nebfld*2-1.or.psd.eq.nebfld*2) then
        do ig=1,ngap
          if(i_gap(ig).ge.xl.and.i_gap(ig).le.xu.and. &
         &   j_gap(ig).ge.yl.and.j_gap(ig).le.yu.and. &
         &   k_gap(ig).ge.zl.and.k_gap(ig).le.zu) then
            ezgap(ig) = 0
            do ieb=0,(nebfld-1)*2,2
              ezgap(ig) = ezgap(ig) &
             &          + eb(EZ,i_gap(ig)-xl,j_gap(ig)-yl,k_gap(ig)-zl,ieb+ps)
            end do
          end if
        end do
      end if


!-------------------- loop for update of the e-field
      if(pmluse) then
        do k=0,zu-zl
        do j=0,yu-yl
        do i=0,xu-xl
          x = i + xl; y = j + yl; z = k + zl
          if(x.ge.pmlxl.and.x.lt.pmlxu.and. &
         &   y.ge.pmlyl.and.y.lt.pmlyu.and. &
         &   z.ge.pmlzl.and.z.lt.pmlzu) then
            eb(EX,i,j,k,psd) = eb(EX,i,j,k,psd) &
           &                 + tcs*( eb(BZ,i,j,k,psd) - eb(BZ,i  ,j-1,k  ,psd) &
           &                       - eb(BY,i,j,k,psd) + eb(BY,i  ,j  ,k-1,psd) ) &
           &                 - dlt*aj(JX,i,j,k,psd)
            eb(EY,i,j,k,psd) = eb(EY,i,j,k,psd) &
           &                 + tcs*( eb(BX,i,j,k,psd) - eb(BX,i  ,j  ,k-1,psd) &
           &                       - eb(BZ,i,j,k,psd) + eb(BZ,i-1,j  ,k  ,psd) ) &
           &                 - dlt*aj(JY,i,j,k,psd)
            eb(EZ,i,j,k,psd) = eb(EZ,i,j,k,psd) &
           &                 + tcs*( eb(BY,i,j,k,psd) - eb(BY,i-1,j  ,k  ,psd) &
           &                       - eb(BX,i,j,k,psd) + eb(BX,i  ,j-1,k  ,psd) ) &
           &                 - dlt*aj(JZ,i,j,k,psd)
          else
            ebsc(EXY,i,j,k,psd) = pmlc(3,y,1)*ebsc(EXY,i,j,k,psd) &
           &                    + pmlc(4,y,1) &
           &                    *( eb(BZ,i,j,k,psd) - eb(BZ,i  ,j-1,k  ,psd) )
            ebsc(EXZ,i,j,k,psd) = pmlc(5,z,1)*ebsc(EXZ,i,j,k,psd) &
           &                    - pmlc(6,z,1) &
           &                    *( eb(BY,i,j,k,psd) - eb(BY,i  ,j  ,k-1,psd) )
            eb(EX,i,j,k,psd) = ebsc(EXY,i,j,k,psd) + ebsc(EXZ,i,j,k,psd) &
           &                 - dlt*aj(JX,i,j,k,psd)
!
            ebsc(EYZ,i,j,k,psd) = pmlc(5,z,1)*ebsc(EYZ,i,j,k,psd) &
           &                    + pmlc(6,z,1) &
           &                    *( eb(BX,i,j,k,psd) - eb(BX,i  ,j  ,k-1,psd) )
            ebsc(EYX,i,j,k,psd) = pmlc(1,x,1)*ebsc(EYX,i,j,k,psd) &
           &                    - pmlc(2,x,1) &
           &                    *( eb(BZ,i,j,k,psd) - eb(BZ,i-1,j  ,k  ,psd) )
            eb(EY,i,j,k,psd) = ebsc(EYZ,i,j,k,psd) + ebsc(EYX,i,j,k,psd) &
           &                 - dlt*aj(JY,i,j,k,psd)
!
            ebsc(EZX,i,j,k,psd) = pmlc(1,x,1)*ebsc(EZX,i,j,k,psd) &
           &                    + pmlc(2,x,1) &
           &                    *( eb(BY,i,j,k,psd) - eb(BY,i-1,j  ,k  ,psd) )
            ebsc(EZY,i,j,k,psd) = pmlc(3,y,1)*ebsc(EZY,i,j,k,psd) &
           &                    - pmlc(4,y,1) &
           &                    *( eb(BX,i,j,k,psd) - eb(BX,i  ,j-1,k  ,psd) )
            eb(EZ,i,j,k,psd) = ebsc(EZX,i,j,k,psd) + ebsc(EZY,i,j,k,psd) &
           &                 - dlt*aj(JZ,i,j,k,psd)
          end if
        end do
        end do
        end do
      else
        do k=0,zu-zl
        do j=0,yu-yl
        do i=0,xu-xl
          eb(EX,i,j,k,psd) = eb(EX,i,j,k,psd) &
         &                 + tcs*( eb(BZ,i,j,k,psd) - eb(BZ,i  ,j-1,k  ,psd) &
         &                       - eb(BY,i,j,k,psd) + eb(BY,i  ,j  ,k-1,psd) ) &
         &                 - dlt*aj(JX,i,j,k,psd)
          eb(EY,i,j,k,psd) = eb(EY,i,j,k,psd) &
         &                 + tcs*( eb(BX,i,j,k,psd) - eb(BX,i  ,j  ,k-1,psd) &
         &                       - eb(BZ,i,j,k,psd) + eb(BZ,i-1,j  ,k  ,psd) ) &
         &                 - dlt*aj(JY,i,j,k,psd)
          eb(EZ,i,j,k,psd) = eb(EZ,i,j,k,psd) &
         &                 + tcs*( eb(BY,i,j,k,psd) - eb(BY,i-1,j  ,k  ,psd) &
         &                       - eb(BX,i,j,k,psd) + eb(BX,i  ,j-1,k  ,psd) ) &
         &                 - dlt*aj(JZ,i,j,k,psd)
        end do
        end do
        end do
      end if


!-------------------- modeling of lumped elements
      if(mode_dipole.eq.3.and.(psd.eq.nebfld*2-1.or.psd.eq.nebfld*2)) then
        do ig=1,ngap
          i = i_gap(ig) - xl
          j = j_gap(ig) - yl
          k = k_gap(ig) - zl
          if(i.ge.0.and.i.le.xu-xl.and. &
         &   j.ge.0.and.j.le.yu-yl.and. &
         &   k.ge.0.and.k.le.zu-zl) then
            capc1 = capc(ig) + 1.0d0
            if(resc(ig).lt.0) then
              eb(EZ,i,j,k,psd) = ezgap(ig) &
             &                + tcs*( eb(BY,i,j,k,psd) - eb(BY,i-1,j,k,psd) &
             &                - eb(BX,i,j,k,psd) + eb(BX,i,j-1,k,psd) )/capc1 &
             &                - dlt*aj(JZ,i,j,k,psd)/capc1
            else
              eb(EZ,i,j,k,psd) = ezgap(ig)*(resc(ig)*capc1 - 1.0d0) &
             &                 /(resc(ig)*capc1 + 1.0d0) &
             &                + tcs*resc(ig) &
             &                 *( eb(BY,i,j,k,psd) - eb(BY,i-1,j,k,psd) &
             &                  - eb(BX,i,j,k,psd) + eb(BX,i,j-1,k,psd) ) &
             &                 /(resc(ig)*capc1 + 1.0d0) &
             &                -dlt*resc(ig)*aj(JZ,i,j,k,psd) &
             &                 /(resc(ig)*capc1 + 1.0d0)
            end if
          end if
        end do
      end if


  return
  end subroutine efield
