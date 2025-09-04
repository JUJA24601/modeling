#ifndef ssurf
          if(pbuf(m)%z.le.zssurf) then
            tfrac = (pbuf(m)%z - zssurf)/zmove
            pbuf(m)%x = pbuf(m)%x - xmove*tfrac
            pbuf(m)%y = pbuf(m)%y - ymove*tfrac
            pbuf(m)%z = zssurf
            pbuf(m)%preside = INCB
          end if
#endif
#if ssurf==1
          if(pbuf(m)%z.le.zssurf) then
            tfrac = (pbuf(m)%z - zssurf)/zmove
            pbuf(m)%x = pbuf(m)%x - xmove*tfrac
            pbuf(m)%y = pbuf(m)%y - ymove*tfrac
            pbuf(m)%z = zssurf
            pbuf(m)%preside = INCB
          end if
#elif ssurf==2
          if(pbuf(m)%z.le.zssurf) then
            xsepa = pbuf(m)%x - xbowlc
            ysepa = pbuf(m)%y - ybowlc
            zsepa = pbuf(m)%z - zbowlc
            tfrac1 = xsepa*xsepa*rbwlsqi(1) &
           &       + ysepa*ysepa*rbwlsqi(2) &
           &       + zsepa*zsepa*rbwlsqi(3)
            if(tfrac1.ge.1.0d0) then
              xsepa = xsepa - xmove
              ysepa = ysepa - ymove
              zsepa = zsepa - zmove
              tfrac2 = xsepa*xsepa*rbwlsqi(1) &
             &       + ysepa*ysepa*rbwlsqi(2) &
             &       + zsepa*zsepa*rbwlsqi(3)
              tfrac1 = (tfrac2 - tfrac1*tfrac2)/(tfrac2 - tfrac1)
              tfrac2 = (pbuf(m)%z - zssurf)/zmove
              tfrac1 = abs(tfrac1 + floor(tfrac1))
              tfrac2 = abs(tfrac2 + floor(tfrac2))
              tfrac = min(tfrac1,tfrac2)
              pbuf(m)%x = pbuf(m)%x - xmove*tfrac
              pbuf(m)%y = pbuf(m)%y - ymove*tfrac
              pbuf(m)%z = pbuf(m)%z - zmove*tfrac
!              if(istep.gt.0.and.abs(tfrac).gt.1.0d0) then
!                print*, "myid,m,too large tfrac", myid, m, tfrac, &
!             &                 tfrac1,tfrac2, &
!             &                 sqrt(xsepa*xsepa*rbwlsqi(1) &
!             &                    + ysepa*ysepa*rbwlsqi(2) &
!             &                    + zsepa*zsepa*rbwlsqi(3))
!              end if
!              if(istep.gt.0.and.mod(m,100).eq.0) then
!             &  print*, "myid,m,rcapturedBA", myid, m, &
!             &  pbuf(m)%z - zssurf, &
!             &  sqrt(xsepa*xsepa*rbwlsqi(1) &
!             &     + ysepa*ysepa*rbwlsqi(2) &
!             &     + zsepa*zsepa*rbwlsqi(3)), &
!             &  sqrt((pbuf(m)%x - xbowlc)*(pbuf(m)%x - xbowlc)*rbwlsqi(1) &
!             &     + (pbuf(m)%y - ybowlc)*(pbuf(m)%y - ybowlc)*rbwlsqi(2) &
!             &     + (pbuf(m)%z - zbowlc)*(pbuf(m)%z - zbowlc)*rbwlsqi(3))
!              end if
              pbuf(m)%preside = INCB
            end if
          end if
!
          xsepa = pbuf(m)%x - xdomec
          ysepa = pbuf(m)%y - ydomec
          zsepa = pbuf(m)%z - zdomec
          if(xsepa*xsepa+ysepa*ysepa+zsepa*zsepa.le.rdomsq) then
            pbuf(m)%preside = INCB
          end if
#elif ssurf==3
          if(pbuf(m)%z.le.zssurf) then
            if(.not. &
           &   ((xlrechole(1).le.pbuf(m)%x.and.pbuf(m)%x.le.xurechole(1).and. &
           &     ylrechole(1).le.pbuf(m)%y.and.pbuf(m)%y.le.yurechole(1).and. &
           &     zlrechole(1).le.pbuf(m)%z.and.pbuf(m)%z.le.zurechole(1)).or. &
           &    (xlrechole(2).le.pbuf(m)%x.and.pbuf(m)%x.le.xurechole(2).and. &
           &     ylrechole(2).le.pbuf(m)%y.and.pbuf(m)%y.le.yurechole(2).and. &
           &     zlrechole(2).le.pbuf(m)%z.and.pbuf(m)%z.le.zurechole(2)))) then
              pbuf(m)%preside = INCB
            end if
          end if
#elif ssurf==4
          xsepa = pbuf(m)%x - xholec
          ysepa = pbuf(m)%y - yholec
          if(pbuf(m)%z.le.lbhole) then
            pbuf(m)%preside = INCB
          else if(pbuf(m)%z.le.ubhole.and. &
         &        xsepa*xsepa+ysepa*ysepa.gt.rholesq) then
            pbuf(m)%preside = INCB
          end if
#elif ssurf==5
          if(pbuf(m)%z.le.lbhole) then
            pbuf(m)%preside = INCB
          else if(pbuf(m)%z.le.ubhole.and. &
         &        (pbuf(m)%x.le.xlwall.or.pbuf(m)%x.ge.xuwall)) then
            pbuf(m)%preside = INCB
          end if
#endif
