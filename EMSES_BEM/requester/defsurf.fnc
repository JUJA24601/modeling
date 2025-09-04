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
            tfrac = (pbuf(m)%z - lbhole)/pbuf(m)%vz
            xpast0 = pbuf(m)%x - pbuf(m)%vx*tfrac
            ypast0 = pbuf(m)%y - pbuf(m)%vy*tfrac
            xpast1 = xpast0 - xholec
            ypast1 = ypast0 - yholec
            if(xpast1*xpast1+ypast1*ypast1.le.rholesq) then
              pbuf(m)%x = xpast0
              pbuf(m)%y = ypast0
              pbuf(m)%z = lbhole
              if(npc.eq.0.or.npc.eq.1) then
                pbuf(m)%preside = INCB
              else if(npc.eq.2) then
                gcount(1)%influx(npc,is) = gcount(1)%influx(npc,is) + 1
                gcount(1)%infhist(pbuf(m)%pid,npc,is) = &
               &  gcount(1)%infhist(pbuf(m)%pid,npc,is) + 1
                pbuf(m)%preside = IPCB
              end if
!
              if(sec_emission.and.pbuf(m)%vz.lt.0.0d0) then
                keprt = pbuf(m)%vx*pbuf(m)%vx &
               &      + pbuf(m)%vy*pbuf(m)%vy &
               &      + pbuf(m)%vz*pbuf(m)%vz
                ikeprt = 1.0d0/keprt
!
                dotprod = -pbuf(m)%vz
                pcosth = max(abs(dotprod*sqrt(ikeprt)),0.5d0)
                ipcosth = 1.0d0/pcosth
!
                rkeprt = keprt*pemaxinvh(is)
                irkeprt= ikeprt*pemaxdble(is)
                irkeprt035 = irkeprt**0.35d0
                rkeprt135 = rkeprt/irkeprt035
                yield = dltamax1114*ipcosth*irkeprt035*(1.0d0-exp(-2.28d0*pcosth*rkeprt135))
                call RANU0(dranu,1,icon)
                nyield(1) = -round_prob(yield,dranu(1))
              end if
            else
              termA = xmove*xmove + ymove*ymove
              termB = xmove*xsepa + ymove*ysepa
              tfrac = (sqrt(termB*termB - &
             &         termA*(xsepa*xsepa + ysepa*ysepa - rholesq)) &
             &         - termB)/termA
              pbuf(m)%x = pbuf(m)%x + xmove*tfrac
              pbuf(m)%y = pbuf(m)%y + ymove*tfrac
              pbuf(m)%z = pbuf(m)%z + zmove*tfrac
              pbuf(m)%preside = INCB
!
              keprt_xy = pbuf(m)%vx*pbuf(m)%vx &
             &         + pbuf(m)%vy*pbuf(m)%vy
              if(sec_emission.and.keprt_xy.gt.0.0d0) then
                keprt = keprt_xy + pbuf(m)%vz*pbuf(m)%vz
                ikeprt = 1.0d0/keprt
!
                xpast1 = pbuf(m)%x - xholec
                ypast1 = pbuf(m)%y - yholec
                dotprod = xpast1*pbuf(m)%vx + ypast1*pbuf(m)%vy
                pcosth = max(abs(dotprod*sqrt(ikeprt)*rholeinv),0.5d0)
                ipcosth = 1.0d0/pcosth
!
                rkeprt = keprt*pemaxinvh(is)
                irkeprt= ikeprt*pemaxdble(is)
                irkeprt035 = irkeprt**0.35d0
                rkeprt135 = rkeprt/irkeprt035
                yield = dltamax1114*ipcosth*irkeprt035*(1.0d0-exp(-2.28d0*pcosth*rkeprt135))
                call RANU0(dranu,1,icon)
                nyield(2) = round_prob(yield,dranu(1))
              end if
            end if
          else if(pbuf(m)%z.le.ubhole.and. &
         &        xsepa*xsepa+ysepa*ysepa.gt.rholesq) then
            tfrac = (pbuf(m)%z - ubhole)/pbuf(m)%vz
            xpast0 = pbuf(m)%x - pbuf(m)%vx*tfrac
            ypast0 = pbuf(m)%y - pbuf(m)%vy*tfrac
            xpast1 = xpast0 - xholec
            ypast1 = ypast0 - yholec
            if(tfrac.gt.0.0d0.and.tfrac.le.gustep.and. &
           &   xpast1*xpast1+ypast1*ypast1.gt.rholesq) then
              pbuf(m)%x = xpast0
              pbuf(m)%y = ypast0
              pbuf(m)%z = ubhole
              if(npc.eq.0) then
                pbuf(m)%preside = INCB
              else if(npc.eq.1) then
                gcount(1)%influx(npc,is) = gcount(1)%influx(npc,is) + 1
                gcount(1)%infhist(pbuf(m)%pid,npc,is) = &
               &  gcount(1)%infhist(pbuf(m)%pid,npc,is) + 1
                pbuf(m)%preside = IPCB
              else if(npc.eq.2) then
                gcount(1)%influx(npc-1,is) = gcount(1)%influx(npc-1,is) + 1
                gcount(1)%infhist(pbuf(m)%pid,npc-1,is) = &
               &  gcount(1)%infhist(pbuf(m)%pid,npc-1,is) + 1
                pbuf(m)%preside = IPCB
              end if
!
              if(sec_emission.and.pbuf(m)%vz.lt.0.0d0) then
                keprt = pbuf(m)%vx*pbuf(m)%vx &
               &      + pbuf(m)%vy*pbuf(m)%vy &
               &      + pbuf(m)%vz*pbuf(m)%vz
                ikeprt = 1.0d0/keprt
!
                dotprod = -pbuf(m)%vz
                pcosth = max(abs(dotprod*sqrt(ikeprt)),0.5d0)
                ipcosth = 1.0d0/pcosth
!
                rkeprt = keprt*pemaxinvh(is)
                irkeprt= ikeprt*pemaxdble(is)
                irkeprt035 = irkeprt**0.35d0
                rkeprt135 = rkeprt/irkeprt035
                yield = dltamax1114*ipcosth*irkeprt035*(1.0d0-exp(-2.28d0*pcosth*rkeprt135))
                call RANU0(dranu,1,icon)
                nyield(1) = round_prob(yield,dranu(1))
              end if
            else
              termA = xmove*xmove + ymove*ymove
              termB = xmove*xsepa + ymove*ysepa
              tfrac = (sqrt(termB*termB - &
             &         termA*(xsepa*xsepa + ysepa*ysepa - rholesq)) &
             &         - termB)/termA
              pbuf(m)%x = pbuf(m)%x + xmove*tfrac
              pbuf(m)%y = pbuf(m)%y + ymove*tfrac
              pbuf(m)%z = pbuf(m)%z + zmove*tfrac
              pbuf(m)%preside = INCB
!
              keprt_xy = pbuf(m)%vx*pbuf(m)%vx &
             &         + pbuf(m)%vy*pbuf(m)%vy
              if(sec_emission.and.keprt_xy.gt.0.0d0) then
                keprt = keprt_xy + pbuf(m)%vz*pbuf(m)%vz
                ikeprt = 1.0d0/keprt
!
                xpast1 = pbuf(m)%x - xholec
                ypast1 = pbuf(m)%y - yholec
                dotprod = xpast1*pbuf(m)%vx + ypast1*pbuf(m)%vy
                pcosth = max(abs(dotprod*sqrt(ikeprt)*rholeinv),0.5d0)
                ipcosth = 1.0d0/pcosth
!
                rkeprt = keprt*pemaxinvh(is)
                irkeprt= ikeprt*pemaxdble(is)
                irkeprt035 = irkeprt**0.35d0
                rkeprt135 = rkeprt/irkeprt035
                yield = dltamax1114*ipcosth*irkeprt035*(1.0d0-exp(-2.28d0*pcosth*rkeprt135))
                call RANU0(dranu,1,icon)
                nyield(2) = round_prob(yield,dranu(1))
              end if
            end if
          end if
#elif ssurf==5
          if(pbuf(m)%z.le.lbhole) then
            tfrac = (pbuf(m)%z - lbhole)/pbuf(m)%vz
            xpast0 = pbuf(m)%x - pbuf(m)%vx*tfrac
            ypast0 = pbuf(m)%y - pbuf(m)%vy*tfrac
            if(xpast0.ge.xlwall.and.xpast0.le.xuwall) then
              pbuf(m)%x = xpast0
              pbuf(m)%y = ypast0
              pbuf(m)%z = lbhole
              if(npc.eq.0.or.npc.eq.1) then
                pbuf(m)%preside = INCB
              else if(npc.eq.2) then
                gcount(1)%influx(npc,is) = gcount(1)%influx(npc,is) + 1
                gcount(1)%infhist(pbuf(m)%pid,npc,is) = &
               &  gcount(1)%infhist(pbuf(m)%pid,npc,is) + 1
                pbuf(m)%preside = IPCB
              end if
!
              if(sec_emission.and.pbuf(m)%vz.lt.0.0d0) then
                keprt = pbuf(m)%vx*pbuf(m)%vx &
               &      + pbuf(m)%vy*pbuf(m)%vy &
               &      + pbuf(m)%vz*pbuf(m)%vz
                ikeprt = 1.0d0/keprt
!
                dotprod = -pbuf(m)%vz
                pcosth = max(abs(dotprod*sqrt(ikeprt)),0.5d0)
                ipcosth = 1.0d0/pcosth
!
                rkeprt = keprt*pemaxinvh(is)
                irkeprt= ikeprt*pemaxdble(is)
                irkeprt035 = irkeprt**0.35d0
                rkeprt135 = rkeprt/irkeprt035
                yield = dltamax1114*ipcosth*irkeprt035*(1.0d0-exp(-2.28d0*pcosth*rkeprt135))
                call RANU0(dranu,1,icon)
                nyield(1) = -round_prob(yield,dranu(1))
              end if
            else if(xpast0.lt.xlwall) then
              tfrac = (pbuf(m)%x - xlwall)/pbuf(m)%vx
              pbuf(m)%x = xlwall
              pbuf(m)%y = pbuf(m)%y - pbuf(m)%vy*tfrac
              pbuf(m)%z = pbuf(m)%z - pbuf(m)%vz*tfrac
              pbuf(m)%preside = INCB
!
              if(sec_emission.and.pbuf(m)%vx.lt.0.0d0) then
                keprt = pbuf(m)%vx*pbuf(m)%vx &
               &      + pbuf(m)%vy*pbuf(m)%vy &
               &      + pbuf(m)%vz*pbuf(m)%vz
                ikeprt = 1.0d0/keprt
!
                dotprod = -pbuf(m)%vx
                pcosth = max(abs(dotprod*sqrt(ikeprt)*rholeinv),0.5d0)
                ipcosth = 1.0d0/pcosth
!
                rkeprt = keprt*pemaxinvh(is)
                irkeprt= ikeprt*pemaxdble(is)
                irkeprt035 = irkeprt**0.35d0
                rkeprt135 = rkeprt/irkeprt035
                yield = dltamax1114*ipcosth*irkeprt035*(1.0d0-exp(-2.28d0*pcosth*rkeprt135))
                call RANU0(dranu,1,icon)
                nyield(2) = round_prob(yield,dranu(1))
              end if
            else if(xpast0.gt.xuwall) then
              tfrac = (pbuf(m)%x - xuwall)/pbuf(m)%vx
              pbuf(m)%x = xuwall
              pbuf(m)%y = pbuf(m)%y - pbuf(m)%vy*tfrac
              pbuf(m)%z = pbuf(m)%z - pbuf(m)%vz*tfrac
              pbuf(m)%preside = INCB
!
              if(sec_emission.and.pbuf(m)%vx.gt.0.0d0) then
                keprt = pbuf(m)%vx*pbuf(m)%vx &
               &      + pbuf(m)%vy*pbuf(m)%vy &
               &      + pbuf(m)%vz*pbuf(m)%vz
                ikeprt = 1.0d0/keprt
!
                dotprod = +pbuf(m)%vx
                pcosth = max(abs(dotprod*sqrt(ikeprt)),0.5d0)
                ipcosth = 1.0d0/pcosth
!
                rkeprt = keprt*pemaxinvh(is)
                irkeprt= ikeprt*pemaxdble(is)
                irkeprt035 = irkeprt**0.35d0
                rkeprt135 = rkeprt/irkeprt035
                yield = dltamax1114*ipcosth*irkeprt035*(1.0d0-exp(-2.28d0*pcosth*rkeprt135))
                call RANU0(dranu,1,icon)
                nyield(2) = round_prob(yield,dranu(1))
              end if
            end if
          else if(pbuf(m)%z.le.ubhole.and. &
         &        (pbuf(m)%x.le.xlwall.or.pbuf(m)%x.ge.xuwall)) then
            tfrac = (pbuf(m)%z - ubhole)/pbuf(m)%vz
            xpast0 = pbuf(m)%x - pbuf(m)%vx*tfrac
            ypast0 = pbuf(m)%y - pbuf(m)%vy*tfrac
            if((tfrac.gt.0.0d0.and.tfrac.le.gustep).and. &
           &   (xpast0.lt.xlwall.or.xpast0.gt.xuwall)) then
              pbuf(m)%x = xpast0
              pbuf(m)%y = ypast0
              pbuf(m)%z = ubhole
              if(npc.eq.0) then
                pbuf(m)%preside = INCB
              else if(npc.eq.1) then
                gcount(1)%influx(npc,is) = gcount(1)%influx(npc,is) + 1
                gcount(1)%infhist(pbuf(m)%pid,npc,is) = &
               &  gcount(1)%infhist(pbuf(m)%pid,npc,is) + 1
                pbuf(m)%preside = IPCB
              else if(npc.eq.2) then
                gcount(1)%influx(npc-1,is) = gcount(1)%influx(npc-1,is) + 1
                gcount(1)%infhist(pbuf(m)%pid,npc-1,is) = &
               &  gcount(1)%infhist(pbuf(m)%pid,npc-1,is) + 1
                pbuf(m)%preside = IPCB
              end if
!
              if(sec_emission.and.pbuf(m)%vz.lt.0.0d0) then
                keprt = pbuf(m)%vx*pbuf(m)%vx &
               &      + pbuf(m)%vy*pbuf(m)%vy &
               &      + pbuf(m)%vz*pbuf(m)%vz
                ikeprt = 1.0d0/keprt
!
                dotprod = -pbuf(m)%vz
                pcosth = max(abs(dotprod*sqrt(ikeprt)),0.5d0)
                ipcosth = 1.0d0/pcosth
!
                rkeprt = keprt*pemaxinvh(is)
                irkeprt= ikeprt*pemaxdble(is)
                irkeprt035 = irkeprt**0.35d0
                rkeprt135 = rkeprt/irkeprt035
                yield = dltamax1114*ipcosth*irkeprt035*(1.0d0-exp(-2.28d0*pcosth*rkeprt135))
                call RANU0(dranu,1,icon)
                nyield(1) = round_prob(yield,dranu(1))
              end if
            else if(pbuf(m)%x.le.xlwall) then
              tfrac = (pbuf(m)%x - xlwall)/pbuf(m)%vx
              pbuf(m)%x = xlwall
              pbuf(m)%y = pbuf(m)%y - pbuf(m)%vy*tfrac
              pbuf(m)%z = pbuf(m)%z - pbuf(m)%vz*tfrac
              pbuf(m)%preside = INCB
!
              if(sec_emission.and.pbuf(m)%vx.lt.0.0d0) then
                keprt = pbuf(m)%vx*pbuf(m)%vx &
               &      + pbuf(m)%vy*pbuf(m)%vy &
               &      + pbuf(m)%vz*pbuf(m)%vz
                ikeprt = 1.0d0/keprt
!
                dotprod = -pbuf(m)%vx
                pcosth = max(abs(dotprod*sqrt(ikeprt)),0.5d0)
                ipcosth = 1.0d0/pcosth
!
                rkeprt = keprt*pemaxinvh(is)
                irkeprt= ikeprt*pemaxdble(is)
                irkeprt035 = irkeprt**0.35d0
                rkeprt135 = rkeprt/irkeprt035
                yield = dltamax1114*ipcosth*irkeprt035*(1.0d0-exp(-2.28d0*pcosth*rkeprt135))
                call RANU0(dranu,1,icon)
                nyield(2) = round_prob(yield,dranu(1))
              end if
            else if(pbuf(m)%x.ge.xuwall) then
              tfrac = (pbuf(m)%x - xuwall)/pbuf(m)%vx
              pbuf(m)%x = xuwall
              pbuf(m)%y = pbuf(m)%y - pbuf(m)%vy*tfrac
              pbuf(m)%z = pbuf(m)%z - pbuf(m)%vz*tfrac
              pbuf(m)%preside = INCB
!
              if(sec_emission.and.pbuf(m)%vx.gt.0.0d0) then
                keprt = pbuf(m)%vx*pbuf(m)%vx &
               &      + pbuf(m)%vy*pbuf(m)%vy &
               &      + pbuf(m)%vz*pbuf(m)%vz
                ikeprt = 1.0d0/keprt
!
                dotprod = +pbuf(m)%vx
                pcosth = max(abs(dotprod*sqrt(ikeprt)*rholeinv),0.5d0)
                ipcosth = 1.0d0/pcosth
!
                rkeprt = keprt*pemaxinvh(is)
                irkeprt= ikeprt*pemaxdble(is)
                irkeprt035 = irkeprt**0.35d0
                rkeprt135 = rkeprt/irkeprt035
                yield = dltamax1114*ipcosth*irkeprt035*(1.0d0-exp(-2.28d0*pcosth*rkeprt135))
                call RANU0(dranu,1,icon)
                nyield(2) = round_prob(yield,dranu(1))
              end if
            end if
          end if
#endif
