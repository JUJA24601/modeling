!          if(line_mode.ne.1.and.pbuf(m)%nid.ne.-1) then
          if(line_mode.ne.1) then
            jpc = 0
            tfrac = 0.0d0
            do ipc=1,npc
              if(pcond(ipc)) then
!                if(pbuf(m)%preside.lt.0) cycle IPCL1
                if(geotype(ipc).eq.0.or.geotype(ipc).eq.1) then
                  if(pbuf(m)%x.ge.xlpc(ipc).and.pbuf(m)%x.le.xupc(ipc).and. &
                 &   pbuf(m)%y.ge.ylpc(ipc).and.pbuf(m)%y.le.yupc(ipc).and. &
                 &   pbuf(m)%z.ge.zlpc(ipc).and.pbuf(m)%z.le.zupc(ipc)) then
                    xmovei = 1.0d0/xmove
                    ymovei = 1.0d0/ymove
                    zmovei = 1.0d0/zmove
                    tfrac1 = (pbuf(m)%x - xlpc(ipc))*xmovei
                    tfrac2 = (pbuf(m)%x - xupc(ipc))*xmovei
                    tfrac3 = (pbuf(m)%y - ylpc(ipc))*ymovei
                    tfrac4 = (pbuf(m)%y - yupc(ipc))*ymovei
                    tfrac5 = (pbuf(m)%z - zlpc(ipc))*zmovei
                    tfrac6 = (pbuf(m)%z - zupc(ipc))*zmovei
                    tfrac1 = abs(tfrac1 + floor(tfrac1))
                    tfrac2 = abs(tfrac2 + floor(tfrac2))
                    tfrac3 = abs(tfrac3 + floor(tfrac3))
                    tfrac4 = abs(tfrac4 + floor(tfrac4))
                    tfrac5 = abs(tfrac5 + floor(tfrac5))
                    tfrac6 = abs(tfrac6 + floor(tfrac6))
                    tfrac = min(tfrac1,tfrac2,tfrac3,tfrac4,tfrac5,tfrac6)
                    pbuf(m)%x = pbuf(m)%x - xmove*tfrac
                    pbuf(m)%y = pbuf(m)%y - ymove*tfrac
                    pbuf(m)%z = pbuf(m)%z - zmove*tfrac
                    gcount(1)%influx(ipc,is) = gcount(1)%influx(ipc,is) + 1
                    gcount(1)%infhist(pbuf(m)%pid,ipc,is) = &
                   &  gcount(1)%infhist(pbuf(m)%pid,ipc,is) + 1
                    pbuf(m)%preside = IPCB
                    exit
                  end if
!
                else if(geotype(ipc).eq.2) then
                  if(cylinder(ipc)%align.eq.1) then
                    disp1 = pbuf(m)%y - cylinder(ipc)%axis(1)
                    disp2 = pbuf(m)%z - cylinder(ipc)%axis(2)
                    tfrac1 = (disp1*disp1 + disp2*disp2)*radsqi(ipc)
                    if(pbuf(m)%x.ge.cylinder(ipc)%edge(1).and. &
                   &   pbuf(m)%x.le.cylinder(ipc)%edge(2).and. &
                   &   tfrac1.lt.1.0d0) then
                      disp1 = disp1 - ymove
                      disp2 = disp2 - zmove
                      tfrac2 = (disp1*disp1 + disp2*disp2)*radsqi(ipc)
                      tfrac1 = (tfrac2 - tfrac1*tfrac2)/(tfrac2 - tfrac1)
                      xmovei = 1.0d0/xmove
                      tfrac2 = (pbuf(m)%x - cylinder(ipc)%edge(1))*xmovei
                      tfrac3 = (pbuf(m)%x - cylinder(ipc)%edge(2))*xmovei
                      tfrac1 = abs(tfrac1 + floor(tfrac1))
                      tfrac2 = abs(tfrac2 + floor(tfrac2))
                      tfrac3 = abs(tfrac3 + floor(tfrac3))
                      tfrac = min(tfrac1,tfrac2,tfrac3)
                      pbuf(m)%x = pbuf(m)%x - xmove*tfrac
                      pbuf(m)%y = pbuf(m)%y - ymove*tfrac
                      pbuf(m)%z = pbuf(m)%z - zmove*tfrac
                      gcount(1)%influx(ipc,is) = gcount(1)%influx(ipc,is) + 1
                      gcount(1)%infhist(pbuf(m)%pid,ipc,is) = &
                     &  gcount(1)%infhist(pbuf(m)%pid,ipc,is) + 1
                      pbuf(m)%preside = IPCB
                      exit
                    end if
!
                  else if(cylinder(ipc)%align.eq.2) then
                    disp1 = pbuf(m)%z - cylinder(ipc)%axis(1)
                    disp2 = pbuf(m)%x - cylinder(ipc)%axis(2)
                    tfrac1 = (disp1*disp1 + disp2*disp2)*radsqi(ipc)
                    if(pbuf(m)%y.ge.cylinder(ipc)%edge(1).and. &
                   &   pbuf(m)%y.le.cylinder(ipc)%edge(2).and. &
                   &   tfrac1.lt.1.0d0) then
                      disp1 = disp1 - zmove
                      disp2 = disp2 - xmove
                      tfrac2 = (disp1*disp1 + disp2*disp2)*radsqi(ipc)
                      tfrac1 = (tfrac2 - tfrac1*tfrac2)/(tfrac2 - tfrac1)
                      ymovei = 1.0d0/ymove
                      tfrac2 = (pbuf(m)%y - cylinder(ipc)%edge(1))*ymovei
                      tfrac3 = (pbuf(m)%y - cylinder(ipc)%edge(2))*ymovei
                      tfrac1 = abs(tfrac1 + floor(tfrac1))
                      tfrac2 = abs(tfrac2 + floor(tfrac2))
                      tfrac3 = abs(tfrac3 + floor(tfrac3))
                      tfrac = min(tfrac1,tfrac2,tfrac3)
                      pbuf(m)%x = pbuf(m)%x - xmove*tfrac
                      pbuf(m)%y = pbuf(m)%y - ymove*tfrac
                      pbuf(m)%z = pbuf(m)%z - zmove*tfrac
                      gcount(1)%influx(ipc,is) = gcount(1)%influx(ipc,is) + 1
                      gcount(1)%infhist(pbuf(m)%pid,ipc,is) = &
                     &  gcount(1)%infhist(pbuf(m)%pid,ipc,is) + 1
                      pbuf(m)%preside = IPCB
                      exit
                    end if
!
                  else if(cylinder(ipc)%align.eq.3) then
                    disp1 = pbuf(m)%x - cylinder(ipc)%axis(1)
                    disp2 = pbuf(m)%y - cylinder(ipc)%axis(2)
                    tfrac1 = (disp1*disp1 + disp2*disp2)*radsqi(ipc)
                    if(pbuf(m)%z.ge.cylinder(ipc)%edge(1).and. &
                   &   pbuf(m)%z.le.cylinder(ipc)%edge(2).and. &
                   &   tfrac1.lt.1.0d0) then
                      disp1 = disp1 - xmove
                      disp2 = disp2 - ymove
                      tfrac2 = (disp1*disp1 + disp2*disp2)*radsqi(ipc)
                      tfrac1 = (tfrac2 - tfrac1*tfrac2)/(tfrac2 - tfrac1)
                      zmovei = 1.0d0/zmove
                      tfrac2 = (pbuf(m)%z - cylinder(ipc)%edge(1))*zmovei
                      tfrac3 = (pbuf(m)%z - cylinder(ipc)%edge(2))*zmovei
                      tfrac1 = abs(tfrac1 + floor(tfrac1))
                      tfrac2 = abs(tfrac2 + floor(tfrac2))
                      tfrac3 = abs(tfrac3 + floor(tfrac3))
                      tfrac = min(tfrac1,tfrac2,tfrac3)
                      pbuf(m)%x = pbuf(m)%x - xmove*tfrac
                      pbuf(m)%y = pbuf(m)%y - ymove*tfrac
                      pbuf(m)%z = pbuf(m)%z - zmove*tfrac
                      gcount(1)%influx(ipc,is) = gcount(1)%influx(ipc,is) + 1
                      gcount(1)%infhist(pbuf(m)%pid,ipc,is) = &
                     &  gcount(1)%infhist(pbuf(m)%pid,ipc,is) + 1
                      pbuf(m)%preside = IPCB
                      exit
                    end if
                  end if
!
                else if(geotype(ipc).eq.3) then
                  disp1 = pbuf(m)%x - sphere(ipc)%center(1)
                  disp2 = pbuf(m)%y - sphere(ipc)%center(2)
                  disp3 = pbuf(m)%z - sphere(ipc)%center(3)
                  tfrac1 = (disp1*disp1 + disp2*disp2 + disp3*disp3)*radsqi(ipc)
                  if(tfrac1.lt.1.0d0) then
                    disp1 = disp1 - xmove
                    disp2 = disp2 - ymove
                    disp3 = disp3 - zmove
                    tfrac2 = (disp1*disp1 + disp2*disp2 + disp3*disp3)*radsqi(ipc)
                    tfrac = (tfrac2 - tfrac1*tfrac2)/(tfrac2 - tfrac1)
                    pbuf(m)%x = pbuf(m)%x - xmove*tfrac
                    pbuf(m)%y = pbuf(m)%y - ymove*tfrac
                    pbuf(m)%z = pbuf(m)%z - zmove*tfrac
                    gcount(1)%influx(ipc,is) = gcount(1)%influx(ipc,is) + 1
                    gcount(1)%infhist(pbuf(m)%pid,ipc,is) = &
                   &  gcount(1)%infhist(pbuf(m)%pid,ipc,is) + 1
                    pbuf(m)%preside = IPCB
                    exit
                  end if
                end if
!
              else !if(ncond(ipc).le.0) then
!                if(pbuf(m)%preside.lt.0) cycle IPCL1
                if(geotype(ipc).eq.0.or.geotype(ipc).eq.1) then
                  if(pbuf(m)%x.ge.xlpc(ipc).and.pbuf(m)%x.le.xupc(ipc).and. &
                 &   pbuf(m)%y.ge.ylpc(ipc).and.pbuf(m)%y.le.yupc(ipc).and. &
                 &   pbuf(m)%z.ge.zlpc(ipc).and.pbuf(m)%z.le.zupc(ipc)) then
                    xmovei = 1.0d0/xmove
                    ymovei = 1.0d0/ymove
                    zmovei = 1.0d0/zmove
                    tfrac1 = (pbuf(m)%x - xlpc(ipc))*xmovei
                    tfrac2 = (pbuf(m)%x - xupc(ipc))*xmovei
                    tfrac3 = (pbuf(m)%y - ylpc(ipc))*ymovei
                    tfrac4 = (pbuf(m)%y - yupc(ipc))*ymovei
                    tfrac5 = (pbuf(m)%z - zlpc(ipc))*zmovei
                    tfrac6 = (pbuf(m)%z - zupc(ipc))*zmovei
                    tfrac1 = abs(tfrac1 + floor(tfrac1))
                    tfrac2 = abs(tfrac2 + floor(tfrac2))
                    tfrac3 = abs(tfrac3 + floor(tfrac3))
                    tfrac4 = abs(tfrac4 + floor(tfrac4))
                    tfrac5 = abs(tfrac5 + floor(tfrac5))
                    tfrac6 = abs(tfrac6 + floor(tfrac6))
                    tfrac1 = min(tfrac1,tfrac2,tfrac3,tfrac4,tfrac5,tfrac6)
                    if(tfrac1.gt.tfrac) then
                      tfrac = tfrac1
                      jpc = ipc
                    end if
                  end if
!
                else if(geotype(ipc).eq.2) then
                  if(cylinder(ipc)%align.eq.1) then
                    disp1 = pbuf(m)%y - cylinder(ipc)%axis(1)
                    disp2 = pbuf(m)%z - cylinder(ipc)%axis(2)
                    tfrac1 = (disp1*disp1 + disp2*disp2)*radsqi(ipc)
                    if(pbuf(m)%x.ge.cylinder(ipc)%edge(1).and. &
                   &   pbuf(m)%x.le.cylinder(ipc)%edge(2).and. &
                   &   tfrac1.lt.1.0d0) then
                      disp1 = disp1 - ymove
                      disp2 = disp2 - zmove
                      tfrac2 = (disp1*disp1 + disp2*disp2)*radsqi(ipc)
                      tfrac1 = (tfrac2 - tfrac1*tfrac2)/(tfrac2 - tfrac1)
                      xmovei = 1.0d0/xmove
                      tfrac2 = (pbuf(m)%x - cylinder(ipc)%edge(1))*xmovei
                      tfrac3 = (pbuf(m)%x - cylinder(ipc)%edge(2))*xmovei
                      tfrac1 = abs(tfrac1 + floor(tfrac1))
                      tfrac2 = abs(tfrac2 + floor(tfrac2))
                      tfrac3 = abs(tfrac3 + floor(tfrac3))
                      tfrac1 = min(tfrac1,tfrac2,tfrac3)
                      if(tfrac1.gt.tfrac) then
                        tfrac = tfrac1
                        jpc = ipc
                      end if
                    end if
!
                  else if(cylinder(ipc)%align.eq.2) then
                    disp1 = pbuf(m)%z - cylinder(ipc)%axis(1)
                    disp2 = pbuf(m)%x - cylinder(ipc)%axis(2)
                    tfrac1 = (disp1*disp1 + disp2*disp2)*radsqi(ipc)
                    if(pbuf(m)%y.ge.cylinder(ipc)%edge(1).and. &
                   &   pbuf(m)%y.le.cylinder(ipc)%edge(2).and. &
                   &   tfrac1.lt.1.0d0) then
                      disp1 = disp1 - zmove
                      disp2 = disp2 - xmove
                      tfrac2 = (disp1*disp1 + disp2*disp2)*radsqi(ipc)
                      tfrac1 = (tfrac2 - tfrac1*tfrac2)/(tfrac2 - tfrac1)
                      ymovei = 1.0d0/ymove
                      tfrac2 = (pbuf(m)%y - cylinder(ipc)%edge(1))*ymovei
                      tfrac3 = (pbuf(m)%y - cylinder(ipc)%edge(2))*ymovei
                      tfrac1 = abs(tfrac1 + floor(tfrac1))
                      tfrac2 = abs(tfrac2 + floor(tfrac2))
                      tfrac3 = abs(tfrac3 + floor(tfrac3))
                      tfrac1 = min(tfrac1,tfrac2,tfrac3)
                      if(tfrac1.gt.tfrac) then
                        tfrac = tfrac1
                        jpc = ipc
                      end if
                    end if
!
                  else if(cylinder(ipc)%align.eq.3) then
                    disp1 = pbuf(m)%x - cylinder(ipc)%axis(1)
                    disp2 = pbuf(m)%y - cylinder(ipc)%axis(2)
                    tfrac1 = (disp1*disp1 + disp2*disp2)*radsqi(ipc)
                    if(pbuf(m)%z.ge.cylinder(ipc)%edge(1).and. &
                   &   pbuf(m)%z.le.cylinder(ipc)%edge(2).and. &
                   &   tfrac1.lt.1.0d0) then
                      disp1 = disp1 - xmove
                      disp2 = disp2 - ymove
                      tfrac2 = (disp1*disp1 + disp2*disp2)*radsqi(ipc)
                      tfrac1 = (tfrac2 - tfrac1*tfrac2)/(tfrac2 - tfrac1)
                      zmovei = 1.0d0/zmove
                      tfrac2 = (pbuf(m)%z - cylinder(ipc)%edge(1))*zmovei
                      tfrac3 = (pbuf(m)%z - cylinder(ipc)%edge(2))*zmovei
                      tfrac1 = abs(tfrac1 + floor(tfrac1))
                      tfrac2 = abs(tfrac2 + floor(tfrac2))
                      tfrac3 = abs(tfrac3 + floor(tfrac3))
                      tfrac1 = min(tfrac1,tfrac2,tfrac3)
                      if(tfrac1.gt.tfrac) then
                        tfrac = tfrac1
                        jpc = ipc
                      end if
                    end if
                  end if
!
                else if(geotype(ipc).eq.3) then
                  disp1 = pbuf(m)%x - sphere(ipc)%center(1)
                  disp2 = pbuf(m)%y - sphere(ipc)%center(2)
                  disp3 = pbuf(m)%z - sphere(ipc)%center(3)
                  tfrac1 = (disp1*disp1 + disp2*disp2 + disp3*disp3)*radsqi(ipc)
                  if(tfrac1.lt.1.0d0) then
                    disp1 = disp1 - xmove
                    disp2 = disp2 - ymove
                    disp3 = disp3 - zmove
                    tfrac2 = (disp1*disp1 + disp2*disp2 + disp3*disp3)*radsqi(ipc)
                    tfrac1 = (tfrac2 - tfrac1*tfrac2)/(tfrac2 - tfrac1)
                    if(tfrac1.gt.tfrac) then
                      tfrac = tfrac1
                      jpc = ipc
                    end if
                  end if
                end if
              end if
!
              if(ipc.eq.npc.and.jpc.gt.0) then
                pbuf(m)%x = pbuf(m)%x - xmove*tfrac
                pbuf(m)%y = pbuf(m)%y - ymove*tfrac
                pbuf(m)%z = pbuf(m)%z - zmove*tfrac
                gcount(1)%influx(ipc,is) = gcount(1)%influx(ipc,is) + 1
                gcount(1)%infhist(pbuf(m)%pid,ipc,is) = &
               &  gcount(1)%infhist(pbuf(m)%pid,ipc,is) + 1
                pbuf(m)%preside = INCB
              end if
            end do
          end if
