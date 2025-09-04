            if(pbuf(m)%preside.eq.INCB.or.pbuf(m)%preside.eq.IPCB) then
              if(nyield(1).ne.0) then
!              if(.false.) then
                call RANU0(dranu,abs(nyield(1))*3,icon)
                iran = 0
                do mm=1,abs(nyield(1))
                  iran = iran + 3
                  pinj(1)%x = pbuf(m)%x
                  pinj(1)%y = pbuf(m)%y
                  pinj(1)%z = pbuf(m)%z
                  addr0 = nprs + npr(isse(is))*dranu(iran-2)
                  addr1 = nprs + npr(isse(is))*dranu(iran-1)
                  addr2 = nprs + npr(isse(is))*dranu(iran)
                  pinj(1)%vx = vtangf(addr0)
                  pinj(1)%vy = vtangf(addr1)
                  pinj(1)%vz = vnormf(addr2)
!
                  pinj(1)%nid = &
                 &  oh3_map_particle_to_subdomain(pinj(1)%x,pinj(1)%y,pinj(1)%z)
                  pinj(1)%spec = isse(is)
                  pinj(1)%pid = 0
                  pinj(1)%preside = -9
                  call oh2_inject_particle(pinj(1))
                  totalp(is,3) = totalp(is,3) + 1
!
                  if(npc.ge.2.and.nyield(1).gt.0) then
                    pinj(1)%preside = RPCB
                    nsecemit(1) = nsecemit(1) + abs(nyield(1))
                  else if(npc.ge.2.and.nyield(1).lt.0) then
                    pinj(1)%preside = RPCB
                    nsecemit(2) = nsecemit(2) + abs(nyield(1))
                  else if(npc.ge.1.and.nyield(1).gt.0) then
                    pinj(1)%preside = RPCB
                    nsecemit(1) = nsecemit(1) + abs(nyield(1))
                  else
                    pinj(1)%preside = RNCB
                  end if
                  call oh2_inject_particle(pinj(1))
                  totalp(is,3) = totalp(is,3) + 1
                end do
              else if(nyield(2).gt.0) then
                call RANU0(dranu,nyield(2)*3,icon)
                iran = 0
                do mm=1,nyield(2)
                  iran = iran + 3
                  pinj(1)%x = pbuf(m)%x
                  pinj(1)%y = pbuf(m)%y
                  pinj(1)%z = pbuf(m)%z
                  addr0 = nprs + npr(isse(is))*dranu(iran-2)
                  addr1 = nprs + npr(isse(is))*dranu(iran-1)
                  addr2 = nprs + npr(isse(is))*dranu(iran)
                  pinj(1)%vx = sign(vnormf(addr0),-pbuf(m)%vx)
                  pinj(1)%vy = vtangf(addr1)
                  pinj(1)%vz = vtangf(addr2)
!
                  pinj(1)%nid = &
                 &  oh3_map_particle_to_subdomain(pinj(1)%x,pinj(1)%y,pinj(1)%z)
                  pinj(1)%spec = isse(is)
                  pinj(1)%pid = 0
                  pinj(1)%preside = -9
                  call oh2_inject_particle(pinj(1))
                  totalp(is,3) = totalp(is,3) + 1
!
                  pinj(1)%preside = RNCB
                  nsecemit(3) = nsecemit(3) + abs(nyield(2))
                  call oh2_inject_particle(pinj(1))
                  totalp(is,3) = totalp(is,3) + 1
                end do
              end if
            end if
