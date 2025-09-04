!
!  setpos.fnc
!
!--------------------- transfer velocity & assign position
          if     (nemd(iepl).eq. 1) then
            pinj%vx = vxf(j1)
            pinj%vy = vyr(j2)
            pinj%vz = vzr(j2)
!            print*, "debug:",vxf(j1),vyr(j2),vzr(j2)
            pinj%x = remf(iepl) + vxf(j1)*injctisw(is)*dranu(iran-2)
            pinj%y = ymine(iepl) + ywdte(iepl)*dranu(iran-1)
            pinj%z = zmine(iepl) + zwdte(iepl)*dranu(iran  )
          else if(nemd(iepl).eq.-1) then
            pinj%vx = -vxf(j1)
            pinj%vy = vyr(j2)
            pinj%vz = vzr(j2)
            pinj%x = remb(iepl) - vxf(j1)*injctisw(is)*dranu(iran-2)
            pinj%y = ymine(iepl) + ywdte(iepl)*dranu(iran-1)
            pinj%z = zmine(iepl) + zwdte(iepl)*dranu(iran  )
          else if(nemd(iepl).eq. 2) then
            pinj%vx = vyr(j2)
            pinj%vy = vxf(j1)
            pinj%vz = vzr(j2)
            pinj%x = xmine(iepl) + xwdte(iepl)*dranu(iran-2)
            pinj%y = remf(iepl) + vxf(j1)*injctisw(is)*dranu(iran-1)
            pinj%z = zmine(iepl) + zwdte(iepl)*dranu(iran  )
          else if(nemd(iepl).eq.-2) then
            pinj%vx = vyr(j2)
            pinj%vy = -vxf(j1)
            pinj%vz = vzr(j2)
            pinj%x = xmine(iepl) + xwdte(iepl)*dranu(iran-2)
            pinj%y = remb(iepl) - vxf(j1)*injctisw(is)*dranu(iran-1)
            pinj%z = zmine(iepl) + zwdte(iepl)*dranu(iran  )
          else if(nemd(iepl).eq. 3) then
            pinj%vx = vyr(j2)
            pinj%vy = vzr(j2)
            pinj%vz = vxf(j1)
            pinj%x = xmine(iepl) + xwdte(iepl)*dranu(iran-2)
            pinj%y = ymine(iepl) + ywdte(iepl)*dranu(iran-1)
            pinj%z = remf(iepl) + vxf(j1)*injctisw(is)*dranu(iran  )
          else if(nemd(iepl).eq.-3) then
            pinj%vx = vyr(j2)
            pinj%vy = vzr(j2)
            pinj%vz = -vxf(j1)
            pinj%x = xmine(iepl) + xwdte(iepl)*dranu(iran-2)
            pinj%y = ymine(iepl) + ywdte(iepl)*dranu(iran-1)
            pinj%z = remb(iepl) - vxf(j1)*injctisw(is)*dranu(iran  )
          end if
!-------------------- 
