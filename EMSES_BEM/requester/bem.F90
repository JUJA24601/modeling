module bem
    implicit none
    integer i, j, k, l, num_of_rectangles, num_of_triangles, info, num_of_faces
    double precision apex1_x, apex1_y, apex1_z, apex2_x, apex2_y, apex2_z, apex3_x, apex3_y, apex3_z,&
                    apex4_x, apex4_y, apex4_z, a_temp, b_temp, p0_x, p0_y, p0_z, n1_x, n1_y, n1_z, n2_x, n2_y, n2_z, n3_x, n3_y, n3_z
    double precision potential_conductor, C, E0(3)
    character(len=100) file_triangle, file_rectangle
    logical triangle_exists, rectangle_exists
    double precision, parameter :: pi = acos(-1.0d0) ! 3.141592653589793238d0
    double precision, parameter :: eps0 = 1.0d0     ! 8.8541878128d-12     ! permittivity
    integer(kind=4), allocatable :: grid_point_temp(:,:)
    integer(kind=4) :: new_rows(8,3), new_row(3), x1, x2, y1, y2, z1, z2
    real(kind=8),dimension(:,:), allocatable :: coord_bem
    logical :: is_duplicate
    integer(kind=4) :: grid_point_rows = 0

    ! type of rectangle
    type Rectangle
        double precision coord(4,3), a, b, p0(3), n1(3), n2(3), n3(3)
    end type Rectangle
    ! type of triangle
    type Triangle
        double precision coord(3,3), a, b, p0(3), n1(3), n2(3), n3(3)
    end type Triangle

    type(Rectangle), allocatable :: rectangles(:)
    type(Triangle), allocatable :: triangles(:)
    double precision, allocatable :: coef(:,:,:,:)
    
contains
    subroutine hello()
        print *, "Hello from bem"
    end subroutine hello

    ! subroutine calc_coef()
    !     do i = -1, nx/nodes(1)+1
    !         do j = -1, ny/nodes(2)+1
    !             do k = -1, nz/nodes(3)+1
    !                 do l = 1, num_of_faces
    !                     coef(i,j,k,l) = potential_elem_tri(                     &
    !                         triangles(l),                                         &
    !                         [                                                     &
    !                             real(i+mod(myid,nodes(1))*nx/nodes(1),8),           &
    !                             real(j+mod(myid/nodes(1),nodes(2))*ny/nodes(2),8),  &
    !                             real(k+(myid/(nodes(1)*nodes(2)))*nz/nodes(3),8)    & 
    !                         ],                                                  &
    !                         1.0d0                                                 &
    !                         )
    !                 end do
    !             end do
    !         end do
    !     end do
    ! end subroutine calc_coef

    subroutine calc_gridpoint()

        allocate(coord_bem(num_of_faces,3))
        ! 電位の知りたい点を計算する
        open (unit=24601,file='./check_1.txt' ,action="write",status="replace")
        do i = 1, num_of_faces
            if (i<=num_of_rectangles) then
                coord_bem(i,:) = centroid_rectangle(rectangles(i))
            else
                coord_bem(i,:) = centroid_triangle(triangles(i-num_of_rectangles))
            end if
            write (24601,*) coord_bem(i,:)
        end do
        close (24601)

        allocate(grid_point_temp(0, 3))

        open (unit=24601,file='./check_esses.txt' ,action="write",status="replace")

        ! 電位の知りたい点が含まる格子点を記録する
        do i = 1, num_of_faces
            x1 = floor(coord_bem(i,1))
            x2 = ceiling(coord_bem(i,1))
            y1 = floor(coord_bem(i,2))
            y2 = ceiling(coord_bem(i,2))
            z1 = floor(coord_bem(i,3))
            z2 = ceiling(coord_bem(i,3))
            if (x1 == x2) x2 = x2 + 1
            if (y1 == y2) y2 = y2 + 1
            if (z1 == z2) z2 = z2 + 1

            new_rows(1,:) = [x1, y1, z1]
            new_rows(2,:) = [x1, y1, z2]
            new_rows(3,:) = [x1, y2, z1]
            new_rows(4,:) = [x1, y2, z2]
            new_rows(5,:) = [x2, y1, z1]
            new_rows(6,:) = [x2, y1, z2]
            new_rows(7,:) = [x2, y2, z1]
            new_rows(8,:) = [x2, y2, z2]

            do k = 1, 8
                new_row = new_rows(k,:)
                is_duplicate = .false.
                do j = 1, grid_point_rows
                    if (all(grid_point_temp(j,:) == new_row)) then
                        is_duplicate = .true.
                        exit
                    end if
                end do

                if (.not. is_duplicate) then
                    write (24601,*) new_row
                    call add_row(grid_point_temp, grid_point_rows, 3, new_row)
                    ! print *, grid_point_rows
                end if
            end do

        end do

        close (24601)

        ! print *, grid_point_rows

        deallocate(coord_bem)
        deallocate(grid_point_temp)
    end subroutine calc_gridpoint

    subroutine add_row(a, n, c, row)
        integer, allocatable, intent(inout) :: a(:,:)
        integer, intent(inout) :: n
        integer, intent(in) :: c
        integer, intent(in) :: row(c)
        integer, allocatable :: temp(:,:)

        allocate(temp(n + 1, c))
        if (n > 0) temp(1:n, :) = a
        temp(n + 1, :) = row

        deallocate(a)
        allocate(a(n + 1, c))
        a = temp

        deallocate(temp)
        n = n + 1
    end subroutine add_row

    subroutine read_triangle(file_triangle)
        character(len=*), intent(in) :: file_triangle
        open (17, file=file_triangle, status="old")
        num_of_triangles = 0
        read (17, "()")
        do
            read (17, *, end=101) apex1_x, apex1_y, apex1_z, apex2_x, apex2_y, apex2_z, apex3_x, apex3_y, apex3_z, &
            a_temp, b_temp, p0_x, p0_y, p0_z, n1_x, n1_y, n1_z, n2_x, n2_y, n2_z, n3_x, n3_y, n3_z
            num_of_triangles = num_of_triangles + 1
        end do
        101 continue
        rewind(17)
        allocate(triangles(num_of_triangles))
        ! print *, "num_of_triangles =", num_of_triangles
        read (17, '()')
        do i = 1, num_of_triangles
            read (17, *) apex1_x, apex1_y, apex1_z, apex2_x, apex2_y, apex2_z, apex3_x, apex3_y, apex3_z, &
            a_temp, b_temp, p0_x, p0_y, p0_z, n1_x, n1_y, n1_z, n2_x, n2_y, n2_z, n3_x, n3_y, n3_z
            triangles(i)%coord(1,:) = [apex1_x, apex1_y, apex1_z]
            triangles(i)%coord(2,:) = [apex2_x, apex2_y, apex2_z]
            triangles(i)%coord(3,:) = [apex3_x, apex3_y, apex3_z]
            triangles(i)%a          = a_temp
            triangles(i)%b          = b_temp
            triangles(i)%p0(:) = [p0_x, p0_y, p0_z]
            triangles(i)%n1(:) = [n1_x, n1_y, n1_z]
            triangles(i)%n2(:) = [n2_x, n2_y, n2_z]
            triangles(i)%n3(:) = [n3_x, n3_y, n3_z]
        end do
        close (17)
    end subroutine read_triangle

    subroutine read_rectangle(file_rectangle)
        character(len=*), intent(in) :: file_rectangle
        open (17, file="data_rectangle.csv", status="old")
        num_of_rectangles = 0
        read (17, "()")
        do
            read (17, *, end=100) apex1_x, apex1_y, apex1_z, apex2_x, apex2_y, apex2_z, apex3_x, apex3_y, apex3_z, &
            apex4_x, apex4_y, apex4_z, a_temp, b_temp, p0_x, p0_y, p0_z, n1_x, n1_y, n1_z, n2_x, n2_y, n2_z, n3_x, n3_y, n3_z
            num_of_rectangles = num_of_rectangles + 1
        end do
        100 continue
        rewind(17)
        allocate(rectangles(num_of_rectangles))
        ! print *, "num_of_rectangles =", num_of_rectangles
        read (17, '()')
        do i = 1, num_of_rectangles
            read (17, *) apex1_x, apex1_y, apex1_z, apex2_x, apex2_y, apex2_z, apex3_x, apex3_y, apex3_z, &
            apex4_x, apex4_y, apex4_z, a_temp, b_temp, p0_x, p0_y, p0_z, n1_x, n1_y, n1_z, n2_x, n2_y, n2_z, n3_x, n3_y, n3_z
            rectangles(i)%coord(1,:) = [apex1_x, apex1_y, apex1_z]
            rectangles(i)%coord(2,:) = [apex2_x, apex2_y, apex2_z]
            rectangles(i)%coord(3,:) = [apex3_x, apex3_y, apex3_z]
            rectangles(i)%coord(4,:) = [apex4_x, apex4_y, apex4_z]
            rectangles(i)%a          = a_temp
            rectangles(i)%b          = b_temp
            rectangles(i)%p0(:) = [p0_x, p0_y, p0_z]
            rectangles(i)%n1(:) = [n1_x, n1_y, n1_z]
            rectangles(i)%n2(:) = [n2_x, n2_y, n2_z]
            rectangles(i)%n3(:) = [n3_x, n3_y, n3_z]
        end do
        close (17)
    end subroutine read_rectangle

    function centroid_rectangle(rect) result(centroid)
        !
        !
        ! 四角形の中心を求める
        ! rect      : 四角形の型
        ! centorid  : 四角形の中心
        !
        !
        implicit none
        type(Rectangle), intent(in) :: rect
        double precision centroid(3), a, b, p0(3), n1(3), n2(3)
        a  = rect%a
        b  = rect%b
        p0 = rect%p0
        n1 = rect%n1
        n2 = rect%n2
        centroid = p0 + 0.5*(a*n1 + b*n2)
    end function

    function centroid_triangle(tri) result(centroid)
        !
        !
        ! 三角形の内心を求める．
        !
        ! Parameters
        ! -----------
        ! triangle : [3,3]の実数配列でそれぞれの行に頂点が保存される
        !            三角形の座標
        !
        ! Returns
        ! -----------
        ! incenter : [1,3]の実数配列
        !            三角形の内心
        !
        !
        implicit none
        type(Triangle), intent(in) :: tri
        double precision tri_coord(3,3), centroid(3), a(3), b(3), c(3), ab, bc, ca
        tri_coord = tri%coord
        a = tri_coord(1,:); b = tri_coord(2,:); c = tri_coord(3,:)
        ab = distance(a,b); bc = distance(b,c); ca = distance(c,a)
        centroid = (bc*a+ca*b+ab*c)/(ab+bc+ca)
    end function


    function distance(a, b) result(dist)
        !
        !
        ! 2点間 a,b の距離を求める．
        !
        !
        implicit none
        double precision a(3), b(3), dist
        dist = sqrt(sum((a-b)*(a-b)))
    end function


    function p_calc_rect(rect, r) result(ans)
        ! 
        !
        ! 電位係数を求める．(4角形がrに作る)
        !
        ! Parameters
        ! -----------
        ! rectangle : 4角形の構造体
        ! r   : 任意の点
        !
        ! Returns
        ! -----------
        ! p_calc : 実数
        !          電位係数
        !
        !
        implicit none
        double precision ans, r(3), a, b, p0(3), n1(3), n2(3), n3(3), xmin, xmax,&
                        & ymin, ymax, p(3), uP, vP, w
        type(Rectangle) rect

        a  = rect%a
        b  = rect%b
        p0 = rect%p0
        n1 = rect%n1
        n2 = rect%n2
        n3 = rect%n3

        p  = r - p0
        uP = sum(p*n1)
        vP = sum(p*n2)
        w  = sum(p*n3)

        xmin = -uP
        xmax = -uP + a
        ymin = -vP
        ymax = -vP + b

        ans = Integral_ln(xmax, ymax, w) - Integral_ln(xmin, ymax, w) - Integral_ln(xmax, ymin, w) + Integral_ln(xmin, ymin, w)
        
    end function

    function Integral_ln(x, y, w) result(ans)
        implicit none
        double precision x, y, w, ans, r, r0, xa, c1, c2, c3
        r  = sqrt(abs(x*x + y*y + w*w))
        r0 = sqrt(abs(y*y + w*w))
        xa = abs(x)
        if (xa < 1.0d-20) then
            c1 = 0.0d0
        else
            c1 = xa*log(abs(y + r) + 1.0d-22)
        end if
        if (abs(y) < 1.0d-22) then
            c2 = 0.0d0
        else
            c2 = y*log(abs((xa + r)/r0) + 1.0d-22)
        end if
        if (abs(w) < 1.0d-22) then
            c3 = 0.0d0
        else
            c3 = w*(atan(xa/w) + atan(y*w/(x*x + w*w + xa*r)) - atan(y/w))
        end if
        ans = c1 + c2 - xa + c3
        if (x<0.0d0) then
            ans = -ans
        end if
    end function

    function p_calc_tri(tri, r) result(ans)
        ! 
        !
        ! 電位係数を求める．(三角形がrに作る)
        !
        ! Parameters
        ! -----------
        ! tri : 三角形の構造体
        ! r   : 任意の点
        !
        ! Returns
        ! -----------
        ! p_calc : 実数
        !          電位係数
        !
        !
        implicit none
        double precision r(3), a, b, p0(3), n1(3), n2(3), n3(3), x1, x2, x3, ans,&
                        & y1, y2, z, a1, a2, b1, b2, u1, u2, p(3), n1dotn2, N2prime(3), n2dotn2prime
        type(Triangle) tri

        ! 三角形の情報を抜き出す
        a  = tri%a
        b  = tri%b
        p0 = tri%p0
        n1 = tri%n1
        n2 = tri%n2
        n3 = tri%n3
        
        p = p0 - r

        n1dotn2 = sum(n1*n2)
        N2prime = (n2 - (n1dotn2 * n1))/sqrt(sum((n2 - (n1dotn2 * n1))*(n2 - (n1dotn2 * n1))))
        n2dotn2prime = sum(n2*N2prime)
        x1 = sum(p*n1)
        y1 = sum(p*N2prime)
        z = -sum(p*n3)

        x2 = x1 + a

        x3 = x1 + b * n1dotn2

        if (z < 0.0d0) then
            y1 = -y1
            y2 = y1 - b * n2dotn2prime
        else
            y2 = y1 + b * n2dotn2prime
        end if

        z = abs(z)

        if (z > 1.0d-24) then
            u1 = y1 / z
            if (u1 == 0.0d0) then
                u1 = 1.0d-28
            end if
            u2 = y2 / z
            if (u2 == 0.0d0) then
                u2 = 1.0d-28
            end if
            a1 = (x1 * y2 - x3 * y1) / (z * (y2 - y1))
            a2 = (x2 * y2 - x3 * y1) / (z * (y2 - y1))
        else
            u1 = y1
            if (u1 == 0.0d0) then
                u1 = 1.0d-28
            end if
            u2 = y2
            if (u2 == 0.0d0) then
                u2 = 1.0d-28
            end if
            a1 = (x1 * y2 - x3 * y1) / (y2 - y1)
            a2 = (x2 * y2 - x3 * y1) / (y2 - y1)
        end if

        b1 = (x3 - x1) / (y2 - y1)
        b2 = (x3 - x2) / (y2 - y1)

        if (z > 1.0d-24) then
            if (abs(b1) < 1.0d-23) then
                ans = z * (I1(a2, b2, u1, u2) - I2(a1, u1, u2))
            else if (abs(b2) < 1.0d-23) then
                ans = z * (I2(a2, u1, u2) - I1(a1, b1, u1, u2))
            else
                ans = z * (I1(a2, b2, u1, u2) - I1(a1, b1, u1, u2))
            end if
        else
            ans = (p_calc_noZ(a2, b2, a1, b1, u2) - p_calc_noZ(a2, b2, a1, b1, u1))
        end if

        ans = abs(ans)
        
    end function

    function p_calc_noZ(a2, b2, a1, b1, y) result(ans)
        implicit none
        double precision a1, a2, b1, b2, y, logArg1, logArg2, ans1, ans2, ans
        logArg2 = (1.0d0 + b2*b2)*y + a2*b2 + sqrt(1.0d0 + b2*b2)*sqrt((1.0d0 + b2*b2)*y*y + 2*a2*b2*y + a2*a2)
        logArg1 = (1.0d0 + b1*b1)*y + a1*b1 + sqrt(1.0d0 + b1*b1)*sqrt((1.0d0 + b1*b1)*y*y + 2*a1*b1*y + a1*a1)
        
        if (logArg2 > 1.0d-24) then
            if (abs(y) > 1.0d-24) then
                ans2 = y * asinh((a2 + b2 * y) / abs(y)) + a2 / sqrt(1.0d0 + b2 * b2) * log(logArg2)
            else
                ans2 = a2 / sqrt(1.0d0 + b2 * b2) * log(logArg2)
            end if
        else
            if (abs(y) > 1.0d-24) then
                ans2 = y * asinh(y * b2 / abs(y))
            else
                ans2 = 0.0d0
            end if
        end if

        if (logArg1 > 1.0d-24) then
            if (abs(y) > 5.0d-24) then
                ans1 = y * asinh((a1 + b1 * y) / abs(y)) + a1 / sqrt(1.0d0 + b1 * b1) * log(logArg1)
            else
                ans1 = a1 / sqrt(1.0d0 + b1 * b1) * log(logArg1)
            end if
        else
            if (abs(y) > 5.0d-24) then
                ans1 = y * asinh(y * b1 / abs(y))
            else
                ans1 = 0.0d0
            end if
        end if

        ans = ans2 - ans1

    end function

    function F1(a, b, u) result(ans)
        implicit none
        double precision a, b, u, ans
        ans = u*asinh((a + b * u) / sqrt(u * u + 1.0d0))
    end function

    function I3(a, b, u1, u2) result(ans)
        implicit none
        double precision a, b, u1, u2, g1, g2, ans
        g1 = (sqrt(b * b + 1.0d0) * sqrt(a*a + 2.0d0*a*b*u1 + (b*b + 1.0d0)*u1*u1 + 1.0d0) + b*(a + b*u1) + u1)
        g2 = (sqrt(b * b + 1.0d0) * sqrt(a*a + 2.0d0*a*b*u2 + (b*b + 1.0d0)*u2*u2 + 1.0d0) + b*(a + b*u2) + u2)

        if (g1 <= 0.0d0) then
            g1 = -(1.0d0 + a*a + b*b) / (2.0d0*(b*b + 1.0d0)*u1)
        end if
        if (g2 <= 0.0d0) then
            g2 = -(1.0d0 + a*a + b*b) / (2.0d0*(b*b + 1.0d0)*u2)
        end if

        if (abs(g1) < 1.0d-22) then
            if (abs(a) < 1.0d-24) then
                g1 = 1.0d-22
            end if
        end if
        if (abs(g2) < 1.0d-22) then
            if (abs(a) < 1.0d-24) then
                g2 = 1.0d-22
            end if
        end if
                
        ans = a / sqrt(b*b + 1.0d0) * log(abs(g2 / g1))
    end function

    function I4_plus(alpha, gamma, q2, prefac, t1, t2) result(ans)
        implicit none
        double precision alpha, gamma, q, q2, prefac, t1, t2, g1, g2, ans
        q = sqrt(q2)
        g1 = sqrt(gamma*t1*t1 + alpha)
        g2 = sqrt(gamma*t2*t2 + alpha)

        if (t1 > 1.0d25 .or. t2 > 1.0d25) then
            if (t2 < 1.0d25) then
                ans = (prefac*1.0d0/q*(atan(g2/q) - pi/2.0d0))
                return
            else if (t1 < 1.0d25) then
                ans = (prefac*1.0d0/q*(pi/2.0d0 - atan(g1/q)))
                return
            else
                ans = 0.0d0
                return
            end if
        end if

        ans = prefac*1.0d0/q*atan(q*(g2 - g1)/(q2 + g1*g2))

    end function

    recursive function I4(a, b, u1, u2) result(ans)
        implicit none
        double precision a, b, u1, u2, alpha, gamma, q2, prefac, t1, t2, sign, ans
        if (abs(u1 - b/a) < 1.0d-24) then
            if (u2 > b/a) then
                ans = I4(a, b, u1 + 1.0d-24, u2)
                return
            else
                ans = I4(a, b, u1 - 1.0d-24, u2)
                return
            end if
        end if

        if (abs(u2 - b/a) < 1.0d-24) then
            if (u1 > b/a) then
                ans = I4(a, b, u1, u2 + 1.0d-24)
                return
            else
                ans = I4(a, b, u1, u2 - 1.0d-24)
                return
            end if
        end if
        
        if (abs(a) < 1.0d-24) then
            if (a > 0.0d0) then
                ans = I4(a + 1.0d-24, b, u1, u2)
                return
            else
                ans = I4(a - 1.0d-24, b, u1, u2)
                return
            end if
        end if

        alpha = 1.0d0 + (a*a)/(b*b)
        gamma = (a*a + b*b)*(a*a + b*b + 1.0d0)/(b*b)
        q2 = (a*a + b*b)*(a*a + b*b)/(b*b)
        prefac = (a*a/b + b)
        if (a*u1 /= b) then
            t1 = (b*u1 + a)/(a*u1 - b)
        else
            t1 = 1.0d25
        end if
        if (a*u2 /= b) then
            t2 = (b*u2 + a)/(a*u2 - b)
        else
            t2 = 1.0d25
        end if

        sign = 1.0d0

        if (a < 0.0d0) then
            sign = - sign
        end if
        if (b < 0.0d0) then
            sign = -sign
        end if
        if (u1 > b/a) then
            sign = -sign
        end if

        if (((u1<b/a) .and. (b/a<u2)) .or. ((u2<b/a) .and. (b/a<u1))) then
            ans = sign * (I4_plus(alpha, gamma, q2, prefac, t1, 1.0d26) + I4_plus(alpha, gamma, q2, prefac, t2, 1.0d26))
            return
        else
            ans = sign * I4_plus(alpha, gamma, q2, prefac, t1, t2)
            return
        end if
    end function I4

    function I1(a, b, u1, u2) result(ans)
        implicit none
        double precision a, b, u1, u2, ans
        ans = F1(a, b, u2) - F1(a, b, u1) + I3(a, b, u1, u2) - I4(a, b, u1, u2)
    end function

    function I6(x, u1, u2) result(ans)
        implicit none
        double precision x, u1, u2, ans
        if (abs(x) < 1.0d-25) then
            ans = 0.0d0
            return
        end if
        ans = x*log((sqrt(u2*u2 + x*x + 1.0d0) + u2)/(sqrt(u1*u1 + x*x + 1.0d0) + u1))
    end function


    function I7(x, u1, u2) result(ans)
        implicit none
        double precision x, u1, u2, t1, t2, g1, g2, ans
        if (abs(u1) > 1.0d-26) then
            t1 = 1.0d0/u1
        else
            t1 = 1.0d26
        end if
        if (abs(u2) > 1.0d-26) then
            t2 = 1.0d0/u2
        else
            t2 = 1.0d26
        end if

        g1 = sqrt(1.0d0 + t1*t1*(1.0d0 + x*x))
        g2 = sqrt(1.0d0 + t2*t2*(1.0d0 + x*x))

        ans = atan(x*(g2 - g1)/(x*x + g2*g1))
    end function

    function I2(x, u1, u2) result(ans)
        implicit none
        double precision x, u1, u2, ans
        if (((u1<0.0d0) .and. (0.0d0<u2)) .or. ((u2<0.0d0) .and. (0.0d0<u1))) then
            if (u1 <= 0.0d0) then
                ans = (F1(x, 0.0d0, u2) - F1(x, 0.0d0, u1)) + I6(x, u1, u2) + I7(x, 0.0d0, abs(u1)) + I7(x, 0.0d0, abs(u2))
            else
                ans = (F1(x, 0.0d0, u2) - F1(x, 0.0d0, u1)) + I6(x, u1, u2) + I7(x, abs(u1), 0.0d0) + I7(x, abs(u2), 0.0d0)
            end if
        else if (u1 <= 0.0d0) then
            ans = (F1(x, 0.0d0, u2) - F1(x, 0.0d0, u1)) + I6(x, u1, u2) + I7(x, u2, u1)
        else
            ans = (F1(x, 0.0d0, u2) - F1(x, 0.0d0, u1)) + I6(x, u1, u2) + I7(x, u1, u2)
        end if
    end function

    function p_calc_o(tri) result(ans)
        !
        !
        ! 三角形要素が自身の内心に作る電位の電位係数
        !
        ! Parameters
        ! -----------
        ! triangle : [3,3]の実数行列．それぞれの行に三角形の頂点．
        !            三角形の座標
        !
        ! Returns
        ! -----------
        ! p_calc_o : 実数
        !            電位係数
        !
        !
        implicit none
        double precision tri_coord(3,3), a(3), b(3), c(3), ab(3), ac(3), ba(3), bc(3), ca(3), cb(3),&
                        & theta_a, theta_b, theta_c, ans
        type(Triangle) tri
        tri_coord = tri%coord
        a = tri_coord(1,:); b = tri_coord(2,:); c = tri_coord(3,:)
        ab = b-a; ac = c-a; ba = a-b; bc = c-b; ca = a-c; cb = b-c
        theta_a = acos(sum(ab*ac)/(distance(a,b)*distance(a,c)))
        theta_b = acos(sum(bc*ba)/(distance(b,c)*distance(b,a)))
        theta_c = acos(sum(ca*cb)/(distance(c,a)*distance(c,b)))
        ans = 2.0d0*(area_tri(tri)*log(((1.0d0+cos(theta_a/2.0d0))*(1.0d0+cos(theta_b/2.0d0))*(1.0d0+cos(theta_c/2.0d0))/&
            & ((1.0d0-cos(theta_a/2.0d0))*(1.0d0-cos(theta_b/2.0d0))*(1.0d0-cos(theta_c/2.0d0))))))/&
            & (distance(a,b)+distance(b,c)+distance(c,a))
    end function

    function potential_tri(r, tris, sigmas) result(ans)
        !
        !
        ! 位置rでの電位を計算する．
        ! 
        ! Parameters
        ! -------------
        ! r         : [1,3]の実数配列
        !             rの座標
        ! tris      : 三角形の構造体を収めた配列
        !       `     三角形の座標
        ! sigmas    : [n,1]の実数配列．nは三角形の個数．
        !             それぞれの三角形の電荷密度
        ! 
        ! Returns
        ! -------------
        ! potential : 実数
        !             rでの電位
        !
        !
        implicit none
        double precision r(3), sigmas(:), ans
        integer i
        type(Triangle) tris(:), tri
        ans = 0.0d0
        do i = 1, num_of_triangles
            tri = tris(i)
            ans = ans + potential_elem_tri(tri,r,sigmas(i))
        end do
    end function

    function potential_elem_tri(tri, r, sigma) result(ans)
        !
        !
        ! 三角形要素triangleによってできる任意の点rでの電位を求める．
        ! 
        ! Parameters
        ! ---------
        ! tri      : [3,3]の実数配列でそれぞれの行に頂点が保存される．
        !            三角形の座標
        ! r        : [1,3]の実数配列
        !            任意の点
        ! sigma    : 実数
        !            三角形の電荷密度
        !
        ! Returns
        ! ---------
        ! potential_elem : 実数
        !                  三角形要素によって点rで生じる電位
        ! 
        !
        implicit none
        double precision r(3), sigma, ans
        type(Triangle) tri
        ans = (p_calc_tri(tri, r)*sigma)/(4.0d0*pi*eps0)
    end function

    function potential_elem_rect(rect, r, sigma) result(ans)
        implicit none
        double precision r(3), sigma, ans
        type(Rectangle) rect
        ans = (p_calc_rect(rect, r)*sigma)/(4.0d0*pi*eps0)
    end function

    function potential_rect(r, rects, sigmas) result(ans)
        implicit none
        double precision r(3), sigmas(:), ans
        integer i
        type(Rectangle) rects(:), rect
        ans = 0.0d0
        do i = 1, num_of_rectangles
            rect = rects(i)
            ans = ans + potential_elem_rect(rect,r,sigmas(i))
        end do
    end function

    ! function cross_product(a,b)
    !     implicit none
    !     double precision a(3), b(3), cross_product(1,3)
    !     cross_product(1,1) = a(2)*b(3) - a(3)*b(2)
    !     cross_product(1,2) = a(3)*b(1) - a(1)*b(3)
    !     cross_product(1,3) = a(1)*b(2) - a(2)*b(1)
    ! end function cross_product

    function area_rect(rect) result(ans)
        implicit none
        double precision a, b, ans
        type(Rectangle) rect
        a  = rect%a
        b  = rect%b
        ans = a*b
    end function

    function area_tri(tri) result(ans)
        implicit none
        double precision tri_coord(3,3), ans, a(3), b(3), o(3)
        type(Triangle) tri
        tri_coord = tri%coord
        o = [0.0d0, 0.0d0, 0.0d0]
        a = tri_coord(2,:) - tri_coord(1,:)
        b = tri_coord(3,:) - tri_coord(1,:)
        ans = 0.5*sqrt((distance(a,o)**2)*(distance(b,o)**2) - sum(a*b)**2)
    end function
    
end module bem