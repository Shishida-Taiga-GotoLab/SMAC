module global
    implicit none
    ! 時間
    real(8), parameter ::dt = 0.001d0
    integer, parameter :: steps = 1000

    ! 格子数
    integer, parameter :: Nx = 10
    integer, parameter :: Ny = 10

    ! 計算領域サイズ
    real(8), parameter :: Lx = 2.0d0 * acos(-1.0d0)
    real(8), parameter :: Ly = 2.0d0 * acos(-1.0d0)
    real(8) :: dx = Lx / real(Nx, kind=8)
    real(8) :: dy = Ly / real(Ny, kind=8)

    ! 粘性係数μ
    real(8), parameter :: mu = 10.0d0 ** (-3)

end module global

module subprogs
    use global
    implicit none
    contains

    ! 対流項の計算
    subroutine convection_cal(u, v, Ax, Ay)
        real(8), intent(in) :: u(0:Nx+1, 0:Ny+1), v(0:Nx+1, 0:Ny+1)
        real(8), intent(out) :: Ax(1:Nx, 1:Ny), Ay(1:Nx, 1:Ny-1)
        real(8) :: duudx, duvdy, dvvdy, duvdx
        real(8) :: re_dx, re_dy
        integer :: i, j
        re_dx = 1.0d0 / dx
        re_dy = 1.0d0 / dy

        ! Axを計算
        do j = 1, Ny
            do i = 1, Nx
                duudx = (- ((u(i-1, j)+u(i, j)) * 0.5d0)**2 + ((u(i, j)+u(i+1, j)) * 0.5d0)**2) * re_dx
                duvdy = (- ((v(i, j-1)+v(i+1, j-1)) * 0.5d0)*((u(i, j-1)+u(i, j)) * 0.5d0) + ((v(i, j)+v(i+1, j)) * 0.5d0)*((u(i, j)+u(i, j+1)) * 0.5d0)) * re_dy
                Ax(i, j) = - (duudx + duvdy)
            enddo
        enddo

        ! Ayを計算
        do j = 1, Ny-1
            do i = 1, Nx
                dvvdy = (- ((v(i, j-1)+v(i, j)) * 0.5d0)**2 + ((v(i, j)+v(i, j+1)) * 0.5d0)**2) * re_dy
                duvdx = (- ((u(i-1, j)+u(i-1, j+1)) * 0.5d0)*((v(i-1, j)+v(i, j)) * 0.5d0) + ((u(i, j)+u(i, j+1)) * 0.5d0)*((v(i, j)+v(i+1, j)) * 0.5d0)) * re_dx
                Ay(i, j) = - (duvdx + dvvdy)
            enddo
        enddo

    end subroutine convection_cal
    

    ! 粘性項の計算
    subroutine viscosity_cal(u, v, Bx, By) ! 粘性項
        real(8), intent(in) :: u(0:Nx+1, 0:Ny+1), v(0:Nx+1, 0:Ny+1)
        real(8), intent(out) :: Bx(1:Nx, 1:Ny), By(1:Nx, 1:Ny-1)
        real(8) :: dudxx, dudyy, dvdxx, dvdyy
        real(8) :: re_dx, re_dy
        integer :: i, j
        re_dx = 1.0d0 / dx**2
        re_dy = 1.0d0 / dy**2

        ! Bxを計算
        do j = 1, Ny
            do i = 1, Nx
                dudxx = (u(i-1, j) - 2.0d0*u(i, j) + u(i+1, j)) * re_dx
                dudyy = (u(i, j-1) - 2.0d0*u(i, j) + u(i, j+1)) * re_dy
                Bx(i, j) = mu * (dudxx + dudyy)
            enddo
        enddo

        ! Byを計算
        do j = 1, Ny-1
            do i = 1, Nx
                dvdxx = (v(i-1, j) - 2.0d0*v(i, j) + v(i+1, j)) * re_dx
                dvdyy = (v(i, j-1) - 2.0d0*v(i, j) + v(i, j+1)) * re_dy
                By(i, j) = mu * (dvdxx + dvdyy)
            enddo
        enddo

    end subroutine viscosity_cal


    ! 予測速度を計算
    subroutine prespeed_cal(Ax, Ay, Bx, By, u, u_pre, v, v_pre, dpdx, dpdy, p)
        real(8), intent(inout) :: Ax(1:Nx, 1:Ny), Ay(1:Nx, 1:Ny-1), Bx(1:Nx, 1:Ny), By(1:Nx, 1:Ny-1)
        real(8), intent(out) :: dpdx(0:Nx+1, 0:Ny+1), dpdy(0:Nx+1, 0:Ny+1)
        real(8), intent(out) :: u_pre(0:Nx+1, 0:Ny+1), v_pre(0:Nx+1, 0:Ny+1)
        real(8), intent(in) :: u(0:Nx+1, 0:Ny+1), v(0:Nx+1, 0:Ny+1), p(0:Nx+1, 0:Ny+1)
        integer :: i, j
        real(8) :: fx = 10.0d0, fy = 0.0d0


        !call convection_cal(u, v, Ax, Ay) ! 対流項を計算
        Ax(:, :) = 0.0d0 ! ポアズイユ流では対流項は0
        Ay(:, :) = 0.0d0
        call viscosity_cal(u, v, Bx, By) ! 粘性項を計算

        do j = 1, Ny
            do i = 1, Nx
                dpdx(i, j) = (-p(i, j) + p(i+1, j)) / dx ! 圧力勾配（スタガード格子上）
                u_pre(i, j) = u(i, j) + dt * (Ax(i, j) + Bx(i, j) - dpdx(i, j))
            enddo
        enddo

        do j = 1, Ny-1
            do i = 1, Nx
                dpdy(i, j) = (-p(i, j) + p(i, j+1)) / dy ! 圧力勾配（スタガード格子上）
                v_pre(i, j) = v(i, j) + dt * (Ay(i, j) + By(i, j) - dpdy(i, j))
            enddo
        enddo

    end subroutine prespeed_cal

    subroutine SOR(BXP, BXM, BYP, BYM, B0, u_pre, v_pre, phi) ! SOR法を用いてポアソン方程式を計算
        real(8), intent(in) :: BXP, BXM, BYP, BYM, B0
        real(8), intent(in) :: u_pre(0:Nx+1, 0:Ny+1), v_pre(0:Nx+1, 0:Ny+1)
        real(8), intent(inout) :: phi(0:Nx+1, 0:Ny+1)
        real(8) ::  er, er0, RHS(1:Nx, 1:Ny), RHS_sum, omg, E(1:Nx, 1:Ny) ! RHS : 右辺　r : 離散残差
        integer :: i, j, itr, itr_max

        itr_max = 10000 ! 最大繰り返し回数
        er = 0.0d0 ! 誤差を初期化
        er0 = 10e-8 ! 収束判定誤差
        omg = 1.7d0 ! 緩和係数
        RHS(:, :) = 0.0d0
        RHS_sum = 0.0d0

        do itr = 1, itr_max
            er = 0.0d0 ! 誤差をリセット

            do j = 1, Ny
                do i = 1, Nx
                    RHS(i, j) = ((-u_pre(i-1,j)+u_pre(i,j))/dx + (-v_pre(i, j-1)+v_pre(i, j))/dy) / dt ! 右辺
                    E(i, j) = BYM*phi(i,j-1) + BXM*phi(i-1,j) - B0*phi(i, j)  + BXP*phi(i+1,j) + BYP*phi(i,j+1) - RHS(i, j) ! E : ∇φ-(RHS)
                    phi(i, j) = phi(i, j) + omg * E(i, j) / B0
                    er = er + (E(i, j)) ** 2
                enddo
            enddo
            ! 固体境界条件
            phi(:, 0) = phi(:, 1)
            phi(:, Ny+1) = phi(:, Ny)

            ! 周期境界条件
            phi(Nx+1, :) = phi(1, :)
            phi(0, :) = phi(Nx, :)
            
            if (er < er0) then
                !write(*, *) itr
                exit
            endif
            
        enddo

    end subroutine SOR

    ! 速度補正
    subroutine speed_cor(u, u_pre, v, v_pre, phi) 
        real(8), intent(in) :: u_pre(0:Nx+1, 0:Ny+1), v_pre(0:Nx+1, 0:Ny+1)
        real(8), intent(in) :: phi(0:Nx+1, 0:Ny+1)
        real(8), intent(inout) :: u(0:Nx+1, 0:Ny+1), v(0:Nx+1, 0:Ny+1)
        integer i, j

        do j = 1, Ny
            do i = 1, Nx
                u(i, j) = u_pre(i, j) - dt*(-phi(i, j) + phi(i+1, j))/dx
            enddo
        enddo

        do j = 1, Ny
            do i = 1, Nx
                v(i, j) = v_pre(i, j) - dt*(-phi(i, j) + phi(i, j+1))/dy
            enddo
        enddo

    end subroutine speed_cor


    ! 圧力更新
    subroutine press_upd(p, phi)
        real(8), intent(in) :: phi(0:Nx+1, 0:Ny+1)
        real(8), intent(inout) :: p(0:Nx+1, 0:Ny+1)
        integer i, j

        do j = 0, Ny+1
            do i = 1, Nx
                p(i, j) = p(i, j) + phi(i, j)
            enddo
        enddo

    end subroutine press_upd


end module subprogs

program main
    use global
    use subprogs
    implicit none
    
    real(8) :: t = 0.0d0
    real(8), allocatable :: x(:, :), y(:, :)
    real(8), allocatable :: x_u(:, :), y_u(:, :) ! u用の座標（x方向に半格子ズレ）
    real(8), allocatable :: x_v(:, :), y_v(:, :) ! v用の座標（y方向に半格子ズレ）
    real(8), allocatable :: u(:, :),  u_pre(:, :), v(:, :), v_pre(:, :)
    real(8), allocatable :: Ax(:, :), Ax_exact(:, :) ! 対流項
    real(8), allocatable :: Ay(:, :), Ay_exact(:, :) 
    real(8), allocatable :: Bx(:, :), Bx_exact(:, :) ! 粘性項
    real(8), allocatable :: By(:, :), By_exact(:, :)
    real(8), allocatable :: p(:, :), dpdx(:, :), dpdy(:, :), phi(:, :),  phi_new(:, :), psi
    real(8) :: BXP, BXM, BYP, BYM, B0 ! 圧力ポアソン方程式　定数項
    integer :: i, j, k
    character(len=32) :: fname ! ファイル名を格納する文字列変数
    integer :: iu, io

    ! 配列割り付け
    allocate(x(0:Nx+1, 0:Ny+1))
    allocate(x_u(0:Nx+1, 0:Ny+1))
    allocate(x_v(0:Nx+1, 0:Ny+1))
    allocate(y(0:Nx+1, 0:Ny+1))
    allocate(y_u(0:Nx+1, 0:Ny+1))
    allocate(y_v(0:Nx+1, 0:Ny+1))

    allocate(u(0:Nx+1, 0:Ny+1))
    allocate(u_pre(0:Nx+1, 0:Ny+1))
    allocate(v(0:Nx+1, 0:Ny+1))
    allocate(v_pre(0:Nx+1, 0:Ny+1))

    allocate(Ax(1:Nx, 1:Ny))
    allocate(Ax_exact(1:Nx, 1:Ny))
    allocate(Ay(1:Nx, 1:Ny-1))
    allocate(Ay_exact(1:Nx, 1:Ny-1))

    allocate(Bx(1:Nx, 1:Ny))
    allocate(Bx_exact(1:Nx, 1:Ny))
    allocate(By(1:Nx, 1:Ny-1))
    allocate(By_exact(1:Nx, 1:Ny-1))

    allocate(p(0:Nx+1, 0:Ny+1))
    allocate(dpdx(0:Nx+1, 0:Ny+1))
    allocate(dpdy(0:Nx+1, 0:Ny+1))
    allocate(phi(0:Nx+1, 0:Ny+1))
    allocate(phi_new(0:Nx+1, 0:Ny+1))

    ! 初期化
    x(:, :) = 0.0d0
    x_u(:, :) = 0.0d0
    x_v(:, :) = 0.0d0
    y(:, :) = 0.0d0
    y_u(:, :) = 0.0d0
    y_v(:, :) = 0.0d0

    u(:, :) = 0.0d0
    u_pre(:, :) = 0.0d0
    v(:, :) = 0.0d0
    v_pre(:, :) = 0.0d0

    Ax(:, :) = 0.0d0
    Ax_exact(:, :) = 0.0d0
    Ay(:, :) = 0.0d0
    Ay_exact(:, :) = 0.0d0

    Bx(:, :) = 0.0d0
    Bx_exact(:, :) = 0.0d0
    By(:, :) = 0.0d0
    By_exact(:, :) = 0.0d0

    p(:, :) = 0.0d0
    dpdx(:, :) = 0.0d0
    dpdy(:, :) = 0.0d0

    phi(:, :) = 0.0d0

    ! 圧力勾配を与える
    p(0, :) = 10.0d0
    do i = 1, Nx
        p(i, :) = 3.0d0
    enddo
    p(nx+1, :) = 0.0d0

    
    ! 格子生成（計算する領域は i:1~Nx, j:1~Ny)
    do j = 0, Ny+1
        do i = 0, Nx+1
            x(i, j) = real(i, kind=8) * dx
            y(i, j) = real(j, kind=8) * dy
            ! x_u(i, j) = real(i+0.5, kind=8) * dx ! u方向の速度計算で用いる座標
            ! y_u(i, j) = real(j, kind=8) * dy
            ! x_v(i, j) = real(i, kind=8) * dx     ! v方向の速度計算で用いる座標
            ! y_v(i, j) = real(j+0.5, kind=8) * dy
        enddo
    enddo

    ! 圧力ポアソン方程式　定数項
    BXP = 1.0d0 / (dx ** 2)
    BXM = 1.0d0 / (dx ** 2)
    BYP = 1.0d0 / (dy ** 2)
    BYM = 1.0d0 / (dy ** 2)
    B0 = 2.0d0 / (dx ** 2) + 2.0d0 / (dy ** 2)

    ! SMAC法　計算
    do k = 1, steps

        call prespeed_cal(Ax, Ay, Bx, By, u, u_pre, v, v_pre, dpdx, dpdy, p) ! 予測速度計算

        ! 周期境界条件
        ! u_pre(0, :) = 1.0d0
        ! u_pre(Nx+1, :) = u_pre(Nx, :)
        u_pre(0, :) = u_pre(Nx, :)
        u_pre(Nx+1, :) = u_pre(1, :)

        ! 滑りなし条件
        u_pre(:, 0) = -u_pre(:, 1)
        v_pre(:, 0) = 0.0d0

        u_pre(:, Ny+1) = -u_pre(:, Ny)
        v_pre(:, Ny) = 0.0d0
        ! write(*,*) 'a'
        call SOR(BXP, BXM, BYP, BYM, B0, u_pre, v_pre, phi) ! 圧力ポアソン方程式を解く
        
        call speed_cor(u, u_pre, v, v_pre, phi) ! 速度補正



        call press_upd(p, phi) ! 圧力更新

        ! 滑りなし条件
        u(:, 0) = -u(:, 1)
        v(:, 0) = 0.0d0

        u(:, Ny+1) = -u(:, Ny)
        v(:, Ny) = 0.0d0

        ! 周期境界条件
        ! u(0, :) = 1.0d0
        ! u(Nx+1, :) = u(Nx, :)
        u(Nx+1, :) = u(1, :)
        v(Nx+1, :) = v(1, :)
        u(0, :) = u(Nx, :)
        v(0, :) = v(Nx, :)

        t = t + dt ! 時間ステップ更新

        ! 10ステップごとに出力
        if (mod(k, 10) == 0) then
            !write(*,*) 'a'
            ! u出力
            write(fname, '(A, I5.5, A)') 'uv_', k, '.d'
            open(newunit=iu, file='/home/shishida/SMAC/.d/'//trim(fname), status='replace', action='write', iostat=io)
            do i = 1, Nx
                do j = 1, Ny
                    write(iu, *) i, j, u(i, j), v(i,j)
                enddo
            enddo
            close(iu)

        endif


    enddo





    end program main

