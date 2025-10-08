program poiseuille_2d_SMAC
    implicit none
    
    ! 変数宣言
    integer :: Nx, Ny, i, j, k, steps, rep, max_rep
    real(8) :: t, dt
    real(8) :: dx, dy, Lx, Ly
    real(8) :: mu =  1.0d0
    real(8), allocatable :: x(:), y(:),  u(:, :), u_pre(:, :), v(:, :), v_pre(:, :),  A(:, :), B(:, :)
    real(8), allocatable :: p(:, :), dpdx(:, :), phi(:, :),  phi_new(:, :), psi
    real(8) :: BXP, BXM, BYP, BYM, B0
    real(8) :: eps, omega, error
    character(len=20) :: fname

    ! 格子数
    Nx = 101
    Ny = 101

    ! 時間
    t = 0.0d0
    dt = 0.001d0
    steps = 10000

    ! 計算領域サイズ
    Lx = 1.0d0
    Ly = 1.0d0
    dx = Lx / (Nx-1)
    dy = Ly / (Ny-1)

    allocate(x(0:Nx+1), y(0:Ny+1))
    allocate(u(0:Nx+1, 0:Ny+1), u_pre(0:Nx+1, 0:Ny+1))
    allocate(v(0:Nx+1, 0:Ny+1), v_pre(0:Nx+1, 0:Ny+1))
    allocate(A(0:Nx+1, 0:Ny+1), B(0:Nx+1, 0:Ny+1))
    allocate(p(0:Nx+1, 0:Ny+1), dpdx(0:Nx+1, 0:Ny+1), phi(0:Nx+1, 0:Ny+1), phi_new(0:Nx+1, 0:Ny+1))


    ! 格子生成
    do i = 1, Nx
        x(i) = (i-1) * dx
    enddo
    do i = 1, Ny
        y(i) = (i-1) * dy
    enddo

    ! 速度場初期化
    u(:, :) = 0.0d0
    u_pre(:, :) = 0.0d0

    ! 圧力場の初期化
    p(1, :) = 100.0d0
    do i = 2, Nx-1
        p(i, :) = 0.0d0
    enddo
    p(Nx, :) = 0.0d0

    ! 圧力ポテンシャルの初期値
    phi(:, :) = 0.0d0

    ! 境界条件(壁面)
    u(:, Ny+1) = -u(:, Ny)
    p(:, Ny+1) = p(:, Ny)

    u(:, 0) = -u(:, 1)
    p(:, 0) = p(:, 1)

    ! 周期境界条件(出入口)


    ! 圧力ポアソン方程式　定数項
    BXP = 1 / (dx ** 2)
    BXM = 1 / (dx ** 2)
    BYP = 1 / (dy ** 2)
    BYM = 1 / (dy ** 2)
    B0 = 2 * ((dx ** 2) + (dy ** 2)) / ((dx ** 2) * (dy ** 2))

    ! 収束条件
    eps = 10d-6
    max_rep = 500
    omega = 1.7d0 ! SORの緩和係数

    


    do k = 1, steps

        ! 予測速度計算
        do i = 2, Nx-1
            do j = 2, Ny-1
                A(i, j) = 0 ! ポアズイユ流では移流項が0
                B(i, j) = mu * ( (u(i+1,j)-2*u(i,j)+u(i-1,j))/dx**2  &
                & + (u(i,j+1)-2*u(i,j)+u(i,j-1))/dy**2 ) ! 粘性項

                dpdx(i, j) = (-p(i, j) + p(i+1, j)) / dx ! 圧力勾配
                u_pre(i, j) = u(i, j) + dt * (A(i, j) + B(i, j) -dpdx(i, j))
            enddo
        enddo

        ! ポアソン方程式を解く→SOR法？
        ! ガウス・ザイデル法
        do rep = 1, max_rep
            do i = 2, Nx-1
                do j = 2, Ny-1
                    psi = ((u_pre(i+1,j)-u_pre(i-1,j))/(2*dx)) / dt
                    phi_new(i, j) = phi(i, j) + &
                    & ( BYM*phi(i, j-1) + BXM*phi(i-1, j) - B0*phi(i, j) &
                    & + BXP*phi(i+1, j) + BYP*phi(i, j+1) - psi) / B0
                    error = error + abs(phi_new(i, j) - phi(i, j))
                enddo
            enddo
            phi(:, :) = phi_new(:, :)
            if (error < eps) exit
            
        enddo

        ! 速度補正
        do i = 2, Nx-1
            do j = 2, Ny-1
                u(i, j) = u_pre(i, j) - dt/dx * (phi(i+1, j)-phi(i-1, j))/2.0d0
            enddo
        enddo

        ! 圧力更新
        do i = 2, Nx-1
            do j = 2, Ny-1
                p(i, j) = p(i, j) + phi(i, j)
            enddo
        enddo

        
        ! 可視化用出力
        if (mod(k, 100) == 0) then
            write(fname,'("output_",I5.5,".dat")') k
            open(20, file=fname, status="replace")
            do j = 1, Ny
                do i = 1, Nx
                    write(20,*) x(i), y(j), u(i,j)
                enddo
                write(20,*)
            enddo
            close(20)
        endif

    enddo


end program poiseuille_2d_SMAC