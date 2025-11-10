module global
    implicit none
    ! 時間
    real(8), parameter ::dt = 0.01d0
    integer, parameter :: steps = 10000

    ! 格子数
    !integer, parameter :: Nx = 100
    !integer, parameter :: Ny = 100

    ! 計算領域サイズ
    real(8), parameter :: Lx = 2.0d0 * acos(-1.0d0)
    real(8), parameter :: Ly = 2.0d0 * acos(-1.0d0)
    real(8) :: dx 
    real(8) :: dy 
    ! 粘性係数μ
    real(8), parameter :: mu = 1.0d0

end module global

module subprogs
    use global
    implicit none
    contains

    subroutine convection_cal(u, v, Ax, Ay, Nx, Ny) ! 対流項
        integer, intent(in) :: Nx, Ny
        real(8), intent(in) :: u(0:Nx, 0:Ny+1), v(0:Nx+1, 0:Ny)
        real(8), intent(out) :: Ax(1:Nx-1, 1:Ny), Ay(1:Nx, 1:Ny-1)
        real(8) :: duudx, duvdy, dvvdy, duvdx
        real(8) :: re_dx, re_dy
        integer :: i, j
        re_dx = 1.0d0 / dx
        re_dy = 1.0d0 / dy

        ! Axを計算
        do j = 1, Ny
            do i = 1, Nx-1
                duudx = (- ((u(i-1, j)+u(i, j)) * 0.5d0)**2 + ((u(i, j)+u(i+1, j)) * 0.5d0)**2) * re_dx
                duvdy = (- ((v(i, j-1)+v(i+1, j-1)) * 0.5d0)*((u(i, j-1) + u(i, j)) * 0.5d0) + ((v(i, j) + v(i+1, j)) * 0.5d0)*((u(i, j) + u(i, j+1)) * 0.5d0)) * re_dy
                Ax(i, j) = duudx + duvdy
            enddo
        enddo

        ! Ayを計算
        do j = 1, Ny-1
            do i = 1, Nx
                dvvdy = (- ((v(i, j-1)+v(i, j)) * 0.5d0)**2 + ((v(i, j)+v(i, j+1)) * 0.5d0)**2) * re_dy
                duvdx = (- ((u(i-1, j)+u(i-1, j+1)) * 0.5d0)*((v(i-1, j)+v(i, j)) * 0.5d0) + ((u(i, j)+u(i, j+1)) * 0.5d0)*((v(i, j)+v(i+1, j)) * 0.5d0)) * re_dx
                Ay(i, j) = duvdx + dvvdy
            enddo
        enddo



    end subroutine convection_cal
    
end module subprogs

program main
    use global
    use subprogs
    implicit none
    real(8) :: t, error_x, error_y
    real(8), allocatable :: x(:, :), y(:, :)
    real(8), allocatable :: x_u(:, :), y_u(:, :) ! u用の座標（x方向に半格子ズレ）
    real(8), allocatable :: x_v(:, :), y_v(:, :) ! v用の座標（y方向に半格子ズレ）
    real(8), allocatable :: u(:, :), v(:, :)
    real(8), allocatable :: Ax(:, :), Ax_exact(:, :)
    real(8), allocatable :: Ay(:, :), Ay_exact(:, :)
    integer :: exp, Nx, Ny, N, i, j

    open(50, file = '/home/shishida/SMAC/.d/err_convection_v2.d') ! 出力ファイルを開く

    do exp = 3, 7 ! Nx, Ny = 8 ~ 128
        Nx = 2 ** exp
        Ny = 2 ** exp
        N = Nx * Ny

        if (allocated(x)) deallocate(x)
        if (allocated(x_u)) deallocate(x_u)
        if (allocated(x_v)) deallocate(x_v)
        if (allocated(y)) deallocate(y)
        if (allocated(y_u)) deallocate(y_u)
        if (allocated(y_v)) deallocate(y_v)
        if (allocated(u)) deallocate(u)
        if (allocated(v)) deallocate(v)
        if (allocated(Ax)) deallocate(Ax)
        if (allocated(Ax_exact)) deallocate(Ax_exact)
        if (allocated(Ay)) deallocate(Ay)
        if (allocated(Ay_exact)) deallocate(Ay_exact)

        ! 配列割り付け
        allocate(x(0:Nx+1, 0:Ny+1))
        allocate(x_u(0:Nx+1, 0:Ny+1))
        allocate(x_v(0:Nx+1, 0:Ny+1))
        allocate(y(0:Nx+1, 0:Ny+1))
        allocate(y_u(0:Nx+1, 0:Ny+1))
        allocate(y_v(0:Nx+1, 0:Ny+1))
        allocate(u(0:Nx, 0:Ny+1))
        allocate(v(0:Nx+1, 0:Ny))
        allocate(Ax(1:Nx-1, 1:Ny))
        allocate(Ax_exact(1:Nx-1, 1:Ny))
        allocate(Ay(1:Nx, 1:Ny-1))
        allocate(Ay_exact(1:Nx, 1:Ny-1))

        ! 初期化
        x(:, :) = 0.0d0
        x_u(:, :) = 0.0d0
        x_v(:, :) = 0.0d0
        y(:, :) = 0.0d0
        y_u(:, :) = 0.0d0
        y_v(:, :) = 0.0d0
        u(:, :) = 0.0d0
        v(:, :) = 0.0d0
        Ax(:, :) = 0.0d0
        Ax_exact(:, :) = 0.0d0
        Ay(:, :) = 0.0d0
        Ay_exact(:, :) = 0.0d0

        dx = Lx / real(Nx-1,kind=8)
        dy = Ly / real(Ny-1,kind=8)


        ! 格子生成（計算する領域は i:1~Nx, j:1~Ny)
        do j = 0, Ny+1
            do i = 0, Nx+1
                x(i, j) = real(i, kind=8) * dx
                y(i, j) = real(j, kind=8) * dy
                x_u(i, j) = real(i+0.5, kind=8) * dx ! u方向の速度計算で用いる座標
                y_u(i, j) = real(j, kind=8) * dy
                x_v(i, j) = real(i, kind=8) * dx     ! v方向の速度計算で用いる座標
                y_v(i, j) = real(j+0.5, kind=8) * dy
            enddo
        enddo

        ! 与える速度(Taylor-green)
        do j = 0, Ny+1
            do i = 0, Nx
                u(i, j) =   sin(x_u(i, j)) * cos(y_u(i, j))
            enddo
        enddo

        do j = 0, Ny
            do i = 0, Nx+1
                v(i, j) = - cos(x_v(i, j)) * sin(y_v(i, j))
            enddo
        enddo

        ! 数値計算
        call convection_cal(u, v, Ax, Ay, Nx, Ny)

        ! 理論解
        do j = 1, Ny
            do i = 1, Nx-1
                Ax_exact(i, j) = 0.5d0 * sin(2.0d0 * x_u(i, j))
            enddo
        enddo

        do j = 1, Ny-1
            do i = 1, Nx
                Ay_exact(i, j) = 0.5d0 * sin(2.0d0 * y_v(i, j))
            enddo
        enddo

        ! RSMEを計算
        error_x = 0.0d0
        error_y = 0.0d0
        do j = 1, Ny
            do i = 1, Nx-1
                error_x = error_x + (Ax(i, j)-Ax_exact(i, j))**2
                !write(*, *) Bx(i, j), Bx_exact(i, j)
            enddo
        enddo

        do j = 1, Ny-1
            do i = 1, Nx
                error_y = error_y + (Ay(i, j)-Ay_exact(i, j))**2
            enddo
        enddo
        error_x = sqrt(error_x / real((Nx-1)*Ny,kind=8))
        error_y = sqrt(error_y / real(Nx*(Ny-1),kind=8))

        write(50, *)  Nx, error_x, error_y


    enddo


end program main
