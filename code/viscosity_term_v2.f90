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

    subroutine viscosity_cal(u, v, Bx, By, Nx, Ny) ! 粘性項
        integer, intent(in) :: Nx, Ny
        real(8), intent(in) :: u(0:Nx, 0:Ny+1), v(0:Nx+1, 0:Ny)
        real(8), intent(out) :: Bx(1:Nx-1, 1:Ny), By(1:Nx, 1:Ny-1)
        real(8) :: dudxx, dudyy, dvdxx, dvdyy
        real(8) :: re_dx, re_dy
        integer :: i, j
        re_dx = 1.0d0 / dx
        re_dy = 1.0d0 / dy

        ! Bxを計算
        do j = 1, Ny
            do i = 1, Nx-1
                dudxx = (u(i-1, j) - 2*u(i, j) + u(i+1, j)) / dx**2
                dudyy = (u(i, j-1) - 2*u(i, j) + u(i, j+1)) / dy**2
                Bx(i, j) = mu * (dudxx + dudyy)
            enddo
        enddo

        ! Byを計算
        do j = 1, Ny-1
            do i = 1, Nx
                dvdxx = (v(i-1, j) - 2*v(i, j) + v(i+1, j)) / dx**2
                dvdyy = (v(i, j-1) - 2*v(i, j) + v(i, j+1)) / dy**2
                By(i, j) = mu * (dvdxx + dvdyy)
            enddo
        enddo



    end subroutine viscosity_cal
    
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
    real(8), allocatable :: Bx(:, :), Bx_exact(:, :)
    real(8), allocatable :: By(:, :), By_exact(:, :)
    integer :: exp, Nx, Ny, N, i, j

    open(50, file = '/home/shishida/SMAC/.d/err_viscosity_v2.d') ! 出力ファイルを開く

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
        if (allocated(Bx)) deallocate(Bx)
        if (allocated(Bx_exact)) deallocate(Bx_exact)
        if (allocated(By)) deallocate(By)
        if (allocated(By_exact)) deallocate(By_exact)

        ! 配列割り付け
        allocate(x(0:Nx+1, 0:Ny+1))
        allocate(x_u(0:Nx+1, 0:Ny+1))
        allocate(x_v(0:Nx+1, 0:Ny+1))
        allocate(y(0:Nx+1, 0:Ny+1))
        allocate(y_u(0:Nx+1, 0:Ny+1))
        allocate(y_v(0:Nx+1, 0:Ny+1))
        allocate(u(0:Nx, 0:Ny+1))
        allocate(v(0:Nx+1, 0:Ny))
        allocate(Bx(1:Nx-1, 1:Ny))
        allocate(Bx_exact(1:Nx-1, 1:Ny))
        allocate(By(1:Nx, 1:Ny-1))
        allocate(By_exact(1:Nx, 1:Ny-1))

        ! 初期化
        x(:, :) = 0.0d0
        x_u(:, :) = 0.0d0
        x_v(:, :) = 0.0d0
        y(:, :) = 0.0d0
        y_u(:, :) = 0.0d0
        y_v(:, :) = 0.0d0
        u(:, :) = 0.0d0
        v(:, :) = 0.0d0
        Bx(:, :) = 0.0d0
        Bx_exact(:, :) = 0.0d0
        By(:, :) = 0.0d0
        By_exact(:, :) = 0.0d0

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
        call viscosity_cal(u, v, Bx, By, Nx, Ny)

        ! 理論解
        do j = 1, Ny
            do i = 1, Nx-1
                Bx_exact(i, j) = -2.0d0 * mu * sin(x_u(i, j)) * cos(y_u(i, j))
            enddo
        enddo

        do j = 1, Ny-1
            do i = 1, Nx
                By_exact(i, j) =  2.0d0 * mu * cos(x_v(i, j)) * sin(y_v(i, j))
            enddo
        enddo

        ! RSMEを計算
        error_x = 0.0d0
        error_y = 0.0d0
        do j = 1, Ny
            do i = 1, Nx-1
                error_x = error_x + (Bx(i, j)-Bx_exact(i, j))**2
                !write(*, *) Bx(i, j), Bx_exact(i, j)
            enddo
        enddo

        do j = 1, Ny-1
            do i = 1, Nx
                error_y = error_y + (By(i, j)-By_exact(i, j))**2
            enddo
        enddo
        error_x = sqrt(error_x / real((Nx-1)*Ny,kind=8))
        error_y = sqrt(error_y / real(Nx*(Ny-1),kind=8))

        write(50, *)  Nx, error_x, error_y


    enddo


end program main
