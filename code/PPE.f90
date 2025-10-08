module global
    implicit none
    ! 時間
    real(8), parameter ::dt = 0.01d0
    integer, parameter :: steps = 10000

    ! 格子数
    integer, parameter :: Nx = 101
    integer, parameter :: Ny = 101

    ! 計算領域サイズ
    real(8), parameter :: Lx = 1.0d0
    real(8), parameter :: Ly = 1.0d0
    real(8), parameter :: dx = Lx / (Nx-1)
    real(8), parameter :: dy = Ly / (Ny-1)
end module global


module subprogs
    implicit none
    contains

    subroutine gauss_seidel(BXP, BXM, BYP, BYM, B0, u_pre, phi, x)
        use global
        real(8), intent(in) :: BXP, BXM, BYP, BYM, B0, u_pre(Nx+1, Ny+1), x(2, Nx+1, Ny+1)
        real(8), intent(inout) :: phi(Nx+1, Ny+1)
        real(8) :: phi_new(Nx+1, Ny+1), er, er0, RHS, phi_pre
        integer :: i, j, itr, itr_max

        itr_max = 10000
        phi_new(:, :) = 0.0d0
        er = 0.0d0
        er0 = 10e-6

        do itr = 1, itr_max
            er = 0.0d0
            do j = 2, Ny-1
                do i = 2, Nx-1

                    RHS = -2 * (x(1, i, j) * (1-x(1, i, j)) + x(2, i, j) * (1-x(2, i, j)))
                    !RHS = ((u_pre(i+1,j)-u_pre(i-1,j))/(2*dx)) / dt →いったん簡単な右辺を考える(20251002)
                    !phi_new(i, j) = phi(i, j) + &
                    !& ( BYM*phi(i, j-1) + BXM*phi(i-1, j) - B0*phi(i, j) &
                    !& + BXP*phi(i+1, j) + BYP*phi(i, j+1) - RHS) / B0 →長いので一旦ボツ(20251002)
                    phi_pre = phi(i,j)
                    phi(i,j) = (BYM*phi(i,j-1) + BXM*phi(i-1,j) + BXP*phi(i+1,j) + BYP*phi(i,j+1) - RHS ) / B0
                    er = er + abs(phi(i, j) - phi_pre)
                enddo
            enddo
            !phi(:, :) = phi_new(:, :)
            write(*, *) 'itr = ', itr, 'err = ', er
            if (er < er0) then
                write(*, *) "exit itr = ", itr
                exit
            endif
            
        enddo
        end subroutine gauss_seidel

        subroutine SOR(BXP, BXM, BYP, BYM, B0, u_pre, phi, x)
        use global
        real(8), intent(in) :: BXP, BXM, BYP, BYM, B0, u_pre(Nx+1, Ny+1), x(2, Nx+1, Ny+1)
        real(8), intent(inout) :: phi(Nx+1, Ny+1)
        real(8) :: phi_new(Nx+1, Ny+1), er, er0, RHS, omg, phi_pre ! RHS : 右辺
        integer :: i, j, itr, itr_max

        itr_max = 10000
        phi_new(:, :) = 0.0d0
        er = 0.0d0
        er0 = 10e-6
        omg = 1.5d0

        do itr = 1, itr_max
            er = 0.0d0 ! 誤差をリセット
            do j = 2, Ny-1
                do i = 2, Nx-1

                    RHS = - 2 * (x(1, i, j) * (1-x(1, i, j)) + x(2, i, j) * (1-x(2, i, j)))
                    !RHS = ((u_pre(i+1,j)-u_pre(i-1,j))/(2*dx)) / dt →いったん簡単な右辺を考える(20251002)
                    phi_pre = phi(i, j)
                    phi(i, j) = phi(i, j) + &
                    & omg * ( BYM*phi(i,j-1) + BXM*phi(i-1,j) - B0*phi(i, j)  + BXP*phi(i+1,j) + BYP*phi(i,j+1) - RHS ) / B0
                    er = er + abs(phi(i, j) - phi_pre)
                enddo
            enddo
            
            !write(*, *) 'itr = ', itr, 'err = ', er
            if (er < er0) then
                write(*, *) "exit itr = " , itr
                exit
            endif
            
        enddo

    end subroutine SOR

end module subprogs

program main
    use global
    use subprogs
    real(8) :: t
    real(8) :: error, diff
    real(8) :: x(2, Nx+1, Ny+1), BXP, BXM, BYP, BYM, B0
    real(8) :: u_pre(Nx+1, Ny+1)
    real(8) :: phi(Nx+1, Ny+1), phi_exact(Nx+1, Ny+1)
    !integer :: itr, itr_max

    u_pre(:, :) = 0.0d0
    phi(:, :) = 0.0d0

    ! 格子生成
    do j = 1, Ny
        do i = 1, Nx
            x(1, i, j) = (i-1) * dx
            x(2, i, j) = (j-1) * dy
        enddo
    enddo

    ! 圧力ポアソン方程式　定数項
    BXP = 1 / (dx ** 2)
    BXM = 1 / (dx ** 2)
    BYP = 1 / (dy ** 2)
    BYM = 1 / (dy ** 2)
    B0 = 2 * ((dx ** 2) + (dy ** 2)) / ((dx ** 2) * (dy ** 2))

    !call gauss_seidel(BXP, BXM, BYP, BYM, B0, u_pre, phi, x)
    call SOR(BXP, BXM, BYP, BYM, B0, u_pre, phi, x)

    error = 0.0d0
    do j = 1, Ny
        do i = 1, Nx
            phi_exact(i, j) = x(1, i, j) * (1.0d0 - x(1, i, j)) + x(2, i, j) * (1.0d0 - x(2, i, j))
            diff = phi(i, j) - phi_exact(i, j)
            error = error + diff ** 2
        enddo
    enddo
    error = sqrt(error / real(Nx*Ny))
    write(*, *) ' error = ', error

    ! open(50, file = 'phi.d') ! 出力ファイルを開く
    ! do j = 1, Ny
    !     do i = 1, Nx
    !         write(50, '(3e12.4)') x(2, 51, j), phi(51, j) ! (x座標)、(y座標)、(ポテンシャル) の順で表示
    !     enddo
    !     write(50, *) ''
    ! enddo



end program main