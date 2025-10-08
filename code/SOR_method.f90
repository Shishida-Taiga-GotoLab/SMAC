module subprogs
    implicit none
    contains

    subroutine set_dbc(phi, x, n1, n2) ! ディリクレ境界条件
        real(8), intent(inout) :: phi(n1, n2), x(2, n1, n2)
        integer, intent(in) :: n1, n2
        integer i, j

        do i = 1, n1
            phi(i, 1) = sin(acos(-1.0d0) * x(1, i, 1))
            phi(i, n2) = 0.0d0
        enddo

        do j = 1, n2
            phi(1, j) = 0.0d0
            phi(n1, j) = 0.0d0
        enddo

    end subroutine set_dbc

    subroutine chk_err(phi, c, d, n1, n2, er) ! 誤差を求めるサブルーチンだが、SORと独立させる意味がよくわからん
        real(8), intent(in) :: phi(n1, n2), c, d
        real(8), intent(out) :: er
        integer, intent(in) :: n1, n2
        real(8) :: rhs
        integer i, j
        er = 0.0d0

        do i = 1, n1
            do j = 1, n2
                rhs = -c * (phi(i-1, j) + phi(i+1, j)) &
                &     -d * (phi(i, j-1) + phi(i, j+1))
                er = er + (rhs - phi(i, j)) **  2
            enddo
        enddo
    end subroutine chk_err

    subroutine SOR(phi, omg, n1, n2, c, d, er) ! SOR法
        real(8), intent(inout) :: phi(n1, n2)
        real(8), intent(out) :: er
        real(8), intent(in) :: omg, c, d
        integer, intent(in) :: n1, n2
        real(8) :: rhs
        integer i, j
        er = 0.0d0

        do j = 2, n2-1
            do i = 2, n1-1
                rhs = -c * (phi(i-1, j) + phi(i+1, j)) &
                &     -d * (phi(i, j-1) + phi(i, j+1))
                phi(i, j) = phi(i, j) + omg * (rhs - phi(i, j))
                er = er + (rhs - phi(i, j)) **  2
            enddo
        enddo

    end subroutine SOR

end module subprogs




program main
    use subprogs

    real(8), allocatable :: phi(:, :), x(:, :, :)
    integer :: n1, n2, N, itrmax
    integer :: i, j, itr
    real(8) :: dx1, dx2, c, d, omg, er, er0

    n1 = 100
    n2 = 100
    N = n1 * n2
    itrmax = 3000
    dx1 = 0.01
    dx2 = 0.01
    c = - (dx2 ** 2) / (2.0d0 * (dx1 ** 2 + dx2 ** 2)) !c,, d :: 連立方程式の係数
    d = (dx1 / dx2) ** 2 * c
    omg = 1.5d0
    er0 = 1e-6

    allocate(phi(n1, n2), x(2, n1, n2))

    phi(:, :) = 0.0d0
    !x(:) = 0.0d0 ! xはどこで使う？
    ! 座標配列の初期化
    do j = 1, n2
        do i = 1, n1
            x(1, i, j) = (i-1) * dx1
            x(2, i, j) = (j-1) * dx2
        enddo
    enddo

    
    call set_dbc(phi, x, n1, n2)
    do itr = 1, itrmax
        call SOR(phi, omg, n1, n2, c, d, er)
        ! call chk_err(phi, c, d, n1, n2, er)
        !write(*, *) 'itr = ', itr, 'err = ', er ! 途中結果
        if (er < er0) then
            write(*, *) "converged at itr = ", itr
            exit
        endif
        phi(:, n2) = phi(:, n2-1) ! ノイマン境界条件
    enddo

    open(50, file = 'phi.d') ! 出力ファイルを開く
    do j = 1, n2
        do i = 1, n1
            write(50, '(3e12.4)') x(:, i, j), phi(i, j) ! (x座標)、(y座標)、(ポテンシャル) の順で表示
        enddo
        write(50, *) ''
    enddo

end program main