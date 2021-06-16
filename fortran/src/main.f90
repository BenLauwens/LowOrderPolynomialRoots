program main
    use, intrinsic :: iso_fortran_env, only: REAL64
    use fcsprp
    implicit none

    real(kind=REAL64), dimension(:), allocatable :: poly
    real(kind=REAL64) :: res
    real(kind=REAL64), dimension(:), allocatable :: times, bench
    real(kind=REAL64) :: start, finish
    integer :: i, j, k, l
    integer, parameter :: N = 10000
    integer, parameter :: S = 100
    real(kind=REAL64), dimension(:), allocatable :: out

    allocate(times(N))

    allocate(bench(S))
    do l = 3, 9
        allocate(poly(l))
        allocate(out(l-1))
        do k = 1, S
            call random_number(poly)
            poly = 2.0D0 * poly - 1.0
            do i = 1, N
                call cpu_time(start)
                do j = 1, 1000
                    call allrealroots(poly, out)
                end do
                call cpu_time(finish)
                times(i) = 1000000 * (finish - start)
            end do
            bench(k) = minval(times)
        end do
        deallocate(out)
        deallocate(poly)
        print *, l, minval(bench), sum(bench) / S, maxval(bench), sqrt((sum(bench**2)-sum(bench)**2/S)/(S-1))
    end do

    do l = 3, 9
        allocate(poly(l))
        do k = 1, S
2           call random_number(poly)
            poly = 2.0D0 * poly - 1.0
            do i = 1, N
                call cpu_time(start)
                do j = 1, 1000
                    call smallestpositiveroot(poly, res)
                end do
                call cpu_time(finish)
                times(i) = 1000000 * (finish - start)
            end do
            bench(k) = minval(times)
            if (bench(k) < epsilon(1.0D0)) then
                go to 2
            end if
        end do
        deallocate(poly)
        print *, l, minval(bench), sum(bench) / S, maxval(bench), sqrt((sum(bench**2)-sum(bench)**2/S)/(S-1))
    end do

end program