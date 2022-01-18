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
    integer, parameter :: S = 1000
    character(len=20) :: ch
    real(kind=REAL64), dimension(:), allocatable :: out

    allocate(times(N))

    allocate(bench(S))
    do l = 3, 9
        write(ch, *) l
        open(150, file='../examples/polynomials'//trim(adjustl(ch))//'.txt', action = 'read')
        allocate(poly(l))
        allocate(out(l-1))
        do k = 1, S
            read(150, *) poly
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
        close(150)
        write(*, *) l, minval(bench), sum(bench) / S, maxval(bench), sqrt((sum(bench**2)-sum(bench)**2/S)/(S-1))
        open(150, file='out.txt', action='write', position='append')
        write(150, *) l, minval(bench), sum(bench) / S, maxval(bench), sqrt((sum(bench**2)-sum(bench)**2/S)/(S-1))
        close(150)
    end do

    do l = 3, 9
        write(ch, *) l
        open(150, file='../examples/polynomials'//trim(adjustl(ch))//'.txt', action = 'read')
        allocate(poly(l))
        do k = 1, S
            read(150, *) poly
            do i = 1, N
                call cpu_time(start)
                do j = 1, 1000
                    call smallestpositiveroot(poly, res)
                end do
                call cpu_time(finish)
                times(i) = 1000000 * (finish - start)
            end do
            bench(k) = minval(times)
        end do
        deallocate(poly)
        close(150)
        write(*, *) l, minval(bench), sum(bench) / S, maxval(bench), sqrt((sum(bench**2)-sum(bench)**2/S)/(S-1))
        open(150, file='out.txt', action='write', position='append')
        write(150, *) l, minval(bench), sum(bench) / S, maxval(bench), sqrt((sum(bench**2)-sum(bench)**2/S)/(S-1))
        close(150)
    end do

end program