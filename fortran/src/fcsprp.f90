module fcsprp
    use, intrinsic :: iso_fortran_env, only: REAL64
    implicit none

    private 
    public allrealrootsquadratic, smallestpositiverootquadratic, allrealroots, smallestpositiveroot

contains

    subroutine smallestpositiverootquadratic (coeffs, res)
        implicit none
        real(kind=REAL64), dimension(3), intent(in) :: coeffs
        real(kind=REAL64), intent(out) :: res
        real(kind=REAL64) :: a1, b, c, delta

        a1 = 1.0_REAL64 / coeffs(1)
        b= -0.5_REAL64 * coeffs(2) * a1
        c = coeffs(3) * a1
        delta = b * b - c
        res = huge(0.0_REAL64)
        if (delta < -4 * epsilon(0.0_REAL64)) then
            
        elseif (delta > 4 * epsilon(0.0_REAL64)) then
            if (b > epsilon(0.0_REAL64)) then
                if (c > epsilon(0.0_REAL64)) then
                    res = c / (b + sqrt(delta))
                else
                    res = b + sqrt(delta)
                end if
            elseif (b < -epsilon(0.0_REAL64)) then
                if (c < epsilon(0.0_REAL64)) then
                    res = c / (b - sqrt(delta))
                end if
            else
                res = sqrt(-c)
            end if
        else
            if (b > -epsilon(0.0_REAL64)) then
                res = b
            end if
        end if
    end subroutine

    subroutine allrealrootsquadratic (coeffs, res)
        implicit none
        real(kind=REAL64), dimension(3), intent(in) :: coeffs
        real(kind=REAL64), dimension(2), intent(out) :: res
        real(kind=REAL64) :: a1, b, c, delta

        a1 = 1.0_REAL64 / coeffs(1)
        b = -0.5_REAL64 * coeffs(2) * a1
        c = coeffs(3) * a1
        delta = b * b - c
        res = huge(0.0_REAL64)
        if (delta < -4 * epsilon(0.0_REAL64)) then
            res = (/huge(0.0_REAL64), huge(0.0_REAL64)/)
        elseif (delta > 4 * epsilon(0.0_REAL64)) then
            if (b > epsilon(0.0_REAL64)) then
                res = (/c / (b + sqrt(delta)), b + sqrt(delta)/)
            elseif (b < -epsilon(0.0_REAL64)) then
                res = (/c / (b - sqrt(delta)), b - sqrt(delta)/)
            else
                res = sqrt(-c)
            end if
        else
            if (b > -epsilon(0.0_REAL64)) then
                res = (/b, b/)
            end if
        end if
    end subroutine

    subroutine horner(coeffs, x, v)
        implicit none
        real(kind=REAL64), dimension(:), intent(in) :: coeffs
        real(kind=REAL64), intent(in) :: x
        real(kind=REAL64), intent(out) :: v
        integer :: i

        v = coeffs(1)
        do i = 2, size(coeffs)
            v = v * x + coeffs(i)
        end do
    end subroutine

    subroutine horner2(coeffs, x, y, v, w)
        implicit none
        real(kind=REAL64), dimension(:), intent(in) :: coeffs
        real(kind=REAL64), intent(in) :: x, y
        real(kind=REAL64), intent(out) :: v, w
        integer :: i

        v = coeffs(1)
        w = coeffs(1)
        do i = 2, size(coeffs)
            v = v * x + coeffs(i)
            w = w * y + coeffs(i)
        end do
    end subroutine

    subroutine hornerd(coeffs, x, v, d)
        implicit none
        real(kind=REAL64), dimension(:), intent(in) :: coeffs
        real(kind=REAL64), intent(in) :: x
        real(kind=REAL64), intent(out) :: v, d
        integer :: i

        v = coeffs(1)
        d = 0.0_REAL64
        do i = 2, size(coeffs)
            d = d * x + v
            v = v * x + coeffs(i)
        end do
    end subroutine

    subroutine posintervalhorner(coeffs, low, high, colow, cohigh)
        implicit none
        real(kind=REAL64), dimension(:), intent(in) :: coeffs
        real(kind=REAL64), intent(in) :: low
        real(kind=REAL64), intent(in) :: high
        real(kind=REAL64), intent(out) :: colow
        real(kind=REAL64), intent(out) :: cohigh
        integer :: i

        colow = coeffs(1)
        cohigh = coeffs(1)
        do i = 2, size(coeffs)
            if (colow > 0.0_REAL64) then
                colow = colow * low + coeffs(i)
                cohigh = cohigh * high + coeffs(i)
            elseif (cohigh < 0.0_REAL64) then
                colow = colow * high + coeffs(i)
                cohigh = cohigh * low + coeffs(i)
            else
                colow = colow * high + coeffs(i)
                cohigh = cohigh * high + coeffs(i)
            end if
        end do
    end subroutine

    subroutine intervalhorner(coeffs, low, high, colow, cohigh)
        implicit none
        real(kind=REAL64), dimension(:), intent(in) :: coeffs
        real(kind=REAL64), intent(in) :: low
        real(kind=REAL64), intent(in) :: high
        real(kind=REAL64), intent(out) :: colow
        real(kind=REAL64), intent(out) :: cohigh
        real(kind=REAL64) :: tmp
        integer :: i

        colow = coeffs(1)
        cohigh = coeffs(1)
        do i = 2, size(coeffs)
            if (colow > 0.0_REAL64) then
                if (low > 0.0_REAL64) then
                    tmp = colow * low + coeffs(i)
                    cohigh = cohigh * high + coeffs(i)
                elseif (high < 0.0_REAL64) then
                    tmp = cohigh * low + coeffs(i)
                    cohigh = colow * high + coeffs(i)
                else
                    tmp = cohigh * low + coeffs(i)
                    cohigh = cohigh * high + coeffs(i)
                end if
            elseif (cohigh < 0.0_REAL64) then
                if (low > 0.0_REAL64) then
                    tmp = colow * high + coeffs(i)
                    cohigh = cohigh * low + coeffs(i)
                elseif (high < 0.0_REAL64) then
                    tmp = cohigh * high + coeffs(i)
                    cohigh = colow * low + coeffs(i)
                else
                    tmp = colow * high + coeffs(i)
                    cohigh = colow * low + coeffs(i)
                end if
            else
                if (low > 0.0_REAL64) then
                    tmp = colow * high + coeffs(i)
                    cohigh = cohigh * high + coeffs(i)
                elseif (high < 0.0_REAL64) then
                    tmp = cohigh * low + coeffs(i)
                    cohigh = colow * low + coeffs(i)
                else
                    tmp = min(cohigh * low + coeffs(i), colow * high + coeffs(i))
                    cohigh = max(colow * low + coeffs(i), cohigh * high + coeffs(i))
                end if
            end if
            colow = tmp
        end do
    end subroutine

    subroutine smallestpositiveroot(coeffs, res)
        implicit none
        real(kind=REAL64), dimension(:), intent(in) :: coeffs
        real(kind=REAL64), intent(out) :: res
        real(REAL64), parameter  :: MAX64 = huge(0.0_REAL64), MIN64 = -huge(0.0_REAL64)
        real(REAL64) :: coeffs1, low, high, mid, colow, comid, cohigh, dcolow, dcomid, dcohigh
        real(REAL64) :: leftlow, lefthigh, rightlow, righthigh, delta
        real(kind=REAL64), dimension(size(coeffs)) :: poly
        real(kind=REAL64), dimension(size(coeffs)-1) :: dpoly
        real(kind=REAL64), dimension(2, 2*size(coeffs)+1) :: list
        integer :: i, N, index

        N = size(coeffs)
        res = MAX64
        if (N == 1) then
            return
        elseif (N == 2) then
            res = - coeffs(2) / coeffs(1)
            if (res < 0.0) then
                res = MAX64
            end if
            return
        end if
        coeffs1 = 1.0_REAL64 / coeffs(1)
        poly(1) = 1.0_REAL64
        high = 1.0_REAL64
        do i = 2, N
            dpoly(i-1) = (N+1-i) * poly(i-1)
            poly(i) = coeffs(i) * coeffs1
            if (poly(i) < high) then
                high = poly(i)
            end if
        end do
        if (high > 0.0_REAL64) then
            return
        end if
        low = 0.0_REAL64
        high = 1.0_REAL64 - high
        index = 0
        do
            mid = 0.5_REAL64 * (low + high)
            call horner(poly, mid, comid)
            call posintervalhorner(dpoly, low, high, dcolow, dcohigh)
            if (dcolow < 0.0_REAL64 .and. dcohigh > 0.0_REAL64) then
                if (comid > 0.0_REAL64) then
                    leftlow = low
                    lefthigh = mid - comid / dcohigh
                    rightlow = mid - comid / dcolow
                    righthigh = high
                elseif (comid < 0.0_REAL64) then
                    leftlow = low
                    lefthigh = mid - comid / dcolow
                    rightlow = mid - comid / dcohigh
                    righthigh = high
                else
                    res = mid
                    return
                end if
                if (leftlow < lefthigh) then
                    if (rightlow < righthigh) then
                        index = index + 1
                        list(1, index) = rightlow
                        list(2, index) = righthigh
                    end if
                    low = leftlow
                    high = lefthigh
                    cycle
                elseif (rightlow < righthigh) then
                    low = rightlow
                    high = righthigh
                    cycle
                end if
            else
                call horner2(poly, low, high, colow, cohigh)
                if (colow * cohigh < 0.0_REAL64) then
                    do
                        call hornerd(poly, mid, comid, dcomid)
                        delta = comid / dcomid
                        res = mid - delta
                        if (abs(delta) < 1.0D-8*mid) then
                            return
                        elseif (low < res .and. res < high) then
                            mid = res
                        else
                            if (comid * colow > 0.0_REAL64) then
                                low = mid
                                colow = comid
                            elseif (comid * cohigh > 0.0_REAL64) then
                                high = mid
                                cohigh = comid
                            else
                                res = mid
                                return
                            end if
                            mid = (low*cohigh - high*colow) / (cohigh - colow)
                        end if
                    end do
                end if
            end if
            if (index == 0) then
                return
            else
                low = list(1, index)
                high = list(2, index)
                index = index - 1
            end if
        end do
    end subroutine

    subroutine allrealroots(coeffs, res)
        implicit none
        real(kind=REAL64), dimension(:), intent(in) :: coeffs
        real(kind=REAL64), dimension(size(coeffs)-1), intent(out) :: res
        real(REAL64), parameter  :: MAX64 = huge(0.0_REAL64), MIN64 = - huge(0.0_REAL64)
        real(REAL64) :: coeffs1, s, low, high, mid, colow, comid, cohigh, dcolow, dcomid, dcohigh
        real(REAL64) :: leftlow, lefthigh, rightlow, righthigh, delta, newmid
        real(kind=REAL64), dimension(size(coeffs)) :: poly
        real(kind=REAL64), dimension(size(coeffs)-1) :: dpoly
        real(kind=REAL64), dimension(2, 2*size(coeffs)+1) :: list
        integer :: i, N, index, counter

        N = size(coeffs)
        if (N == 1) then
            return
        elseif (N == 2) then
            res(1) = - coeffs(2) / coeffs(1)
            return
        end if
        coeffs1 = 1.0_REAL64 / coeffs(1)
        poly(1) = 1.0_REAL64
        low = 0.0_REAL64
        high = 0.0_REAL64
        s = -1.0_REAL64
        do i = 2, N
            dpoly(i-1) = (N+1-i) * poly(i-1)
            poly(i) = coeffs(i) * coeffs1
            if (poly(i) < high) then
                high = poly(i)
            end if
            if (s * poly(i) < low) then
                low = s * poly(i)
            end if
            s = -s
        end do
        index = 0
        if (low < 0.0_REAL64) then
            if (high < 0.0_REAL64) then
                index = 1
                list(1, index) = 0.0_REAL64
                list(2, index) = 1.0_REAL64 - high
            end if
            low = low - 1.0_REAL64
            high = 0.0_REAL64
        else
            low = 0.0_REAL64
            high = 1.0_REAL64 - high
        end if
        counter = 0
        do
            !print *, low, high
            mid = 0.5_REAL64 * (low + high)
            call horner(poly, mid, comid)
            call intervalhorner(dpoly, low, high, dcolow, dcohigh)
            !print *, mid, comid, dcolow, dcohigh
            if (dcolow < 0.0_REAL64 .and. dcohigh > 0.0_REAL64) then
                if (comid > 0.0_REAL64) then
                    leftlow = low
                    lefthigh = mid - comid / dcohigh
                    rightlow = mid - comid / dcolow
                    righthigh = high
                elseif (comid < 0.0_REAL64) then
                    leftlow = low
                    lefthigh = mid - comid / dcolow
                    rightlow = mid - comid / dcohigh
                    righthigh = high
                else
                    counter = counter + 1
                    res(counter) = mid
                    leftlow = 0.0_REAL64
                    lefthigh = 0.0_REAL64
                    rightlow = 0.0_REAL64
                    righthigh = 0.0_REAL64
                end if
                if (leftlow < lefthigh) then
                    if (rightlow < righthigh) then
                        index = index + 1
                        list(1, index) = rightlow
                        list(2, index) = righthigh
                    end if
                    low = leftlow
                    high = lefthigh
                    cycle
                elseif (rightlow < righthigh) then
                    low = rightlow
                    high = righthigh
                    cycle
                end if
            else
                call horner2(poly, low, high, colow, cohigh)
                !print *, colow, cohigh
                if (colow * cohigh < 0.0) then
                    do
                        call hornerd(poly, mid, comid, dcomid)
                        delta = comid / dcomid
                        newmid = mid - delta
                        if (abs(delta) < 1.0D-8 * abs(mid)) then
                            counter = counter + 1
                            res(counter) = newmid
                            exit
                        elseif (low < newmid .and. newmid < high) then
                            mid = newmid
                        else
                            if (comid * colow > 0.0_REAL64) then
                                low = mid
                                colow = comid
                            elseif (comid * cohigh > 0.0_REAL64) then
                                high = mid
                                cohigh = comid
                            else
                                counter = counter + 1
                                res(counter) = mid
                                exit
                            end if
                            mid = (low*cohigh - high*colow) / (cohigh - colow)
                        end if
                    end do
                end if
            end if
            if (index == 0) then
                return
            else
                low = list(1, index)
                high = list(2, index)
                index = index - 1
            end if
        end do
    end subroutine
end module
