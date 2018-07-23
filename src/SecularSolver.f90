module SecularSolver
    implicit none
    real, parameter      :: OFFSET = 1.0E-5
    real, parameter      :: TOL = 1.0E-6

contains
    function SecularEq(bm,N,d,xi,x) result(total)
        implicit none
        real                 :: total
        real, intent(in)     :: bm, d(:), x, xi(:)
        integer              :: i
        integer, intent(in)  :: N

        total = 0.0
        do i = 1, N
            total = total + xi(i)**2/((d(i)-x))
            total = bm*total + 1.0
        end do
    end function SecularEq

    subroutine SolveSecularEq(bm,N,d,xi,lambda);
        implicit none
        real                 :: xl, xr, yl, yr, xm, ym
        real, intent(in)     :: bm, d(:), xi(:)
        real, intent(inout)  :: lambda(:)
        integer              :: i
        integer, intent(in)  :: N

        do i = 1, N-1
            xl = d(i) + OFFSET
            yl = SecularEq(bm,N,d,xi,xl)
            xr = d(i+1) - OFFSET
            yr = SecularEq(bm,N,d,xi,xr)
            xm = 0.5*(xl+xr)
            ym = SecularEq(bm,N,d,xi,xm)

            if (yl*yr > 0) then
                lambda(i) = xl
                cycle
            end if

            do while (abs(ym) > TOL)
                if (yl*ym < 0) then
                    xr = xm
                else
                    xl = xm
                end if
                xm = 0.5*(xl+xr)
                ym = SecularEq(bm,N,d,xi,xm)
            end do
            lambda(i) = xm
        end do

        xl = d(N) + OFFSET
        yl = SecularEq(bm,N,d,xi,xl)
        xr = 2*d(N)
        yr = SecularEq(bm,N,d,xi,xr)
        xm = 0.5*(xl+xr)
        ym = SecularEq(bm,N,d,xi,xm)

        if (yl*yr > 0) then
            lambda(N) = xl
        else
            do while (abs(ym) > TOL)
                if (yl*ym < 0) then
                    xr = xm
                else
                    xl = xm
                end if
                xm = 0.5*(xl+xr)
                ym = SecularEq(bm,N,d,xi,xm)
            end do
            lambda(N) = xm
        end if
    end
end module SecularSolver
