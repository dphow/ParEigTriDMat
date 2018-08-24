module SecularSolver
    implicit none
    real (KIND (0.0D0)), parameter      :: OFFSET = 1.0E-6
    real (KIND (0.0D0)), parameter      :: TOL = 1.0E-8

contains
    function SecularEq(bm,N,d,xi,x) result(total)
        implicit none
        real (KIND (0.0D0))                 :: total
        real (KIND (0.0D0)), intent(in)     :: bm, d(:), x, xi(:)
        integer              :: i
        integer, intent(in)  :: N

        total = 0.0d0
        do i = 1, N
            total = total + xi(i)**2/(d(i)-x)
        end do
        total = bm*total + 1.0d0
    end function SecularEq

    subroutine SolveSecularEq(bm,N,d,xi,lambda, my_rank);
        implicit none
        real (KIND (0.0D0))                 :: xl, xr, yl, yr, xm, ym
        real (KIND (0.0D0)), intent(in)     :: bm, d(:), xi(:)
        real (KIND (0.0D0)), intent(out)  :: lambda(:)
        integer              :: i
        integer, intent(in)  :: N, my_rank

        do i = 1, N-1
            xl = d(i) + OFFSET
            yl = SecularEq(bm,N,d,xi,xl)
            xr = d(i+1) - OFFSET
            yr = SecularEq(bm,N,d,xi,xr)
            xm = 0.5d0*(xl + xr)
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
                xm = 0.5d0*(xl+xr)
                ym = SecularEq(bm,N,d,xi,xm)
            end do
            lambda(i) = xm
        end do

        xl = d(N) + OFFSET
        yl = SecularEq(bm,N,d,xi,xl)
        xr = 2*d(N)
        yr = SecularEq(bm,N,d,xi,xr)
        xm = 0.5d0*(xl+xr)
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
                xm = 0.5d0*(xl+xr)
                ym = SecularEq(bm,N,d,xi,xm)
            end do
            lambda(N) = xm
        end if
    end
end module SecularSolver
