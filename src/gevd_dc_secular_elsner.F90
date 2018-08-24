!> @brief                   Compute the root of g(\\lambda) = 0, which is smaller than d_1 accurately,
!!                          where g(x) := \\sum_{k = 1}^{n} \\frac{z_k^2}{d_k - x} - \\frac{1}{\\rho - x}.
!! @note                    \\lambda is implicitly expressed as y := d(1) - \\lambda > 0.
!! @note                    This subroutine does not check the argument validity.
!! @todo                    Use error code instead of 'stop'.
!!
!! @param n                 The dimension of the vectors d and z.
!! @param d                 The vector d. dimension(n).
!! @param z                 The vector z. dimension(n)..
!! @param rho               The scalar \\rho.
!! @param y                 The root of the function g.
!! @param s_numerator       Workspace for computing s. dimension(n).
!! @param r_numerator       Workspace for computing r (will not change in the main loop). dimension(n).
!! @param p_numerator       Workspace for computing p. dimension(n).
!! @param invdenominator    Workspace for computing s, r, p. dimension(n).
!! @param term              Workspace for computing s, r, p. dimension(n).
!! @param info              Error code.
subroutine gevd_dc_secular_elsner_leftmost(n, d, z, rho, y, s_numerator, r_numerator, p_numerator, &
                                           invdenominator, term, info)
  implicit none
  integer, intent(in) :: n
  real(kind(0.0d0)), intent(in) :: d(*), z(*), rho
  real(kind(0.0d0)), intent(out) :: y
  real(kind(0.0d0)), intent(out) :: s_numerator(*), r_numerator(*), p_numerator(*)
  real(kind(0.0d0)), intent(out) :: invdenominator(*), term(*)
  integer, intent(out) :: info
  !
  integer :: k, h
  integer, parameter :: k_max = 50
  real(kind(0.0d0)) :: s, r, p 
  real(kind(0.0d0)) :: beta ! An auxiliary variable for computing coef(:)
  real(kind(0.0d0)) :: coef(4) ! Coefficients of the cubic equation.
  real(kind(0.0d0)) :: lower, upper, one_minus_zzt, bound_diff
  integer :: m
  real(kind(0.0d0)) :: root(3), rwork(3 * 3)
  integer :: m_in_interval
  real(kind(0.0d0)) :: y_candidate

  !> Compute the lower and upper bounds of the root and set an initial approximation.
  one_minus_zzt = 1.0d0
  do k = 1, n
    one_minus_zzt = one_minus_zzt - z(k) ** 2.0d0
  enddo
  bound_diff = 0.0d0
  do k = 1, n
    bound_diff = bound_diff + (z(k) ** 2.0d0) * dabs(d(k) - rho)
  enddo
  bound_diff = dble(n) * bound_diff / one_minus_zzt
  lower = bound_diff
  upper = 0.0d0
  ! An initial approximation
  y = 0.5d0 * lower

  !> Iteratively compute the root.
  info = 100 ! Zero will be set to info if converged.
  r_numerator(1 : n) = (z(1 : n) ** 2.0d0) * (d(1 : n) - d(1))
  do k = 1, k_max
    ! Compute s, r, p
    s_numerator(1 : n) = (z(1 : n) ** 2.0d0) * (y ** 3.0d0)
    p_numerator(1 : n) = (z(1 : n) ** 2.0d0) * (d(1 : n) - d(1)) * (d(1 : n) - 2.0d0 * d(1) + 3.0d0 * y)
    invdenominator(1 : n) = ((d(1 : n) - d(1)) + y) ** (- 3.0d0)
    !
    term(1 : n - 1) = s_numerator(2 : n) * invdenominator(2 : n)
    call accurate_sum(n - 1, term, s)
    s = z(1) ** 2.0d0 + s
    !
    term(1 : n - 1) = r_numerator(2 : n) * invdenominator(2 : n)
    call accurate_sum(n - 1, term, r)
    !
    term(1 : n - 1) = p_numerator(2 : n) * invdenominator(2 : n)
    call accurate_sum(n - 1, term, p)
    ! Check
    if (r .le. 0.0d0 .or. s .le. 0.0d0) then
      write (*, *) '[leftmost] Cannot happen. Bad s or r. k, s, r = ', k, s, r
      write (*, *) '[leftmost] Cannot happen. d = ', d(1 : n) 
      write (*, *) '[leftmost] Cannot happen. z = ', z(1 : n) 
      stop
    endif
    ! Compute new coefficients from s, r, p
    beta = rho - d(1)
    coef(1) = s * beta 
    coef(2) = (p + r * d(1)) * beta + s - 1
    coef(3) = (p + r * d(1)) - r * beta
    coef(4) = - r
    ! Solve a cubic equation
    call solve_real_algebraic_equation(3, coef, maxval(dabs(coef(1 : 4))) * 1.0e-18, m, root, rwork)
    ! Select a solution in the target interval as a new approximation.
    if (m .lt. 1 .or. m .gt. 3) then
      write (*, *) 'Cannot happen. Bad m. m = ', m
      stop
    else
      m_in_interval = 0
      y_candidate = - huge(0.0d0) ! To suppress an uninitialized warning.
      do h = 1, m
        if (root(h) .gt. upper) then
          if (m_in_interval .ne. 0) then
            write (*, *) '[Warning] Multiple roots exist in the target interval.'
            if (dabs(y_candidate - y) .gt. dabs(root(h) - y)) then
              ! Select the nearest root from lambda.
              m_in_interval = m_in_interval + 1
              y_candidate = root(h)
            endif
          else
            m_in_interval = 1
            y_candidate = root(h)
          endif
        endif
      enddo
      if (m_in_interval .le. 0) then
        write (*, *) '[Leftmost] Cannot happen. Bad m_in_interval. m, m_in_interval = ', m, m_in_interval
        stop
      endif
      if (dabs(y - y_candidate) <= 0.0d0) then
        info = 0
        exit
      endif
      y = y_candidate
    endif
  enddo

  return
end

!> @brief    Compute the root of g(\\lambda) = 0, which is larger than d_n accurately.
!!           This subroutine is essentially same as gevd_dc_secular_elsner_leftmost(...).
!! @sa       gevd_dc_secular_elsner_leftmost.
subroutine gevd_dc_secular_elsner_rightmost(n, d, z, rho, y,  s_numerator, r_numerator, p_numerator, &
                                            invdenominator, term, info)
  implicit none
  integer, intent(in) :: n
  real(kind(0.0d0)), intent(in) :: d(*), z(*), rho
  real(kind(0.0d0)), intent(out) :: y ! y := d(n) - \\lambda < 0
  real(kind(0.0d0)), intent(out) :: s_numerator(*), r_numerator(*), p_numerator(*)
  real(kind(0.0d0)), intent(out) :: invdenominator(*), term(*)
  integer, intent(out) :: info
  !
  integer :: k, h
  integer, parameter :: k_max = 50
  real(kind(0.0d0)) :: s, r, p 
  real(kind(0.0d0)) :: beta ! Auxiliary variables for computing coef(:)
  real(kind(0.0d0)) :: coef(4) ! Coefficients of the cubic equation.
  real(kind(0.0d0)) :: lower, upper, one_minus_zzt, bound_diff
  integer :: m
  real(kind(0.0d0)) :: root(3), rwork(3 * 3)
  integer :: m_in_interval
  real(kind(0.0d0)) :: y_candidate

  !> Compute the lower and upper bounds of the root and set an initial approximation.
  one_minus_zzt = 1.0d0
  do k = 1, n
    one_minus_zzt = one_minus_zzt - z(k) ** 2.0d0
  enddo
  bound_diff = 0.0d0
  do k = 1, n
    bound_diff = bound_diff + (z(k) ** 2.0d0) * dabs(d(k) - rho)
  enddo
  bound_diff = dble(n) * bound_diff / one_minus_zzt
  lower = 0.0d0
  upper = - bound_diff
  ! An initial approximation
  y = 0.5d0 * upper

  !> Iteratively compute the root.
  info = 100 ! Zero will be set to info if converged.
  r_numerator(1 : n) = (z(1 : n) ** 2.0d0) * (d(1 : n) - d(n))
  do k = 1, k_max
    ! Compute s, r, p
    s_numerator(1 : n) = (z(1 : n) ** 2.0d0) * (y ** 3.0d0)
    p_numerator(1 : n) = (z(1 : n) ** 2.0d0) * (d(1 : n) - d(n)) * (d(1 : n) - 2.0d0 * d(n) + 3.0d0 * y)
    invdenominator(1 : n) = (d(1 : n) - d(n) + y) ** (- 3.0d0)
    term(1 : n - 1) = s_numerator(1 : n - 1) * invdenominator(1 : n - 1)
    call accurate_sum(n - 1, term, s)
    s = z(n) ** 2.0d0 + s 
    !
    term(1 : n - 1) = r_numerator(1 : n - 1) * invdenominator(1 : n - 1)
    call accurate_sum(n - 1, term, r)
    !
    term(1 : n - 1) = p_numerator(1 : n - 1) * invdenominator(1 : n - 1)
    call accurate_sum(n - 1, term, p)
    ! Check
    if (r .le. 0.0d0 .or. s .le. 0.0d0) then
      write (*, *) '[rightmost] Cannot happen. Bad s or r. s, r = ', s, r
      stop
    endif
    ! Compute new coefficients from s, r, p
    beta = rho - d(n)
    coef(1) = s * beta
    coef(2) = (p + r * d(n)) * beta + s - 1.0d0 
    coef(3) = (p + r * d(n)) - r * beta
    coef(4) = - r
    ! Solve a cubic equation
    call solve_real_algebraic_equation(3, coef, maxval(dabs(coef(1 : 4))) * 1.0e-18, m, root, rwork)
    ! Select a solution in the target interval as a new approximation.
    if (m .lt. 1 .or. m .gt. 3) then
      write (*, *) 'Cannot happen. Bad m. m = ', m
      stop
    else
      m_in_interval = 0
      y_candidate = huge(0.0d0) ! To suppress an uninitialized warning.
      do h = 1, m
        if (root(h) .lt. lower) then
          if (m_in_interval .ne. 0) then
            write (*, *) '[Warning] Multiple roots exist in the target interval.'
            if (dabs(y_candidate - y) .gt. dabs(root(h) - y)) then
              ! Select the nearest root from lambda.
              m_in_interval = m_in_interval + 1
              y_candidate = root(h)
            endif
          else
            m_in_interval = 1
            y_candidate = root(h)
          endif
        endif
      enddo
      if (m_in_interval .le. 0) then
        write (*, *) '[Leftmost] Cannot happen. Bad m_in_interval. m, m_in_interval = ', m, m_in_interval
        stop
      endif
      if (dabs(y - y_candidate) <= 0.0d0) then
        info = 0
        exit
      endif
      y = y_candidate
    endif
  enddo

  return
end

!> Compute the root of g(\\lambda) = 0 in (d_{j-1}, d_j),
!! whose nearest poll is d_j, accurately
!! where g(x) := \\sum_{k = 1}^{n} \\frac{z_k^2}{d_k - x} - \\frac{1}{\\rho - x}.
!!
!! \\lambda is implicitly expressed as y := d(j) - \\lambda > 0.
!! This subroutine does not check the argument validity.
!!
!! @todo                      Use error code. Do not use 'stop'.
!!
!! @param n                   The dimension of the vectors d and z.
!! @param j                   The interval number.
!! @param d                   The vector d. dimension(n).
!! @param z                   The vector z. dimension(n).
!! @param rho                 The scalar \\rho.
!! @param y                   The root of the function g.
!! @param s_numerator         Workspace for computing s (will not change in the main loop). dimension(n).
!! @param r_numerator         Workspace for computing r (will not change in the main loop). dimension(n).
!! @param p_numerator         Workspace for computing p (will not change in the main loop). dimension(n).
!! @param invdenominator      Workspace for computing s, r, p. dimension(n).
!! @param term                Workspace for computing s, r, p. dimension(n).
!! @param dmdj                Workspace for d(1 : n) - d(j). dimension(n).
!! @param info                Error code.
subroutine gevd_dc_secular_elsner_middle_right(n, j, d, z, rho, y, s_numerator, r_numerator, p_numerator, &
                                               invdenominator, term, dmdj, info)
  implicit none
  integer, intent(in) :: n, j
  real(kind(0.0d0)), intent(in) :: d(*), z(*), rho
  real(kind(0.0d0)), intent(out) :: y
  real(kind(0.0d0)), intent(out) :: s_numerator(*), r_numerator(*), p_numerator(*)
  real(kind(0.0d0)), intent(out) :: invdenominator(*),  term(*), dmdj(*)
  integer, intent(out) :: info
  !
  integer :: k, h
  integer, parameter :: k_max = 50
  real(kind(0.0d0)) :: s, r, p 
  real(kind(0.0d0)) :: alpha, beta ! Auxiliary variables for computing coef(:)
  real(kind(0.0d0)) :: coef(4) ! Coefficients of the cubic equation.
  real(kind(0.0d0)) :: lower, upper
  integer :: m
  real(kind(0.0d0)) :: root(3), rwork(3 * 3)
  integer :: m_in_interval
  real(kind(0.0d0)) :: y_candidate

  !> Compute the lower and upper bounds and set an initial approximation.
  lower = d(j) - d(j - 1)
  upper = 0.0d0
  y = 0.25d0 * lower ! An initial approximation

  !> Iteratively compute the root.
  info = 100 ! Zero will be set to info if converged.
  s_numerator(1 : n) = (z(1 : n) ** 2.0d0) * (d(1 : n) - d(j - 1))
  r_numerator(1 : n) = (z(1 : n) ** 2.0d0) * (d(1 : n) - d(j))
  p_numerator(1 : n) = (z(1 : n) ** 2.0d0) * (d(1 : n) - d(j - 1)) * (d(1 : n) - d(j))
  dmdj(1 : n) = d(1 : n) - d(j) 
  do k = 1, k_max
    ! Compute s, r, p
    invdenominator(1 : n) = (dmdj(1 : n) + y) ** (- 3.0d0)
    !
    term(1 : j - 2) = s_numerator(1 : j - 2) * invdenominator(1 : j - 2)
    term(j - 1 : n - 2) = s_numerator(j + 1 : n) * invdenominator(j + 1 : n)
    call accurate_sum(n - 2, term, s)
    s = z(j) ** 2.0d0 + s * (y ** 3.0d0) / (d(j) - d(j - 1))
    !
    term(1 : j - 2) = r_numerator(1 : j - 2) * invdenominator(1 : j - 2)
    term(j - 1 : n - 2) = r_numerator(j + 1 : n) * invdenominator(j + 1 : n)
    call accurate_sum(n - 2, term, r)
    r = z(j - 1) ** 2.0d0 - r * ((d(j - 1) - d(j) + y) ** 3.0d0) / (d(j) - d(j - 1))
    !
    term(1 : j - 2) = p_numerator(1 : j - 2) * invdenominator(1 : j - 2)
    term(j - 1 : n - 2) = p_numerator(j + 1 : n) * invdenominator(j + 1 : n)
    call accurate_sum(n - 2, term, p)
    ! Check
    if (r .le. 0.0d0 .or. s .le. 0.0d0) then
      write (*, *) 'Cannot happen. Bad s or r. s, r = ', s, r
      stop
    endif
    ! Compute new coefficients from s, r, p
    alpha = d(j - 1) - d(j)
    beta = rho - d(j)
    coef(1) = s * alpha * beta
    coef(2) = p * alpha * beta + s * (alpha + beta) + r * beta - alpha
    coef(3) = p * (alpha + beta) + r + s - 1.0d0
    coef(4) = p
    ! Solve a cubic equation
    call solve_real_algebraic_equation(3, coef, maxval(dabs(coef(1 : 4))) * 1.0e-18, m, root, rwork)
    ! Select a solution in the target interval as a new approximation.
    if (m .lt. 1 .or. m .gt. 3) then
      write (*, *) 'Cannot happen. Bad m. m = ', m
      stop
    else
      m_in_interval = 0
      y_candidate = - huge(0.0d0) ! To suppress an uninitialized warning.
      do h = 1, m
        if (root(h) .gt. upper .and. root(h) .lt. lower) then
          if (m_in_interval .ne. 0) then
            write (*, *) '[Warning] Multiple roots exist in the target interval.'
            if (dabs(y_candidate - y) .gt. dabs(root(h) - y)) then
              ! Select the nearest root from lambda.
              m_in_interval = m_in_interval + 1
              y_candidate = root(h)
            endif
          else
            m_in_interval = 1
            y_candidate = root(h)
          endif
        endif
      enddo
      if (m_in_interval .le. 0) then
        write (*, *) '[right] Cannot happen. Bad m_in_interval. j, m, m_in_interval = ', j, m, m_in_interval
        write (*, *) '[right] Cannot happen. Bad m_in_interval. root(1 : m) = ', root(1 : m)
        write (*, *) '[right] Cannot happen. Bad m_in_interval. d(j) - root(1 : m) = ', d(j) - root(1 : m)
        write (*, *) '[right] Cannot happen. Bad m_in_interval. d(j - 1), d(j) = ', d(j - 1), d(j)
        write (*, *) '[right] Cannot happen. Bad m_in_interval. coef(1 : 4) = ', coef(1 : 4)
        write (*, *) '[right] Cannot happen. Bad m_in_interval.', coef(1) / coef(4) / root(2) / root(3)
        stop
      endif
      if (dabs(y - y_candidate) <= 0.0d0) then
        info = 0
        exit
      endif
      y = y_candidate
    endif
  enddo

  return
end

!> @brief    This subroutine is essentially same as gevd_dc_secular_elsner_middle_left(...).
!! @sa       gevd_dc_secular_elsner_middle_left.
subroutine gevd_dc_secular_elsner_middle_left(n, j, d, z, rho, y, s_numerator, r_numerator, p_numerator, &
                                              invdenominator, term, dmdjm1, info)
  implicit none
  integer, intent(in) :: n, j ! Interval number, (d_{j-1}, d_j)
  real(kind(0.0d0)), intent(in) :: d(*), z(*), rho
  real(kind(0.0d0)), intent(out) :: y ! d(j - 1) - lambda < 0
  real(kind(0.0d0)), intent(out) :: s_numerator(*), r_numerator(*), p_numerator(*)
  real(kind(0.0d0)), intent(out) :: invdenominator(*),  term(*), dmdjm1(*)
  integer, intent(out) :: info
  !
  integer :: k, h
  integer, parameter :: k_max = 50
  real(kind(0.0d0)) :: s, r, p 
  real(kind(0.0d0)) :: alpha, beta ! Auxiliary variables for computing coef(:)
  real(kind(0.0d0)) :: coef(4) ! Coefficients of the cubic equation.
  real(kind(0.0d0)) :: lower, upper
  integer :: m
  real(kind(0.0d0)) :: root(3), rwork(3 * 3)
  integer :: m_in_interval
  real(kind(0.0d0)) :: y_candidate

  !> Compute the lower and upper bounds and set an initial approximation.
  lower = 0.0d0
  upper = d(j - 1) - d(j)
  y = 0.25d0 * upper

  !> Iteratively compute the root.
  info = 100 ! Zero will be set to info if converged.
  s_numerator(1 : n) = (z(1 : n) ** 2.0d0) * (d(1 : n) - d(j - 1))
  r_numerator(1 : n) = (z(1 : n) ** 2.0d0) * (d(1 : n) - d(j))
  p_numerator(1 : n) = (z(1 : n) ** 2.0d0) * (d(1 : n) - d(j - 1)) * (d(1 : n) - d(j))
  dmdjm1(1 : n) = d(1 : n) - d(j - 1)
  do k = 1, k_max
    ! Compute s, r, p
    invdenominator(1 : n) = (dmdjm1(1 : n) + y) ** (-3.0d0)
    term(1 : j - 2) = s_numerator(1 : j - 2) * invdenominator(1 : j - 2)
    term(j - 1 : n - 2) = s_numerator(j + 1 : n) * invdenominator(j + 1 : n)
    call accurate_sum(n - 2, term, s)
    s = z(j) ** 2.0d0 + s * (((d(j) - d(j- 1)) + y) ** 3.0d0) / (d(j) - d(j - 1))
    !write (*, *) 'y', y
    !
    term(1 : j - 2) = r_numerator(1 : j - 2) * invdenominator(1 : j - 2)
    term(j - 1 : n - 2) = r_numerator(j + 1 : n) * invdenominator(j + 1 : n)
    call accurate_sum(n - 2, term, r)
    r = z(j - 1) ** 2.0d0 - r * (y ** 3.0d0) / (d(j) - d(j - 1))
    !
    term(1 : j - 2) = p_numerator(1 : j - 2) * invdenominator(1 : j - 2)
    term(j - 1 : n - 2) = p_numerator(j + 1 : n) * invdenominator(j + 1 : n)
    call accurate_sum(n - 2, term, p)
    ! Check
    if (s .le. 0.0d0 .or. r .le. 0.0d0) then
      write (*, *) '[left] Cannot happen. Bad s or r. s, r = ', s, r
      write (*, *) '[left] Cannot happen. Bad s or r. n, j = ', n, j
      write (*, *) '[left] Cannot happen. Bad s or y.', y
      write (*, *) '[left] Cannot happen. Bad s or r.', z(j) ** 2.0d0
      stop
    endif
    ! Compute new coefficients from s, r, p
    alpha = d(j) - d(j - 1)
    beta = rho - d(j - 1)
    coef(1) = r * alpha * beta
    coef(2) = p * alpha * beta + r * (alpha + beta) + s * beta - alpha
    coef(3) = p * (alpha + beta) + r + s - 1.0d0
    coef(4) = p
    ! Solve a cubic equation
    call solve_real_algebraic_equation(3, coef, maxval(dabs(coef(1 : 4))) * 1.0e-18, m, root, rwork)
    ! Select a solution in the target interval as a new approximation.
    if (m .lt. 1 .or. m .gt. 3) then
      write (*, *) 'Cannot happen. Bad m. m = ', m
      stop
    else
      m_in_interval = 0
      y_candidate = huge(0.0d0) ! To suppress an uninitialized warning.
      do h = 1, m
        if (root(h) .gt. upper .and. root(h) .lt. lower) then
          if (m_in_interval .ne. 0) then
            write (*, *) '[Warning] Multiple roots exist in the target interval.'
            if (dabs(y_candidate - y) .gt. dabs(root(h) - y)) then
              ! Select the nearest root from lambda.
              m_in_interval = m_in_interval + 1
              y_candidate = root(h)
            endif
          else
            m_in_interval = 1
            y_candidate = root(h)
          endif
        endif
      enddo
      if (m_in_interval .le. 0) then
        write (*, *) '[right] Cannot happen. Bad m_in_interval. j, m, m_in_interval = ', j, m, m_in_interval
        stop
      endif
      if (dabs(y - y_candidate) <= 0.0d0) then
        info = 0
        exit
      endif
      y = y_candidate
    endif
  enddo

  return
end

!> @brief    Compute the root of g(\\lambda) = 0 in (d_{j-1}, d_j).
!!           where g(x) := \\sum_{k = 1}^{n} \\frac{z_k^2}{d_k - x} - \\frac{1}{\\rho - x}.
!! @note     \\lambda is implicitly expressed as y := d(nearest_root_index) - \\lambda.
!! @note     This subroutine does not check the argument validity.
!!
!! @param n                   The dimension of the vectors d and z.
!! @param j                   The interval number, (d_{j-1}, d_j).
!! @param d                   The vector d. dimension(n). d(1) < d(2) < ... < d(n).
!! @param z                   The vector z. dimension(n). z(i) .ne. 0 (i = 1, ..., n) .
!! @param rho                 The scalar \\rho.
!! @param y                   The root of the function g.
!! @param nearest_root_index  d(nearest_root_index) is the nearest poll from the root.
!! @param rwork               Temporary array. dimension(6 * n).
!! @param info                Error code.
subroutine gevd_dc_secular_elsner_middle(n, j, d, z, rho, y, nearest_root_index, rwork, info)
  implicit none
  integer, intent(in) :: n, j
  real(kind(0.0d0)), intent(in) :: d(*), z(*), rho
  real(kind(0.0d0)), intent(out) :: y 
  integer, intent(out) :: nearest_root_index
  real(kind(0.0d0)), intent(out) :: rwork(*)
  integer, intent(out) :: info
  !
  real(kind(0.0d0)) :: g

  ! Find the nearest poll.
  call gevd_dc_secular_elsner_evaluation(n, d, z, rho, (d(j - 1) + d(j)) * 0.5d0, g, info)
  if (g > 0.0d0) then
    nearest_root_index = j - 1
    call gevd_dc_secular_elsner_middle_left(n, j, d, z, rho, y, rwork(1), rwork(n + 1), rwork(n * 2 + 1), &
                                            rwork(n * 3 + 1), rwork(n * 4 + 1), rwork(n * 5 + 1), info)
  else
    nearest_root_index = j
    call gevd_dc_secular_elsner_middle_right(n, j, d, z, rho, y, rwork(1), rwork(n + 1), rwork(n * 2 + 1), &
                                             rwork(n * 3 + 1), rwork(n * 4 + 1), rwork(n * 5 + 1), info)
  endif

  return
end

!> @brief    Compute the root of i-th smallest root \\lambda of the function
!!           g(x) := \\sum_{k = 1}^{n} \\frac{z_k^2}{d_k - x} - \\frac{1}{\\rho - x}
!!           accurately (with less cancellations)
!!           (i.e. compute i-th smallest eigenvalues of D - \rho z z' - \lambda (I - z z')).
!! @note     \\lambda is implicitly expressed as y := d(nearest_root_index) - \\lambda.
!!           This subroutine does not check the argument validity.
!!
!! @param n                   The dimension of the vectors d and z.
!! @param i                   The root number.
!! @param d                   The vector d. dimension(n). d(1) < d(2) < ... < d(n).
!! @param z                   The vector z. dimension(n).z(i) .ne. 0 (i = 1, ..., n).
!! @param rho                 The scalar \\rho.
!! @param y                   The root of the function g.
!! @param nearest_root_index  d(nearest_root_index) is the nearest poll from the root.
!! @param rwork               Temporary array. dimension(6 * n).
!! @param info                Error code.
subroutine gevd_dc_secular_elsner(n, i, d, z, rho, y, nearest_root_index, rwork, info)
  implicit none
  integer, intent(in) :: n, i
  real(kind(0.0d0)), intent(in) :: d(*), z(*), rho
  real(kind(0.0d0)), intent(out) :: y
  integer, intent(out) :: nearest_root_index
  real(kind(0.0d0)), intent(out) :: rwork(*)
  integer, intent(out) :: info
  !
  integer :: j ! Current interval (d_{j-1}, d_j)

  !> Find the interval which includes the root.
  if (d(i) .lt. rho) then
    ! The current interval is left of rho
    j = i
  else
    j = i + 1
  endif

  !> Call one of the three root finder subroutines (leftmost, rightmost and middle).
  if (j < 1 .or. j > n + 1) then
    write (*, *) '[Error] Cannot happen. Bad j. n, j = ', n, j
    info = 10
    return 
  else if (j == 1) then
    ! Leftmost 
    nearest_root_index = 1
    call gevd_dc_secular_elsner_leftmost(n, d, z, rho, y, rwork(1), rwork(n + 1), rwork(n * 2 + 1), &
                                         rwork(n * 3 + 1), rwork(n * 4 + 1), info)
  else if (j == (n + 1)) then
    ! Rightmost
    nearest_root_index = n
    call gevd_dc_secular_elsner_rightmost(n, d, z, rho, y, rwork(1), rwork(n + 1), rwork(n * 2 + 1), &
                                          rwork(n * 3 + 1), rwork(n * 4 + 1), info)
  else
    ! Middle
    call gevd_dc_secular_elsner_middle(n, j, d, z, rho, y, nearest_root_index, rwork, info)
  endif
  if (info /= 0) then
    info = 20 ! Most probably the iteration in gevd_dc_secular_elsner_* did not converge.
    return 
  endif
  call assert_y_is_valid(j, nearest_root_index, d, y) ! Check the validity of the solution.
  
  return
contains
  subroutine assert_y_is_valid(j, nearest_root_index, d, y)
    implicit none
    integer, intent(in) :: j 
    integer, intent(in) :: nearest_root_index
    real(kind(0.0d0)), intent(in) :: d(*)
    real(kind(0.0d0)), intent(in) :: y

    if (j == 1) then
      ! Error check (leftmost).
      if (y <= 0.0d0) then
        write (*, *) '[Error] Condition y > 0 is not satisfied. y = ', y
        call abort()
      endif
    else if (j == n + 1) then
      ! Error check (rightmost).
      if (y >= 0.0d0) then
        write (*, *) '[Error] Condition y < 0 is not satisfied. y = ', y
        call abort()
      endif
    else if (nearest_root_index == j - 1) then
      if (y >= 0.0 .or. y <= d(j - 1) - d(j)) then
        write (*, *) '[Error] Bad y. j - 1, d(j - 1) - d(j), y = ', j - 1, d(j - 1) - d(j), y
        call abort()
      endif
    else if (nearest_root_index == j) then
      if (y <= 0.0 .or. y >= d(j) - d(j - 1)) then
        write (*, *) '[Error] Bad y. j, d(j) - d(j - 1), y = ', j, d(j) - d(j - 1), y
        call abort()
      endif
    endif
  end subroutine
end

!> @brief    Evaluate g(x) = \\sum_{i = 1}^n \\frac{w_i}{d_i - \\lambda} - \\frac{1}{\\rho - \\lambda}
!!           The result is not very accurate. This subroutine is used for finding nearest poll
!!           from the root in Elsner's iterative method.
!! @note     This subroutine does not check the argument validity.
!!
!! @param n       The dimension of the vectors d and z.
!! @param j       The interval number.
!! @param d       The vector d. dimension(n). d(1) < d(2) < ... < d(n).
!! @param z       The vector z. dimension(n). z(i) .ne. 0 (i = 1, ..., n).
!! @param rho     The scalar \rho.
!! @param x       The argument of the function.
!! @param g       The function value.
!! @param info    Error code. This subroutine always returns info = 0.
subroutine gevd_dc_secular_elsner_evaluation(n, d, z, rho, x, g, info)
  implicit none
  integer, intent(in) :: n
  real(kind(0.0d0)), intent(in) :: d(*), z(*), rho, x
  real(kind(0.0d0)), intent(out) :: g
  integer, intent(out) :: info
  !
  integer :: k
  real(kind(0.0d0)) :: plus, minus, this_term

  plus = 0.0d0
  minus = 0.0d0
  do k = 1, n
    this_term = z(k) * z(k) / (d(k) - x)
    if (this_term > 0.0d0) then
      plus = plus + this_term
    else
      minus = minus + this_term
    endif
  enddo
  this_term = - 1.0d0 / (rho - x)
  if (this_term > 0.0d0) then
    plus = plus + this_term
  else
    minus = minus + this_term
  endif
  g = plus + minus

  info = 0
  return
end

