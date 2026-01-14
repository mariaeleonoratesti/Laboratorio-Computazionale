! ================================================================================
! XY MODEL - COMPLETE FORTRAN MODULE FOR PROJECT 17.25
! ================================================================================
! Gould & Tobochnik style with proper f2py interfaces
!
! Compile with:
!   python -m numpy.f2py -c xy_fortran.f90 -m xy_fortran --opt='-O3 -march=native'
! ================================================================================

module xy_model
  implicit none
  
  ! Constants
  real(8), parameter :: PI = 3.14159265358979323846d0
  real(8), parameter :: TWO_PI = 6.28318530717958647692d0
  
  ! Lookup tables
  integer, parameter :: EXP_TABLE_SIZE = 8001
  real(8) :: exp_table(0:EXP_TABLE_SIZE)
  real(8) :: current_beta = -1.0d0
  
contains

  subroutine init_exp_table(beta)
    real(8), intent(in) :: beta
    integer :: j
    real(8) :: dE
    
    if (abs(beta - current_beta) < 1.0d-10) return
    
    do j = 0, EXP_TABLE_SIZE
      dE = real(j, 8) * 0.001d0
      exp_table(j) = exp(-beta * dE)
    end do
    current_beta = beta
  end subroutine init_exp_table
  
  subroutine init_random_config(L, theta)
    !f2py intent(in) :: L
    !f2py intent(out) :: theta
    !f2py depend(L) :: theta
    integer, intent(in) :: L
    real(8), intent(out) :: theta(L, L)
    
    call random_number(theta)
    theta = theta * TWO_PI
  end subroutine init_random_config
  
  subroutine init_cold_config(L, theta)
    !f2py intent(in) :: L
    !f2py intent(out) :: theta
    !f2py depend(L) :: theta
    integer, intent(in) :: L
    real(8), intent(out) :: theta(L, L)
    
    theta = 0.0d0
  end subroutine init_cold_config
  
  subroutine set_random_seed(seed)
    !f2py intent(in) :: seed
    integer, intent(in) :: seed
    integer :: n, i
    integer, allocatable :: seed_array(:)
    
    call random_seed(size=n)
    allocate(seed_array(n))
    do i = 1, n
      seed_array(i) = seed + i * 37
    end do
    call random_seed(put=seed_array)
    deallocate(seed_array)
  end subroutine set_random_seed

  pure function boltzmann_factor(dE) result(prob)
    real(8), intent(in) :: dE
    real(8) :: prob
    integer :: j
    
    if (dE <= 0.0d0) then
      prob = 1.0d0
    else
      j = int(dE * 1000.0d0)
      if (j > EXP_TABLE_SIZE) then
        prob = 0.0d0
      else
        prob = exp_table(j)
      end if
    end if
  end function boltzmann_factor
  
  subroutine compute_energy(L, theta, E)
    !f2py intent(in) :: L, theta
    !f2py intent(out) :: E
    !f2py depend(L) :: theta
    integer, intent(in) :: L
    real(8), intent(in) :: theta(L, L)
    real(8), intent(out) :: E
    integer :: x, y, xp, yp
    
    E = 0.0d0
    do y = 1, L
      yp = mod(y, L) + 1
      do x = 1, L
        xp = mod(x, L) + 1
        E = E - cos(theta(x,y) - theta(xp,y))
        E = E - cos(theta(x,y) - theta(x,yp))
      end do
    end do
  end subroutine compute_energy
  
  subroutine metropolis_sweep(L, theta, beta, delta, E, n_accept)
    !f2py intent(in) :: L, beta, delta
    !f2py intent(inout) :: theta, E
    !f2py intent(out) :: n_accept
    !f2py depend(L) :: theta
    integer, intent(in) :: L
    real(8), intent(inout) :: theta(L, L)
    real(8), intent(in) :: beta, delta
    real(8), intent(inout) :: E
    integer, intent(out) :: n_accept
    
    integer :: x, y, xm, xp, ym, yp
    real(8) :: theta_old, theta_new, E_old, E_new, dE, rnd
    
    call init_exp_table(beta)
    n_accept = 0
    
    do y = 1, L
      ym = merge(L, y-1, y == 1)
      yp = merge(1, y+1, y == L)
      
      do x = 1, L
        xm = merge(L, x-1, x == 1)
        xp = merge(1, x+1, x == L)
        
        theta_old = theta(x, y)
        
        call random_number(rnd)
        theta_new = theta_old + delta * (rnd - 0.5d0)
        
        E_old = -cos(theta_old - theta(xm,y)) - cos(theta_old - theta(xp,y)) &
               -cos(theta_old - theta(x,ym)) - cos(theta_old - theta(x,yp))
        
        E_new = -cos(theta_new - theta(xm,y)) - cos(theta_new - theta(xp,y)) &
               -cos(theta_new - theta(x,ym)) - cos(theta_new - theta(x,yp))
        
        dE = E_new - E_old
        
        call random_number(rnd)
        if (rnd < boltzmann_factor(dE)) then
          theta(x, y) = theta_new
          E = E + dE
          n_accept = n_accept + 1
        end if
      end do
    end do
  end subroutine metropolis_sweep

  subroutine compute_magnetization(L, theta, mx, my, M_abs)
    !f2py intent(in) :: L, theta
    !f2py intent(out) :: mx, my, M_abs
    !f2py depend(L) :: theta
    integer, intent(in) :: L
    real(8), intent(in) :: theta(L, L)
    real(8), intent(out) :: mx, my, M_abs
    integer :: x, y
    
    mx = 0.0d0
    my = 0.0d0
    do y = 1, L
      do x = 1, L
        mx = mx + cos(theta(x,y))
        my = my + sin(theta(x,y))
      end do
    end do
    M_abs = sqrt(mx*mx + my*my) / real(L*L, 8)
  end subroutine compute_magnetization
  
  subroutine compute_theta2(L, theta, theta2)
    !f2py intent(in) :: L, theta
    !f2py intent(out) :: theta2
    !f2py depend(L) :: theta
    integer, intent(in) :: L
    real(8), intent(in) :: theta(L, L)
    real(8), intent(out) :: theta2
    integer :: x, y
    real(8) :: mx, my, M_mag, sx, sy, dot_prod, angle
    
    mx = 0.0d0
    my = 0.0d0
    do y = 1, L
      do x = 1, L
        mx = mx + cos(theta(x,y))
        my = my + sin(theta(x,y))
      end do
    end do
    M_mag = sqrt(mx*mx + my*my)
    if (M_mag < 1.0d-15) M_mag = 1.0d-15
    
    theta2 = 0.0d0
    do y = 1, L
      do x = 1, L
        sx = cos(theta(x,y))
        sy = sin(theta(x,y))
        dot_prod = (sx*mx + sy*my) / M_mag
        dot_prod = max(-1.0d0, min(1.0d0, dot_prod))
        angle = acos(dot_prod)
        theta2 = theta2 + angle * angle
      end do
    end do
    theta2 = theta2 / real(L*L, 8)
  end subroutine compute_theta2

  subroutine find_vortices(L, theta, max_vortices, n_positive, n_negative, &
                           n_vortices, vortex_x, vortex_y, vortex_charge)
    !f2py intent(in) :: L, theta, max_vortices
    !f2py intent(out) :: n_positive, n_negative, n_vortices
    !f2py intent(out) :: vortex_x, vortex_y, vortex_charge
    !f2py depend(L) :: theta
    !f2py depend(max_vortices) :: vortex_x, vortex_y, vortex_charge
    integer, intent(in) :: L, max_vortices
    real(8), intent(in) :: theta(L, L)
    integer, intent(out) :: n_positive, n_negative, n_vortices
    real(8), intent(out) :: vortex_x(max_vortices), vortex_y(max_vortices)
    integer, intent(out) :: vortex_charge(max_vortices)
    
    integer :: x, y, xp, yp, i
    real(8) :: angles(5), delta_theta, total_rotation
    
    n_positive = 0
    n_negative = 0
    n_vortices = 0
    vortex_x = 0.0d0
    vortex_y = 0.0d0
    vortex_charge = 0
    
    do y = 1, L
      yp = mod(y, L) + 1
      do x = 1, L
        xp = mod(x, L) + 1
        
        angles(1) = theta(x, y)
        angles(2) = theta(xp, y)
        angles(3) = theta(xp, yp)
        angles(4) = theta(x, yp)
        angles(5) = theta(x, y)
        
        total_rotation = 0.0d0
        do i = 1, 4
          delta_theta = angles(i+1) - angles(i)
          do while (delta_theta > PI)
            delta_theta = delta_theta - TWO_PI
          end do
          do while (delta_theta < -PI)
            delta_theta = delta_theta + TWO_PI
          end do
          total_rotation = total_rotation + delta_theta
        end do
        
        if (total_rotation > PI) then
          n_positive = n_positive + 1
          n_vortices = n_vortices + 1
          if (n_vortices <= max_vortices) then
            vortex_x(n_vortices) = real(x, 8) + 0.5d0
            vortex_y(n_vortices) = real(y, 8) + 0.5d0
            vortex_charge(n_vortices) = 1
          end if
        else if (total_rotation < -PI) then
          n_negative = n_negative + 1
          n_vortices = n_vortices + 1
          if (n_vortices <= max_vortices) then
            vortex_x(n_vortices) = real(x, 8) + 0.5d0
            vortex_y(n_vortices) = real(y, 8) + 0.5d0
            vortex_charge(n_vortices) = -1
          end if
        end if
      end do
    end do
  end subroutine find_vortices

  subroutine compute_correlation(L, theta, max_r, C_r)
    !f2py intent(in) :: L, theta, max_r
    !f2py intent(out) :: C_r
    !f2py depend(L) :: theta
    !f2py depend(max_r) :: C_r
    integer, intent(in) :: L, max_r
    real(8), intent(in) :: theta(L, L)
    real(8), intent(out) :: C_r(0:max_r)
    
    integer :: x, y, dx, dy, r, x2, y2
    integer :: count(0:max_r)
    real(8) :: corr
    
    C_r = 0.0d0
    count = 0
    
    do y = 1, L
      do x = 1, L
        do dy = 0, L/2
          do dx = 0, L/2
            r = nint(sqrt(real(dx*dx + dy*dy, 8)))
            if (r <= max_r) then
              x2 = mod(x + dx - 1, L) + 1
              y2 = mod(y + dy - 1, L) + 1
              corr = cos(theta(x,y) - theta(x2,y2))
              C_r(r) = C_r(r) + corr
              count(r) = count(r) + 1
            end if
          end do
        end do
      end do
    end do
    
    do r = 0, max_r
      if (count(r) > 0) then
        C_r(r) = C_r(r) / real(count(r), 8)
      end if
    end do
  end subroutine compute_correlation

  subroutine run_equilibration(L, theta, beta, delta_in, n_sweeps, target_accept, &
                               delta_out, accept_out)
    !f2py intent(in) :: L, beta, delta_in, n_sweeps, target_accept
    !f2py intent(inout) :: theta
    !f2py intent(out) :: delta_out, accept_out
    !f2py depend(L) :: theta
    integer, intent(in) :: L, n_sweeps
    real(8), intent(inout) :: theta(L, L)
    real(8), intent(in) :: beta, delta_in, target_accept
    real(8), intent(out) :: delta_out, accept_out
    
    real(8) :: E, delta, accept_rate
    integer :: sweep, n_accept, total_accept, adjust_interval
    
    call compute_energy(L, theta, E)
    delta = delta_in
    adjust_interval = 10
    total_accept = 0
    
    do sweep = 1, n_sweeps
      call metropolis_sweep(L, theta, beta, delta, E, n_accept)
      total_accept = total_accept + n_accept
      
      if (mod(sweep, adjust_interval) == 0) then
        accept_rate = real(total_accept, 8) / real(adjust_interval * L * L, 8)
        
        if (accept_rate > target_accept + 0.05d0) then
          delta = delta * 1.1d0
        else if (accept_rate < target_accept - 0.05d0) then
          delta = delta * 0.9d0
        end if
        
        delta = max(0.1d0, min(delta, TWO_PI))
        total_accept = 0
      end if
    end do
    
    delta_out = delta
    
    total_accept = 0
    do sweep = 1, 100
      call metropolis_sweep(L, theta, beta, delta, E, n_accept)
      total_accept = total_accept + n_accept
    end do
    accept_out = real(total_accept, 8) / real(100 * L * L, 8)
  end subroutine run_equilibration

  subroutine run_measurement(L, theta, beta, delta, n_sweeps, &
                             E_series, M_series, M2_series, theta2_series, &
                             vort_series, E2_series)
    !f2py intent(in) :: L, beta, delta, n_sweeps
    !f2py intent(inout) :: theta
    !f2py intent(out) :: E_series, M_series, M2_series, theta2_series, vort_series, E2_series
    !f2py depend(L) :: theta
    !f2py depend(n_sweeps) :: E_series, M_series, M2_series, theta2_series, vort_series, E2_series
    integer, intent(in) :: L, n_sweeps
    real(8), intent(inout) :: theta(L, L)
    real(8), intent(in) :: beta, delta
    real(8), intent(out) :: E_series(n_sweeps), M_series(n_sweeps)
    real(8), intent(out) :: M2_series(n_sweeps), theta2_series(n_sweeps)
    real(8), intent(out) :: vort_series(n_sweeps), E2_series(n_sweeps)
    
    real(8) :: E, mx, my, M_abs, t2
    integer :: sweep, n_accept, n_pos, n_neg, n_vort, max_v
    real(8), allocatable :: vx(:), vy(:)
    integer, allocatable :: vc(:)
    
    max_v = L * L / 2
    allocate(vx(max_v), vy(max_v), vc(max_v))
    
    call compute_energy(L, theta, E)
    
    do sweep = 1, n_sweeps
      call metropolis_sweep(L, theta, beta, delta, E, n_accept)
      
      E_series(sweep) = E / real(L*L, 8)
      E2_series(sweep) = (E / real(L*L, 8))**2
      
      call compute_magnetization(L, theta, mx, my, M_abs)
      M_series(sweep) = M_abs
      M2_series(sweep) = M_abs * M_abs
      
      call compute_theta2(L, theta, t2)
      theta2_series(sweep) = t2
      
      call find_vortices(L, theta, max_v, n_pos, n_neg, n_vort, vx, vy, vc)
      vort_series(sweep) = real(n_vort, 8)
    end do
    
    deallocate(vx, vy, vc)
  end subroutine run_measurement
  
  subroutine study_theta2_vs_n(n_L, L_values, T, n_equil, n_meas, &
                               theta2_means, theta2_errors, delta_out, accept_out)
    !f2py intent(in) :: n_L, L_values, T, n_equil, n_meas
    !f2py intent(out) :: theta2_means, theta2_errors, delta_out, accept_out
    !f2py depend(n_L) :: L_values, theta2_means, theta2_errors, delta_out, accept_out
    integer, intent(in) :: n_L, n_equil, n_meas
    integer, intent(in) :: L_values(n_L)
    real(8), intent(in) :: T
    real(8), intent(out) :: theta2_means(n_L), theta2_errors(n_L), delta_out(n_L), accept_out(n_L)
    
    real(8), allocatable :: theta(:,:), theta2_series(:)
    real(8) :: beta, delta, final_delta, final_accept, E, t2
    integer :: i, L, sweep, n_accept
    
    beta = 1.0d0 / T
    
    do i = 1, n_L
      L = L_values(i)
      allocate(theta(L, L), theta2_series(n_meas))
      
      call init_random_config(L, theta)
      
      delta = 1.0d0
      call run_equilibration(L, theta, beta, delta, n_equil, 0.4d0, &
                            final_delta, final_accept)
      delta_out(i) = final_delta
      accept_out(i) = final_accept
      
      call compute_energy(L, theta, E)
      do sweep = 1, n_meas
        call metropolis_sweep(L, theta, beta, final_delta, E, n_accept)
        call compute_theta2(L, theta, t2)
        theta2_series(sweep) = t2
      end do
      
      theta2_means(i) = sum(theta2_series) / real(n_meas, 8)
      theta2_errors(i) = sqrt(sum((theta2_series - theta2_means(i))**2) / &
                         real(n_meas * (n_meas - 1), 8))
      
      deallocate(theta, theta2_series)
    end do
  end subroutine study_theta2_vs_n

  subroutine temperature_quench(L, T_final, n_sweeps, theta, &
                                n_pos, n_neg, delta_out, accept_out)
    !f2py intent(in) :: L, T_final, n_sweeps
    !f2py intent(out) :: theta, n_pos, n_neg, delta_out, accept_out
    !f2py depend(L) :: theta
    integer, intent(in) :: L, n_sweeps
    real(8), intent(in) :: T_final
    real(8), intent(out) :: theta(L, L)
    integer, intent(out) :: n_pos, n_neg
    real(8), intent(out) :: delta_out, accept_out
    
    real(8) :: beta, delta, final_delta, final_accept
    integer :: n_vort, max_v
    real(8), allocatable :: vx(:), vy(:)
    integer, allocatable :: vc(:)
    
    max_v = L * L / 2
    allocate(vx(max_v), vy(max_v), vc(max_v))
    
    call init_random_config(L, theta)
    
    beta = 1.0d0 / T_final
    delta = 1.0d0
    
    call run_equilibration(L, theta, beta, delta, n_sweeps, 0.4d0, &
                          final_delta, final_accept)
    
    call find_vortices(L, theta, max_v, n_pos, n_neg, n_vort, vx, vy, vc)
    
    delta_out = final_delta
    accept_out = final_accept
    
    deallocate(vx, vy, vc)
  end subroutine temperature_quench

  subroutine temperature_scan(L, n_T, T_values, n_equil, n_meas, &
                              E_means, E_errors, M_means, M_errors, &
                              M2_means, Chi_values, Cv_values, &
                              vort_means, vort_errors, &
                              theta_final, accept_rates, deltas)
    !f2py intent(in) :: L, n_T, T_values, n_equil, n_meas
    !f2py intent(out) :: E_means, E_errors, M_means, M_errors
    !f2py intent(out) :: M2_means, Chi_values, Cv_values
    !f2py intent(out) :: vort_means, vort_errors
    !f2py intent(out) :: theta_final, accept_rates, deltas
    !f2py depend(n_T) :: T_values, E_means, E_errors, M_means, M_errors
    !f2py depend(n_T) :: M2_means, Chi_values, Cv_values, vort_means, vort_errors
    !f2py depend(n_T) :: accept_rates, deltas
    !f2py depend(L, n_T) :: theta_final
    integer, intent(in) :: L, n_T, n_equil, n_meas
    real(8), intent(in) :: T_values(n_T)
    real(8), intent(out) :: E_means(n_T), E_errors(n_T)
    real(8), intent(out) :: M_means(n_T), M_errors(n_T)
    real(8), intent(out) :: M2_means(n_T)
    real(8), intent(out) :: Chi_values(n_T), Cv_values(n_T)
    real(8), intent(out) :: vort_means(n_T), vort_errors(n_T)
    real(8), intent(out) :: theta_final(L, L, n_T)
    real(8), intent(out) :: accept_rates(n_T), deltas(n_T)
    
    real(8), allocatable :: theta(:,:)
    real(8), allocatable :: E_series(:), M_series(:), M2_series(:)
    real(8), allocatable :: theta2_series(:), vort_series(:), E2_series(:)
    real(8) :: T, beta, delta, final_delta, final_accept
    real(8) :: E_mean, E2_mean, M_mean, M2_mean, vort_mean
    integer :: i
    
    allocate(theta(L, L))
    allocate(E_series(n_meas), M_series(n_meas), M2_series(n_meas))
    allocate(theta2_series(n_meas), vort_series(n_meas), E2_series(n_meas))
    
    call init_cold_config(L, theta)
    delta = 1.0d0
    
    do i = 1, n_T
      T = T_values(i)
      beta = 1.0d0 / T
      
      call run_equilibration(L, theta, beta, delta, n_equil, 0.4d0, &
                            final_delta, final_accept)
      delta = final_delta
      deltas(i) = final_delta
      accept_rates(i) = final_accept
      
      call run_measurement(L, theta, beta, delta, n_meas, &
                          E_series, M_series, M2_series, theta2_series, &
                          vort_series, E2_series)
      
      E_mean = sum(E_series) / real(n_meas, 8)
      E2_mean = sum(E2_series) / real(n_meas, 8)
      M_mean = sum(M_series) / real(n_meas, 8)
      M2_mean = sum(M2_series) / real(n_meas, 8)
      vort_mean = sum(vort_series) / real(n_meas, 8)
      
      E_means(i) = E_mean
      M_means(i) = M_mean
      M2_means(i) = M2_mean
      vort_means(i) = vort_mean / real(L*L, 8)
      
      E_errors(i) = sqrt(sum((E_series - E_mean)**2) / real(n_meas*(n_meas-1), 8))
      M_errors(i) = sqrt(sum((M_series - M_mean)**2) / real(n_meas*(n_meas-1), 8))
      vort_errors(i) = sqrt(sum((vort_series/real(L*L,8) - vort_means(i))**2) / &
                       real(n_meas*(n_meas-1), 8))
      
      Chi_values(i) = real(L*L, 8) * M2_mean / T
      Cv_values(i) = real(L*L, 8) * (E2_mean - E_mean**2) / (T*T)
      
      theta_final(:,:,i) = theta
    end do
    
    deallocate(theta, E_series, M_series, M2_series)
    deallocate(theta2_series, vort_series, E2_series)
  end subroutine temperature_scan

  subroutine fit_kt_transition(n_T, T_values, Chi_values, T_min, T_max, nu, n_trials, &
                             TKT_best, A_best, b_best, sigma_A, sigma_b, variance_best)
    !f2py intent(in) :: n_T, T_values, Chi_values, T_min, T_max, nu, n_trials
    !f2py intent(out) :: TKT_best, A_best, b_best, sigma_A, sigma_b, variance_best
    !f2py depend(n_T) :: T_values, Chi_values

    integer, intent(in) :: n_T, n_trials
    real(8), intent(in) :: T_values(n_T), Chi_values(n_T)
    real(8), intent(in) :: T_min, T_max, nu
    real(8), intent(out) :: TKT_best, A_best, b_best
    real(8), intent(out) :: sigma_A, sigma_b, variance_best

    real(8) :: TKT, TKT_min, TKT_max, dTKT
    real(8) :: epsilon, x, y, sum_x, sum_y, sum_xx, sum_xy
    real(8) :: slope, intercept, variance, y_pred, denom
    integer :: i, j, count

    ! --- Range di TKT candidate ---
    TKT_min = 0.7d0
    TKT_max = 1.0d0
    dTKT = (TKT_max - TKT_min) / real(n_trials - 1, 8)

    ! --- Inizializzazione ---
    variance_best = 1.0d30
    TKT_best = 0.89d0
    A_best = 1.0d0
    b_best = 1.0d0
    sigma_A = 0.0d0
    sigma_b = 0.0d0

    ! --- Scansione dei TKT ---
    do j = 1, n_trials
        TKT = TKT_min + real(j-1, 8) * dTKT
        sum_x = 0.0d0
        sum_y = 0.0d0
        sum_xx = 0.0d0
        sum_xy = 0.0d0
        count = 0

        ! --- Selezione punti per fit ---
        do i = 1, n_T
            if (T_values(i) >= T_min .and. T_values(i) <= T_max .and. &
                T_values(i) > TKT .and. Chi_values(i) > 0.0d0) then
                epsilon = (T_values(i) - TKT) / TKT
                if (epsilon > 0.01d0) then
                    x = 1.0d0 / (epsilon ** nu)
                    y = log(Chi_values(i))
                    sum_x = sum_x + x
                    sum_y = sum_y + y
                    sum_xx = sum_xx + x*x
                    sum_xy = sum_xy + x*y
                    count = count + 1
                end if
            end if
        end do

        ! Se meno di 3 punti, salta
        if (count < 3) cycle

        ! --- Fit lineare ---
        denom = count*sum_xx - sum_x*sum_x
        if (denom == 0.0d0) cycle  ! evita divisione per zero
        slope = (count*sum_xy - sum_x*sum_y) / denom
        intercept = (sum_y - slope*sum_x) / real(count,8)

        ! --- Varianza residui ---
        variance = 0.0d0
        do i = 1, n_T
            if (T_values(i) >= T_min .and. T_values(i) <= T_max .and. &
                T_values(i) > TKT .and. Chi_values(i) > 0.0d0) then
                epsilon = (T_values(i) - TKT) / TKT
                if (epsilon > 0.01d0) then
                    x = 1.0d0 / (epsilon ** nu)
                    y = log(Chi_values(i))
                    y_pred = intercept + slope*x
                    variance = variance + (y - y_pred)**2
                end if
            end if
        end do

        if (count > 2) then
            variance = variance / real(count - 2, 8)
        else
            variance = 1.0d30
        end if

        ! --- Salvataggio parametri migliori ---
        if (variance < variance_best) then
            variance_best = variance
            TKT_best = TKT
            A_best = exp(intercept)
            b_best = slope

            ! --- Errori parametri ---
            sigma_b = sqrt(variance / (sum_xx - sum_x*sum_x/count))
            sigma_A = A_best * sqrt(variance * sum_xx / (count*(sum_xx - sum_x*sum_x/count)))
        end if
    end do
  end subroutine fit_kt_transition


  subroutine binning_error(n_data, data, max_bins, n_bins_out, bin_sizes, bin_errors)
    !f2py intent(in) :: n_data, data, max_bins
    !f2py intent(out) :: n_bins_out, bin_sizes, bin_errors
    !f2py depend(n_data) :: data
    !f2py depend(max_bins) :: bin_sizes, bin_errors
    integer, intent(in) :: n_data, max_bins
    real(8), intent(in) :: data(n_data)
    integer, intent(out) :: n_bins_out
    integer, intent(out) :: bin_sizes(max_bins)
    real(8), intent(out) :: bin_errors(max_bins)
    
    integer :: bin_size, n_bins, i
    real(8) :: bin_mean, bin_var
    real(8), allocatable :: bin_means(:)
    
    n_bins_out = 0
    bin_sizes = 0
    bin_errors = 0.0d0
    
    do bin_size = 1, min(max_bins, n_data/2)
      n_bins = n_data / bin_size
      if (n_bins < 2) exit
      
      allocate(bin_means(n_bins))
      
      do i = 1, n_bins
        bin_means(i) = sum(data((i-1)*bin_size+1 : i*bin_size)) / real(bin_size, 8)
      end do
      
      bin_mean = sum(bin_means) / real(n_bins, 8)
      bin_var = sum((bin_means - bin_mean)**2) / real(n_bins - 1, 8)
      
      n_bins_out = n_bins_out + 1
      bin_sizes(n_bins_out) = bin_size
      bin_errors(n_bins_out) = sqrt(bin_var / real(n_bins, 8))
      
      deallocate(bin_means)
    end do
  end subroutine binning_error
  
  subroutine autocorrelation(n_data, data, max_lag, autocorr)
    !f2py intent(in) :: n_data, data, max_lag
    !f2py intent(out) :: autocorr
    !f2py depend(n_data) :: data
    !f2py depend(max_lag) :: autocorr
    integer, intent(in) :: n_data, max_lag
    real(8), intent(in) :: data(n_data)
    real(8), intent(out) :: autocorr(0:max_lag)
    
    integer :: t, tau
    real(8) :: mean_val, var_val, sum_prod
    
    mean_val = sum(data) / real(n_data, 8)
    var_val = sum((data - mean_val)**2) / real(n_data, 8)
    
    autocorr = 0.0d0
    
    if (var_val < 1.0d-15) then
      autocorr(0) = 1.0d0
      return
    end if
    
    do tau = 0, max_lag
      sum_prod = 0.0d0
      do t = 1, n_data - tau
        sum_prod = sum_prod + (data(t) - mean_val) * (data(t+tau) - mean_val)
      end do
      autocorr(tau) = sum_prod / real(n_data - tau, 8) / var_val
    end do
  end subroutine autocorrelation

  subroutine compute_vorticity(theta, vorticity, n)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: theta(n,n)
    integer, intent(inout) :: vorticity(n,n)
    !f2py intent(in) :: theta
    !f2py intent(in,out) :: vorticity
    !f2py intent(hide), depend(theta) :: n = shape(theta,0)
    
    integer :: i, j, ip, jp
    real(8) :: dtheta1, dtheta2, dtheta3, dtheta4
    real(8), parameter :: pi = 3.141592653589793d0

    do i = 1, n
        ip = mod(i, n) + 1
        do j = 1, n
            jp = mod(j, n) + 1

            dtheta1 = theta(ip,j)   - theta(i,j)
            dtheta2 = theta(ip,jp)  - theta(ip,j)
            dtheta3 = theta(i,jp)   - theta(ip,jp)
            dtheta4 = theta(i,j)    - theta(i,jp)

            dtheta1 = dtheta1 - 2*pi*nint(dtheta1/(2*pi))
            dtheta2 = dtheta2 - 2*pi*nint(dtheta2/(2*pi))
            dtheta3 = dtheta3 - 2*pi*nint(dtheta3/(2*pi))
            dtheta4 = dtheta4 - 2*pi*nint(dtheta4/(2*pi))

            vorticity(i,j) = 0
            if (dtheta1+dtheta2+dtheta3+dtheta4 > pi) vorticity(i,j) = 1
            if (dtheta1+dtheta2+dtheta3+dtheta4 < -pi) vorticity(i,j) = -1
        end do
    end do
  end subroutine compute_vorticity

  subroutine temperature_scan_2(L, n_T, T_values, n_equil, n_meas, &
                            E_means, E_errors_bin, E_errors_boot, &
                            M_means, M_errors_bin, M_errors_boot, &
                            M2_means, Chi_values, Cv_values, &
                            vort_means, vort_errors_bin, vort_errors_boot, &
                            theta_final, accept_rates, deltas)
    !f2py intent(in) :: L, n_T, T_values, n_equil, n_meas
    !f2py intent(out) :: E_means, E_errors_bin, E_errors_boot
    !f2py intent(out) :: M_means, M_errors_bin, M_errors_boot
    !f2py intent(out) :: M2_means, Chi_values, Cv_values
    !f2py intent(out) :: vort_means, vort_errors_bin, vort_errors_boot
    !f2py intent(out) :: theta_final, accept_rates, deltas
    !f2py depend(n_T) :: T_values, E_means, E_errors_bin, E_errors_boot
    !f2py depend(n_T) :: M_means, M_errors_bin, M_errors_boot
    !f2py depend(n_T) :: M2_means, Chi_values, Cv_values
    !f2py depend(n_T) :: vort_means, vort_errors_bin, vort_errors_boot
    !f2py depend(L, n_T) :: theta_final

    integer, intent(in) :: L, n_T, n_equil, n_meas
    real(8), intent(in) :: T_values(n_T)
    real(8), intent(out) :: E_means(n_T), E_errors_bin(n_T), E_errors_boot(n_T)
    real(8), intent(out) :: M_means(n_T), M_errors_bin(n_T), M_errors_boot(n_T)
    real(8), intent(out) :: M2_means(n_T)
    real(8), intent(out) :: Chi_values(n_T), Cv_values(n_T)
    real(8), intent(out) :: vort_means(n_T), vort_errors_bin(n_T), vort_errors_boot(n_T)
    real(8), intent(out) :: theta_final(L,L,n_T)
    real(8), intent(out) :: accept_rates(n_T), deltas(n_T)

    ! Serie temporanee
    real(8), allocatable :: theta(:,:)
    real(8), allocatable :: E_series(:), M_series(:), M2_series(:)
    real(8), allocatable :: theta2_series(:), vort_series(:), E2_series(:)
    real(8) :: T, beta, delta, final_delta, final_accept
    real(8) :: E_mean, E2_mean, M_mean, M2_mean, vort_mean
    integer :: i, max_bins, n_bins_out
    integer, allocatable :: bin_sizes(:)
    real(8), allocatable :: bin_errors(:)

    max_bins = 100
    allocate(theta(L,L))
    allocate(E_series(n_meas), M_series(n_meas), M2_series(n_meas))
    allocate(theta2_series(n_meas), vort_series(n_meas), E2_series(n_meas))
    allocate(bin_sizes(max_bins), bin_errors(max_bins))

    call init_cold_config(L, theta)
    delta = 1.0d0

    do i = 1, n_T
        T = T_values(i)
        beta = 1.0d0 / T

        ! Equilibriazione
        call run_equilibration(L, theta, beta, delta, n_equil, 0.4d0, &
                               final_delta, final_accept)
        deltas(i) = final_delta
        accept_rates(i) = final_accept
        delta = final_delta

        ! Misurazioni
        call run_measurement(L, theta, beta, delta, n_meas, &
                             E_series, M_series, M2_series, theta2_series, &
                             vort_series, E2_series)

        ! Medie
        E_mean = sum(E_series) / real(n_meas,8)
        E2_mean = sum(E2_series) / real(n_meas,8)
        M_mean = sum(M_series) / real(n_meas,8)
        M2_mean = sum(M2_series) / real(n_meas,8)
        vort_mean = sum(vort_series) / real(n_meas,8)

        E_means(i) = E_mean
        M_means(i) = M_mean
        M2_means(i) = M2_mean
        vort_means(i) = vort_mean / real(L*L,8)

        ! Calcolo errori con binning
        call binning_error(n_meas, E_series, max_bins, n_bins_out, bin_sizes, bin_errors)
        E_errors_bin(i) = bin_errors(n_bins_out)

        call binning_error(n_meas, M_series, max_bins, n_bins_out, bin_sizes, bin_errors)
        M_errors_bin(i) = bin_errors(n_bins_out)

        call binning_error(n_meas, vort_series/real(L*L,8), max_bins, n_bins_out, bin_sizes, bin_errors)
        vort_errors_bin(i) = bin_errors(n_bins_out)

        ! Bootstrap (semplice, in Fortran)
        E_errors_boot(i) = bootstrap_std(E_series, n_meas, 1000)
        M_errors_boot(i) = bootstrap_std(M_series, n_meas, 1000)
        vort_errors_boot(i) = bootstrap_std(vort_series/real(L*L,8), n_meas, 1000)

        ! Suscettibilità e capacità termica
        Chi_values(i) = real(L*L,8) * M2_mean / T
        Cv_values(i) = real(L*L,8) * (E2_mean - E_mean**2) / (T*T)

        ! Salva configurazione finale
        theta_final(:,:,i) = theta
    end do

    deallocate(theta, E_series, M_series, M2_series)
    deallocate(theta2_series, vort_series, E2_series, bin_sizes, bin_errors)
  end subroutine temperature_scan_2

  ! Funzione interna per bootstrap error
  real(8) function bootstrap_std(data, n_data, n_boot)
      integer, intent(in) :: n_data, n_boot
      real(8), intent(in) :: data(n_data)
      integer :: i, j
      real(8) :: sample_mean
      real(8), allocatable :: boot_means(:)
      allocate(boot_means(n_boot))
      call random_seed()  ! inizializza RNG

      do i = 1, n_boot
          sample_mean = 0.0d0
          do j = 1, n_data
              sample_mean = sample_mean + data(1 + mod(int(n_data*rand()), n_data))
          end do
          boot_means(i) = sample_mean / real(n_data,8)
      end do
      bootstrap_std = sqrt(sum((boot_means - sum(boot_means)/real(n_boot,8))**2) / real(n_boot-1,8))
      deallocate(boot_means)
  end function bootstrap_std

end module xy_model