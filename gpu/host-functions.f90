module host_functions

  use parameters

contains

  ! compute the number of energy levels
  subroutine n_levels(nl)

    implicit none
    integer :: nl, E

    nl = 0
    E = -2*N

    do while (E <= 0)
       nl = nl + 1
       if (E == -2*N) then
          E = E + 12
       else
          E = E + 4
       end if
    end do

  end subroutine n_levels

  ! compute the boundaries for each energy window
  subroutine prepare_windows(E_min, E_max, nw)

    implicit none

    integer :: i, nw
    integer :: E_min(:), E_max(:)

    do i = 1, nw
       E_min(i) = -2*N + (i - 1)*(deltaE/2)
       E_max(i) = -2*N + deltaE + (i - 1)*(deltaE/2)
    end do

  end subroutine prepare_windows

  ! find the order of sublattice for the spin with indexes (i, j)
  subroutine find_sublattice(i, j, sl)

    implicit none
    integer :: i, j, sl

    if (mod(i, 3) == 0 .and. mod(j, 3) == 0) then
       sl = 1
    end if
    if (mod(i, 3) == 1 .and. mod(j, 3) == 1) then
       sl = 1
    end if
    if (mod(i, 3) == 2 .and. mod(j, 3) == 2) then
       sl = 1
    end if

    if (mod(i, 3) == 2 .and. mod(j, 3) == 1) then
       sl = 2
    end if
    if (mod(i, 3) == 0 .and. mod(j, 3) == 2) then
       sl = 2
    end if
    if (mod(i, 3) == 1 .and. mod(j, 3) == 0) then
       sl = 2
    end if

    if (mod(i, 3) == 0 .and. mod(j, 3) == 1) then
       sl = 3
    end if
    if (mod(i, 3) == 1 .and. mod(j, 3) == 2) then
       sl = 3
    end if
    if (mod(i, 3) == 2 .and. mod(j, 3) == 0) then
       sl = 3
    end if

  end subroutine find_sublattice

  ! compute the energy and magnetization
  subroutine calc_energy_mag(latt)

    use parameters
    implicit none
    type(lattice_init) :: latt
    integer :: i, j, pos
    integer :: a, b, c, d, sl

    latt%E = 0
    latt%M1 = 0
    latt%M2 = 0
    latt%M3 = 0

    do i = 1, L
       do j = 1, L
          pos = (i - 1)*L + j

          call find_sublattice(i, j, sl)

          if (sl == 1) then
             latt%M1 = latt%M1 + latt%s(pos)
          end if
          if (sl == 2) then
             latt%M2 = latt%M2 + latt%s(pos)
          end if
          if (sl == 3) then
             latt%M3 = latt%M3 + latt%s(pos)
          end if

          if (i == 1) then
             a = L
          else
             a = i - 1
          end if
          if (i == L) then
             b = 1
          else
             b = i + 1
          end if
          if (j == 1) then
             c = L
          else
             c = j - 1
          end if
          if (j == L) then
             d = 1
          else
             d = j + 1
          end if

          latt%E = latt%E - latt%s((i - 1)*L + j)* &
               (latt%s((a - 1)*L + j)*latt%s((i - 1)*L + c) + &
               latt%s((i - 1)*L + c)*latt%s((b - 1)*L + c) + &
               latt%s((b - 1)*L + c)*latt%s((b - 1)*L + j) + &
               latt%s((b - 1)*L + j)*latt%s((i - 1)*L + d) + &
               latt%s((i - 1)*L + d)*latt%s((a - 1)*L + d) + &
               latt%s((a - 1)*L + d)*latt%s((a - 1)*L + j))
       end do
    end do

    latt%E = latt%E/3

  end subroutine calc_energy_mag


  ! test if energy value is acceptable
  subroutine verify_energy(energy, test)

    use parameters
    implicit none

    integer :: test, energy, i

    test = 0
    if (energy == -2*N) test = 1
    if (energy == -2*N + 12) test = 1
    do i = 1, N/2 - 3
       if (energy == -2*N + 12 + i*4) then
          test = 1
       end if
    end do

  end subroutine verify_energy


  ! one step for the Wang-Landau algorithm
  subroutine wl_step(latt, wlp)

    use parameters
    implicit none
    type(lattice_init) :: latt
    type(wl_parameters_init) :: wlp

    integer :: index_old, index_new
    integer :: energy_old, energy_new, e_old, e_new
    integer :: a, b, c, d, s_old, s_new, i, j, step, pos, mag_new
    integer :: m1_old, m2_old, m3_old, m1_new, m2_new, m3_new, sl
    double precision :: rn, factor

    call random_number(rn)
    pos = int(N*rn)
    i = mod(pos, L) + 1
    j = pos/L + 1

    s_old = latt%s((i - 1)*L + j)
    s_new = -s_old

    call find_sublattice(i, j, sl)

    m1_old = latt%M1
    m2_old = latt%M2
    m3_old = latt%M3

    if (sl == 1) then
       m1_new = latt%M1 + s_new - s_old
       m2_new = latt%M2
       m3_new = latt%M3
    end if
    if (sl == 2) then
       m1_new = latt%M1
       m2_new = latt%M2 + s_new - s_old
       m3_new = latt%M3
    end if
    if (sl == 3) then
       m1_new = latt%M1
       m2_new = latt%M2
       m3_new = latt%M3 + s_new - s_old
    end if

    if (i == 1) then
       a = L
    else
       a = i - 1
    end if
    if (i == L) then
       b = 1
    else
       b = i + 1
    end if
    if (j == 1) then
       c = L
    else
       c = j - 1
    end if
    if (j == L) then
       d = 1
    else
       d = j + 1
    end if

    e_old = -s_old*(latt%s((a - 1)*L + j)*latt%s((i - 1)*L + c) + &
         latt%s((i - 1)*L + c)*latt%s((b - 1)*L + c) + &
         latt%s((b - 1)*L + c)*latt%s((b - 1)*L + j) + &
         latt%s((b - 1)*L + j)*latt%s((i - 1)*L + d) + &
         latt%s((i - 1)*L + d)*latt%s((a - 1)*L + d) + &
         latt%s((a - 1)*L + d)*latt%s((a - 1)*L + j))
    e_new = -s_new*(latt%s((a - 1)*L + j)*latt%s((i - 1)*L + c) + &
         latt%s((i - 1)*L + c)*latt%s((b - 1)*L + c) + &
         latt%s((b - 1)*L + c)*latt%s((b - 1)*L + j) + &
         latt%s((b - 1)*L + j)*latt%s((i - 1)*L + d) + &
         latt%s((i - 1)*L + d)*latt%s((a - 1)*L + d) + &
         latt%s((a - 1)*L + d)*latt%s((a - 1)*L + j))

    energy_old = latt%E
    energy_new = latt%E + (e_new - e_old)

    if (energy_old == -2*N) then
       index_old = 1
    end if
    if (energy_old == -2*N + 12) then
       index_old = 2
    end if
    if (energy_old > -2*N + 12) then
       index_old = (energy_old + 2*N - 12)/4 + 2
    end if

    if (energy_new == -2*N) then
       index_new = 1
    end if
    if (energy_new == -2*N + 12) then
       index_new = 2
    end if
    if (energy_new > -2*N + 12) then
       index_new = (energy_new + 2*N - 12)/4 + 2
    end if

    factor = 0
    if (energy_new <= 0) then
       factor = exp(wlp%lng(index_old) - wlp%lng(index_new))
    end if

    call random_number(rn)

    if (factor > rn) then
       wlp%lng(index_new) = wlp%lng(index_new) + wlp%lnf
       wlp%H(index_new) = wlp%H(index_new) + 1
       latt%E = energy_new
       latt%s((i - 1)*L + j) = s_new
       latt%M1 = m1_new
       latt%M2 = m2_new
       latt%M3 = m3_new
    else
       wlp%lng(index_old) = wlp%lng(index_old) + wlp%lnf
       wlp%H(index_old) = wlp%H(index_old) + 1
    end if

  end subroutine wl_step

  ! Histogram flatness test
  subroutine flatness_test(H, test, nw, E_min)

    use parameters
    implicit none

    integer :: H(:), E_min(:)
    integer :: test, i, j, nw, e_test, energy, index
    double precision :: Hmin, Hmax

    test = 1
    do i = 1, nw*NR
       Hmin = H((i - 1)*(deltaE/4 + 1) + 1)
       Hmax = H((i - 1)*(deltaE/4 + 1) + j)
       do j = 1, deltaE/4 + 1

          energy = E_min((i - 1)/NR + 1) + 4*(j - 1)

          index = (i - 1)*(deltaE/4 + 1) + j
          e_test = 1

          if (i <= NR .and. j == 2) then
             e_test = 0
          end if
          if (i <= NR .and. j == 3) then
             e_test = 0
          end if

          if (e_test == 1 .and. H((i - 1)*(deltaE/4 + 1) + j) < Hmin) then
             Hmin = H((i - 1)*(deltaE/4 + 1) + j)
          end if
          if (e_test == 1 .and. H((i - 1)*(deltaE/4 + 1) + j) > Hmax) then
             Hmax = H((i - 1)*(deltaE/4 + 1) + j)
          end if
       end do

       do j = 1, deltaE/4 + 1
          energy = E_min((i - 1)/NR + 1) + 4*(j - 1)

          index = (i - 1)*(deltaE/4 + 1) + j
          e_test = 1

          if (i <= NR .and. j == 2) then
             e_test = 0
          end if
          if (i <= NR .and. j == 3) then
             e_test = 0
          end if

          if (e_test == 1 .and. (Hmax - Hmin)/Hmax > 0.2) then
             test = 0
          end if
       end do
    end do
  end subroutine flatness_test


  ! Calculate average values for the logarithm of g(E) estimators
  subroutine average_lng(lng, lng_av, nw, av)

    use parameters
    implicit none

    double precision :: lng(:), lng_av(:)
    integer :: nw, i, j, k, av

    do i = 1, nw
       do j = 1, deltaE/4 + 1
          lng_av((i - 1)*(deltaE/4 + 1) + j) = 0
          do k = 1, NR
             lng_av((i - 1)*(deltaE/4 + 1) + j) = &
                  lng_av((i - 1)*(deltaE/4 + 1) + j) + &
                  lng((i - 1)*NR*(deltaE/4 + 1) + (k - 1)*(deltaE/4 + 1) + j)
          end do
          lng_av((i - 1)*(deltaE/4 + 1) + j) = lng_av((i - 1)*(deltaE/4 + 1) + j)/NR
       end do
    end do

    if (av == 1) then
       do i = 1, nw
          do j = 1, deltaE/4 + 1
             do k = 1, NR
                lng((i - 1)*NR*(deltaE/4 + 1) + (k - 1)*(deltaE/4 + 1) + j) = &
                     lng_av((i - 1)*(deltaE/4 + 1) + j)
             end do
          end do
       end do
    end if

  end subroutine average_lng


  ! Compute the final values for the logarithm of g(E) estimators
  subroutine calc_lng(lng_av, lng_final, nw, E_min, nl)

    use parameters
    implicit none

    integer :: nw, i, j, E, index, nl

    double precision :: lng_av(:), lng_final(:)
    integer :: E_min(:)
    double precision :: delta_g, lng0

    do i = 1, nw - 1
       delta_g = lng_av(i*(deltaE/4 + 1) + 1) - &
            lng_av((i - 1)*(deltaE/4 + 1) + deltaE/8 + 1)
       do j = 1, deltaE/4 + 1
          lng_av(i*(deltaE/4 + 1) + j) = lng_av(i*(deltaE/4 + 1) + j) - delta_g
       end do
    end do

    do i = 1, nw - 1
       do j = 1, deltaE/8
          E = E_min(i) + 4*(j - 1)
          index = 0
          if (E == -2*N) then
             index = 1
          end if
          if (E == -2*N + 12) then
             index = 2
          end if
          if (E > -2*N + 12 .and. mod(E + 2*N, 4) == 0) then
             index = (E + 2*N - 12)/4 + 2
          end if

          if (index /= 0) then
             lng_final(index) = lng_av((i - 1)*(deltaE/4 + 1) + j)
          end if
       end do
    end do

    do j = 1, deltaE/4 + 1
       E = E_min(nw) + 4*j - 1
       index = (E + 2*N - 12)/4 + 2

       lng_final(index) = lng_av((nw - 1)*(deltaE/4 + 1) + j)

    end do

    if (mod(L, 3) == 0) then
       lng0 = lng_final(1) - log(4.0)
    else
       lng0 = lng_final(1)
    end if

    do i = 1, nl
       lng_final(i) = lng_final(i) - lng0
    end do

  end subroutine calc_lng

  ! Compute thermodynamic functions
  subroutine calc_parameters(lng, Emean, E2mean, CV, UE, T, nl)

    use parameters
    implicit none

    double precision :: lng(:)
    double precision :: Emean, E2mean, E4mean, CV, T, sum1, sum2, sum4, sum, lng_max
    double precision :: e_per_spin, UE
    integer :: test, i, E, nl

    lng_max = lng(1)

    do i = 1, nl
       if (i == 1) then
          E = -2*N
       end if
       if (i == 2) then
          E = -2*N + 12
       end if
       if (i > 2) then
          E = -2*N + 12 + (i - 2)*4
       end if

       if (lng(i) - dble(E)/T > lng_max) then
          lng_max = lng(i) - dble(E)/T
       end if

    end do

    sum = 0
    sum1 = 0
    sum2 = 0
    sum4 = 0

    do i = 1, nl

       if (i == 1) then
          E = -2*N
       end if
       if (i == 2) then
          E = -2*N + 12
       end if
       if (i > 2) then
          E = -2*N + 12 + (i - 2)*4
       end if

       e_per_spin = dble(E)/N
       sum = sum + exp(lng(i) - dble(E)/T - lng_max)
       sum1 = sum1 + e_per_spin*exp(lng(i) - dble(E)/T - lng_max)
       sum2 = sum2 + e_per_spin*e_per_spin*exp(lng(i) - dble(E)/T - lng_max)
       sum4 = sum4 + e_per_spin*e_per_spin*e_per_spin*e_per_spin* &
            exp(lng(i) - dble(E)/T - lng_max)

    end do

    Emean = sum1/sum
    E2mean = sum2/sum
    E4mean = sum4/sum

    CV = N*(E2mean - Emean*Emean)/(T*T)

    UE = (E4mean)/(E2mean*E2mean)

  end subroutine calc_parameters

  ! Compute critical temperature and maximum value of heat capacity
  subroutine calc_Tc_CV(lng, Tc_CV, CV_max, nl)

    use parameters
    implicit none

    double precision :: lng(:)
    double precision :: Emean, E2mean, UE, T, CV, Tc_CV, CV_max, UE_max
    integer :: nl

    T = T_min

    CV_max = 0
    Tc_CV = T
    UE_max = 0

    do while (T < T_max)

       call calc_parameters(lng, Emean, E2mean, CV, UE, T, nl)

       if (CV > CV_max) then
          CV_max = CV
          Tc_CV = T
       end if

       T = T + 5E-6
    end do

  end subroutine calc_Tc_CV

  subroutine replace_spins(s, idx1, idx2)

    implicit none
    integer :: s(:)

    integer :: i, s1, s2, temp, idx1, idx2

    do i = 1, N
       s1 = s((idx1 - 1)*N + i)
       s2 = s((idx2 - 1)*N + i)
       s((idx1 - 1)*N + i) = s2
       s((idx2 - 1)*N + i) = s1
    end do

  end subroutine replace_spins

  ! Replica exchange phase
  ! We compute the number of replica exchange attempts
  ! and the number of accepted exchanges
  subroutine replica_exchange(latt, lng, E_min, E_max, nw, count_tot, count_accept)

    implicit none

    type(lattice) :: latt
    integer :: idx1, idx2, i, j, energy1, energy2, nw, m1, m2
    integer :: m1_1, m1_2, m2_1, m2_2, m3_1, m3_2
    integer :: E_min(:), E_max(:)
    double precision lng(:)
    double precision :: lng1_old, lng1_new, lng2_old, lng2_new, rn
    integer :: index1_old, index1_new, index2_old, index2_new
    integer :: count_tot, count_accept

    do i = 1, nw - 1
       do j = 1, NR
          idx1 = (i - 1)*NR + j
          idx2 = idx1 + NR
          energy1 = latt%E(idx1)
          energy2 = latt%E(idx2)
          m1_1 = latt%M1(idx1)
          m1_2 = latt%M1(idx2)
          m2_1 = latt%M2(idx1)
          m2_2 = latt%M2(idx2)
          m3_1 = latt%M3(idx1)
          m3_2 = latt%M3(idx2)

          count_tot = count_tot + 1

          if (energy1 >= E_min(i + 1) .and. energy2 <= E_max(i)) then
             index1_old = (energy1 - E_min(i))/4 + 1 + (idx1 - 1)*(deltaE/4 + 1)
             index2_old = (energy2 - E_min(i + 1))/4 + 1 + (idx2 - 1)*(deltaE/4 + 1)
             index1_new = (energy2 - E_min(i))/4 + 1 + (idx1 - 1)*(deltaE/4 + 1)
             index2_new = (energy1 - E_min(i + 1))/4 + 1 + (idx2 - 1)*(deltaE/4 + 1)
             lng1_old = lng(index1_old)
             lng2_old = lng(index2_old)
             lng1_new = lng(index1_new)
             lng2_new = lng(index2_new)

             call random_number(rn)

             if (exp(lng1_old + lng2_old - lng1_new - lng2_new) > rn) then

                call replace_spins(latt%s, idx1, idx2)
                latt%E(idx1) = energy2
                latt%E(idx2) = energy1
                latt%M1(idx1) = m1_2
                latt%M1(idx2) = m1_1
                latt%M2(idx1) = m2_2
                latt%M2(idx2) = m2_1
                latt%M3(idx1) = m3_2
                latt%M3(idx2) = m3_1
                count_accept = count_accept + 1
             end if
          end if

       end do

    end do

  end subroutine replica_exchange

end module host_functions
