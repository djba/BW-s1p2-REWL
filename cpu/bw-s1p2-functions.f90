module bw_routines
   implicit none
contains

   subroutine n_levels(nl)
      use parameters
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

   subroutine calc_energy_mag(latt)
      use parameters
      implicit none
      type(lattice) :: latt
      integer :: i, j, sl
      integer :: a, b, c, d

      latt%E = 0
      latt%M1 = 0
      latt%M2 = 0
      latt%M3 = 0
      do i = 1, L
         do j = 1, L

            call find_sublattice(i, j, sl)

            if (sl == 1) then
               latt%M1 = latt%M1 + latt%s(i, j)
            end if
            if (sl == 2) then
               latt%M2 = latt%M2 + latt%s(i, j)
            end if
            if (sl == 3) then
               latt%M3 = latt%M3 + latt%s(i, j)
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

            latt%E = latt%E - latt%s(i, j)* &
                     (latt%s(a, j)*latt%s(i, c) + latt%s(i, c)*latt%s(b, c) + &
                      latt%s(b, c)*latt%s(b, j) + latt%s(b, j)*latt%s(i, d) + &
                      latt%s(i, d)*latt%s(a, d) + latt%s(a, d)*latt%s(a, j))
         end do
      end do

      latt%E = latt%E/3

   end subroutine calc_energy_mag

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

   subroutine wl_step(latt, wlp)
      use parameters
      implicit none

      type(lattice) :: latt
      type(wl_parameters) :: wlp
      integer :: energy_old, energy_new, e_old, e_new
      integer :: m1_old, m2_old, m3_old, m1_new, m2_new, m3_new
      integer :: index_old, index_new
      integer :: a, b, c, d, s_old, s_new, i, j, step, pos, sl
      double precision :: rn, factor, mag_old, mag_new

      do step = 1, N

         call random_number(rn)
         pos = int(N*rn)
         i = mod(pos, L) + 1
         j = pos/L + 1

         call find_sublattice(i, j, sl)

         s_old = latt%s(i, j)

         s_new = -s_old

         m1_old = latt%M1
         m2_old = latt%M2
         m3_old = latt%M3

         mag_old = sqrt(3.0*(m1_old**2 + m2_old**2 + m3_old**2))/N

         if (sl == 1) then
            m1_new = m1_old + s_new - s_old
            m2_new = m2_old
            m3_new = m3_old
         end if
         if (sl == 2) then
            m1_new = m1_old
            m2_new = m2_old + s_new - s_old
            m3_new = m3_old
         end if
         if (sl == 3) then
            m1_new = m1_old
            m2_new = m2_old
            m3_new = m3_old + s_new - s_old
         end if

         mag_new = sqrt(3.0*(m1_new**2 + m2_new**2 + m3_new**2))/N

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

         e_old = -s_old*(latt%s(a, j)*latt%s(i, c) + latt%s(i, c)*latt%s(b, c) + &
                         latt%s(b, c)*latt%s(b, j) + latt%s(b, j)*latt%s(i, d) + &
                         latt%s(i, d)*latt%s(a, d) + latt%s(a, d)*latt%s(a, j))
         e_new = -s_new*(latt%s(a, j)*latt%s(i, c) + latt%s(i, c)*latt%s(b, c) + &
                         latt%s(b, c)*latt%s(b, j) + latt%s(b, j)*latt%s(i, d) + &
                         latt%s(i, d)*latt%s(a, d) + latt%s(a, d)*latt%s(a, j))

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
            latt%s(i, j) = s_new
            latt%M1 = m1_new
            latt%M2 = m2_new
            latt%M3 = m3_new
            wlp%M_av(index_new) = wlp%M_av(index_new) + mag_new
            wlp%M2_av(index_new) = wlp%M2_av(index_new) + mag_new*mag_new
         else
            wlp%lng(index_old) = wlp%lng(index_old) + wlp%lnf
            wlp%H(index_old) = wlp%H(index_old) + 1
            wlp%M_av(index_old) = wlp%M_av(index_old) + mag_old
            wlp%M2_av(index_old) = wlp%M2_av(index_old) + mag_old*mag_old
         end if

      end do

   end subroutine wl_step

   subroutine flatness_test(flat, H, nl)

      use parameters
      implicit none
      integer :: flat, test, i, nl
      integer :: H(:)
      double precision :: Hmin, Hmax

      Hmin = H(1)
      Hmax = H(1)

      flat = 1
      do i = 1, nl
         if (H(i) < Hmin) then
            Hmin = H(i)
         end if
         if (H(i) > Hmax) then
            Hmax = H(i)
         end if
      end do

      do i = 1, nl
         if ((Hmax - Hmin)/Hmax > 0.2) then
            flat = 0
         end if
      end do

   end subroutine flatness_test

   subroutine calc_parameters(lng, Emean, E2mean, CV, UE, T, nl)

      use parameters
      implicit none

      double precision :: lng(:)
      double precision :: Emean, E2mean, E4mean, CV, T, sum1, sum2, sum4, sum, lng_max
      double precision :: e_per_spin, UE
      integer :: test, i, nl, E

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

   subroutine calc_Tc_CV(lng, Tc_CV, CV_max, nl)

      use parameters
      implicit none

      double precision :: lng(:)
      double precision :: Emean, E2mean, CV, UE, T, CV_max, UE_max, Tc_CV
      integer :: nl

      T = T_min

      CV_max = 0
      UE_max = 0
      Tc_CV = T

      do while (T < T_max)

         call calc_parameters(lng, Emean, E2mean, CV, UE, T, nl)

         if (CV > CV_max) then
            CV_max = CV
            Tc_CV = T
         end if

         T = T + 5E-6
      end do

   end subroutine calc_Tc_CV

end module bw_routines
