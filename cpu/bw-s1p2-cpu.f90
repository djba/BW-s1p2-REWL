program bw
   use parameters
   use bw_routines

   implicit none
   type(lattice) :: latt
   type(wl_parameters) :: wlp
   double precision :: rn, lng0
   double precision :: Emean, E2mean, CV, UE, CV_max, UE_max, T, Tc_CV, delta_Tc, Tc_0
   double precision :: CV_max_0, eps

   integer :: i, j, step, test, flat, stage, finished, finish_test, prepare, nl, energy
   integer :: index
   integer, allocatable :: H_tot(:)

   integer r_time(8)
   integer seed(2)

   real :: t1, t2, etime, time
   real :: elapsed(2)

   call date_and_time(values=r_time)
   seed(1) = r_time(4)*(360000*r_time(5) + 6000*r_time(6) + 100*r_time(7) + r_time(8))
   call random_seed(put=seed)

   allocate (latt%s(L, L))

   do i = 1, L
      do j = 1, L
         latt%s(i, j) = 1
      end do
   end do

   call calc_energy_mag(latt)

   write (unit, '(a,i3)') "L = ", L
   write (unit, '(a,e10.2)') "eps_min = ", eps_min
   write (unit, '(a,f4.2)') "T_min = ", T_min
   write (unit, '(a,f4.2)') "T_max = ", T_max
   call flush (unit)

   call n_levels(nl)

   write (unit, '(a,i)') "nl = ", nl

   prepare = 0

   do while (prepare == 0)
      do i = 1, L
         do j = 1, L
            call random_number(rn)

            if (rn > 0.5) then
               latt%s(i, j) = 1
            else
               latt%s(i, j) = -1
            end if
         end do
      end do

      call calc_energy_mag(latt)

      if (latt%E <= 0) then
         prepare = 1
      end if
   end do

   allocate (wlp%lng(nl))
   allocate (wlp%H(nl))
   allocate (wlp%M_av(nl))
   allocate (wlp%M2_av(nl))
   allocate (H_tot(nl))

   do i = 1, nl
      wlp%lng(i) = 0
      wlp%H(i) = 0
      wlp%M_av(i) = 0
      wlp%M2_av(i) = 0
      H_tot(i) = 0
   end do

   t1 = etime(elapsed)

   wlp%lnf = 1.0

   stage = 1

   finished = 0

   do while (finished == 0)

      flat = 0
      step = 0
      finish_test = 1
      do while (flat == 0)
         call wl_step(latt, wlp)
         step = step + 1

         if (stage >= stage0 .and. step == 1) then
            call calc_Tc_CV(wlp%lng, Tc_CV, CV_max, nl)
            Tc_0 = Tc_CV
            CV_max_0 = CV_max
            write (unit, '(a,f8.6)') "Tc_0 = ", Tc_0
            write (unit, '(a,f12.6)') "CV_max_0 = ", CV_max_0
         end if

         if (mod(step, 10000) == 0) then
            call flatness_test(flat, wlp%H, nl)
            if (stage >= stage0) then
               call calc_Tc_CV(wlp%lng, Tc_CV, CV_max, nl)
               delta_Tc = abs(Tc_CV - Tc_0)
               eps = abs(CV_max - CV_max_0)/CV_max_0
               if (eps >= eps_min) then
                  finish_test = 0
               end if
               write (unit, '(a,f8.6,a,f8.6)') "Tc_CV = ", &
                  Tc_CV, "    delta_Tc = ", delta_Tc
               write (unit, '(a,f12.6,a,f8.6)') "CV_max = ", CV_max, "    eps = ", eps
               call flush (unit)
            end if
         end if
      end do

      write (unit, '(a,i3,a,g)') "Finished stage = ", stage, "    lnf = ", wlp%lnf

      wlp%lnf = wlp%lnf/2.0
      do i = 1, nl
         H_tot(i) = H_tot(i) + wlp%H(i)
         wlp%H(i) = 0
      end do

      write (unit, '(a,i)') "step = ", step
      call flush (unit)

      if (stage >= stage0 + 1 .and. finish_test == 1) then
         finished = 1
      end if

      stage = stage + 1
   end do

   t2 = etime(elapsed)

   time = t2 - t1

   write (unit, '(a,f)') "Time = ", time
   call flush (unit)

   if (mod(L, 3) == 0) then
      lng0 = wlp%lng(1) - log(4.0)
   else
      lng0 = wlp%lng(1)
   end if

   do i = 1, nl
      wlp%lng(i) = wlp%lng(i) - lng0
      wlp%M_av(i) = wlp%M_av(i)/H_tot(i)
      wlp%M2_av(i) = wlp%M2_av(i)/H_tot(i)
   end do

   do i = 1, nl

      if (i == 1) then
         energy = -2*N
      end if
      if (i == 2) then
         energy = -2*N + 12
      end if
      if (i > 2) then
         energy = -2*N + 12 + (i - 2)*4
      end if

      write (unit, '(i,g,g,g)') energy, wlp%lng(i), wlp%M_av(i), wlp%M2_av(i)

   end do

   call calc_Tc_CV(wlp%lng, Tc_CV, CV_max, nl)

   write (unit, '(a,f8.6)') "Tc_CV = ", Tc_CV
   call flush (unit)

   deallocate (latt%s)
   deallocate (wlp%lng)
   deallocate (wlp%H)
   deallocate (H_tot)
   deallocate (wlp%M_av)
   deallocate (wlp%M2_av)

end program bw
