module parameters

   implicit none

   integer, parameter :: L = 36
   integer, parameter :: N = L*L
   integer, parameter :: stage0 = 20
   integer, parameter :: unit = 6
   double precision, parameter :: eps_min = 5D-3
   double precision, parameter :: T_min = 2.2
   double precision, parameter :: T_max = 2.3

   type lattice
      integer, allocatable :: s(:, :); 
      integer :: E; 
      integer :: M1; 
      integer :: M2; 
      integer :: M3; 
   end type lattice

   type wl_parameters
      integer, allocatable :: H(:); 
      double precision, allocatable :: lng(:); 
      double precision, allocatable :: M_av(:); 
      double precision, allocatable :: M2_av(:); 
      double precision :: lnf; 
   end type wl_parameters

end module parameters
