C

C     prop_magnus_ord4_step.F
C
C     Steps dens mat forward dt in time (in MO basis) using fourth order
C     (single term) Magnus propagator.
C
C     P'(t+dt) = e^W P'(t) (e^W)^+ 
C              = e^W P'(t) e^-W  (assuming W is Hermitian, which is not true if we have a complex absorbing potential)
C
C     where W = -i (Dt/6(F'(t) + 4 F'(t+dt/2)  F'(t+dt)) + h^2/12[F'(t),F'(t+dt])
C
C     To exponentiate, we can diagonalize, use a power series, or use
C     the Baker-Campbell-Hausdorff expansion to evaluate the whole
C     expression via commutators (note only for Hermitian W).  P', F'
C     and W are all ns_mo x ns_mo.  Note that dt can be different from
C     the simulation time step!
C
      subroutine prop_magnus_ord4_step (params, dt,
     $     i_coef, g_zfock, g_zdens)
c      subroutine prop_magnus_ord4_step (params, dt,
c     $     i_coef, g_zfock_0, g_zfock_half, g_zfock_end,
c     $      g_zdens_half, g_zdens_end)
      implicit none

#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "matutils.fh"
#include "rtdb.fh"
#include "rt_tddft.fh"
C     == Inputs ==
      type (rt_params_t), intent(in) :: params
      double precision, intent(in)   :: dt              ! time step dt
!      integer, intent(in) :: g_zfock_0        !  F'(t_n)
!      integer, intent(in) :: g_zfock_half     !  F'(t_n+dt/2)
      integer, intent(in) :: g_zfock(3) !  F'(t_n+dt)
      integer, intent(in) :: g_zdens(3)      !  F'(t_n+dt)      
      double complex, intent(in) :: i_coef(3,2) !  Coefficients for first integral

C     == Parameters ==
      character(*), parameter   :: pname = "prop_magnus_ord4_step: "

C     == Arrays ==
      integer g_zscr(3)
      
C     == Variables ==
      integer g_w               ! W = Omega_1 = -i F'(t+dt/2)
      integer g_expw            ! e^W
      logical ok
      integer me
      double complex zdt, zc1,zc2
      double precision spec_range, max_elem


C     (junk)
      integer ivals, lvals
      integer g_vecs,docom
      double complex zval, comm_coef, sixth

      me = ga_nodeid ()

C
C     Allocation.
C
      ok = .true.
      ok = ok .and. ga_duplicate(g_zfock(1), g_w, "ord4w")
      ok = ok .and. ga_duplicate(g_zfock(1), g_expw, "ord4expw")
      ok = ok .and. ga_duplicate(g_zfock(1), g_zscr(1), "ord4scr1")
      ok = ok .and. ga_duplicate(g_zfock(1), g_zscr(2), "ord4scr2")
      ok = ok .and. ga_duplicate(g_zfock(1), g_zscr(3), "ord4scr3")
      if (.not. ok) call errquit (pname//"allocation failed",0,GA_ERR)

      zdt = dcmplx (dt, 0.0d0)
      zc1 = dcmplx(dt*dt/12.0d0,0.0d0)
      zc2 = dcmplx(dt*dt/48.0d0,0.0d0)      
      !  Get the linear combination for the integral of the first term
!      call ga_add(i_coef(1,1), g_zfock(1),
!     $     i_coef(2,1), g_zfock(2), g_zscr(1))
      call ga_add(i_coef(1,1), g_zfock(1),
     $     i_coef(2,1), g_zfock(2), g_zscr(1))
      call ga_add(z1, g_zscr(1),i_coef(3,1),g_zfock(3),g_zscr(1))
      
      call ga_add(i_coef(1,2), g_zfock(1),
     $     i_coef(2,2), g_zfock(2), g_zscr(2))
      call ga_add(z1, g_zscr(2),i_coef(3,2), g_zfock(3),g_zscr(2))

      docom=1
      if (docom .eq. 1) then
      !  Compute the commutator terms
      !F1F0
         call ga_zgemm ("N", "N",
     $        params%ns_mo, params%ns_mo, params%ns_mo,
     $        z1, g_zfock(3), g_zfock(1), z0, g_zscr(3))

c        F0F1 - F1F0
         call ga_zgemm ("N", "N",
     $        params%ns_mo, params%ns_mo, params%ns_mo,
     $        z1, g_zfock(1), g_zfock(3), zn1, g_zscr(3))
         call ga_add(zni*zdt, g_zscr(2),zc1, g_zscr(3), g_zscr(2))

!         F0Fhalf - FhalfF0
c$$$         call ga_zgemm ("N", "N",
c$$$     $        params%ns_mo, params%ns_mo, params%ns_mo,
c$$$     $        z1, g_zfock(2), g_zfock(1), z0, g_zscr(3))
c$$$
c$$$         call ga_zgemm ("N", "N",
c$$$     $     params%ns_mo, params%ns_mo, params%ns_mo,
c$$$     $     z1, g_zfock(1), g_zfock(2), zn1, g_zscr(3))

      call ga_add(zni*zdt,g_zscr(1),zc2, g_zscr(3),g_zscr(1))      
      else
         call ga_scale (g_zscr(1), zni*zdt)
         call ga_scale (g_zscr(2), zni*zdt)         
      end if
      

!     
      call exp_pseries (params, g_zscr(2), g_expw)

C
C     Compute eigenvalues of e^W (for diagnostics only)
C
      if (params%checklvl .ge. 2) then
          ok = .true.
          ok = ok .and. ga_duplicate(g_zdens(3), g_vecs, "vecs")
          ok = ok .and. ma_push_get(mt_dcpl, params%ns_mo,
     $        "evals", lvals, ivals)
          if (.not. ok) call errquit (pname//"alloc failed", 0, 0)

          call zmat_diag (g_expw, g_vecs, dcpl_mb(ivals),
     $         .false., .true.)
          call zmat_diag_check (g_expw, g_vecs,
     $            dcpl_mb(ivals), params%tol_zero)
          call zmat_diag_print_eigeninfo (g_expw, g_vecs,
     $         dcpl_mb(ivals))
          call zmat_diag_print_eigenvals (g_expw, g_vecs,
     $         dcpl_mb(ivals))

          ok = .true.
          ok = ok .and. ga_destroy (g_vecs)
          ok = ok .and. ma_pop_stack (lvals)
          if (.not. ok) call errquit (pname//"free failed", 0, 0)

C
C     Since F is Hermitian, W = -iF is skew Hermitian, and thus e^W must
C     be unitary.
C
         if (.not. mat_is_hermitian (g_zfock(1), params%tol_zero))
     $        call errquit (pname//"fock0 not hermitian", 0, 0)
         if (.not. mat_is_hermitian (g_zfock(2), params%tol_zero))
     $        call errquit (pname//"fockhalf not hermitian", 0, 0)
         if (.not. mat_is_hermitian (g_zfock(3), params%tol_zero))
     $        call errquit (pname//"fockend not hermitian", 0, 0)
         if (.not. mat_is_unitary (g_expw, 1d-4))
     $        call errquit (pname//"e^W not unitary", 0, 0)
      endif


C
C     P'(t+dt) = e^W P'(t) (e^W)^+.  Recycle g_w as scratch space to
C     hold (e^W)^+.
C
C     g_zscr(2) = P'(t) (e^W)^+   Note the "C" mean conjugate transpose
      if (.not. mat_is_unitary (g_expw, 1d-10))
     $        call errquit (pname//"e^W not unitary", 0, 0)
      
      call ga_zero (g_zscr(2))
      call ga_zgemm ("N", "C", params%ns_mo, params%ns_mo,
     $     params%ns_mo, z1, g_zdens(3), g_expw, z0, g_zscr(2))

C     zdens_end =  (e^W)g_zscr(2)   
      call ga_zgemm ("N", "N", params%ns_mo, params%ns_mo,
     $     params%ns_mo, z1, g_expw, g_zscr(2), z0, g_zdens(3))
  
      call ga_zero (g_expw)      
      call exp_pseries (params, g_zscr(1), g_expw)
      if (.not. mat_is_unitary (g_expw, 1d-10))
     $     call errquit (pname//"e^W not unitary", 0, 0)

      call ga_zero (g_zscr(1))      
      call ga_zgemm ("N", "C", params%ns_mo, params%ns_mo,
     $     params%ns_mo, z1, g_zdens(2), g_expw, z0, g_zscr(1))

C     zdens_end =  (e^W)g_zscr(2)   
      call ga_zgemm ("N", "N", params%ns_mo, params%ns_mo,
     $     params%ns_mo, z1, g_expw, g_zscr(1), z0, g_zdens(2))

C
C     Clean up.
C
      ok = .true.
      ok = ok .and. ga_destroy(g_expw)
      if (.not. ok)
     $     call errquit (pname//"failed to destroy GAs 1", 0, GA_ERR)
      
      ok = ok .and. ga_destroy(g_w)
      if (.not. ok)
     $     call errquit (pname//"failed to destroy GAs 2", 0, GA_ERR)
      
      ok = ok .and. ga_destroy(g_zscr(1))
            if (.not. ok)
     $     call errquit (pname//"failed to destroy GAs 3", 0, GA_ERR)

            ok = ok .and. ga_destroy(g_zscr(2))
      if (.not. ok)
     $     call errquit (pname//"failed to destroy GAs 4", 0, GA_ERR)
            
      ok = ok .and. ga_destroy(g_zscr(3))

      if (.not. ok)
     $     call errquit (pname//"failed to destroy GAs 5", 0, GA_ERR)

      end subroutine

