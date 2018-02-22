C

C     prop_magnus_ord2_step.F
C
C     Steps dens mat forward dt in time (in MO basis) using second order
C     (single term) Magnus propagator.
C
C     P'(t+dt) = e^W P'(t) (e^W)^+ 
C              = e^W P'(t) e^-W  (assuming W is Hermitian, which is not true if we have a complex absorbing potential)
C
C     where W = -i F'(t+dt/2).
C
C     To exponentiate, we can diagonalize, use a power series, or use
C     the Baker-Campbell-Hausdorff expansion to evaluate the whole
C     expression via commutators (note only for Hermitian W).  P', F'
C     and W are all ns_mo x ns_mo.  Noote that dt can be different from
C     the simulation time step!
C
      subroutine prop_magnus_ord4_step (params, dt,
     $    integral_coefs, g_zfock_0, g_zfock_half, g_zfock_end, g_zdens)
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
      integer, intent(in) :: g_zfock_0, g_zfock_half, g_zfock_end
      double complex, intent(in) :: integral_coefs(3)


C     == In/out ==
      integer, intent(in) :: g_zdens           ! in = P'(t), out = P'(t+dt)


C     == Parameters ==
      character(*), parameter   :: pname = "prop_magnus_ord4_step: "


C     == Variables ==
      integer g_zscr1, g_zscr2, g_zscr3  ! scratch
      integer g_w               ! W = Omega_1 = -i F'(t+dt/2)
      integer g_expw            ! e^W
      logical ok
      integer me
      double complex zdt
      double precision spec_range, max_elem


C     (junk)
      integer ivals, lvals
      integer g_vecs
      double complex zval, comm_coef, sixth

      me = ga_nodeid ()

C
C     Allocation.
C
      ok = .true.
      ok = ok .and. ga_duplicate(g_zfock_0, g_w, "ord4w")
      ok = ok .and. ga_duplicate(g_zfock_0, g_expw, "ord4expw")
      ok = ok .and. ga_duplicate(g_zfock_0, g_zscr1, "ord4scr1")
      ok = ok .and. ga_duplicate(g_zfock_0, g_zscr2, "ord4scr2")
      ok = ok .and. ga_duplicate(g_zfock_0, g_zscr3, "ord4scr3")
      if (.not. ok) call errquit (pname//"allocation failed",0,GA_ERR)

      zdt = dcmplx (dt, 0d0)

      call ga_copy(g_zfock_0, g_zscr2)
      call ga_scale(g_zscr2, integral_coefs(1))

      call ga_add(integral_coefs(2), g_zfock_half,
     $     integral_coefs(3), g_zfock_end, g_zscr2)

      if (params%dt == dt) then
      !F1F0
      call ga_zgemm ("N", "N",
     $     params%ns_mo, params%ns_mo, params%ns_mo,
     $     z1, g_zfock_end, g_zfock_0, z0, g_zscr1)

      !F0F1 - F1F0
      call ga_zgemm ("N", "N",
     $     params%ns_mo, params%ns_mo, params%ns_mo,
     $     z1, g_zfock_0, g_zfock_end, zn1, g_zscr1)
      else
         call ga_zgemm ("N", "N",
     $        params%ns_mo, params%ns_mo, params%ns_mo,
     $        z1, g_zfock_half, g_zfock_0, z0, g_zscr1)

!     F0F1 - F1F0
         call ga_zgemm ("N", "N",
     $        params%ns_mo, params%ns_mo, params%ns_mo,
     $        z1, g_zfock_0, g_zfock_half, zn1, g_zscr1)
      endif

      call ga_add (zni*zdt, g_zscr2, -(zdt**2)/12.d0, g_zscr1, g_w)

      call exp_pseries (params, g_w, g_expw)

C
C     Compute eigenvalues of e^W (for diagnostics only)
C
      if (params%checklvl .ge. 2) then
          ok = .true.
          ok = ok .and. ga_duplicate(g_zdens, g_vecs, "vecs")
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
         if (.not. mat_is_hermitian (g_zfock_0, params%tol_zero))
     $        call errquit (pname//"fock0 not hermitian", 0, 0)
         if (.not. mat_is_hermitian (g_zfock_half, params%tol_zero))
     $        call errquit (pname//"fockhalf not hermitian", 0, 0)
         if (.not. mat_is_hermitian (g_zfock_end, params%tol_zero))
     $        call errquit (pname//"fockend not hermitian", 0, 0)
         if (.not. mat_is_unitary (g_expw, 1d-4))
     $        call errquit (pname//"e^W not unitary", 0, 0)
      endif


C
C     P'(t+dt) = e^W P'(t) (e^W)^+.  Recycle g_w as scratch space to
C     hold (e^W)^+.
C
      call ga_zgemm ("N", "C", params%ns_mo, params%ns_mo,
     $     params%ns_mo, z1, g_zdens, g_expw, z0, g_zscr1)

      call ga_zgemm ("N", "N", params%ns_mo, params%ns_mo,
     $     params%ns_mo, z1, g_expw, g_zscr1, z0, g_zdens)

C
C     Clean up.
C
      ok = .true.
      ok = ok .and. ga_destroy(g_expw)
      ok = ok .and. ga_destroy(g_w)
      ok = ok .and. ga_destroy(g_zscr1)
      ok = ok .and. ga_destroy(g_zscr2)
      ok = ok .and. ga_destroy(g_zscr3)

      if (.not. ok)
     $     call errquit (pname//"failed to destroy GAs", 0, GA_ERR)


      end subroutine
c $Id: prop_magnus_ord2_step.F 24344 2013-06-22 06:15:21Z d3y133 $
