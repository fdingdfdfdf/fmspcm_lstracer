module lscale_cond_mod

use simple_sat_vapor_pres_mod, only: escomp, descomp

use constants_mod, only: grav, rdgas, rvgas, kappa, latent_heat, &
        cp_air, cp_vapor, cp_water

implicit none
private

public :: lscale_cond

real(8), public, parameter :: d622 = rdgas/rvgas
real(8), public, parameter :: d378 = 1.-d622
real(8), public, parameter :: d608 = d378/d622
!##==============================================
contains

subroutine compute_geopotential(t, q, p_half, p_full, z_full, z_half)

real(8), intent(in), dimension(:) :: t, q, p_full
real(8), intent(in), dimension(:) :: p_half
real(8), intent(out), dimension(:) :: z_full
real(8), intent(out), dimension(:) :: z_half

integer :: n, k
real(8) :: virtual_t

n = size(t)
z_half(n+1) = 0.
    !virtual_t = t * (1. + d608 * q)
    !for k in range(len(z_full) - 1, 0, -1):
    do k = n, 2, -1
        virtual_t = t(k) * (1. + d608 * q(k)) 
        z_half(k) = z_half(k+1) + rdgas* virtual_t *&
                    (log(p_half(k+1) / p_half(k)))

        z_full(k) = z_half(k+1) + rdgas* virtual_t *&
                    (log(p_half(k+1) / p_full(k)))
    end do

    virtual_t = t(1) * (1. + d608 * q(1)) 
    z_full(1) = z_half(2) + rdgas* virtual_t *&
                (log(p_half(2) / p_full(1)))
    z_half(1) = 0.
    !return z_full, z_half
end subroutine compute_geopotential
!##==============================================
!hc = 1.
!do_evap = False

subroutine qsat(T,p, qs)
 
real(8), intent(in) :: T, p
real(8), intent(out) :: qs
real(8) :: esat
    call escomp(T, esat)
    qs = 1.
    if (esat .le. p) qs = d622 * esat / (p - d378 * esat)
    !return qs
end subroutine qsat

subroutine dqsatdT(T,p,dqs)

real(8), intent(in) :: T, p
real(8), intent(out):: dqs
real(8) :: esat, desat
    call escomp(T, esat)
    call descomp(T, desat)
    dqs = 1.e6
    if ( esat .le. p) dqs = d622 * p / (p - d378 * esat)**2 * desat
end subroutine dqsatdT

!#considering liquid in the air
subroutine rsat(T,p, rs)

real(8), intent(in) :: T, p
real(8), intent(out):: rs
real(8) :: esat
    call escomp(T, esat)
    if (esat .le. p) rs = d622 * esat / (p -  esat)
end subroutine rsat

subroutine drsatdT(T,p, drs)

real(8), intent(in) :: T, p
real(8), intent(out):: drs
real(8) :: esat, desat
    call escomp(T, esat)
    call descomp(T, desat)
    if ( esat .le. p) drs = d622 * p / (p - esat)**2 * desat
end subroutine drsatdT

subroutine ff(rt, T, p, x, val, val_deriv)

real(8), intent(in) :: rt, T, p, x
real(8), intent(out):: val, val_deriv
real(8) :: cp, T1, l, l1, rs, drs
cp = cp_air + cp_water * rt
T1 = T + x
call latent_heat(T1, l1)
call latent_heat(T, l)
call rsat(T1, p, rs)
call drsatdT(T1, p, drs)
val_deriv = cp  + (cp_vapor - cp_water)*rs + l1* drs
val = cp * x + l1*rs - l*rt
end subroutine ff

subroutine newtSolve_ff(rt, T, p, xi, solu)

real(8), intent(in) :: rt, T, p, xi
real(8), intent(out):: solu
real(8) :: x, val, x1, val1, x2, val2, deriv
real(8) :: epss = 1.e-6
integer :: nmax = 100, n 

    x = xi
    val = 10.
    n = 0
    do while (abs(val) .gt. 1.e-6)
       call ff(rt, T, p, x, val, deriv)
       
       !x1 = x + epss
       !call ff(rt, T, p, x1, val1, solu)
       !x2 = x - epss
       !call ff(rt, T, p, x2, val2, solu) 
       !deriv = (val1 - val2) / 2./ epss
       
       x = x - val / deriv  
       
       !print *, val, deriv, x
       
       n = n+1
       if (n > nmax) then
          print *, 'does not converge'
          exit
       end if  
    end do
    solu = x
    !print *, val, deriv, x
end subroutine newtSolve_ff


subroutine compute_k(tin,qin,phalf, k)

real(8), intent(in), dimension(:) :: tin, qin
real(8), intent(in), dimension(:) :: phalf
real(8), intent(out) :: k
real(8) :: kpum
integer :: i
    k = 0.
    do i = 1, size(tin)
       kpum = (cp_air*(1 - qin(i)) + cp_vapor*qin(i))*tin(i)
       k = k + kpum * (phalf(i+1) - phalf(i)) / grav  
    end do
    !return sum(k* (phalf[1:] - phalf[:-1])/cons.grav)
end subroutine compute_k

!#initial profile
subroutine dlnpdlnT_initial(lnT,lnp, slope)

real(8), intent(in) :: lnT, lnp
real(8), intent(out):: slope
real(8) :: T, p, qs, D, rs, pa, es, A, B, C, ll

    T = exp(lnT)
    p = exp(lnp)
    
    call qsat(T,p, qs)
    call latent_heat(T, ll)
    D = ll / (rvgas * T)

    if (qs .eq. 1) then
        slope = D
    else
        rs = qs / (1. - qs)
        pa = p / (rs /d622 + 1.)
        es = p - pa
        A = ll * rs / rdgas / T
        B = rs * (cp_vapor / cp_air + &
              (ll / rvgas / T-1.) * ll / cp_air / T)
        C = 1. / (kappa * (1. + A) / (1. + B)) !#dlnpa / dlnT
        
        slope = (pa * C + es * D) / p
    end if
end subroutine dlnpdlnT_initial


subroutine rk4_f(x0, y0, dx, nx, x, y)

real(8), intent(in) :: x0, y0, dx
integer, intent(in) :: nx
real(8), intent(out), dimension(nx+1) :: x, y
integer :: n
real(8) :: k1, k2, k3, k4

x(1) = x0
y(1) = y0

do n = 1, nx
   call dlnpdlnT_initial(x(n), y(n), k1)
   call dlnpdlnT_initial(x(n) + dx/2., y(n) + dx/2.*k1, k2)
   call dlnpdlnT_initial(x(n) + dx/2., y(n) + dx/2.*k2, k3)
   call dlnpdlnT_initial(x(n) + dx, y(n) + dx *k3, k4)
   x(n+1) = x(n) + dx
   y(n+1) = y(n) + dx/6.*(k1+k4 + 2.*(k2+k3))
end do
end subroutine rk4_f

subroutine warm_start(ts, pd, nx, psg, tt, lnpp, qs)
real(8), intent(in) :: ts, pd
integer, intent(in) :: nx
real(8), intent(out):: psg
real(8), intent(out), dimension(nx+1) :: tt, lnpp, qs
integer :: i
real(8) :: psc, x1, x2, dx
real(8), dimension(nx+1) :: lntlist, lnplist, lntt, pp

    call escomp(ts, psc)
    psg = pd + psc !sat.escomp(ts)

    x1 = log(ts)
    x2 = log(20.)
    !nx = 400
    dx = (x2 - x1) /nx
    call rk4_f(x1, log(psg), dx, nx, lntlist, lnplist)
    do i = 1, nx+1
       lntt(i) = lntlist(nx+2-i)
       lnpp(i) = lnplist(nx+2-i)
    end do
    tt = exp(lntt)
    pp = exp(lnpp)
    do i = 1, nx+1
       call qsat(tt(i), pp(i), qs(i))
    end do
    !lntlist = [x1]
    !lnplist = [np.log(psg)]
    !int_g = integrator(dlnpdlnT_initial, x1, np.log(psg), dx)
    !for i in range(nx):
        !lnt, lnp = int_g.next()
        !lntlist.append(lnt)
        !lnplist.append(lnp)
    !lntlist.reverse()
    !lnplist.reverse()
    !lntt = np.array(lntlist)
    !lnpp = np.array(lnplist)
    !tt = np.exp(lntt)
    !pp = np.exp(lnpp)
    !qs = tt *1.
    !for i in range(len(tt)):
    !    qs[i] = cond.qsat(tt[i], pp[i])

    !return psg, tt, lnpp, qs
end subroutine warm_start

!#==================================================================    
subroutine dlnpdlnT(T,p, slope)

real(8), intent(in) :: T, p
real(8), intent(out):: slope
real(8) :: qs, lh, D, rs, pa, es, A, B, C

    call qsat(T,p, qs)
    call latent_heat(T, lh)
    D = lh / (rvgas * T)

    if (qs .eq. 1) then
       slope = D
    else
        rs = qs / (1. - qs)
        pa = p / (rs /d622 + 1.)
        es = p - pa
        A = lh * rs / rdgas / T
        B = rs * (cp_vapor / cp_air + &
              (lh / rvgas / T-1.) * lh / cp_air / T)
        C = 1. / (kappa * (1. + A) / (1. + B)) !#dlnpa / dlnT
        slope =  (pa * C + es * D) / p
    end if
end subroutine dlnpdlnT

subroutine convec(T1,P1,P2,PH, T2)

real(8), intent(in) :: T1, P1, P2, PH
real(8), intent(out):: T2
real(8) :: lr1, TH, lr2
    call dlnpdlnT(T1,P1, lr1)

    TH = T1 * (PH/P1)**(1./lr1)

    call dlnpdlnT(TH,PH, lr2)

    T2 = T1 * (P2 / P1)**(1./lr2)
end subroutine convec
    !return T2
!##==============================================
subroutine lscale_cond(tin, qin, p_full_in, p_half_in, &
                rain, Tref, qref, p_full, p_half, Ep, rain_profile)

real(8), intent(in), dimension(:) :: tin, qin, p_full_in
real(8), intent(in), dimension(:) :: p_half_in
real(8), intent(out) :: rain, Ep
real(8), intent(out), dimension(:) :: Tref, qref, p_full, rain_profile
real(8), intent(out), dimension(:) :: p_half

!real(8) :: k1, kk1, k2, k3
!real(8) :: T1, P1, P2, PH, T2
!real(8) :: q1, q2, deltap1, deltap2, fq
real(8) :: cp1, cp2, cp3, cp4, deltat
real(8) :: rin, dp, qs, rs, rl, lmass
real(8) :: dpn, kappa2
real(8) :: lh, t1

real(8), allocatable, dimension(:) :: pfulln, kappa1, zf
real(8), allocatable, dimension(:) :: phalfn, zh

integer :: iter1, k, i, n, j, nc

    n = size(tin)
    allocate(pfulln(n))
    allocate(kappa1(n))
    allocate(zf(n))
    allocate(phalfn(n+1))
    allocate(zh(n+1))
    Ep = 0.
    rain = 0.
!    Ep_2 = 0. !check energy conservation

    Tref = tin *1.
    qref = qin *1.
    p_full = p_full_in *1.
    p_half = p_half_in *1.
    rain_profile = 0.
        !for iter1 in range(5): #20):   
            !for k in range(len(p_full)-1, 0, -1):
!        do iter1 = 1, 5 !20 !5
           do k = size(p_full), 2, -1
!                T1 = Tref(k)
!                P1 = p_full(k)
!                P2 = p_full(k-1)
!                PH = p_half(k)
!                call convec(T1, P1, P2, PH, T2)
                !call compute_k(Tref, qref, p_half, k1)
!#do convective adjustment
!                if (Tref(k-1) .lt. T2) then
!                    call qsat(T1,P1, q1)
!                    call qsat(T2,P2, q2)
!                    deltap1 = p_half(k+1) - p_half(k)
!                    deltap2 = p_half(k) - p_half(k-1)
!                    fq = (qref(k) * deltap1 + qref(k-1)*deltap2) &
!                         / (q1*deltap1 + q2 * deltap2)
                    !if (fq .le. 1.) then
                    !#final specific humidity while fixing the full level pressure
!                        q1 = q1 * fq
!                        q2 = q2 * fq
!                    !else
!                    !    q1 = qref(k) + (qref(k-1) - q2)*deltap2 / deltap1
!##q2 is the saturation humidity
!                    !end if
!#energy conservation
!                    cp2 = cp_air * deltap2 * (1 - qref(k-1)) + &
!                          cp_vapor * deltap2 * qref(k-1)
!                    cp1 = cp_air * deltap1 * (1 - qref(k)) + &
!                          cp_vapor * deltap1 * qref(k)
!                    cp4 = cp_air * deltap2 * (1 - q2) + &
!                          cp_vapor * deltap2 * q2
!                    cp3 = cp_air * deltap1 * (1 - q1) + &
!                          cp_vapor * deltap1 * q1
!                    deltat = -(cp3*T1 - cp1*Tref(k) + cp4*T2 - cp2*Tref(k-1)) &
!                             / (cp4 + cp3)                    
!                    Tref(k) = T1 + deltat
!                    Tref(k-1) = T2 + deltat
!                    qref(k) = q1
!                    qref(k-1) = q2
                    !call compute_k(Tref, qref, p_half, kk1)
                    !print *, kk1 - k1
!                end if

                pfulln = p_full *1.
                phalfn = p_half *1.
!#condensation in the kth layer
                i = k
                rin = qref(i) / (1. - qref(i))
                dp = p_half(i+1) - p_half(i)
                
                kappa1 = (rdgas * (1 - qref) + rvgas * qref) &
                         / (cp_air*(1-qref) + cp_vapor* qref)
                   
                call qsat(Tref(i), p_full(i), qs)
                if ( ((qref(i) - qs) .gt. 0) .and. (qs .lt. 1.)) then
                    !#no liquid
                    rs = qs / (1. - qs)
                    rl = rin - rs
                    !roots = newtSolve(ff)
                    !constants = Dummy()
                    !constants.rt = rin
                    !constants.T = Tref[i]
                    !constants.p = p_full[i]
                    !roots.setParams(constants)
                    !deltat = roots(0.01)
                    call  newtSolve_ff(rin, Tref(i), p_full(i), 0.01d+0, deltat)

                    !print *, deltat, rin, Tref(i), p_full(i), i

                    Tref(i) = Tref(i) + deltat
                    call rsat(Tref(i),p_full(i), rs)

                    t1 = Tref(i)
                    !##now there is liquid in the air
                    qs = rs / (1. + rs)
                    
                    rl = rin - rs
                    lmass = dp * (1 - qref(i)) * rl
                    !call compute_k(Tref, qref, p_half, k2)
                    !do j = 1, n
                    !   if (j .eq. i) then
                    !      zf(j) = (cp_air*(1 - qref(j)) + &
                    !                 cp_vapor*(1 - qref(j))*rs + &
                    !                 cp_water*(1 - qref(j))*rl)*Tref(j)
                    !   else
                    !      zf(j) = (cp_air*(1 - qref(j)) + &
                    !                 cp_vapor*qref(j))*Tref(j)
                    !   end if
                    !end do
                    !k2 = sum(zf(1:n)*(p_half(2:n+1) - p_half(1:n))/grav)
!                    #compute geopotential after condensation
                    call compute_geopotential(Tref, qref, p_half, p_full, zf, zh)
                    
                    zh(i) = zh(i+1) + rdgas* Tref(i) / (1. + rin)* &
                            (1. + rs/d622)* &
                            (log(p_half(i+1) / p_half(i)))

                    if (i .eq. 0) then
                        zh(i) = zh(i+1) + 2.*rdgas* Tref(i) / (1. + rin)* &
                                (1. + rs/d622)* &
                                (log(p_half(i+1) / p_full(i)))
                    end if

                    !#new pressure levels             
                    dpn = dp *(1. - qref(i)) * (1. + rs)
                    !p_half[i+1:] = p_half[i+1:] - dp + dpn
                    p_half(i+1:n+1) = p_half(i+1:n+1) - dp + dpn
                    pfulln(i) = p_half(i) + (p_full(i) - p_half(i)) &
                               *(1 - qref(i)) * (1 + rs)

                    kappa2 = (rdgas * (1 - qs) + rvgas * qs) &
                             / (cp_air*(1-qs) + cp_vapor* qs)
                    Tref(i) = Tref(i) * (pfulln(i) / p_full(i))**kappa2
                    if (i .lt. n) then
                        pfulln(i+1:n) = p_full(i+1:n) - dp + dpn
                        Tref(i+1:n) = Tref(i+1:n) * &
                                     (pfulln(i+1:n) / p_full(i+1:n))**kappa1(i+1:n)
                    end if 
                    p_full = pfulln *1.

                    qref(i) = qs
                    rain = rain + lmass / grav

                    rain_profile(i) = rain_profile(i) + lmass / grav

                    Ep = Ep + cp_water * t1 * lmass / grav &
                          + lmass * zf(i) / grav ! (zh(i) + zh(i+1))/2. / grav

                    !call latent_heat(t1, lh)
                    !Ep_2 = Ep_2 + (lh + (cp_water - cp_vapor)*t1) &
                    !        * lmass/ grav
                    !call compute_k(Tref, qref, p_half, k3)
                    !print *, (k3 - k1),  (lh - cp_vapor * t1 - zf(i)) * lmass /grav 
                    !print *, (k3 - k2), (cp_water * t1 * lmass / grav + lmass * &
                    !        (zh(i) + zh(i+1))/2. /grav)
                    !print *, (k2 - kk1), (lh + (cp_water - cp_vapor)*t1) * lmass/ grav, deltat
                    nc = i

                end if
           end do
!        end do

!        if (nc .gt. 1) then
!           !call cold_trap(Tref, qref, p_half, p_full, nc)
!        end if

end subroutine lscale_cond

subroutine cold_trap(t, q, phalf, pfull, nc)
real(8), intent(inout), dimension(:) :: t, q
real(8), intent(in), dimension(:) :: phalf
real(8), intent(in), dimension(:) :: pfull
integer, intent(in) :: nc

real(8) :: q1, q2
integer :: i

q1 = q(nc)
do i = nc-1, 1, -1
   call qsat(t(i), pfull(i), q2)
   if (q2 .gt. q1) exit
   q1 = q2
end do

q(1:i) = q(i+1)
end subroutine cold_trap



end module lscale_cond_mod

