module r_lambda

    implicit none
    
    real(8), allocatable, save :: omega(:), A(:,:)
    real(8), save :: l, lin, k11, k12, k21, k22, over_k11, over_2k12, over_2k21, over_k22, &
                     pi, pi2
    integer, save :: nosc, nosc1, nosc2, nosc1_plus1
    
end module r_lambda







program main

    use r_lambda

    implicit none

    external F
    integer, parameter :: rk = kind ( 1.0D+00 )
    integer flag, i, j, l_steps, n_pts
    real ( kind = rk ), allocatable :: y(:), yp(:), r1_t(:), r2_t(:), &
                                       r_time(:)
    real ( kind = rk ) t, relerr, abserr, aux, r1, r2, r, &
                       w, d, freq, prob_con, t_1, t_2, mean_omega, &
                       l_ini, l_fin, soma_r, soma2_r, sigma_r, mean_r
    !real ( kind = rk ) phi, psi

    call semente(1)
        
    open ( unit = 1, file = 'r_lambda_vez1.dat', status = 'unknown')
    write (1, *) "l ", "lin ", "r1 ", "r2 ", "r ", "mean_r ", "sigma_r"

    pi = acos(-1.0)
    pi2 = 2.0*pi
    abserr = sqrt ( epsilon ( abserr ) )
    relerr = sqrt ( epsilon ( relerr ) )

    nosc1 = 200
    nosc2 = 100
    nosc = nosc1 + nosc2

    allocate(omega(nosc), y(nosc), yp(nosc), A(nosc,nosc))
    allocate(r1_t(40), r2_t(40))
    
    n_pts = 40
    allocate(r_time(n_pts))
    
    d = 1.0
    w = 1.5
    
    l_ini = 0.0
    l_fin = 3.0
    l_steps = 50

    lin = 1.0
        
    prob_con = 0.5

    do j=1,l_steps
    
        l = l_ini + j*(l_fin - l_ini)/float(l_steps)
    
        call matrix(prob_con)

        do i=1,nosc
            call random_number(aux)
            y(i) = aux*pi2
        end do
    
        do i=1,nosc
            call gaussian(d,freq)
            omega(i) = freq
        end do
        mean_omega = sum(omega)/nosc
        omega = omega - mean_omega
        do i=1,nosc1
            omega(i) = omega(i) + w
        end do
    
        t = 0.0
        
        r1_t = 0.0
        r2_t = 0.0
        
        do while ((minval(r1_t)<0.9).or.(minval(r2_t)<0.9))
            
            do i=0,39
            
                flag = 1
                
                1 t_1 = t + 0.2*i
                t_2 = t + 0.2*(i+1)
    
                call rkf45( F, nosc, y, yp, t_1, t_2, relerr, abserr, flag)
                
                call para_ord(y,r1,r2,r) !para_ord(y,r1,r2,r,phi,psi)
                
                r1_t(i+1) = r1
                r2_t(i+1) = r2
                
                if (flag.eq.4) go to 1
                
            end do
            
            t = t_2
            
            lin = lin + 0.1
            
        end do
        
        lin = lin - 0.1
        
        do i=0,n_pts-1
            
            flag = 1
    
            2 t_1 = t + 0.5*i
            t_2 = t + 0.5*(i+1)
        
            call rkf45( F, nosc, y, yp, t_1, t_2, relerr, abserr, flag)
            
            call para_ord(y,r1,r2,r) !para_ord(y,r1,r2,r,phi,psi)
            
            r_time(i+1) = r
            
            if (flag.eq.4) go to 2
    
        end do
        
        soma_r = sum(r_time)
        mean_r = soma_r/n_pts
        soma2_r = dot_product(r_time,r_time)
        sigma_r = sqrt((n_pts*soma2_r-soma_r*soma_r)/(n_pts*(n_pts-1)))
        
        write (1, *) l, lin, r1, r2, r, mean_r, sigma_r
        
        print*, l, lin
        
        lin = lin - 0.7
        
    end do

    close (1)

    call semente(2)

end program main







subroutine F(t,y,yp)

    use r_lambda

    implicit none

    integer, parameter :: rk = kind ( 1.0D+00 )
    real ( kind = rk ) t, y(nosc), yp(nosc), soma
    integer i

    call faz_nada(t)

    yp = omega
        
    do i=1,nosc1

        soma = sum( sin(y(1:nosc1) - y(i)), mask = A(i,1:nosc1) == 1.0_rk )

        yp(i) = yp(i) + lin * soma * over_k11

        soma = sum( sin(y(nosc1+1:nosc) - y(i)), mask = A(i,nosc1+1:nosc) == 1.0_rk )
        
        yp(i) = yp(i) + l * soma * over_2k12

    end do

    do i=nosc1+1,nosc

        soma = sum( sin(y(1:nosc1) - y(i)), mask = A(i,1:nosc1) == 1.0_rk )

        yp(i) = yp(i) + l * soma * over_2k21

        soma = sum( sin(y(nosc1+1:nosc) - y(i)), mask = A(i,nosc1+1:nosc) == 1.0_rk )
        
        yp(i) = yp(i) + lin * soma * over_k22
        
    end do

    return

end subroutine







subroutine gaussian(larg,freq)
    
    use r_lambda
    
    implicit none

    integer, parameter :: rk = kind ( 1.0D+00 )
    real ( kind = rk ) larg,aleat,freq,error,freq0,delta_freq,F,rho,anorm,b

    call random_number(aleat)

    error=0.001
    b = 1.0/(2.0*larg*larg)
    anorm = sqrt(b/pi)
    freq0 = 0.0
    F = 0.5*(1+erf(freq0*sqrt(b))) - aleat
    rho = anorm*exp(-b*freq0*freq0)

    do while(abs(F) > error)
        delta_freq = - F/rho
        freq = freq0 + delta_freq
        F = 0.5*(1+erf(freq*sqrt(b))) - aleat
        rho = anorm*exp(-b*freq*freq)
        freq0 = freq
    end do

    return

end






subroutine faz_nada(t)

    implicit none

    integer, parameter :: rk = kind ( 1.0D+00 )
    real ( kind = rk ) t

    if (t.ne.t) then
        write (*,*) 't is NAN'
    end if

    return

end subroutine







subroutine para_ord(y,r1,r2,r) !para_ord(y,r1,r2,r,phi,psi)

    use r_lambda

    implicit none

    integer, parameter :: rk = kind ( 1.0D+00 )
    integer i
    !real ( kind = rk ) aux
    real ( kind = rk ), intent(in) :: y(nosc)
    real ( kind = rk ), intent(out) :: r1, r2, r
    !real ( kind = rk ), intent(out) :: phi, psi
    complex ( kind = rk ) imaginario, z1cplx(nosc1), z2cplx(nosc2), &
                          zcplx(nosc), z1, z2, z

    imaginario = (0.0,1.0)

    do i=1,nosc1
        z1cplx(i) = exp(imaginario*y(i))
        zcplx(i) = z1cplx(i)
    end do
    do i=1,nosc2
        z2cplx(i) = exp(imaginario*y(nosc1+i))
        zcplx(nosc1+i) = z2cplx(i)
    end do

    z1 = sum(z1cplx)/nosc1
    z2 = sum(z2cplx)/nosc2
    z = sum(zcplx)/nosc

    r1 = abs(z1)
    r2 = abs(z2)
    r = abs(z)
    
    !psi = atan2(aimag(z),real(z))
    !aux = atan2(aimag(z1),real(z1)) - atan2(aimag(z2),real(z2))
    !if (aux.gt.pi) then
    !  phi = aux - pi2
    !else if (aux.lt.-pi) then
    !  phi = aux + pi2
    !else
    !  phi = aux
    !end if

    return

end subroutine







subroutine matrix(prob_con)

    use r_lambda

    implicit none

    integer, parameter :: rk = kind ( 1.0D+00 )
    real ( kind = rk ), intent(in) :: prob_con
    real ( kind = rk ) aleat, C11, C12, C22
    integer i, j

    C11 = 0.0
    C12 = 0.0
    C22 = 0.0

    do i=1,nosc
        do j=i,nosc
            call random_number(aleat)
            if ( i == j ) then
                A(i,j) = 0
            end if
            if ( aleat.le.prob_con ) then
                A(i,j) = 1
                A(j,i) = 1
                if ( j.le.nosc1 ) then
                    C11 = C11 + 1.0
                else if ( i.gt.nosc1 ) then
                    C22 = C22 + 1.0
                else
                    C12 = C12 + 1.0
                end if
            else
                A(i,j) = 0
                A(j,i) = 0
            end if
        end do
    end do

    C12 = C12 / 2.0
    
    k11 = C11/nosc1
    k12 = C12/nosc1
    k21 = C12/nosc2
    k22 = C22/nosc2
  
    over_k11 = 1.0_rk / (k11)
    over_2k12 = 1.0_rk / (2*k12)
    over_2k21 = 1.0_rk / (2*k21)
    over_k22 = 1.0_rk / (k22)

    return

end subroutine







subroutine semente(i)

    implicit none

    integer i, iseed(33)

    if (i.eq.1) then
        open(unit=50,file='seed.in',status='OLD')
            read(50,*) iseed(1),iseed(2)
        close(50)
        call random_seed(put=iseed)
    elseif (i.eq.2) then
        call random_seed(get=iseed)
        open(unit=50,file='seed.in',status='OLD',position='REWIND')
            write (50,*) iseed(1),iseed(2)
        close(50)
    end if

end subroutine