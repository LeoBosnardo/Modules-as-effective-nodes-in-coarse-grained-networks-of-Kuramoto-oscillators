module r_lambda

    implicit none
    
    real(8), save :: omega(5)
    real(8), save :: l, pi, pi2
    
end module r_lambda







program main

    use r_lambda

    implicit none

    external F
    integer, parameter :: rk = kind ( 1.0D+00 )
    integer flag, i, j, l_steps, n_pts
    real ( kind = rk ) :: y(5), yp(5)
    real ( kind = rk ), allocatable :: r_t(:)
    real ( kind = rk ) t, t_1, t_2, &
                       relerr, abserr, aux, r, &
                       l_ini, l_fin, &
                       mean_r, soma_r, soma2_r, sigmar, &
                       N1, N2, N3, N4, N5, Nt, &
                       rho1, rho2, rho3, rho4, rho5, &
                       r1, r2, r3, r4, r5, &
                       aux0, p12, p13, p14, p15, p23, p24, p25, p34, p35, p45
      
    open ( unit = 1, file = 'r_lambda.dat', status = 'unknown')
    write (1, *) "l ",  "r ", "meanr ", "sigmar"

    !open ( unit = 2, file = 'rs_lin_tempo.dat', status = 'unknown')
    !write (2, *) "t ", "r"

    pi = acos(-1.0)
    pi2 = 2.0*pi
    abserr = sqrt ( epsilon ( abserr ) )
    relerr = sqrt ( epsilon ( relerr ) )

    N1 = 120
    N2 = 83
    N3 = 34
    N4 = 7
    N5 = 4
    Nt = N1 + N2 + N3 + N4 + N5
    rho1 = N1 / Nt
    rho2 = N2 / Nt
    rho3 = N3 / Nt
    rho4 = N4 / Nt
    rho5 = N5 / Nt
    r1 = 0.9
    r2 = 0.9
    r3 = 0.9
    r4 = 0.9
    r5 = 0.9
    aux0 = (rho1*r1)**2 + (rho2*r2)**2 + (rho3*r3)**2 + (rho4*r4)**2 + (rho5*r5)**2
    p12 = 2*rho1*rho2*r1*r2
    p13 = 2*rho1*rho3*r1*r3
    p14 = 2*rho1*rho4*r1*r4
    p15 = 2*rho1*rho5*r1*r5
    p23 = 2*rho2*rho3*r2*r3
    p24 = 2*rho2*rho4*r2*r4
    p25 = 2*rho2*rho5*r2*r5
    p34 = 2*rho3*rho4*r3*r4
    p35 = 2*rho3*rho5*r3*r5
    p45 = 2*rho4*rho5*r4*r5
    
    allocate(r_t(1000))
    n_pts = 1000
    
    omega(1) = 0.0
    omega(2) = 1.0
    omega(3) = 1.0
    omega(4) = -2.0
    omega(5) = 7.0
    
    l_ini = 0.0
    l_fin = 10.0
    l_steps = 100

    do j=0,l_steps
    
        l = l_ini + j*(l_fin - l_ini)/float(l_steps)

        do i=1,5
            call random_number(aux)
            y(i) = aux*pi2
        end do
        
        r_t = 0.0
            
        do i=0,999
            
            flag = 1
                
            1 t_1 = 0.1*i
            t_2 = 0.1*(i+1)
                
            call rkf45( F, 5, y, yp, t_1, t_2, relerr, abserr, flag)
                
            r = sqrt(aux0+p12*cos(y(1)-y(2))+&
                          p13*cos(y(1)-y(3))+&
                          p14*cos(y(1)-y(4))+&
                          p15*cos(y(1)-y(5))+&
                          p23*cos(y(2)-y(3))+&
                          p24*cos(y(2)-y(4))+&
                          p25*cos(y(2)-y(5))+&
                          p34*cos(y(3)-y(4))+&
                          p35*cos(y(3)-y(5))+&
                          p45*cos(y(4)-y(5)))
                
            if (flag.eq.4) go to 1
                
            !write (2, *) t_2, r
                
        end do
        
        t = t_2
        
        do i=1,1000
            
            flag = 1
    
            2 t_1 = t + 0.1*i
            t_2 = t + 0.1*(i+1)
        
            call rkf45( F, 5, y, yp, t_1, t_2, relerr, abserr, flag)
            
            r = sqrt(aux0+p12*cos(y(1)-y(2))+&
                          p13*cos(y(1)-y(3))+&
                          p14*cos(y(1)-y(4))+&
                          p15*cos(y(1)-y(5))+&
                          p23*cos(y(2)-y(3))+&
                          p24*cos(y(2)-y(4))+&
                          p25*cos(y(2)-y(5))+&
                          p34*cos(y(3)-y(4))+&
                          p35*cos(y(3)-y(5))+&
                          p45*cos(y(4)-y(5)))

            r_t(i) = r
            
            if (flag.eq.4) go to 2
            
            !write (2, *) t_2, r
    
        end do
            
        soma_r = sum(r_t)
        mean_r = soma_r/n_pts
        soma2_r = dot_product(r_t,r_t)
        sigmar = sqrt((n_pts*soma2_r-soma_r*soma_r)/(n_pts*(n_pts-1)))
        
        write (1, *) l, r, mean_r, sigmar
        
    end do

    close (1)
    !close (2)

end program main







subroutine F(t,y,yp)

    use r_lambda

    implicit none

    integer, parameter :: rk = kind ( 1.0D+00 )
    real ( kind = rk ) t, y(5), yp(5)
    integer i

    call faz_nada(t)
        
    do i=1,5

        yp(i) = omega(i) + l * sum( sin(y(1:5) - y(i)) ) / 5

    end do

    return

end subroutine







subroutine faz_nada(t)

    implicit none

    integer, parameter :: rk = kind ( 1.0D+00 )
    real ( kind = rk ) t

    if (t.ne.t) then
        write (*,*) 't is NAN'
    end if

    return

end subroutine