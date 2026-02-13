module r_lambda

    implicit none
    
    real(8), save :: omega(3)
    real(8), save :: l, pi, pi2
    
end module r_lambda







program main

    use r_lambda

    implicit none

    external F
    integer, parameter :: rk = kind ( 1.0D+00 )
    integer flag, i, j, l_steps, n_pts
    real ( kind = rk ) :: y(3), yp(3)
    real ( kind = rk ), allocatable :: r_t(:)
    real ( kind = rk ) t, t_1, t_2, &
                       relerr, abserr, aux, r, &
                       l_ini, l_fin, &
                       mean_r, soma_r, soma2_r, sigmar, &
                       N1, N2, N3, Nt, &
                       rho1, rho2, rho3, &
                       r1, r2, r3, &
                       aux0, aux12, aux13, aux23
      
    open ( unit = 1, file = 'r_lambda.dat', status = 'unknown')
    write (1, *) "l ",  "r ", "meanr ", "sigmar"

    !open ( unit = 2, file = 'rs_lin_tempo.dat', status = 'unknown')
    !write (2, *) "t ", "r"

    pi = acos(-1.0)
    pi2 = 2.0*pi
    abserr = sqrt ( epsilon ( abserr ) )
    relerr = sqrt ( epsilon ( relerr ) )

    N1 = 130
    N2 = 77
    N3 = 41
    Nt = N1 + N2 + N3
    rho1 = N1 / Nt
    rho2 = N2 / Nt
    rho3 = N3 / Nt
    r1 = 1.0
    r2 = 1.0
    r3 = 1.0
    aux0 = (rho1*r1)**2 + (rho2*r2)**2 + (rho3*r3)**2
    aux12 = 2*rho1*rho2*r1*r2
    aux13 = 2*rho1*rho3*r1*r3
    aux23 = 2*rho2*rho3*r2*r3
    
    allocate(r_t(1000))
    n_pts = 1000
    
    omega(1) = 0.0
    omega(2) = 1.0
    omega(3) = 6.0
    
    l_ini = 0.0
    l_fin = 10.0
    l_steps = 100

    do j=0,l_steps
    
        l = l_ini + j*(l_fin - l_ini)/float(l_steps)

        do i=1,3
            call random_number(aux)
            y(i) = aux*pi2
        end do
        
        r_t = 0.0
            
        do i=0,999
            
            flag = 1
                
            1 t_1 = 0.1*i
            t_2 = 0.1*(i+1)
                
            call rkf45( F, 3, y, yp, t_1, t_2, relerr, abserr, flag)
                
            r = sqrt(aux0+aux12*cos(y(1)-y(2))+aux13*cos(y(1)-y(3))+aux23*cos(y(2)-y(3)))
                
            if (flag.eq.4) go to 1
                
            !write (2, *) t_2, r
                
        end do
        
        t = t_2
        
        do i=1,1000
            
            flag = 1
    
            2 t_1 = t + 0.1*i
            t_2 = t + 0.1*(i+1)
        
            call rkf45( F, 3, y, yp, t_1, t_2, relerr, abserr, flag)
            
            r = sqrt(aux0+aux12*cos(y(1)-y(2))+aux13*cos(y(1)-y(3))+aux23*cos(y(2)-y(3)))

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
    real ( kind = rk ) t, y(3), yp(3)
    integer i

    call faz_nada(t)
        
    do i=1,3

        yp(i) = omega(i) + l * sum( sin(y(1:3) - y(i)) ) / 3

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