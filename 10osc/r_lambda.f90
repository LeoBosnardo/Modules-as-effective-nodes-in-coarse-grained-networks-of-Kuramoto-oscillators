module r_lambda

    implicit none
    
    real(8), save :: omega(10)
    real(8), save :: l, pi, pi2
    
end module r_lambda







program main

    use r_lambda

    implicit none

    external F
    integer, parameter :: rk = kind ( 1.0D+00 )
    integer flag, i, j, l_steps, n_pts
    real ( kind = rk ) :: y(10), yp(10)
    real ( kind = rk ), allocatable :: r_t(:)
    real ( kind = rk ) t, t_1, t_2, &
                       relerr, abserr, aux, r, &
                       l_ini, l_fin, &
                       mean_r, soma_r, soma2_r, sigmar, &
                       N1, N2, N3, N4, N5, N6, N7, N8, N9, N10, Nt, &
                       rho1, rho2, rho3, rho4, rho5, rho6, rho7, rho8, rho9, rho10, &
                       r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, &
                       aux0, p12, p13, p14, p15, p16, p17, p18, p19, p110, &
                       p23, p24, p25, p26, p27, p28, p29, p210, &
                       p34, p35, p36, p37, p38, p39, p310, &
                       p45, p46, p47, p48, p49, p410, &
                       p56, p57, p58, p59, p510, &
                       p67, p68, p69, p610, &
                       p78, p79, p710, &
                       p89, p810, &
                       p910
      
    open ( unit = 1, file = 'r_lambda.dat', status = 'unknown')
    write (1, *) "l ",  "r ", "meanr ", "sigmar"

    !open ( unit = 2, file = 'rs_lin_tempo.dat', status = 'unknown')
    !write (2, *) "t ", "r"

    pi = acos(-1.0)
    pi2 = 2.0*pi
    abserr = sqrt ( epsilon ( abserr ) )
    relerr = sqrt ( epsilon ( relerr ) )

    N1 = 68
    N2 = 33
    N3 = 76
    N4 = 11
    N5 = 13
    N6 = 7
    N7 = 6
    N8 = 8
    N9 = 8
    N10 = 18
    Nt = N1 + N2 + N3 + N4 + N5 + N6 + N7 + N8 + N9 + N10
    rho1 = N1 / Nt
    rho2 = N2 / Nt
    rho3 = N3 / Nt
    rho4 = N4 / Nt
    rho5 = N5 / Nt
    rho6 = N6 / Nt
    rho7 = N7 / Nt
    rho8 = N8 / Nt
    rho9 = N9 / Nt
    rho10 = N10 / Nt
    r1 = 0.9
    r2 = 0.9
    r3 = 0.9
    r4 = 0.9
    r5 = 0.9
    r6 = 0.9
    r7 = 0.9
    r8 = 0.9
    r9 = 0.9
    r10 = 0.9
    aux0 = (rho1*r1)**2 + (rho2*r2)**2 + (rho3*r3)**2 + (rho4*r4)**2 + (rho5*r5)**2 +&
           (rho6*r6)**2 + (rho7*r7)**2 + (rho8*r8)**2 + (rho9*r9)**2 + (rho10*r10)**2
    p12 = 2*rho1*rho2*r1*r2
    p13 = 2*rho1*rho3*r1*r3
    p14 = 2*rho1*rho4*r1*r4
    p15 = 2*rho1*rho5*r1*r5
    p16 = 2*rho1*rho6*r1*r6
    p17 = 2*rho1*rho7*r1*r7
    p18 = 2*rho1*rho8*r1*r8
    p19 = 2*rho1*rho9*r1*r9
    p110 = 2*rho1*rho10*r1*r10
    p23 = 2*rho2*rho3*r2*r3
    p24 = 2*rho2*rho4*r2*r4
    p25 = 2*rho2*rho5*r2*r5
    p26 = 2*rho2*rho6*r2*r6
    p27 = 2*rho2*rho7*r2*r7
    p28 = 2*rho2*rho8*r2*r8
    p29 = 2*rho2*rho9*r2*r9
    p210 = 2*rho2*rho10*r2*r10
    p34 = 2*rho3*rho4*r3*r4
    p35 = 2*rho3*rho5*r3*r5
    p36 = 2*rho3*rho6*r3*r6
    p37 = 2*rho3*rho7*r3*r7
    p38 = 2*rho3*rho8*r3*r8
    p39 = 2*rho3*rho9*r3*r9
    p310 = 2*rho3*rho10*r3*r10
    p45 = 2*rho4*rho5*r4*r5
    p46 = 2*rho4*rho6*r4*r6
    p47 = 2*rho4*rho7*r4*r7
    p48 = 2*rho4*rho8*r4*r8
    p49 = 2*rho4*rho9*r4*r9
    p410 = 2*rho4*rho10*r4*r10
    p56 = 2*rho5*rho6*r5*r6
    p57 = 2*rho5*rho7*r5*r7
    p58 = 2*rho5*rho8*r5*r8
    p59 = 2*rho5*rho9*r5*r9
    p510 = 2*rho5*rho10*r5*r10
    p67 = 2*rho6*rho7*r6*r7
    p68 = 2*rho6*rho8*r6*r8
    p69 = 2*rho6*rho9*r6*r9
    p610 = 2*rho6*rho10*r6*r10
    p78 = 2*rho7*rho8*r7*r8
    p79 = 2*rho7*rho9*r7*r9
    p710 = 2*rho7*rho10*r7*r10
    p89 = 2*rho8*rho9*r8*r9
    p810 = 2*rho8*rho10*r8*r10
    p910 = 2*rho9*rho10*r9*r10
    
    allocate(r_t(1000))
    n_pts = 1000
    
    omega(1) = 0.0
    omega(2) = 1.0
    omega(3) = -3.5
    omega(4) = -2.0
    omega(5) = 2.0
    omega(6) = 0.6
    omega(7) = 0.4
    omega(8) = 1.2
    omega(9) = 1.5
    omega(10) = -1.5
    
    l_ini = 0.0
    l_fin = 10.0
    l_steps = 100

    do j=0,l_steps
    
        l = l_ini + j*(l_fin - l_ini)/float(l_steps)

        do i=1,10
            call random_number(aux)
            y(i) = aux*pi2
        end do
        
        r_t = 0.0
            
        do i=0,999
            
            flag = 1
                
            1 t_1 = 0.1*i
            t_2 = 0.1*(i+1)
                
            call rkf45( F, 10, y, yp, t_1, t_2, relerr, abserr, flag)
                
            r = sqrt(aux0+p12*cos(y(1)-y(2))+&
                          p13*cos(y(1)-y(3))+&
                          p14*cos(y(1)-y(4))+&
                          p15*cos(y(1)-y(5))+&
                          p16*cos(y(1)-y(6))+&
                          p17*cos(y(1)-y(7))+&
                          p18*cos(y(1)-y(8))+&
                          p19*cos(y(1)-y(9))+&
                          p110*cos(y(1)-y(10))+&
                          p23*cos(y(2)-y(3))+&
                          p24*cos(y(2)-y(4))+&
                          p25*cos(y(2)-y(5))+&
                          p26*cos(y(2)-y(6))+&
                          p27*cos(y(2)-y(7))+&
                          p28*cos(y(2)-y(8))+&
                          p29*cos(y(2)-y(9))+&
                          p210*cos(y(2)-y(10))+&
                          p34*cos(y(3)-y(4))+&
                          p35*cos(y(3)-y(5))+&
                          p36*cos(y(3)-y(6))+&
                          p37*cos(y(3)-y(7))+&
                          p38*cos(y(3)-y(8))+&
                          p39*cos(y(3)-y(9))+&
                          p310*cos(y(3)-y(10))+&
                          p45*cos(y(4)-y(5))+&
                          p46*cos(y(4)-y(6))+&
                          p47*cos(y(4)-y(7))+&
                          p48*cos(y(4)-y(8))+&
                          p49*cos(y(4)-y(9))+&
                          p410*cos(y(4)-y(10))+&
                          p56*cos(y(5)-y(6))+&
                          p57*cos(y(5)-y(7))+&
                          p58*cos(y(5)-y(8))+&
                          p59*cos(y(5)-y(9))+&
                          p510*cos(y(5)-y(10))+&
                          p67*cos(y(6)-y(7))+&
                          p68*cos(y(6)-y(8))+&
                          p69*cos(y(6)-y(9))+&
                          p610*cos(y(6)-y(10))+&
                          p78*cos(y(7)-y(8))+&
                          p79*cos(y(7)-y(9))+&
                          p710*cos(y(7)-y(10))+&
                          p89*cos(y(8)-y(9))+&
                          p810*cos(y(8)-y(10))+&
                          p910*cos(y(9)-y(10)))
                
            if (flag.eq.4) go to 1
                
            !write (2, *) t_2, r
                
        end do
        
        t = t_2
        
        do i=1,1000
            
            flag = 1
    
            2 t_1 = t + 0.1*i
            t_2 = t + 0.1*(i+1)
        
            call rkf45( F, 10, y, yp, t_1, t_2, relerr, abserr, flag)
            
            r = sqrt(aux0+p12*cos(y(1)-y(2))+&
                          p13*cos(y(1)-y(3))+&
                          p14*cos(y(1)-y(4))+&
                          p15*cos(y(1)-y(5))+&
                          p16*cos(y(1)-y(6))+&
                          p17*cos(y(1)-y(7))+&
                          p18*cos(y(1)-y(8))+&
                          p19*cos(y(1)-y(9))+&
                          p110*cos(y(1)-y(10))+&
                          p23*cos(y(2)-y(3))+&
                          p24*cos(y(2)-y(4))+&
                          p25*cos(y(2)-y(5))+&
                          p26*cos(y(2)-y(6))+&
                          p27*cos(y(2)-y(7))+&
                          p28*cos(y(2)-y(8))+&
                          p29*cos(y(2)-y(9))+&
                          p210*cos(y(2)-y(10))+&
                          p34*cos(y(3)-y(4))+&
                          p35*cos(y(3)-y(5))+&
                          p36*cos(y(3)-y(6))+&
                          p37*cos(y(3)-y(7))+&
                          p38*cos(y(3)-y(8))+&
                          p39*cos(y(3)-y(9))+&
                          p310*cos(y(3)-y(10))+&
                          p45*cos(y(4)-y(5))+&
                          p46*cos(y(4)-y(6))+&
                          p47*cos(y(4)-y(7))+&
                          p48*cos(y(4)-y(8))+&
                          p49*cos(y(4)-y(9))+&
                          p410*cos(y(4)-y(10))+&
                          p56*cos(y(5)-y(6))+&
                          p57*cos(y(5)-y(7))+&
                          p58*cos(y(5)-y(8))+&
                          p59*cos(y(5)-y(9))+&
                          p510*cos(y(5)-y(10))+&
                          p67*cos(y(6)-y(7))+&
                          p68*cos(y(6)-y(8))+&
                          p69*cos(y(6)-y(9))+&
                          p610*cos(y(6)-y(10))+&
                          p78*cos(y(7)-y(8))+&
                          p79*cos(y(7)-y(9))+&
                          p710*cos(y(7)-y(10))+&
                          p89*cos(y(8)-y(9))+&
                          p810*cos(y(8)-y(10))+&
                          p910*cos(y(9)-y(10)))

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
    real ( kind = rk ) t, y(10), yp(10)
    integer i

    call faz_nada(t)
        
    do i=1,10

        yp(i) = omega(i) + l * sum( sin(y(1:10) - y(i)) ) / 10

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