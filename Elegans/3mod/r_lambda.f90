module r_lambda

    implicit none
    
    real(8), save :: omega(248)
    integer, save :: m1(130), m2(77), m3(41), modules(248)
    real(8), save :: l, lin, pi, pi2
    integer, save :: nosc, nosc1, nosc2, nosc3, nmod
    real, dimension(248,248) :: A
    real, dimension(3,3) :: K
    
end module r_lambda







program main

    use r_lambda

    implicit none

    external F
    integer, parameter :: rk = kind ( 1.0D+00 )
    integer flag, i, j, l_steps, n_pts
    real ( kind = rk ) :: y(248), yp(248)
    real ( kind = rk ), allocatable :: r_t(:)
    real ( kind = rk ) t, t_1, t_2, &
                       relerr, abserr, aux, &
                       w1, w2, w3, &
                       l_ini, l_fin, &
                       r, mean_r, soma_r, soma2_r, sigmar
    
    nmod = 3
    nosc1 = 130
    nosc2 = 77
    nosc3 = 41
    nosc = 248

    open(unit=10, file="A.txt", status="old", action="read")                 
    do i = 1, 248
        read(10,*) (A(j, i), j=1, 248)
    end do                
    close(10)

    open(unit=10, file="K.txt", status="old", action="read")                 
    do i = 1, nmod
        read(10,*) (K(j, i), j=1, nmod)
    end do                
    close(10)

    open(unit=10, file="m1.txt", status="old", action="read")                 
    do i = 1, nosc1
        read(10,*) m1(i)
    end do                
    close(10)
      
    open(unit=10, file="m2.txt", status="old", action="read")                 
    do i = 1, nosc2
        read(10,*) m2(i)
    end do                
    close(10)
      
    open(unit=10, file="m3.txt", status="old", action="read")                 
    do i = 1, nosc3
        read(10,*) m3(i)
    end do                
    close(10)
    
    open(unit=10, file="module.txt", status="old", action="read")                 
    do i = 1, nosc
        read(10,*) modules(i)
    end do                
    close(10)

    open ( unit = 1, file = 'r_lambda.dat', status = 'unknown')
    write (1, *) "l ", "lin ", "r ", "mean_r ", "sigma_r "

    !open ( unit = 2, file = 'rs_lin_tempo.dat', status = 'unknown')
    !write (2, *) "t ", "r "

    pi = acos(-1.0)
    pi2 = 2.0*pi
    abserr = sqrt ( epsilon ( abserr ) )
    relerr = sqrt ( epsilon ( relerr ) )
    
    allocate(r_t(5000))
    n_pts = 5000
    
    w1 = 0.0
    w2 = 1.0
    w3 = 6.0

    do i=1,nosc1
        omega(m1(i)) = w1
    end do
    do i=1,nosc2
        omega(m2(i)) = w2
    end do
    do i=1,nosc3
        omega(m3(i)) = w3
    end do

    l_ini = 0.0
    l_fin = 10.0
    l_steps = 50

    lin = 20.0
    
    do j=0,l_steps
    
        l = l_ini + j*(l_fin - l_ini)/float(l_steps)

        do i=1,nosc
            call random_number(aux)
            y(i) = aux*pi2
        end do
    
        t = 0.0
            
        do i=1,5000

            flag = 1
    
            1 t_1 = t + 0.2*(i-1)
            t_2 = t + 0.2*i
        
            call rkf45( F, nosc, y, yp, t_1, t_2, relerr, abserr, flag)
            if (flag.eq.4) go to 1
            
            call para_ord(y,r)
            
            r_t(i) = r
            
            !write (2, *) t_2, r
    
        end do
        
        soma_r = sum(r_t)
        mean_r = soma_r/n_pts
        soma2_r = dot_product(r_t,r_t)
        sigmar = sqrt((n_pts*soma2_r-soma_r*soma_r)/(n_pts*(n_pts-1)))
        
        write (1, *) l, lin, r, mean_r, sigmar
        
        print*, l, mean_r
        
    end do

    close (1)
    !close (2)

end program main







subroutine F(t,y,yp)

    use r_lambda

    implicit none

    integer, parameter :: rk = kind ( 1.0D+00 )
    real ( kind = rk ) t, y(248), yp(248), soma
    integer i, j, mod_i, mod_j

    call faz_nada(t)

    yp = omega

    do i=1,nosc

        mod_i = modules(i)

        soma = 0.0

        do j=1,nosc

            if ( A(i,j) /= 0.0 ) then

                mod_j = modules(j)

                if ( mod_i == mod_j ) then

                    soma = soma + sin(y(j) - y(i)) * A(i,j) * lin / K(mod_j,mod_i)

                else

                    soma = soma + sin(y(j) - y(i)) * A(i,j) * l / (2 * K(mod_j,mod_i))

                end if

            end if

        end do

        yp(i) = yp(i) + soma

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







subroutine para_ord(y,r)

    use r_lambda

    implicit none

    integer, parameter :: rk = kind ( 1.0D+00 )
    integer i
    real ( kind = rk ), intent(in) :: y(248)
    real ( kind = rk ), intent(out) :: r
    complex ( kind = rk ) imaginario, zcplx(nosc), z

    imaginario = (0.0,1.0)

    do i=1,nosc
        zcplx(i) = exp(imaginario*y(i))
    end do

    z = sum(zcplx)/nosc
    r = abs(z)

    return

end subroutine