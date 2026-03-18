module r_lambda

    implicit none
    
    real(8), save :: omega(130)
    real(8), save :: lin, pi, pi2
    integer, save :: nosc1, nosc2, nosc3, nmod, focal
    integer, dimension(130,130) :: A11
    integer, dimension(77,77) :: A22
    integer, dimension(41,41) :: A33
    real, dimension(3,3) :: K
    
end module r_lambda







program main

    use r_lambda

    implicit none

    external F
    integer, parameter :: rk = kind ( 1.0D+00 )
    integer flag, i, j, s, lin_steps, eps_steps
    real ( kind = rk ) :: y(130), yp(130)
    real ( kind = rk ) y0m1(130), y0m2(77), y0m3(41)
    real ( kind = rk ) t, t_1, t_2, &
                       relerr, abserr, aux, &
                       w, &
                       lin_ini, lin_fin, &
                       r, eps, epsc, &
                       eps_ini, eps_fin
    
    nmod = 3
    focal = 1
    nosc1 = 130
    nosc2 = 77
    nosc3 = 41

    open(unit=10, file="A11.txt", status="old", action="read")                 
    do i = 1, nosc1
        read(10,*) (A11(j, i), j=1, nosc1)
    end do                
    close(10)

    open(unit=10, file="A22.txt", status="old", action="read")                 
    do i = 1, nosc2
        read(10,*) (A22(j, i), j=1, nosc2)
    end do                
    close(10)

    open(unit=10, file="A33.txt", status="old", action="read")                 
    do i = 1, nosc3
        read(10,*) (A33(j, i), j=1, nosc3)
    end do                
    close(10)

    open(unit=10, file="K.txt", status="old", action="read")                 
    do i = 1, nmod
        read(10,*) (K(j, i), j=1, nmod)
    end do                
    close(10)

    open(unit=10, file="cond_init_m1.dat", status="old", action="read")                 
    do i = 1, nosc1
        read(10,*) y0m1(i)
    end do                
    close(10)

    open(unit=10, file="cond_init_m2.dat", status="old", action="read")                 
    do i = 1, nosc2
        read(10,*) y0m2(i)
    end do                
    close(10)
    
    open(unit=10, file="cond_init_m3.dat", status="old", action="read")                 
    do i = 1, nosc3
        read(10,*) y0m3(i)
    end do                
    close(10)

    !open ( unit = 1, file = 'epsilon_lambda.dat', status = 'unknown')
    !write (1, *) "lin ", "ep "

    open ( unit = 2, file = 'rs_lin_tempo.dat', status = 'unknown')
    !write (2, *) "t ", "r "

    pi = acos(-1.0)
    pi2 = 2.0*pi
    abserr = sqrt ( epsilon ( abserr ) )
    relerr = sqrt ( epsilon ( relerr ) )
    
    w = 0.0

    do i=1,nosc1
        omega(i) = w
    end do

    lin_ini = 1.0
    lin_fin = 1.0
    lin_steps = 1

    eps_ini = pi!0.0
    eps_fin = pi!2.0*pi
    eps_steps = 1
    
    do j=1,lin_steps
    
        lin = lin_ini + j*(lin_fin - lin_ini)/float(lin_steps)

        do s=1,eps_steps

            eps = eps_ini + s*(eps_fin - eps_ini)/float(eps_steps)
            !epsc = 2.0*pi

            do i=1,nosc1
                !call random_number(aux)
                !y(i) = eps*aux
                !y(i) = 0.0
                y(i) = y0m1(i)
            end do
        
            t = 0.0
                
            do i=1,500

                flag = 1
        
                1 t_1 = t + 0.01*(i-1)
                t_2 = t + 0.01*i
            
                call rkf45( F, nosc1, y, yp, t_1, t_2, relerr, abserr, flag)
                if (flag.eq.4) go to 1
                
                call para_ord(y,r)
                write (2, *) t_2, r
        
            end do

            call para_ord(y,r)

            if (r.le.0.9) then
                epsc = eps
                exit
            end if

        end do

        !write (1, *) lin, epsc
            
        print*, lin, epsc
        
    end do

    close (1)
    !close (2)

end program main







subroutine F(t,y,yp)

    use r_lambda

    implicit none

    integer, parameter :: rk = kind ( 1.0D+00 )
    real ( kind = rk ) t, y(nosc1), yp(nosc1), soma
    integer i, j

    call faz_nada(t)

    yp = omega

    do i=1,nosc1

        soma = 0.0

        do j=1,nosc1

            if ( A11(i,j) /= 0 ) then

                soma = soma + sin(y(j) - y(i)) * lin / K(focal,focal)

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
    real ( kind = rk ), intent(in) :: y(nosc1)
    real ( kind = rk ), intent(out) :: r
    complex ( kind = rk ) imaginario, zcplx(nosc1), z

    imaginario = (0.0,1.0)

    do i=1,nosc1
        zcplx(i) = exp(imaginario*y(i))
    end do

    z = sum(zcplx)/nosc1
    r = abs(z)

    return

end subroutine