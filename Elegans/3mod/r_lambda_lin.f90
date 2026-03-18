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
    integer flag, i, j, lin_steps
    real ( kind = rk ) :: y(248), yp(248)
    real ( kind = rk ) t, t_1, t_2, eps, &
                       relerr, abserr, aux, &
                       w1, w2, w3, &
                       lin_ini, lin_fin, &
                       r1, r2, r3, &
                       r
    
    nmod = 3
    nosc1 = 130
    nosc2 = 77
    nosc3 = 41
    nosc = 248

    open(unit=10, file="A.txt", status="old", action="read")                 
    do i = 1, nosc
        read(10,*) (A(j, i), j=1, nosc)
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

    !open ( unit = 1, file = 'r_lambda.dat', status = 'unknown')
    !write (1, *) "l ", "lin ", "r ", "mean_r ", "sigma_r "

    open ( unit = 2, file = 'rs_lin_tempo.dat', status = 'unknown')
    !write (2, *) "t ", "r "

    open ( unit = 3, file = 'cond_init_m1.dat', status = 'unknown')
    open ( unit = 4, file = 'cond_init_m2.dat', status = 'unknown')
    open ( unit = 5, file = 'cond_init_m3.dat', status = 'unknown')

    pi = acos(-1.0)
    pi2 = 2.0*pi
    abserr = sqrt ( epsilon ( abserr ) )
    relerr = sqrt ( epsilon ( relerr ) )
    
    w1 = 0.0
    w2 = 1.0
    w3 = 6.0

    do i=1,nosc1
        omega(m1(i)+1) = w1
    end do
    do i=1,nosc2
        omega(m2(i)+1) = w2
    end do
    do i=1,nosc3
        omega(m3(i)+1) = w3
    end do

    lin_ini = 1.0
    lin_fin = 20.0
    lin_steps = 19

    l = 2.0

    eps = pi
    
    do i=1,nosc
        call random_number(aux)
        y(i) = eps*aux
        !y(i) = 0.0
    end do

    do i=1,nosc1
        write(3, *) y(m1(i)+1)
    end do
    do i=1,nosc2
        write(4, *) y(m2(i)+1)
    end do
    do i=1,nosc3
        write(5, *) y(m3(i)+1)
    end do
    close(3)
    close(4)
    close(5)

    t = 0.0

    do j=0,lin_steps
    
        lin = lin_ini + j*(lin_fin - lin_ini)/float(lin_steps)
            
        do i=1,400

            flag = 1
    
            1 t_1 = t + 0.1*(i-1)
            t_2 = t + 0.1*i
        
            call rkf45( F, nosc, y, yp, t_1, t_2, relerr, abserr, flag)
            if (flag.eq.4) go to 1
            
            !call para_ord(y,r)
            call para_ord_full(y,r,r1,r2,r3)
            
            !r_t(i) = r
            
            write (2, *) t_2, r, r1, r2, r3, lin
            !write (2, *) t_2, r
    
        end do

        t = t_2
        
        print*, lin
        
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

                    soma = soma + sin(y(j) - y(i)) * lin / K(mod_i,mod_i)

                else

                    soma = soma + sin(y(j) - y(i)) * l / (3.0 * K(mod_j,mod_i))

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
subroutine para_ord_full(y,r,r1,r2,r3)

    use r_lambda

    implicit none

    integer, parameter :: rk = kind ( 1.0D+00 )
    integer i
    real ( kind = rk ), intent(in) :: y(248)
    real ( kind = rk ), intent(out) :: r1, r2, r3, r
    complex ( kind = rk ) imaginario, zcplx(nosc), z1, z2, z3, z, &
                          z1cplx(nosc1), z2cplx(nosc2), z3cplx(nosc3)

    imaginario = (0.0,1.0)

    do i=1,nosc1
        z1cplx(i) = exp(imaginario*y(m1(i)+1))
    end do
    do i=1,nosc2
        z2cplx(i) = exp(imaginario*y(m2(i)+1))
    end do
    do i=1,nosc3
        z3cplx(i) = exp(imaginario*y(m3(i)+1))
    end do
    do i=1,nosc
        zcplx(i) = exp(imaginario*y(i))
    end do

    z1 = sum(z1cplx)/nosc1
    r1 = abs(z1)
    z2 = sum(z2cplx)/nosc2
    r2 = abs(z2)
    z3 = sum(z3cplx)/nosc3
    r3 = abs(z3)
    z = sum(zcplx)/nosc
    r = abs(z)

    return

end subroutine