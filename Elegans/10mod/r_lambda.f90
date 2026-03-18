module r_lambda

    implicit none
    
    real(8), save :: omega(248)
    integer, save :: modules(248)
    real(8), save :: l, lin, pi, pi2
    integer, save :: nosc
    real, dimension(248,248) :: A
    real, dimension(10,10) :: K
    
end module r_lambda







program main

    use r_lambda

    implicit none

    external F
    integer, parameter :: rk = kind ( 1.0D+00 )
    integer flag, i, j, l_steps, n_pts, mod_i
    real ( kind = rk ) :: y(248), yp(248)
    real ( kind = rk ), allocatable :: r_t(:)
    real ( kind = rk ) t, t_1, t_2, &
                       relerr, abserr, aux, &
                       w1, w2, w3, w4, w5, &
                       w6, w7, w8, w9, w10, &
                       l_ini, l_fin, &
                       r, mean_r, soma_r, soma2_r, sigmar
    
    open(unit=10, file="A.txt", status="old", action="read")                 
    do i = 1, 248
        read(10,*) (A(j, i), j=1, 248)
    end do                
    close(10)

    open(unit=10, file="K.txt", status="old", action="read")                 
    do i = 1, 10
        read(10,*) (K(j, i), j=1, 10)
    end do                
    close(10)

    open(unit=10, file="module.txt", status="old", action="read")                 
    do i = 1, 248
        read(10,*) modules(i)
    end do                
    close(10)
      
    open ( unit = 1, file = 'r_lambda.dat', status = 'unknown')
    !write (1, *) "l ", "lin ", "r ", "mean_r ", "sigma_r "

    !open ( unit = 2, file = 'rs_lin_tempo.dat', status = 'unknown')
    !write (2, *) "t ", "r "

    pi = acos(-1.0)
    pi2 = 2.0*pi
    abserr = sqrt ( epsilon ( abserr ) )
    relerr = sqrt ( epsilon ( relerr ) )

    nosc = 248
    
    allocate(r_t(5000))
    n_pts = 5000
    
    w1 = 0.0
    w2 = 1.0
    w3 = -3.5
    w4 = -2.0
    w5 = 2.0
    w6 = 0.6
    w7 = 0.4
    w8 = 1.2
    w9 = 1.5
    w10 = -1.5

    do i=1,nosc
        mod_i = modules(i)
        if (mod_i == 1) then
            omega(i) = w1
        else if (mod_i == 2) then
            omega(i) = w2
        else if (mod_i == 3) then
            omega(i) = w3
        else if (mod_i == 4) then
            omega(i) = w4
        else if (mod_i == 5) then
            omega(i) = w5
        else if (mod_i == 6) then
            omega(i) = w6
        else if (mod_i == 7) then
            omega(i) = w7
        else if (mod_i == 8) then
            omega(i) = w8
        else if (mod_i == 9) then
            omega(i) = w9
        else if (mod_i == 10) then
            omega(i) = w10
        end if
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
    
            1 t_1 = t + 0.01*(i-1)
            t_2 = t + 0.01*i
        
            call rkf45( F, nosc, y, yp, t_1, t_2, relerr, abserr, flag)
            if (flag.eq.4) go to 1
    
        end do

        t = t_2

        do i=1,5000

            flag = 1
    
            2 t_1 = t + 0.01*(i-1)
            t_2 = t + 0.01*i
        
            call rkf45( F, nosc, y, yp, t_1, t_2, relerr, abserr, flag)
            if (flag.eq.4) go to 2
    
            call para_ord(y,r)

            r_t(i) = r

        end do
        
        soma_r = sum(r_t)
        mean_r = soma_r/n_pts
        soma2_r = dot_product(r_t,r_t)
        sigmar = sqrt((n_pts*soma2_r-soma_r*soma_r)/(n_pts*(n_pts-1)))
        
        write (1, *) l, lin, r, mean_r, sigmar
        
        print*, l
        
    end do

    close (1)
    !close (2)

    call semente(2)

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

        soma = 0.0

        mod_i = modules(i)

        do j=1,nosc

            if ( A(i,j) /= 0.0 ) then

                mod_j = modules(j)

                if ( mod_i == mod_j ) then

                    soma = soma + sin(y(j) - y(i)) * A(i,j) * lin / K(mod_i,mod_i)

                else

                    soma = soma + sin(y(j) - y(i)) * A(i,j) * l / (10.0 * K(mod_j,mod_i))

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