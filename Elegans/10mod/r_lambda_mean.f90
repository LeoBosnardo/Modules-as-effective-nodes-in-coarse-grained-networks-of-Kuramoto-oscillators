module r_lambda

    implicit none
    
    real(8), save :: omega(248)
    integer, save :: m1(68), m2(33), m3(76), m4(11), m5(13), &
                     m6(7), m7(6), m8(8), m9(8), m10(18), modules(248)
    real(8), save :: l, lin, pi, pi2
    integer, save :: nosc, nosc1, nosc2, nosc3, nosc4, nosc5, &
                     nosc6, nosc7, nosc8, nosc9, nosc10, nmod
    real, dimension(248,248) :: A
    real, dimension(10,10) :: K
    
end module r_lambda







program main

    use r_lambda

    implicit none

    external F
    integer, parameter :: rk = kind ( 1.0D+00 )
    integer flag, i, j, u, lin_steps, l_steps, n_pts
    real ( kind = rk ) :: y(248), yp(248)
    real ( kind = rk ) :: r1_t(4000), r2_t(4000), r3_t(4000), r4_t(4000), r5_t(4000), &
                          r6_t(4000), r7_t(4000), r8_t(4000), r9_t(4000), r10_t(4000)
    real ( kind = rk ) t, t_1, t_2, eps, &
                       relerr, abserr, aux, &
                       w1, w2, w3, w4, w5, &
                       w6, w7, w8, w9, w10, &
                       lin_ini, lin_fin, &
                       r1, r2, r3, r4, r5, &
                       r6, r7, r8, r9, r10, &
                       r, l_ini, l_fin, &
                       mean_r1, mean_r2, mean_r3, mean_r4, mean_r5, &
                       mean_r6, mean_r7, mean_r8, mean_r9, mean_r10
    
    nmod = 10
    nosc1 = 68
    nosc2 = 33
    nosc3 = 76
    nosc4 = 11
    nosc5 = 13
    nosc6 = 7
    nosc7 = 6
    nosc8 = 8
    nosc9 = 8
    nosc10 = 18
    nosc = 248

    n_pts = 4000

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
    
    open(unit=10, file="m4.txt", status="old", action="read")                 
    do i = 1, nosc4
        read(10,*) m4(i)
    end do                
    close(10)
      
    open(unit=10, file="m5.txt", status="old", action="read")                 
    do i = 1, nosc5
        read(10,*) m5(i)
    end do                
    close(10)

    open(unit=10, file="m6.txt", status="old", action="read")                 
    do i = 1, nosc6
        read(10,*) m6(i)
    end do                
    close(10)
      
    open(unit=10, file="m7.txt", status="old", action="read")                 
    do i = 1, nosc7
        read(10,*) m7(i)
    end do                
    close(10)
      
    open(unit=10, file="m8.txt", status="old", action="read")                 
    do i = 1, nosc8
        read(10,*) m8(i)
    end do                
    close(10)
    
    open(unit=10, file="m9.txt", status="old", action="read")                 
    do i = 1, nosc9
        read(10,*) m9(i)
    end do                
    close(10)
      
    open(unit=10, file="m10.txt", status="old", action="read")                 
    do i = 1, nosc10
        read(10,*) m10(i)
    end do                
    close(10)

    open(unit=10, file="module.txt", status="old", action="read")                 
    do i = 1, nosc
        read(10,*) modules(i)
    end do                
    close(10)

    open ( unit = 1, file = 'r_lambda_mean.dat', status = 'unknown')

    pi = acos(-1.0)
    pi2 = 2.0*pi
    abserr = sqrt ( epsilon ( abserr ) )
    relerr = sqrt ( epsilon ( relerr ) )
    
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

    do i=1,nosc1
        omega(m1(i)+1) = w1
    end do
    do i=1,nosc2
        omega(m2(i)+1) = w2
    end do
    do i=1,nosc3
        omega(m3(i)+1) = w3
    end do
    do i=1,nosc4
        omega(m4(i)+1) = w4
    end do
    do i=1,nosc5
        omega(m5(i)+1) = w5
    end do
    do i=1,nosc6
        omega(m6(i)+1) = w6
    end do
    do i=1,nosc7
        omega(m7(i)+1) = w7
    end do
    do i=1,nosc8
        omega(m8(i)+1) = w8
    end do
    do i=1,nosc9
        omega(m9(i)+1) = w9
    end do
    do i=1,nosc10
        omega(m10(i)+1) = w10
    end do

    lin_ini = 20.0
    lin_fin = 20.0
    lin_steps = 1

    l_ini = 0.0
    l_fin = 10.0
    l_steps = 50

    eps = pi

    do u=0,l_steps

        l = l_ini + u*(l_fin - l_ini)/float(l_steps)

        t = 0.0

        do i=1,nosc
            call random_number(aux)
            y(i) = eps*aux
            !y(i) = 0.0
        end do

        do j=1,lin_steps
        
            lin = lin_ini + j*(lin_fin - lin_ini)/float(lin_steps)
                
            do i=1,400

                flag = 1
        
                1 t_1 = t + 0.1*(i-1)
                t_2 = t + 0.1*i
            
                call rkf45( F, nosc, y, yp, t_1, t_2, relerr, abserr, flag)
                if (flag.eq.4) go to 1
        
            end do

            t = t_2
            
        end do

        do i=1,4000

            flag = 1
    
            2 t_1 = t + 0.01*(i-1)
            t_2 = t + 0.01*i
        
            call rkf45( F, nosc, y, yp, t_1, t_2, relerr, abserr, flag)
            if (flag.eq.4) go to 2
            
            call para_ord_full(y,r,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10)
            
            r1_t(i) = r1
            r2_t(i) = r2
            r3_t(i) = r3
            r4_t(i) = r4
            r5_t(i) = r5
            r6_t(i) = r6
            r7_t(i) = r7
            r8_t(i) = r8
            r9_t(i) = r9
            r10_t(i) = r10
    
        end do

        mean_r1 = sum(r1_t)/n_pts
        mean_r2 = sum(r2_t)/n_pts
        mean_r3 = sum(r3_t)/n_pts
        mean_r4 = sum(r4_t)/n_pts
        mean_r5 = sum(r5_t)/n_pts
        mean_r6 = sum(r6_t)/n_pts
        mean_r7 = sum(r7_t)/n_pts
        mean_r8 = sum(r8_t)/n_pts
        mean_r9 = sum(r9_t)/n_pts
        mean_r10 = sum(r10_t)/n_pts

        write (1, *) l, mean_r1, mean_r2, mean_r3, mean_r4, mean_r5, mean_r6, mean_r7, mean_r8, mean_r9, mean_r10 
        
        print*, l, mean_r1, mean_r2, mean_r3, mean_r4, mean_r5, mean_r6, mean_r7, mean_r8, mean_r9, mean_r10 

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







subroutine para_ord_full(y,r,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10)

    use r_lambda

    implicit none

    integer, parameter :: rk = kind ( 1.0D+00 )
    integer i
    real ( kind = rk ), intent(in) :: y(248)
    real ( kind = rk ), intent(out) :: r1, r2, r3, r4, r5, &
                                       r6, r7, r8, r9, r10, r
    complex ( kind = rk ) imaginario, zcplx(nosc), z1, z2, z3, z4, z5, z, &
                          z6, z7, z8, z9, z10, &
                          z1cplx(nosc1), z2cplx(nosc2), z3cplx(nosc3), z4cplx(nosc4), z5cplx(nosc5), &
                          z6cplx(nosc6), z7cplx(nosc7), z8cplx(nosc8), z9cplx(nosc9), z10cplx(nosc10)

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
    do i=1,nosc4
        z4cplx(i) = exp(imaginario*y(m4(i)+1))
    end do
    do i=1,nosc5
        z5cplx(i) = exp(imaginario*y(m5(i)+1))
    end do
    do i=1,nosc6
        z6cplx(i) = exp(imaginario*y(m6(i)+1))
    end do
    do i=1,nosc7
        z7cplx(i) = exp(imaginario*y(m7(i)+1))
    end do
    do i=1,nosc8
        z8cplx(i) = exp(imaginario*y(m8(i)+1))
    end do
    do i=1,nosc9
        z9cplx(i) = exp(imaginario*y(m9(i)+1))
    end do
    do i=1,nosc10
        z10cplx(i) = exp(imaginario*y(m10(i)+1))
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
    z4 = sum(z4cplx)/nosc4
    r4 = abs(z4)
    z5 = sum(z5cplx)/nosc5
    r5 = abs(z5)
    z6 = sum(z6cplx)/nosc6
    r6 = abs(z6)
    z7 = sum(z7cplx)/nosc7
    r7 = abs(z7)
    z8 = sum(z8cplx)/nosc8
    r8 = abs(z8)
    z9 = sum(z9cplx)/nosc9
    r9 = abs(z9)
    z10 = sum(z10cplx)/nosc10
    r10 = abs(z10)
    z = sum(zcplx)/nosc
    r = abs(z)

    return

end subroutine