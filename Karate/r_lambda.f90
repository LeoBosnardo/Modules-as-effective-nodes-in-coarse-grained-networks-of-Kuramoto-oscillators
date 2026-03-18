module r_lambda

    implicit none
    
    real(8), save :: omega(34)
    integer, save :: mr_hi(17), officer(17)
    real(8), save :: k(4)
    real(8), save :: l, lin, &
                     pi, pi2, &
                     over_k11, over_2k12, over_2k21, over_k22
    integer, save :: nosc = 34, nosc1 = 17, nosc2 = 17
    integer, dimension(34,34) :: A
    
end module r_lambda







program main

    use r_lambda

    implicit none

    external F
    integer flag, i, j, l_steps, n_pts
    real (8) :: y(34), yp(34)
    real (8), allocatable :: r1_t(:), r2_t(:)
    real (8), allocatable :: r1_t_sd(:), r2_t_sd(:), r_t_sd(:)
    real (8) t, t_1, t_2, t_stop, &
                       relerr, abserr, aux, &
                       w, &
                       k11, k12, k21, k22, &
                       l_ini, l_fin, &
                       r, r1, r2, &
                       psi1, psi2, psi, phi, &
                       mean_r, mean_r1, mean_r2, &
                       soma_r, soma_r1, soma_r2, &
                       soma2_r, soma2_r1, soma2_r2, &
                       sigmar, sigmar1, sigmar2

    pi = acos(-1.0)
    pi2 = 2.0*pi
    abserr = sqrt ( epsilon ( abserr ) )
    relerr = sqrt ( epsilon ( relerr ) )

    open(unit=10, file="A.txt", status="old", action="read")                 
    do i = 1, 34
        read(10,*) (A(j, i), j=1, 34)
    end do                
    close(10)

    open(unit=10, file="k.txt", status="old", action="read")                 
    do i = 1, 4
        read(10,*) k(i)
    end do                
    close(10)
    k11 = k(1)
    k12 = k(2)
    k21 = k(3)
    k22 = k(4)
    over_k11 = 1.0 / k11
    over_2k12 = 1.0 / (2*k12)
    over_2k21 = 1.0 / (2*k21)
    over_k22 = 1.0 / k22

    open(unit=10, file="mr_hi.txt", status="old", action="read")                 
    do i = 1, nosc1
        read(10,*) mr_hi(i)
    end do                
    close(10)

    open(unit=10, file="officer.txt", status="old", action="read")                 
    do i = 1, nosc2
        read(10,*) officer(i)
    end do                
    close(10)

    do i=1,17
        mr_hi(i) = mr_hi(i) + 1
        officer(i) = officer(i) + 1
    end do

    open ( unit = 1, file = 'r_lambda.dat', status = 'unknown')
    !write (1, *) "l ", "lin ", "meanr1 ", "meanr2 ", "meanr ", "sigmar1 ", "sigmar2 ", "sigmar ", "r1 ", "r2 ", "r ", "psi1 ", &
    !"psi2 ", "psi ", "phi " 

    !open ( unit = 2, file = 'rs_lin_tempo.dat', status = 'unknown')
    !write (2, *) "t ", "r1 ", "r2 ", "r ", "lin"
    
    allocate(r1_t(500), r2_t(500))
    allocate(r1_t_sd(10000), r2_t_sd(10000), r_t_sd(10000))
    n_pts = 10000
    
    w = 1.5

    omega = 0.0
    do i=1,nosc1
        omega(mr_hi(i)) = 1.5
    end do
    
    l_ini = 0.0
    l_fin = 3.0
    l_steps = 100

    t_stop = 10000

    do j=0,l_steps
    
        l = l_ini + j*(l_fin - l_ini)/float(l_steps)

        !lin = 0.0
        lin = 20.0

        do i=1,nosc
            call random_number(aux)
            y(i) = aux*pi
        end do
    
        t = 0.0
        
        r1_t = 0.0
        r2_t = 0.0
        r1_t_sd = 0.0
        r2_t_sd = 0.0
        r_t_sd = 0.0
        
        !do while ((minval(r1_t)<0.9).or.(minval(r2_t)<0.9))
            
            do i=0,499
            
                flag = 1
                
                1 t_1 = t + 0.1*i
                t_2 = t + 0.1*(i+1)
                
                call rkf45( F, nosc, y, yp, t_1, t_2, relerr, abserr, flag)
                
        !        call para_ord(y,r1,r2,r,psi1,psi2,psi,phi)
                
        !        r1_t(i+1) = r1
        !        r2_t(i+1) = r2
                
                if (flag.eq.4) go to 1
                
                !write (2, *) t_2, r1, r2, r, lin
                
            end do
            
            t = t_2
            
        !    if (t_2>t_stop) go to 3
            
        !    lin = lin + 0.1
            
        !end do
        
        !lin = lin - 0.1
        
        do i=1,10000
        
            flag = 1
    
            2 t_1 = t + 0.01*i
            t_2 = t + 0.01*(i+1)
        
            call rkf45( F, nosc, y, yp, t_1, t_2, relerr, abserr, flag)
            
            call para_ord(y,r1,r2,r,psi1,psi2,psi,phi)
            
            r1_t_sd(i) = r1
            r2_t_sd(i) = r2
            r_t_sd(i) = r
            
            if (flag.eq.4) go to 2
            
            !write (2, *) t_2, r1, r2, r, lin
    
        end do
        
        soma_r1 = sum(r1_t_sd)
        mean_r1 = soma_r1/n_pts
        soma2_r1 = dot_product(r1_t_sd,r1_t_sd)
        sigmar1 = sqrt((n_pts*soma2_r1-soma_r1*soma_r1)/(n_pts*(n_pts-1)))
        
        soma_r2 = sum(r2_t_sd)
        mean_r2 = soma_r2/n_pts
        soma2_r2 = dot_product(r2_t_sd,r2_t_sd)
        sigmar2 = sqrt((n_pts*soma2_r2-soma_r2*soma_r2)/(n_pts*(n_pts-1)))
            
        soma_r = sum(r_t_sd)
        mean_r = soma_r/n_pts
        soma2_r = dot_product(r_t_sd,r_t_sd)
        sigmar = sqrt((n_pts*soma2_r-soma_r*soma_r)/(n_pts*(n_pts-1)))
        
        write (1, *) l, lin, mean_r1, mean_r2, mean_r, sigmar1, sigmar2, sigmar, r1, r2, r, psi1, psi2, psi, phi
        
        print*, l, lin, r1, r2, r
        
    end do

    close (1)
    !close (2)

    call semente(2)

end program main







subroutine F(t, y, yp)

    use r_lambda
    
    implicit none

    real (8) t, y(34), yp(34), soma
    integer i, j

    call faz_nada(t)

    yp = omega

    do i = 1, nosc1

        soma = 0.0
        
        do j = 1, nosc1
            soma = soma + sin(y(mr_hi(j)) - y(mr_hi(i))) * A(mr_hi(i), mr_hi(j))
        end do

        yp(mr_hi(i)) = yp(mr_hi(i)) + lin * soma * over_k11

        soma = 0.0

        do j = 1, nosc2
            soma = soma + sin(y(officer(j)) - y(mr_hi(i))) * A(mr_hi(i), officer(j))
        end do

        yp(mr_hi(i)) = yp(mr_hi(i)) + l * soma * over_2k12

    end do

    do i = 1, nosc2

        soma = 0.0
        
        do j = 1, nosc1
            soma = soma + sin(y(mr_hi(j)) - y(officer(i))) * A(officer(i), mr_hi(j))
        end do
        
        yp(officer(i)) = yp(officer(i)) + l * soma * over_2k21

        soma = 0.0
        
        do j = 1, nosc2
            soma = soma + sin(y(officer(j)) - y(officer(i))) * A(officer(i), officer(j))
        end do
        
        yp(officer(i)) = yp(officer(i)) + lin * soma * over_k22
    
    end do

    return

end subroutine







subroutine faz_nada(t)

    implicit none

    real (8) t

    if (t.ne.t) then
        write (*,*) 't is NAN'
    end if

    return

end subroutine






subroutine para_ord(y, r1, r2, r, psi1, psi2, psi, phi)
    
    use r_lambda

    implicit none

    integer i
    real (8), intent(in) :: y(34)
    real (8), intent(out) :: r1, r2, r, psi1, psi2, psi, phi
    complex (8) imaginario, z1cplx(nosc1), z2cplx(nosc2), &
                          zcplx(nosc), z1, z2, z

    imaginario = (0.0,1.0)
    zcplx = (0.0, 0.0)

    do i = 1, size(mr_hi)
        z1cplx(i) = exp(imaginario * y(mr_hi(i)))
        zcplx(mr_hi(i)) = z1cplx(i)
    end do

    do i = 1, size(officer)
        z2cplx(i) = exp(imaginario * y(officer(i)))
        zcplx(officer(i)) = z2cplx(i)
    end do

    z1 = sum(z1cplx) / nosc1
    z2 = sum(z2cplx) / nosc2
    z  = sum(zcplx)  / nosc

    r1 = abs(z1)
    r2 = abs(z2)
    r  = abs(z)
    psi1 = atan2(aimag(z1), real(z1))
    psi2 = atan2(aimag(z2), real(z2))
    psi  = atan2(aimag(z),  real(z))
    phi  = atan2(sin(psi1 - psi2), cos(psi1 - psi2))

    return

end subroutine








subroutine semente(i)

    implicit none

    integer i, iseed(8)

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