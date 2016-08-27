program heatsolve
    implicit none
    double precision :: k,Q,L,A,B,C,D,E,F
    integer :: n
    double precision, allocatable, dimension(:) :: T
    print *, "Enter K, Q, L, n, A, B, C, D, E, F"
    read (*,*), K, Q, L, n, A, B, C, D, E, F
    allocate(T(n))
    call solve_heat_eqn(K, Q, L, n, A, B, C, D, E, F, T)
    call print_array(T,n)
end program heatsolve

subroutine solve_heat_eqn(K, Q, L, n, A, B, C, D, E, F, T)
    implicit none
    double precision :: delta,k,Q,L,A,B,C,D,E,F
    integer :: i,n
    double precision :: temp, T(n)
    double precision, allocatable, dimension(:) :: sub_d, main_d, super_d, b_cons
    allocate(sub_d(n))
    allocate(main_d(n))
    allocate(super_d(n))
    allocate(b_cons(n))
    delta = L/n
    do i=1,n
        if (i==1) then
            sub_d(i)=0.0
            main_d(i)=A-(B/delta)
            super_d(i) = B/delta
            b_cons(i)=C
        else if (i==n) then
            sub_d(i) =-E/delta
            main_d(i) =D+(E/delta)
            b_cons(i) = F
        else
            sub_d(i) = 1/(delta*delta)
            main_d(i) = -2/(delta*delta)
            super_d(i) = 1/(delta*delta)
            b_cons(i) = -Q/k
        end if
    end do
    ! Forward Elimination
    do i=2,n
        temp=sub_d(i)/main_d(i-1)
        main_d(i) = main_d(i)-temp*super_d(i-1)
        b_cons(i)= b_cons(i) - temp*b_cons(i-1)
    end do
    ! Backwards Substitution
    T(n) = b_cons(n)/main_d(n)
    do i=(n-1),1,-1
        T(i)=(b_cons(i)-super_d(i)*T(i+1))/main_d(i)
    end do

end subroutine solve_heat_eqn

subroutine print_array(x,n)
    implicit none
    double precision, intent(in) :: x(n)
    integer, intent(in) :: n
    integer :: i
    write(*, "(a)", advance='no') "["
    do i = 1, n
        write(*, "(f8.4 )", advance = 'no') x(i)
        if (i .eq. n) write(*, "(a)", advance='no') " ]"
        if (modulo(i,5) .eq. 0) write(*, "(a)") " "
        if (modulo(i,5) .eq. 0) write(*, "(a)", advance='no') " "
    end do
end subroutine print_array


