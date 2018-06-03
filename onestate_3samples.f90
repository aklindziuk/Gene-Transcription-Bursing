program twostate_ssa

implicit none

integer :: trials, n, m, i, j, k, dim
real :: time, t_max, random1, random2, a, b, alpha_0, alpha_1, alpha_2, alpha_3, mean, var, sum, f_on, tau, t1, t0, t2, tt
real, dimension(0:2) :: dist
real, dimension(0:2) :: dist_norm
real, dimension(0:2) :: histogram
real :: start, finish


call cpu_time(start)

! I
a = 2
b = 10

trials = 100
t_max = 100000
dim = 2

call random_seed

dist_norm = 0

dist = 0
sum = 0
time = 0.0
t0 = 0
t1 = 0
t2 = 0

	n = 1
	time = 0.0

		do while (time < t_max) 
		
		call random_number (random1)
		call random_number (random2)
		
		if ( n == 0) then

			!update time
			tau = 1/a*log(1/random1)
			time = time + tau
			if (time > t_max*.20) then
				t0 = t0 + tau
			end if
			
			!update n
			n = n + 1

		else if ( n == 1 ) then 

			alpha_0 = a + b

			!update time
			tau = (1/alpha_0)*log(1/random1)
			time = time + tau
			if ( time > t_max*.20 ) then
				t1 = t1 + tau
			end if

			!update n
			if ( random2 < a/alpha_0) then
				n = n + 1
			else
				n = n - 1
			end if

		else 

			!update time
			tau = (1/(b*n))*log(1/random1)
			time = time + tau
			if ( time > t_max*.20 ) then
				t2 = t2 + tau
			end if

			!update n
			n = n - 1

		end if

		if (k == trials) then
			write(2,*) time, n
		end if

		end do




tt = t0 + t1 + t2
print*, t0/(tt)
print*, t1/(tt)
print*, t2/(tt)

print*, "simulation avg", t1/tt*1 + t2/tt*2



call cpu_time(finish)
!print *,"simul time =", finish-start


end program twostate_ssa


