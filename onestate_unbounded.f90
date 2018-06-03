program twostate_ssa

implicit none

integer :: trials, n, m, i, j, k, dim
real :: time, t_max, random1, random2, a, b, alpha_0, alpha_1, alpha_2, alpha_3, mean, var, sum, f_on, tau, t1, t0, t2, tt
real, dimension(0:250) :: dist
real, dimension(0:250) :: dist_norm
real :: start, finish


call cpu_time(start)

! I
a = .9
b = 4.2

t_max = 10000
dim = 250

call random_seed

dist_norm = 0
dist = 0

	n = 0
	time = 0.0

	open(2, file = "timeseries.dat", status = "replace")

		do while (time < t_max) 
		
		call random_number (random1)
		call random_number (random2)
		
		if ( n == 0) then

			!update time
			tau = 1/a*log(1/random1)
			if (time > t_max*.20) then
				dist(n) = dist(n) + tau
			end if
			time = time + tau
			
			!update n
			n = n + 1

		else 

			alpha_0 = a + b*n

			!update time
			tau = (1/alpha_0)*log(1/random1)
			if ( time > t_max*.20 ) then
				dist(n) = dist(n) + tau
			end if

			time = time + tau

			!update n
			if ( random2 < a/alpha_0) then
				n = n + 1
			else
				n = n - 1
			end if

		end if

		write(2,*) time, n

		end do


sum= 0
mean = 0


 do i = 0, dim
 	!print*, i, ":", dist(i)
 	sum = sum + dist(i)
 end do
 
 print*, "total time =", sum

 do i = 0, dim
 	dist_norm(i) = real(dist(i))/sum
 	!print*, i, ":", dist_norm(i)
 	mean = mean + i*dist_norm(i)
 end do
print*, "sumulations mean =", mean


print*, "calculated avg", a/b


call cpu_time(finish)
!print *,"simul time =", finish-start

open(3, file = "dist.dat", status = "replace")
do i = 0,dim
		write(3, *) i, dist(i)
end do
close(3)

open(1, file = "dist_norm.dat", status = "replace")
do i = 0,dim
		write(1, *) i, dist_norm(i)
end do
close(3)


end program twostate_ssa


