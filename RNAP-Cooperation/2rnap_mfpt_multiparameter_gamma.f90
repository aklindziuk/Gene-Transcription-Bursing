program twostate_ssa

implicit none

integer :: n, m, i, j, k, dim, counter, c, h, nn, initail_m, t, rounds
real :: time, t_max, g, e, x,random1, random2, a, b, alpha_0, alpha_1, alpha_2, mean, mean2, var, var2, sum, tau, fon, dnt
real :: y, phi,mean_theory, fano_theory, n2_epsilon, gamma, epsilon, var_epsilon, f2_epsilon, part2, n_epsilon
real :: fano, fano_yone, fano_yinf, f2_yone, f2, q1, q2, r1, r2, n2_yone,variance_yone, n_yone, fano_yone_simplified
real :: target, fpmean, O, Q, fact1, fact2
real, dimension(1:2, 0:10000) :: dist
real, dimension(1:2, 0:10000) :: dist_norm
real, dimension(0:10000) :: histogram
real, dimension (:,:), allocatable :: fp_time 
real, dimension(1:4) :: kns
real, dimension(1:4) :: kfs
real :: start, finish
character(len=10) :: file_id
character(len=50) :: file_name
character(20) filename


call cpu_time(start)

! gamma = Inf, 100, 10, 1, .1, .01, 0


!===set 4======

kfs(1) = .1
kfs(2) = 1
kfs(3) = 10
kfs(4) = 100

kns(1) = .1
kns(2) = 1
kns(3) = 10
kns(4) = 100

!===set 3======

!kfs(1) = 1      ! Inf
!kfs(2) = 10     ! 100
!kfs(3) = 10     ! 10
!kfs(4) = 10
!kfs(5) = 10
!kfs(6) = 10
!kfs(7) = 0

!kns(1) = 0
!kns(2) = .1
!kns(3) = 1
!kns(4) = 10
!kns(5) = 100
!kns(6) = 1000
!kns(7) = 1

!====set 2=====

!kfs(1) = 1          ! Inf
!kfs(2) = 10000      ! 100
!kfs(3) = 1000       ! 10
!kfs(4) = 1000       ! 1
!kfs(5) = 100        ! 0.1
!kfs(6) = 100        ! 0.01
!kfs(7) = 0          ! 0

!kns(1) = 0
!kns(2) = 100
!kns(3) = 100
!kns(4) = 1000
!kns(5) = 1000
!kns(6) = 10000
!kns(7) = 1

!====set 1=======

!kfs(1) =  0        ! 0
!kfs(2) = .0001     ! 0.01
!kfs(3) = .001      ! 0.1
!kfs(4) = .001      ! 1
!kfs(5) = .01       ! 10
!kfs(6) = .01       ! 100
!kfs(7) = 1         ! Inf

!kns(1) = 1
!kns(2) = .01
!kns(3) = .01
!kns(4) = .001
!kns(5) = .001
!kns(6) = .0001
!kns(7) = 0

rounds = 10

a = 10
b = 1
y = 1.1
x = a/b

allocate ( fp_time(1:rounds,1:size(kns)) )  

t_max = 300000
dim = 10000

fp_time = 0

call random_seed


do t = 5, 5 !================iterate through dirrerent targets

	target = real(t)

	write(filename,'(a,i0,a)') 'target', t,'.txt'
	open(3, file = filename , status = "replace")

	write(filename,'(a,i0,a)') 'fp_time', t,'.txt'
	open(5, file = filename , status = "replace")
	write(5,*), "ys",kfs(1)/kns(1),kfs(2)/kns(2),kfs(3)/kns(3),kfs(4)/kns(4)

do h = 1, size(kfs) !========iterate through different gammas 

do c = 1, rounds !=========== iterate through 100 rounds 
	
	n = 0
	m = 1
	if (c .ge. int(rounds/2)) then 
		m = 2
	end if

	!Special cases
	if (kfs(h)/kns(h) < .001) then 
		m = 2
	end if
	if ( kfs(h)/kns(h) > 1000) then
		m = 1
	end if

	dist_norm = 0
	dist = 0

	time = 0.0
	counter = 1

	do while (time < t_max) 
		
		call random_number (random1)
		call random_number (random2)
		
		if ( n == 0) then
			if ( m == 1) then

				!!! case 1 ---------->

				phi = kns(h) + a/y
				tau = (1/phi)*log(1/random1)
				if (time > t_max*.20) then
					dist(m,n) = dist(m,n) + tau
				end if
				time = time + tau
				
				if ( random2 < (a/y)/phi) then
					n = n + 1
				else
					m = m + 1
				end if
			else  

				!!! case 2 ------------>

			    phi = kfs(h) + a + a/y
				tau = (1/phi)*log(1/random1)
				if (time > t_max*.20) then
					dist(m,n) = dist(m,n) + tau
				end if
				time = time + tau
				if (random2 < (a + a/y)/phi) then
					n = n + 1
				else
					m = m - 1
				end if
			end if
		else
			if ( m == 1) then

				!!! case 3 ---------->

				phi = kns(h) + a/(y**(n+1)) + n*b
				tau = (1/phi)*log(1/random1)
				if (time > t_max*.20) then
					dist(m,n) = dist(m,n) + tau
				end if
				time = time + tau

				if (random2 < kns(h)/phi) then
					m = m + 1 
				else if (kns(h)/phi < random2 .and. random2 < (kns(h)+a/(y**(n+1)))/phi) then
					n = n + 1
				else 
					n = n - 1 
				end if
			else

				!!! case 4 --------------->

				phi = kfs(h) + a + a/(y**(n+1)) + n*b
				tau = (1/phi)*log(1/random1)
				if (time > t_max*.20) then
					dist(m,n) = dist(m,n) + tau
				end if
				time = time + tau

				if ( random2 < kfs(h)/phi) then
					m = m - 1
				else if (kfs(h)/phi < random2 .and. random2 < (kfs(h) + a + a/(y**(n+1)))/phi) then 
					n = n + 1
				else
					n = n - 1
				end if 

			end if
		end if

	if (n == target) then
		if (counter == 1) then
			fp_time(c, h) = time
		end if
		counter = counter + 1
	end if


	if (n .eq. dim) then
		print*, "n is out of bounds"
	end if


	end do !========= one round over ============

if (counter == 1) then 
	print*, "target not reached in ", t_max, "time units"
	fp_time(c,h) = t_max
end if 

end do !======== all 100 rounds over ====================

fpmean = 0
do i = 1,rounds
	fpmean = fpmean + fp_time(i,h)
end do
fpmean = fpmean/real(rounds)
!print*, kfs(h)/kns(h), fpmean, target
write(3,*), target,",", kfs(h)/kns(h),",", fpmean

end do ! =======================finished with all kfs-===========

close(3)

do i = 1, rounds
	write(5, *), i,",",fp_time(i,1),",",fp_time(i,2),",",fp_time(i,3),",",fp_time(i,4)!,",", &
	!fp_time(i,5),",",fp_time(i,6),",",fp_time(i,7),",",fp_time(i,8),",",fp_time(i,9),",",fp_time(i,10)
end do

end do !==================finished all targets

close(5)
print*, "done"





call cpu_time(finish)




end program twostate_ssa