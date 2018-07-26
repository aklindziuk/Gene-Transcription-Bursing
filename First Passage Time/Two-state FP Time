!------------------------------------------------------------------
! Building on the one-state Gillespie simulation for transcription, this 
! code determines the first passage time to reach a given target. The first
! passage time (fp) is the time it takes to reach a target number of mRNAs, 
! designated by n, for the first time. a, b are the kinetic parameters
! of the system. target  is the target number of mRNAs. trails is the number
! of times the system is launeched (aka how many times fp time is collected)
! The output files of interest are fptime.dat, which is a record of the fp
! time recorded during each trial, and the fphist.dat, which gives a histogram 
! of the recorded fp values (there might be a small error (slight shift) in the 
! histogram. The command line output gives mean and variance of the the fp time - 
! these calculations re accurate. 
!-------------------------------------------------------------------


program twostate_ssa

implicit none

integer :: trials, n, m, i, j, k, dim
real :: time, t_max, random1, random2, a, b, alpha_0, alpha_1, alpha_2, alpha_3, mean, var, sum, f_on, tau, t1, t0, t2, tt
real, dimension(:), allocatable :: dist
real, dimension(:), allocatable :: dist_norm
real :: start, finish

real, dimension(:),allocatable :: fp_time
real, dimension(1:100) :: fp_hist
real :: fpmean, fpvar, fpsum, min, max, nbin, binsize
integer :: target, counter

call cpu_time(start)

! --------------select parameters -----------------------------
a = 100
b = 1

t_max = 1000
target = 130

dim = 250 ! set to the size of dist (max poss n of system)
trials = 500 ! how many times fp time is sampled

! ----------------------------------------------------------------

open(2, file = "timeseries.dat", status = "replace")

allocate(dist(0:dim))
allocate(dist_norm(0:dim))
allocate(fp_time(1:trials))
fp_time = 0
fp_hist = 0


call random_seed

do i = 1, size(fp_time)

dist_norm = 0
dist = 0

	n = 0
	time = 0.0
	counter = 1

		do while (time < t_max) 
		
		call random_number (random1)
		call random_number (random2)

		if (random1 == 0) then
			call random_number (random1)
		end if

		if ( n == 0) then

			!update time
			tau = - 1/a*log(random1)
			if (time > t_max*.20) then
				dist(n) = dist(n) + tau
			end if
			time = time + tau
			
			!update n
			n = n + 1

		else 

			alpha_0 = a + b*n

			!update time
			tau = -(1/alpha_0)*log(random1)
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

		if ( n == target ) then
			if (counter == 1) then
				!print*, "target found after", i, time, "time units"
				fp_time(i) = time
			end if
			counter = counter + 1
		end if

	if (i == 1) then
		write(2,*), time, n
	end if

	end do
	close(2)

	if (counter == 1) then
		print*, "target not reached in ", t_max, "time units" 
		fp_time(i) = t_max
	end if


end do !---------------end loop-------------------

sum = 0
mean = 0

open(3, file = "dist.dat", status = "replace")
do i = 0,dim
		write(3, *), i, dist(i)
end do
close(3)

open(4, file = "fptime.dat", status = "replace")
do i = 1,size(fp_time)
		write(4, *), i, fp_time(i)
end do
close(4)

!-------------fisrt passage mean ------------------------------

fpmean = 0
do i = 1, size(fp_time)
	fpmean = fpmean + fp_time(i)
end do
fpmean = fpmean/size(fp_time)
print*, 'fp mean', fpmean

!--------------calculate histogram------------------------------

min = minval(fp_time)
max = maxval(fp_time)
!print*, "min", min
!print*, "max", max
binsize = (max - min)/real(size(fp_hist))
!print*, 'binsize', binsize

do i = 1,size(fp_time)
do j = 1,size(fp_hist)
	if ((fp_time(i) .ge. (min + binsize*(j-1))) .and. (fp_time(i) .le. (min + j*binsize))) then
		fp_hist(j) = fp_hist(j) + 1
	end if
end do
end do

fpsum = 0
do i = 1,size(fp_hist)
	fpsum = fpsum + fp_hist(i)
end do
!print*, 'histogram sum', fpsum

open(5, file = "fphist.dat", status = "replace")
do i = 1,size(fp_hist)
		fp_hist(i) = fp_hist(i)/fpsum
		write(5, *), min + (i + .5)*binsize, fp_hist(i)
end do
close(5)

fpsum = 0
do i = 1,size(fp_hist)
	fpsum = fpsum + fp_hist(i)
end do
!print*, 'histogram sum', fpsum

!----------------calculate variance ----------------------------------

fpmean = 0
do i = 1,size(fp_hist)
		fpmean = fpmean + (min + (i + .5)*binsize)*fp_hist(i)
end do
!print*, 'fp mean', fpmean

fpvar = 0
do i = 1, size(fp_hist)
	fpvar = fpvar + fp_hist(i)*(fpmean - (min + (i + .5)*binsize))**2
end do
print*, 'fp var', fpvar
print*, " "


!----------------calculate mean -------------------------------------


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
close(1)


end program twostate_ssa
