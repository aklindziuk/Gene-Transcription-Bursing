program twostate_ssa

implicit none

integer :: n, s, p, i, j, h, k, dimJ, dimN, values(1:8), jmax
real :: time, t_max, mean, mean2, var, var2, sum, tau, fon, dnt, record_after
real :: nn, ff, on, phi, n_inf, part2
real :: P1_simultion, P2_simultion, amean, jmean
real :: a, b, kn, kf, rn, rf, y, gamma, omega, gamma0, gammaj, n1, n2, Aterm, Bterm, phi0, phij
real :: P1, P2, Pon1, Pon2, P01, P02, Pj1, Pj2, Px1, Px2, sumPj1, sumPj2, term1, term2
real :: Pjm1_plane2, Pj_plane2, Pon_plane2, P0_plane2, Pterm_plane2, fterm_plane2, f0_plane2, fj_plane2, fjm1_plane2
real :: n_plane2, n2_plane2, fon_plane2, factor, pterm1, pterm2, den11, den22
real :: Pjm1_plane1, Pj_plane1, Pon_plane1, P0_plane1, Pterm_plane1, fterm_plane1, f01, fj_plane1, fjm1_plane1
real :: n_plane1, n2_plane1, fon1, fon12, fon2, f02, num1, num2, den1, den2, num11, num22
real :: prod1, prod2, prodm1, prodm2, sum1, sum2, fon1_den, fon1_num, fon2_num1, fon2_num2, fon2_num3, fon2_den
real, dimension(:), allocatable :: prob_dist
real, dimension(:), allocatable :: jprob_dist
real, dimension(:, :), allocatable :: nhistogram
real, dimension(:, :), allocatable :: jhistogram
real, dimension(:, :, :), allocatable :: TimeSpent
real, dimension(:, :, :), allocatable :: dist_norm
real, dimension(:, :), allocatable :: ahistogram
real, dimension(:), allocatable :: Pjhist1
real, dimension(:), allocatable :: Pjhist2
real, dimension(:), allocatable :: fjhist1
real, dimension(:), allocatable :: fjhist2
real, dimension(1:7) :: rns
real, dimension(1:7) :: rfs

character(len=10) :: file_id
character(len=50) :: file_name
character(20) filename
integer, dimension(:), allocatable :: seed
real(8) :: r1, r2
real :: start, finish

call cpu_time(start)
call date_and_time(values=values) 
call random_seed(size=i)
allocate(seed(1:i))
seed(:) = values(8)
call random_seed(put=seed)

!!! =======chose parameters ============

a = 15
b = 1
kn =10
kf =10
y = 100

!v2 
rfs(1) = 0          
rfs(2) = 50
rfs(3) = 50
rfs(4) = 50       
rfs(5) = 50     
rfs(6) = 50      
rfs(7) = 1          

rns(1) = 1                
rns(2) = 5000
rns(3) = 500
rns(4) = 50   
rns(5) = 5      
rns(6) = .5      
rns(7) = 0   

!v3
!rfs(1) = 0          
!rfs(2) = .5
!rfs(3) = 5
!rfs(4) = 50       
!rfs(5) = 500     
!rfs(6) = 5000      
!rfs(7) = 1          

!rns(1) = 1                
!rns(2) = 50
!rns(3) = 50
!rns(4) = 50   
!rns(5) = 50      
!rns(6) = 50      
!rns(7) = 0  

t_max = 100000
jmax = 1000
dimN = 1000
dimJ = 1000
record_after = 0.2

allocate(TimeSpent(0:1, 0:dimJ, 0:dimN))
allocate(dist_norm(0:1, 0:dimJ, 0:dimN))
allocate(nhistogram(0:dimJ, 0:dimN))
allocate(prob_dist(0:dimN))
allocate(jprob_dist(0:dimJ))
allocate(jhistogram(0:dimJ, 0:dimN))
allocate(ahistogram(0:1, 0:dimJ))
allocate(Pjhist1(0:jmax))
allocate(Pjhist2(0:jmax))
allocate(fjhist1(0:jmax))
allocate(fjhist2(0:jmax))

open(5, file = "results_o_v2_s_y100.dat", status = "replace")
write(5, *) "y, mean, var, Fano (simulation, analytical)"

do h = 1, size(rns)

rn = rns(h)
rf = rfs(h)

if (rn .eq. 0) then
    omega = 1000000
else
    omega = rf/rn
end if
gamma = kf/kn

nhistogram = 0
dist_norm = 0
TimeSpent = 0
ahistogram = 0
prob_dist = 0
jhistogram = 0
P1_simultion = 0
P2_simultion = 0

n = 0 !(0 to dimN)
s = 0 !(0 to dimJ; s=0 is ON, s=1 is j=0)
p = 0 !(0 or 1)

time = 0.0

!write(filename,'(a,i0,a)') 'timeseries', 0,'.dat'
!open(2, file = filename , status = "replace")

do while (time < t_max) 
    
    call random_number (r1)
    call random_number (r2)

    if (p == 0) then
        if ( s == 0) then
            if ( n == 0) then
                !!! 1 RNAP plane, ON state, #RNAs = 0
                phi = a + kf + rn 
                tau = (1/phi)*log(1/r1)
                if (time > t_max*record_after) then
                    TimeSpent(p,s,n) = TimeSpent(p,s,n) + tau
                end if
                time = time + tau
                if ( r2 < a/phi) then
                    n = n + 1
                else if ( r2 > a/phi .and. r2 < (a + kf)/phi ) then
                    s = 1
                else
                    p = 1
                end if
            else  
                !!! 1 RNAP plate, ON state, #RNAs > 0
                phi = a + kf + n*b + rn
                tau = (1/phi)*log(1/r1)
                if (time > t_max*record_after) then
                    TimeSpent(p,s,n) = TimeSpent(p,s,n) + tau
                end if
                time = time + tau
                if ( r2 < a/phi ) then 
                    n = n + 1
                else if (r2 > a/phi .and. r2 < (a+kf)/phi) then
                    s = 1
                else if ( r2 > (a+kf)/phi .and. r2 < (a+kf+n*b)/phi ) then
                    n = n - 1
                else 
                    p = 1
                end if 
            end if 
        else
            if ( n ==0 ) then
                !!! 1 RNAP plane, j=s-1 state, #RNAs = 0
                phi = a/(y**s) + kn + rn
                tau = (1/phi)*log(1/r1)
                if (time > t_max*record_after) then
                    TimeSpent(p,s,n) = TimeSpent(p,s,n) + tau
                end if
                time = time + tau
                if ( r2 < a/(y**s)/phi ) then
                    n = n + 1
                    s = s + 1
                else if ( r2 > a/(y**s)/phi .and. r2 < (a/(y**s) + kn)/phi) then
                    s = 0
                else 
                    p = 1
                end if 
            else 
                !!! 1 RNAP plane, j=s-1 state, #RNAs > 0
                phi = a/(y**s) + kn + n*b + rn
                tau = (1/phi)*log(1/r1)
                if (time > t_max*record_after) then
                    TimeSpent(p,s,n) = TimeSpent(p,s,n) + tau
                end if
                time = time + tau
                if ( r2 < a/(y**s)/phi ) then
                    n = n + 1
                    s = s + 1
                else if ( r2 > a/(y**s)/phi .and. r2 < (a/(y**s)+kn)/phi ) then
                    s = 0
                else if ( r2 > (a/(y**s)+kn)/phi .and. r2 < (a/(y**s)+kn+n*b)/phi )then
                    n = n - 1
                else 
                    p = 1
                end if 
            end if 
        end if 
    else
        if ( s == 0 ) then
            if ( n == 0 ) then
                !!! 2 RNAP plane, j=0 state, # RNAs = 0
                phi = 2*a + kf + rf
                tau = (1/phi)*log(1/r1)
                if (time > t_max*record_after) then
                    TimeSpent(p,s,n) = TimeSpent(p,s,n) + tau
                end if
                time = time + tau
                if ( r2 < 2*a/phi ) then
                    n = n + 1
                else if ( r2 > 2*a/phi .and. r2 < (2*a + kf)/phi ) then
                    s = 1
                else 
                    p = 0
                end if 
            else
                !!! 2 RNAP plane, j=0 state, # RNAs > 0
                phi = 2*a + kf + n*b + rf
                tau = (1/phi)*log(1/r1)
                if (time > t_max*record_after) then
                    TimeSpent(p,s,n) = TimeSpent(p,s,n) + tau
                end if
                time = time + tau
                if ( r2 < 2*a/phi ) then
                    n = n + 1
                else if ( r2 > 2*a/phi .and. r2 < (2*a+kf)/phi ) then
                    s = 1
                else if ( r2 > (2*a+kf)/phi .and. r2 < (2*a+kf+n*b)/phi ) then
                    n = n - 1
                else
                    p = 0
                end if 
            end if 
        else
            if ( n == 0 ) then 
                !!! 2 RNAP plane, j=s-1 state, # RNAs = 0
                phi = a + a/(y**s) + kn + rf
                tau = (1/phi)*log(1/r1)
                if (time > t_max*record_after) then
                    TimeSpent(p,s,n) = TimeSpent(p,s,n) + tau
                end if
                time = time + tau
                if ( r2 < (a+a/(y**s))/phi ) then
                    n = n + 1
                    s = s + 1
                else if ( r2 > (a+a/(y**s))/phi .and. r2 < (a+a/(y**s)+kn)/phi ) then
                    s = 0
                else
                    p = 0
                end if 
            else 
                !!! 2 RNAP plane, j=s-1 state, # RNAs = 0
                phi = a + a/(y**s) + kn + n*b + rf 
                tau = (1/phi)*log(1/r1)
                if (time > t_max*record_after) then
                    TimeSpent(p,s,n) = TimeSpent(p,s,n) + tau
                end if
                time = time + tau
                if ( r2 < (a+a/(y**s))/phi ) then
                    n = n + 1
                    s = s + 1
                else if ( r2 > (a+a/(y**s))/phi .and. r2 < (a+a/(y**s)+kn)/phi ) then
                    s = 0
                else if ( r2 > (a+a/(y**s)+kn)/phi .and. r2 < (a+a/(y**s)+kn+n*b)/phi ) then
                    n = n - 1
                else
                    p = 0
                end if 
            end if 
        end if 
    end if

    if (n .eq. dimN) then
        print*, "n is out of bounds"
    end if
    if (j .eq. dimJ) then
        print*, "n is out of bounds"
    end if
    !write (2, *) time, ",", p,",", s,",", n 
    !write (2, *) time, p, s, n 
end do
!close(2)

sum= 0
mean = 0
mean2 = 0
var2 = 0
fon = 0
dnt = 0
jmean = 0
amean = 0

! total time
!open(3, file = "time_spent.dat", status = "replace")
do k = 0, 1
do j = 0, dimJ
 do i = 0, dimN
    sum = sum + TimeSpent(k,j,i)
    !write(3, *) k, ",",j,",", i, ",", TimeSpent(k,j,i)
 end do
end do
end do
!close(3)

! normalize distribution, find mean
!open(1, file = "dist_norm.dat", status = "replace")
do k = 0,1
do j = 0,dimJ
 do i = 0, dimN
    dist_norm(k,j,i) = real(TimeSpent(k,j,i))/sum
    mean = mean + i*dist_norm(k,j,i)
    !write(1,*) k, j, i, dist_norm(k, j, i)
 end do
end do
end do 
!close(1)

! find variance, check normalized time sum
do k = 0,1
 do j = 0, dimJ
    do i = 0, dimN
        var2 = var2 + dist_norm(k,j,i)*(real(i)-mean)**2
        dnt = dnt + dist_norm(k,j, i)
    end do 
 end do
end do
!print*, "this should be = 1", dnt

do i = 0, dimN
    jhistogram(0,:) = jhistogram(0,:) + dist_norm(0,:,i)
    jhistogram(1,:) = jhistogram(1,:) + dist_norm(1,:,i)
end do

do i = 1, dimN ! start from 1 because 0 is ON state, which is excluded in P1_simultion, P2_simultion
    P1_simultion = P1_simultion + jhistogram(0, i)
    P2_simultion = P2_simultion + jhistogram(1, i)
end do 

jprob_dist(:) = jhistogram(0,:)+jhistogram(1,:)
do i = 1, dimJ
    jmean = jmean + jprob_dist(i)*i
end do 

do j = 0, dimJ
    nhistogram(j,:) = dist_norm(0,j,:) + dist_norm(1,j,:)
    prob_dist(:) = prob_dist(:) + nhistogram(j,:) 
end do 

!find mean, find fraction on
do i = 0, dimJ
    mean2 = mean2 + prob_dist(i)*i
    !fon = fon + dist_norm(1, i)
end do 

!open(4, file = "prob_dist.dat", status = "replace")
!do i = 0,dimJ
!        write(4, *) i, ",", prob_dist(i)
!end do
!close(4)

!open(6, file= "jhistogram.dat", status = "replace")
!do i = 0, dimN
!    write(6,*), i, ",", jhistogram(0,i), jhistogram(1, i)
!end do
!close(6)

!close(5) !===============finished all gammas=============



write(5,*), omega,",", mean,",", var2,",", var2/mean
print*, "done", omega

end do

close(5) 


end program twostate_ssa