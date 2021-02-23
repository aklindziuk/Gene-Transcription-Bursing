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
real, dimension(1:7) :: kns
real, dimension(1:7) :: kfs

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
rn = 50
rf = 50
y = 1.0

kfs(1) = 0          ! 0
kfs(2) = .1        ! 0.01
kfs(3) = 1        ! 0.1
kfs(4) = 10       ! 1
kfs(5) = 100      ! 10
kfs(6) = 1000      ! 100
kfs(7) = 1          ! Inf

kns(1) = 1
kns(2) = 10
kns(3) = 10
kns(4) = 10
kns(5) = 10
kns(6) = 10
kns(7) = 0

t_max = 1000000
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

open(5, file = "results.dat", status = "replace")
write(5, *) "y, mean, var, Fano (simulation, analytical)"

do h = 2, 6

kn = kns(h)
kf = kfs(h)

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


!!!===============Analytical Calculations====================

Pjhist1 = 0
Pjhist2 = 0
fjhist1 = 0
fjhist2 = 0

sum1 = 0
sum2 = 0
prod1 = 0
prod2 = 0

omega = rf/rn
if (h .eq. 7) then
    gamma = 10000000.0
else 
    gamma = kf/kn
end if 

Pon1 = omega/(1+gamma)/(1+omega)
Pon2 = 1/(1+gamma)/(1+omega)
gamma0 = (kn+a/y)*(kn+a+a/y)+rn*(a+a/y+kn)+rf*(a/y+kn)
P01 = kf*(kn+rf+a+a/y)*Pon1/gamma0 + rf*kf*Pon2/gamma0
P02 = kf*(kn+rn+a/y)*Pon2/gamma0 + rn*kf*Pon1/gamma0

Pjhist1(0) = P01
Pjhist2(0) = P02
!Evaluating Pjs
do j = 1, jmax 
    gammaj = (kn+a/y**(j+1))*(kn+a+a/y**(j+1))+rf*(a/y**(j+1)+kn)+rn*(kn+a+a/y**(j+1))
    Pj1 = a/y**j*(a+a/y**(j+1)+kn+rf)*Pjhist1(j-1)/gammaj+rf*(a+a/y**j)*Pjhist2(j-1)/gammaj
    Pj2 = (a+a/y**j)*(a/y**(j+1)+kn+rn)*Pjhist2(j-1)/gammaj+rn*a/y**j*Pjhist1(j-1)/gammaj
    Pjhist1(j) = Pj1
    Pjhist2(j) = Pj2
end do
!Summing over Pjs for mean
sum1 = 0
sum2 = 0
do j = 0, jmax
    sum1 = sum1+ Pjhist1(j)/(y**(j+1))
    sum2 = sum2+ Pjhist2(j)*(1.0+1.0/(y**(j+1)))
end do 
Pterm_plane1 = sum1
Pterm_plane2 = sum2 

n_plane2 = a/b*(2.0*Pon2 + Pterm_plane2)
n_plane1 = a/b*(Pon1 + Pterm_plane1)
n1 = a/b*(Pon1+2.0*Pon2+Pterm_plane1+Pterm_plane2)

!====================================================================================================

sum1 = 0
sum2 = 0
prod1 = 0
prod2 = 0

!fon12 = (kn*n1+a*Pon1+2.0*a*Pon2)/(b+kf+kn)

!fon2 = 1/(1+rf/rn)*fon12
!fon1 = rf/rn/(1+rf/rn)*fon12

!fon1_num = -((-b*kn-kf*kn-kn**2)*(-kn*(kn-rf)+a*Pon1*(rf+rn))-(kn*a*Pon1-kn*(-kn-2*a*Pon2))*(-(-kf-kn)*kn+rf*(rf+rn)))
!fon1_den =-(-kn*(kn-rn)-kn*(b+kf+rn))*(-(-kf-kn)*kn+rf*(rf+rn))+(-b*kn-kf*kn-kn**2)**(-(b+kf+rn)*(rf+rn)-kn*(-kf-kn+rf+rn))

!fon2_num1 = -b*kn**2+kf*kn*a*Pon1+kn**2*a*Pon1+kf*kn*2*a*Pon2+kn**2*2*a*Pon2-b*2*a*Pon2*rf-kf*2*a*Pon2*rf-kn*2*a*Pon2*rf-b*kn*rn
!fon2_num2 = -kf*kn*rn-kn**2*rn-b*2*a*Pon2*rn-kf*2*a*Pon2*rn-kn*2*a*Pon2*rn-kn*rf*rn-a*Pon1*rf*rn-2*a*Pon2*rf*rn-kn*rn**2
!fon2_num3 = -a*Pon1*rn**2-2*a*Pon2*rn**2
!fon2_den = (b+kf+kn)*(rf+rn)*(b+kf+kn+rf+rn)

!fon2 = (fon2_num1+fon2_num2+fon2_num3)/fon2_den
!fon1 = fon1_num/fon1_den
num1 = 0
num11 = 0
den1 = 0
num2 = 0
num22 = 0
den2 = 0
pterm1 = a*Pon1
pterm2 = 2*a*Pon2
num1 = (kn*pterm1 - kn*(-kn*n1 - pterm2))* rf*(rf + rn)
num11 = -((kn*rf - kn*(b + kf + kn + rf))*(kn*n1*rf + pterm1*(rf + rn)))
den1 = rf*(rf + rn)*(-kn*(kn - rn) - kn*(b + kf + rn))
den11 =  -(kn*rf - kn*(b + kf + kn + rf))*(-kn*(rf+rn)-(b+kf+rn)*(rf+rn))
num2 = (-b*kn-kf*kn-kn**2)*(kn*n1*rf+pterm1*(rf + rn))
num22 = -(kn*pterm1 - kn*(-kn*n1 - pterm2))*(-kn*(rf + rn)-(b+kf+rn)*(rf+rn))
den2 = (-b*kn-kf*kn-kn**2)*rf*(rf+rn)
den22 = -(kn*rf-kn*(b+kf+kn+rf))*(-kn*(rf+rn)-(b+kf+rn)*(rf+rn))
fon1 = -(num1+num11)/(den1+den11)
fon2 = -(num2+num22)/(den2+den22)

!print*, fon1, fon2

phi0 = (b+kn+rn+a/y)*(b+kn+rf+a+a/y)-rf*rn
!f01 = (b+kn+rf+a+a/y)*kf*fon1/phi0+rf*kf*fon2/phi0
!f02 = (b+kn+rn+a/y)*kf*fon2/phi0+rn*kf*fon1/phi0

f01 = -(kf*rf*fon2+kf*fon1*(b+kn+rf+a*(1+1/y)))/(rf*rn-(b+kn+rf+a*(1+1/y))*(b+kn+rn+a/y))
f02 = -((kf*rn*fon1+kf*fon2*(b+kn+rn+a/y))/(rf*rn-(b+kn+rf+a*(1+1/y))*(b+kn+rn+a/y)))

fjhist1(0) = f01
fjhist2(0) = f02

num1 = 0
den1 = 0
num2 = 0
den2 = 0

!Evaluating fjs
do j = 1, jmax
    Aterm = (b+kn+rn+a/y**(j+1))
    Bterm = (b+kn+rf+a+a/y**(j+1))
    phij = Aterm*Bterm-rn*rf

    !sum1 = Bterm*(a/y**j)*Pjhist1(j-1)/phij+rf*(a+a/y**j)*Pjhist2(j-1)/phij
    !sum2 = Aterm*(a+a/y**j)*Pjhist2(j-1)/phij+rn*(a/y**j)*Pjhist1(j-1)/phij

    !prod1 = Bterm*(a/y**j)*fjhist1(j-1)/phij+rf*(a+a/y**j)*fjhist2(j-1)/phij
    !prod2 = Aterm*(a+a/y**j)*fjhist2(j-1)/phij+rn*(a/y**j)*fjhist1(j-1)/phij

    !fjhist1(j) = sum1 + prod1
    !fjhist2(j) = sum2 + prod2 
    !print*, fjhist1(j)
    !print*, fjhist2(j)

    num1= -(-a*fjhist1(j-1)*y**(-j)-a*Pjhist1(j-1)*y**(-j))*(b+kn+rf+a*(1+y**(-1-j)))
    num11= -rf*(-a*fjhist2(j-1)*(1+y**(-j))-a*Pjhist2(j-1)*(1+y**(-j)))
    den1 = (rf*rn-(b+kn+rn+a*y**(-1-j))*(b+kn+rf+a*(1+y**(-1-j))))

    num2=-rn*(-a*fjhist1(j-1)*y**(-j)-a*Pjhist1(j-1)*y**(-j))
    num22 =-(b+kn+rn+a*y**(-1-j))*(-a*fjhist2(j-1)*(1+y**(-j))-a*Pjhist2(j-1)*(1+y**(-j)))
    den2 =(rf*rn-(b+kn+rn+a*y**(-1-j))*(b+kn+rf+a*(1+y**(-1-j))))

    fjhist1(j) = -(num1+num11)/den1
    fjhist2(j) = -(num2+num22)/den2


end do

sum1 = 0 
sum2 = 0
prod1 = 0
prod2 = 0
do j = 0, jmax
    sum1 = sum1 + fjhist1(j)/y**(j+1)
    sum2 = sum2 + fjhist2(j)*(1.0+1.0/y**(j+1))
    prod1 = prod1 + fjhist1(j)
    prod2 = prod2 + fjhist2(j)
end do
fterm_plane1 = sum1
fterm_plane2 = sum2 
!print*, "sum fon1+fon2", fon1+fon2, (kn*n1+a*Pon1+2*a*Pon2)/(b+kf+kn)
!print*, "sum fj 1+2", prod1+ prod2, "vs", n-(kn*n1+a*Pon1+2*a*Pon2)/(b+kf+kn)
!print*, "this should be n", prod1+prod2+fon1+fon2
!print*, "fterms", sum1, sum2

n2_plane2 = n_plane2 + a/b*(2.0*fon2+ fterm_plane2)
n2_plane1 = n_plane1 + a/b*(fon1 +fterm_plane1)
n2 = n1 + a/b*(fon1 + 2.0*fon2 + fterm_plane1 +fterm_plane2)



write(5,*), gamma, ",", mean,",", n1,",", var2,",", n2-n1**2, ",",var2/mean, ",",(n2-n1**2)/n1

print*, "done", gamma

end do

close(5) 


end program twostate_ssa