! File: Driver_LinAL.f90
! Author: Adrian Garcia
! Date: 03/18/2023
Program Driver_LinAl
  implicit none
  external :: dgesvd
  integer,parameter :: M = 1920,N = 1279,LDA = M,LDU = M,LDVT = N
  character(len = 100) :: myFileName
  character(len = 100),dimension(8) :: filename
  integer :: i,j,l,p,INFO,LWORK
  real(kind = kind(0.d0)) :: fnorm,error
  real(kind = kind(0.d0)),dimension(8) :: k
  real(kind = kind(0.d0)),dimension(N) :: S
  real(kind = kind(0.d0)),dimension(:),allocatable :: WORK
  real(kind = kind(0.d0)),dimension(N,M) :: mat
  real(kind = kind(0.d0)),dimension(M,N) :: SIGMA,AApprox,E
  real(kind = kind(0.d0)),dimension(LDA,N) :: A
  real(kind = kind(0.d0)),dimension(LDU,M) :: U
  real(kind = kind(0.d0)),dimension(LDVT,N) :: VT
  ! Read data file
  myFileName = 'dog_bw_data.dat'
  open(10,file = myFileName)
  do i = 1,N
    read(10,*) (mat(i,j),j = 1,M)
  enddo
  close(10)
  ! Initialize matrix A
  A = transpose(mat)
  ! Query the optimal workspace
  LWORK = -1
  allocate(WORK(1))
  call dgesvd('All','All',M,N,A,LDA,S,U,LDU,VT,LDVT,WORK,LWORK,INFO)
  LWORK = nint(WORK(1))
  deallocate(WORK)
  allocate(WORK(LWORK))
  ! Compute SVD
  call dgesvd('All','All',M,N,A,LDA,S,U,LDU,VT,LDVT,WORK,LWORK,INFO)
  ! Check for convergence
  if(INFO.GT.0) then
    write(*,*) 'The algorithm computing SVD failed to converge.'
    stop
  end if
  ! Initialize desired low-rank approximation values
  k = (/10,20,40,80,160,320,640,1279/)
  filename = (/'Image_appn_100010.dat','Image_appn_100020.dat','Image_appn_100040.dat','Image_appn_100080.dat',& 
             & 'Image_appn_100160.dat','Image_appn_100320.dat','Image_appn_100640.dat','Image_appn_101279.dat'/)
  do i = 1,size(k)
    ! Initialize approximated k-SIGMA
    SIGMA = 0.0
    do j = 1,int(k(i))
      SIGMA(j,j) = S(j)
    end do
    ! Compute the k-SVD approximation
    AApprox = matmul(U,matmul(SIGMA,VT))
    ! Save to appropriate data file
    open(11,file = filename(i))
    do l = 1,M
      write(11,*) (AApprox(l,p),p = 1,N)
    enddo
    close(11)
    ! Calculate error matrix
    E = transpose(mat)-AApprox
    ! Initiallize the Frobenious norm as 0
    fnorm = 0.0
    ! Calculate the Frobenious norm
    do l = 1,M
      do p = 1,N
        fnorm = fnorm+E(l,p)**2
      enddo
    enddo
    fnorm = sqrt(fnorm)
    ! Calculate average Frobenious norm
    error = fnorm/(M*N)
    ! Save to error data file
    open(12,file = 'error.dat',action = 'write',position = 'append')
    write(12,*) k(i), error
    close(12)
  end do
  write(*,*) '**********************************************************'
  do i = 1,size(k)
    if (int(k(i)) .eq. 10) then
      do j = 1,10
        write(*,*) j, S(j)
      end do
    else
      write(*,*) int(k(i)), S(int(k(i)))
    end if
  end do
  write(*,*) '**********************************************************'
  deallocate(WORK)
End Program Driver_LinAl