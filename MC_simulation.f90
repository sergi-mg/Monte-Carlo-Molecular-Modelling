!Program: Monte Carlo simulation
!Author: Julia Carrillo Alonso and Sergi Martinez Galindo

program P2_Monte_Carlo
    implicit none
    integer :: L,Nmcs,Nmeas,N,z
    character*32 :: sL,sT,sNMCS,sNMEAS
    character*128 :: file_to_read
    double precision :: E,M,T

    integer, dimension(:), allocatable :: S
    integer, dimension(:,:), allocatable :: nbr

    !compilation: gfortran simulation.f90 r1279.f90 ran2.f -o file.exe
    !execution ".\file.exe L T NMCS NMEAS "file_name""

    !input
    call getarg(1 , sL) !lattice size: LxL (int)
    call getarg(2 , sT) !Temperature (double precision)
    call getarg(3 , sNMCS) !MC steps (int)
    call getarg(4 , sNMEAS) !values saved every NMEAS steps, also random seed (int)
    call getarg(5, file_to_read) !file containing the initial configuration


    READ (sL,*) L
    READ (sT,*) T
    READ (sNMCS,*) Nmcs
    READ (sNMEAS,*) Nmeas

    write(*,*) "Values read."

    !constructing nbr
    N=L**2
    z=4 !square 2D lattice
    allocate(S(1:N),nbr(1:N,1:z))
    !initial configuration
    call read_spins(file_to_read,L,S)
    !nbr matrix
    call square_nbr(N,L,nbr)
    write(*,*) "S and nbr created"

    !Montecarlo simulation
    !call random seed Nmeas
    call setr1279(Nmeas)
    write(*,*) "Random seed initialized."
    !Monte Carlo simulation
    call metropolis(S,L,nbr,z,T,Nmcs,Nmeas)

    write(*,*)"Final of the program reached."
end program P2_Monte_Carlo

double precision function ENERG(S,N,nbr,z)
    implicit none
    integer :: N,z,i,k
    integer :: S(1:N),nbr(1:N,1:z) 
    !system's energy
    ENERG=0.d0 
    do i=1,N
        do k=1,z
            ENERG=ENERG-0.5d0*dble(S(i)*S(nbr(i,k)))
        enddo
    enddo
    return
end function

double precision function MAGNE(S,N)
    implicit none
    integer :: N,i
    integer :: S(1:N) 
    !system's magnetization
    MAGNE=0.d0 
    do i=1,N
        MAGNE=MAGNE+dble(S(i))
    enddo
    return
end function

subroutine measures(S,N,nbr,z,E,M)
    implicit none
    integer :: N,z,i,k
    integer :: S(1:N),nbr(1:N,1:z)
    double precision :: ENERG,MAGNE,E,M
    E=ENERG(S,N,nbr,z)
    M=MAGNE(S,N,z)
    return
end subroutine

subroutine PBC_matrix(L,PBC)
    implicit none
    integer :: L,i,PBC(0:1,1:L)
    do i = 1, L
        PBC(0,i)=i-1
        PBC(1,i)=i+1
    end do
    PBC(0,1)=L
    PBC(1,L)=1
    return
end subroutine

subroutine square_nbr(N,L,nbr)
    implicit none
    integer :: N,L,nbr(1:N,1:4),PBC(0:1,1:L),i,x,y
    call PBC_matrix(L,PBC) 
    i=0
    do y=1,L
        do x=1,L
            i=i+1
            nbr(i,1)=PBC(1,x)+L*(y-1)
            nbr(i,2)=PBC(0,x)+L*(y-1)
            nbr(i,3)=x+L*(PBC(1,y)-1)
            nbr(i,4)=x+L*(PBC(0,y)-1)
        end do
    end do
    return
end subroutine

subroutine read_spins(file_name,L,S)
    implicit none
    integer :: L, S(1:L**2), i,node,spin,iostat
    character*128 :: file_name
    open(1,file=file_name)
    do i=1,2*L**2
        read(1,*,iostat=iostat)node,spin
        if (iostat.lt.0) then
            write(*,*)"Exit"
            exit
        else
            if (node.le.L**2) then
                S(node)=spin
            else
                write(*,*)"There are more spins that expected."
            endif
        endif
    enddo
    close(1)
    return
end subroutine

subroutine metropolis(S,L,nbr,z,T,MCTOT,dMC)
    implicit none
    !inputs
    integer :: z,L,S(1:L**2),nbr(1:L**2,1:z),MCTOT,dMC
    double precision :: T
    !variables used 
    integer :: sr,i,IMC,IPAS,dE, N,dM
    double precision :: ENE, MAG, exponent,q,w(0:4*z),t0,tf
    REAL :: r1279
    character*128 :: file_name
    !creating the file to save the values needed
    WRITE(file_name, '(A11,F4.1, A1, I0,A4)') "results_MC_", T,"_" ,L ,".txt"  
    open(3,file=file_name)
    write(*,*)"Montecarlo: Filed created."
    !total number of spins
    N=L*L
    !Initial state
    call measures(S,N,nbr,z,ENE,MAG)
    write(3,"(2(F14.8,3X))") ENE/dble(N),MAG/dble(N)
    write(*,*)"Initial values measured"

    !Vector containing the possible Boltzmann factors for a single spin flip
    do i=0,4*z 
        exponent=-dble(i-2*z)/T
        w(i)=exp(exponent)
    enddo
    write(*,*)"Boltzmann Factor vector created"

    call cpu_time(t0)
    write(*,*)"Cpu time"

    do IMC=1,MCTOT
        !MC step
        do IPAS=1,N
            !spin selection 
            sr=mod(int(N*r1279()),N)+1
            !energy variation
            dE=0
            do i=1,z
                dE=dE+S(nbr(sr,i))
            enddo
            dE=dE*2*S(sr)
            !magnetization variation
            dM=-2*S(sr)
            !do we accept the flip?
            if (dE.le.0) then
                S(sr)=-S(sr)
                ENE=ENE+dble(dE)
                MAG=MAG+dble(dM)
            else 
                q=dble(r1279())
                if (q.lt.w(dE+2*z)) then
                    S(sr)=-S(sr)
                    ENE=ENE+dble(dE)
                    MAG=MAG+dble(dM)
                endif
            endif
        enddo
        !writing the results
        if (mod(IMC,dMC).eq.0) then
            write(3,"(2(F14.8,3X))") ENE/dble(N),MAG/dble(N)
        endif
    enddo
    call cpu_time(tf)
    write(*,*)dble(MCTOT*IPAS)/(tf-t0), "spins/second"
    close(3)
    return
end subroutine



