!Program: random configuration generator, 2D square lattice
!Author: Julia Carrillo Alonso and Sergi Martinez Galindo

program configuration
    implicit none
    integer :: L,N,i,seed
    double precision :: x
    character*32 :: sL
    character*128 :: file_to_write,r_seed
    REAL :: r1279

    !compilation: gfortran configuration.f90 r1279.f90 ran2.f -o file.exe
    !execution: .\file.exe L "file_name" random_seed
    !L and random_seed are integers
    !random seed used: 123

    !input variables
    call getarg(1 , sL)
    call getarg(2, file_to_write)
    call getarg(3, r_seed)

    READ (sL,*) L
    READ (r_seed,*) seed

    !random seed
    call setr1279(seed)
    !dimensions
    N=L**2
    !writing the random configuration
    open(1,file=file_to_write)
    do i=1,N
        write(1,*)i,2*mod(int(2*r1279()),2)-1
    enddo
    close(1)
    write(*,*)"End of the program reached."
end program


