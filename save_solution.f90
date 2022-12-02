    subroutine save_solution
    
    include 'com.txt'
    
    open(unit = 1,file = 'Q1.txt')
    open(unit = 2,file = 'Q2.txt')
    open(unit = 3,file = 'Q3.txt')
    open(unit = 4,file = 'Q4.txt')
    open(unit = 5,file = 'Q5.txt')
    open(unit = 6,file = 'Q6.txt')
    
    uG = 0
    
    do i = 1,Nx
        do j = 1,Ny
            do n = 1,NumEq
                do d = 1,dimPk
                    uG(i,j,:,:,n) = uG(i,j,:,:,n) + uh(i,j,d,n)*phiG(:,:,d)
                end do
            end do
        end do
    end do
    
    do i = 1,Nx
        do i1 = 1,NumGLP
            do j = 1,Ny
                do j1 = 1,NumGLP
                    write(1,*) uG(i,j,i1,j1,1)
                    write(2,*) uG(i,j,i1,j1,2)
                    write(3,*) uG(i,j,i1,j1,3)
                    write(4,*) uG(i,j,i1,j1,4)
                    write(5,*) uG(i,j,i1,j1,5)
                    write(6,*) uG(i,j,i1,j1,6)
                end do
            end do
        end do
    end do
    
    close(1)
    close(2)
    close(3)
    close(4)
    close(5)
    close(6)
    
    end subroutine save_solution