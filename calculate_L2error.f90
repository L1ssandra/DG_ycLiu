    subroutine calculate_L2error
    
    include 'com.txt'
    
    L2 = 0
    uG = 0
    
    ! The value of num solution on the GL points
    do i = 1,Nx
        do j = 1,Ny
            do n = 1,NumEq
                do d = 1,dimPk
                    uG(i,j,:,:,n) = uG(i,j,:,:,n) + uh(i,j,d,n)*phiG(:,:,d)
                end do
            end do
        end do
    end do
    
    uE = (uG - ureal)**2
    
    do n = 1,NumEq
        do i = 1,Nx
            do j = 1,Ny
                do i1 = 1,NumGLP
                    do j1 = 1,NumGLP
                        L2(n) = L2(n) + weight(i1)*weight(j1)*uE(i,j,i1,j1,n)
                    end do
                end do
            end do
        end do
    end do
    
    L2 = (L2*hx1*hy1)**0.5d0
    
    print *,L2
    
    end subroutine calculate_L2error