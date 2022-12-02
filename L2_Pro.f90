    subroutine L2_Pro
    
    include 'com.txt'
    
    do i = 1,Nx
        do j = 1,Ny
            do n = 1,numEq
                uh(i,j,:,n) = 0
                do d = 1,dimPk
                    do i1 = 1,NumGLP
                        do j1 = 1,NumGLP
                            uh(i,j,d,n) = uh(i,j,d,n) + weight(i1)*weight(j1)*phiG(i1,j1,d)*ureal(i,j,i1,j1,n)
                        end do
                    end do
                end do
            end do
        end do
    end do
    
    uh = 0.25*uh
    
    do d = 1,dimPk
        uh(:,:,d,:) = uh(:,:,d,:)/mm(d)
    end do
    
    end subroutine L2_Pro