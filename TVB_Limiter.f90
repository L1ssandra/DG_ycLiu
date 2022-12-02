    subroutine TVB_Limiter
    
    include 'com.txt'
    
    uhmod = uh
    
    do i = 1,Nx
        do j = 1,Ny
            
            change = 0
            
            ! x-direction
            DeltaUR(:,1) = uh(i + 1,j,1,:) - uh(i,j,1,:)
            DeltaUL(:,1) = uh(i,j,1,:) - uh(i - 1,j,1,:)
            DeltaU(:,1) = uh(i,j,2,:)
            
            call dimtodim1(DeltaUR,DeltaUR1)
            call dimtodim1(DeltaUL,DeltaUL1)
            call dimtodim1(DeltaU,DeltaU1)
            
            direction = 1
            
            call minmod
            
            call dim1todim(DeltaUmod1,DeltaUxmod)
            
            do d = 1,NumEq
                if (DeltaUxmod(d,1) /= DeltaU(d,1)) then
                    change(d) = 1
                end if
            end do
            
            ! y-direction
            DeltaUR(:,1) = uh(i,j + 1,1,:) - uh(i,j,1,:)
            DeltaUL(:,1) = uh(i,j,1,:) - uh(i,j - 1,1,:)
            DeltaU(:,1) = uh(i,j,3,:)
            
            call dimtodim1(DeltaUR,DeltaUR1)
            call dimtodim1(DeltaUL,DeltaUL1)
            call dimtodim1(DeltaU,DeltaU1)
            
            direction = 2
            
            call minmod
            
            call dim1todim(DeltaUmod1,DeltaUymod)
            
            do d = 1,NumEq
                if (DeltaUymod(d,1) /= DeltaU(d,1)) then
                    change(d) = 1
                end if
            end do
            
            do d = 1,NumEq
                if (change(d) == 1) then
                    uhmod(i,j,4:dimPk,d) = 0
                    uhmod(i,j,2,d) = DeltaUxmod(d,1)
                    uhmod(i,j,3,d) = DeltaUymod(d,1)
                end if
            end do
            
        end do
    end do
    
    uh = uhmod
    
    end subroutine TVB_Limiter
    
    