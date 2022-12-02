    subroutine get_basis
    
    include 'com.txt'
    
    do i = 1,NumGLP
        do j = 1,NumGLP
            phiG(i,j,1) = 1
            phixG(i,j,1) = 0
            phiyG(i,j,1) = 0
            mm(1) = 1
            if (k > 0) then
                phiG(i,j,2) = lambda(i)
                phixG(i,j,2) = 1d0/hx1
                phiyG(i,j,2) = 0
                mm(2) = 1d0/3d0
                
                phiG(i,j,3) = lambda(j)
                phixG(i,j,3) = 0
                phiyG(i,j,3) = 1d0/hy1
                mm(3) = 1d0/3d0
                if (k > 1) then
                    phiG(i,j,4) = lambda(i)**2 - 1d0/3d0
                    phixG(i,j,4) = 2d0*lambda(i)/hx1
                    phiyG(i,j,4) = 0
                    mm(4) = 4d0/45d0
                    
                    phiG(i,j,5) = lambda(i)*lambda(j)
                    phixG(i,j,5) = lambda(j)/hx1
                    phiyG(i,j,5) = lambda(i)/hy1
                    mm(5) = 1d0/9d0
                    
                    phiG(i,j,6) = lambda(j)**2 - 1d0/3d0
                    phixG(i,j,6) = 0
                    phiyG(i,j,6) = 2d0*lambda(j)/hy1
                    mm(6) = 4d0/45d0
                end if
            end if
        end do
    end do
    
    end subroutine get_basis