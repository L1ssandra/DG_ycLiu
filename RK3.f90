    subroutine RK3
    
    include 'com.txt'
    
    t = 0
    count = 0
    
    call set_bc
    
    call TVB_Limiter
    
    do while (t < tend)
        
        call calculate_dt
    
        if (t + dt >= tend) then
            dt = tend - t
            t = tend
        else
            t = t + dt
        end if
    
        call set_bc
    
        call Lh
        
        uI = uh + dt*du
        
        uh0 = uh
        
        uh = uI
        
        call set_bc
        
        call TVB_Limiter
        
        if (RKorder == 3) then
        
            call set_bc
        
            call Lh
        
            uII = (3d0/4d0)*uh0 + (1d0/4d0)*uI + (1d0/4d0)*dt*du
        
            uh = uII
        
            call set_bc
        
            call TVB_Limiter
        
            call set_bc
        
            call Lh
        
            uh = (1d0/3d0)*uh0 + (2d0/3d0)*uII + (2d0/3d0)*dt*du
            
            call set_bc
        
            call TVB_Limiter
            
        end if
        
        count = count + 1
        
        call calculate_umax
        
        print *,count,"  ",t,"  ",umax
        
    end do
    
    end subroutine RK3