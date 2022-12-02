    program main
    
    include 'com.txt'
    
    call get_GLP
    
    call init_data
    
    call get_basis
    
    call L2_Pro
    
    call RK3
    
    call save_XY
    
    !call real_solution
    
    call save_solution
    
    call calculate_L2error

    end program main

