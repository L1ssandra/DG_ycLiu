    subroutine get_GLP
    
    include 'com.txt'
    
    if (NumGLP == 3) then
        lambda(1) = -1
        lambda(2) = 0
        lambda(3) = 1
        
        weight(1) = 1d0/3d0
        weight(2) = 4d0/3d0
        weight(3) = 1d0/3d0
    else if (NumGLP == 4) then
        lambda(1) = -1
        lambda(2) = -0.447213595499957939282
        lambda(3) = 0.447213595499957939282
        lambda(4) = 1
        
        weight(1) = 1d0/6d0
        weight(2) = 5d0/6d0
        weight(3) = 5d0/6d0
        weight(4) = 1d0/6d0
    else if (NumGLP == 5) then
        lambda(1) = -1
        lambda(2) = -0.6546536707079771437983
        lambda(3) = 0
        lambda(4) = 0.6546536707079771437983
        lambda(5) = 1
        
        weight(1) = 0.1
        weight(2) = 49d0/90d0
        weight(3) = 32d0/45d0
        weight(4) = 49d0/90d0
        weight(5) = 0.1
    end if
    
    end subroutine get_GLP