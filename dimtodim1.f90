    subroutine dimtodim1(U,U1)
    
    real U(6),U1(8)
    
    U1 = 0
    
    U1(1:3) = U(1:3)
    U1(5:7) = U(4:6)
    
    end subroutine dimtodim1
    
    
    