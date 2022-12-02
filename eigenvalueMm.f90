    subroutine eigenvalueMm(Amax,Amin,rho,rhou,rhov,E,B1,B2,n1,n2)
    
    include 'com.txt'
    real u,v,p,c,BP,Bn,un,cf
    
    u = rhou/rho
    v = rhov/rho
    
    BP = B1**2 + B2**2
    Bn = B1*n1 + B2*n2
    un = u*n1 + v*n2
    
    p = gamma1*(E - 0.5d0*rho*(u**2 + v**2) - 0.5d0*BP)
    
    c = sqrt(abs(gamma*p/rho))
    
    cf = sqrt(abs( 0.5d0*(c**2 + BP/rho + sqrt((c**2 + BP/rho)**2 - 4*c**2*Bn**2/rho) ) ))
    
    Amax = un + cf
    Amin = un - cf
    
    end subroutine eigenvalueMm
    
    