    subroutine HLLC_Flux
    
    include 'com.txt'
    
    real Sstar,rhoR,rhoL,uRbot,uLbot,PR,PL,ER,EL,u1R,u1L,u2R,u2L,B1R,B1L,B2R,B2L
    real rhoRstar,rhoLstar,B1star,B2star,u1star,u2star,u1Rstar,u1Lstar,u2Rstar,u2Lstar
    real BRbot,BLbot,Bstarbot,Pstar,ERstar,ELstar
    
    if (SR < 0) then
        Fhat1 = FR1
    else if (SL > 0) then
        Fhat1 = FL1
    else
        Ustar = ( SR*UR1 - SL*UL1 + FL1 - FR1 )/(SR - SL)
        
        
        rhoR = UR1(1)
        rhoL = UL1(1)
        
        ER = UR1(4)
        u1R = UR1(2)/rhoR
        u2R = UR1(3)/rhoR
        B1R = UR1(5)
        B2R = UR1(6)
        
        EL = UL1(4)
        u1L = UL1(2)/rhoL
        u2L = UL1(3)/rhoL
        B1L = UL1(5)
        B2L = UL1(6)
        
        if (direction == 1) then
            uRbot = u1R
            uLbot = u1L
        else if (direction == 2) then
            uRbot = u2R
            uLbot = u2L
        end if
        
        PR = gamma1*(ER - 0.5d0*rhoR*(u1R**2 + u2R**2) - 0.5d0*(B1R**2 + B2R**2))
        PL = gamma1*(EL - 0.5d0*rhoL*(u1L**2 + u2L**2) - 0.5d0*(B1L**2 + B2L**2))
        
        Sstar = ( (SR - uRbot)*rhoR*uRbot - (SL - uLbot)*rhoL*uLbot - (PR - PL) )/( (SR - uRbot)*rhoR - (SL - uLbot)*rhoL )
        
        rhoRstar = rhoR*(SR - uRbot)/(SR - Sstar)
        rhoLstar = rhoL*(SL - uLbot)/(SL - Sstar)
        
        B1star = Ustar(5)
        B2star = Ustar(6)
        u1star = Ustar(2)/Ustar(1)
        u2star = Ustar(3)/Ustar(1)
        
        if (direction == 1) then
            u1Rstar = Sstar
            u1Lstar = Sstar
            u2Rstar = u2R + (B1R*B2R - B1star*B2star)/(rhoR*(SR - uRbot))
            u2Lstar = u2L + (B1L*B2L - B1star*B2star)/(rhoL*(SL - uLbot))
            BRbot = B1R
            BLbot = B1L
            Bstarbot = B1star
        else if (direction == 2) then
            u1Rstar = u1R + (B1R*B2R - B1star*B2star)/(rhoR*(SR - uRbot))
            u1Lstar = u1L + (B1L*B2L - B1star*B2star)/(rhoL*(SL - uLbot))
            u2Rstar = Sstar
            u2Lstar = Sstar
            BRbot = B2R
            BLbot = B2L
            Bstarbot = B2star
        end if
        
        Pstar = PL + rhoL*(SL - uLbot)*(Sstar - uLbot)
        
        ERstar = ( (SR - uRbot)*ER - PR*uRbot + Pstar*Sstar + BRbot*(u1R*B1R + u2R*B2R) - Bstarbot*(u1star*B1star + u2star*B2star) )/(SR - Sstar)
        ELstar = ( (SL - uLbot)*EL - PL*uLbot + Pstar*Sstar + BLbot*(u1L*B1L + u2L*B2L) - Bstarbot*(u1star*B1star + u2star*B2star) )/(SL - Sstar)
        
        URstar(1) = rhoRstar
        URstar(2) = rhoRstar*u1Rstar
        URstar(3) = rhoRstar*u2Rstar
        URstar(4) = ERstar
        URstar(5) = B1star
        URstar(6) = B2star
        
        ULstar(1) = rhoLstar
        ULstar(2) = rhoLstar*u1Lstar
        ULstar(3) = rhoLstar*u2Lstar
        ULstar(4) = ELstar
        ULstar(5) = B1star
        ULstar(6) = B2star
        if (Sstar > 0) then
            Fhat1 = FL1 + SL*(ULstar - UL1)
        else
            Fhat1 = FR1 + SR*(URstar - UR1)
        end if
        
    end if
    
    end subroutine HLLC_Flux