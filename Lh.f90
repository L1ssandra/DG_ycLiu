    subroutine Lh
    
    include 'com.txt'
    
    du = 0
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
    
    rhoM = uG(:,:,:,:,1)
    uM = uG(:,:,:,:,2)/rhoM
    vM = uG(:,:,:,:,3)/rhoM
    EM = uG(:,:,:,:,4)
    B1M = uG(:,:,:,:,5)
    B2M = uG(:,:,:,:,6)
    
    pM = gamma1*(EM - 0.5d0*rhoM*(uM**2 + vM**2) - 0.5d0*(B1M**2 + B2M**2))
    
    SM = pM + 0.5d0*(B1M**2 + B2M**2)
    TM = EM + SM
    KM = uM*B1M + vM*B2M
    
    Fx(:,:,:,:,1) = uG(:,:,:,:,2)
    Fx(:,:,:,:,2) = rhoM*uM**2 + SM - B1M**2
    Fx(:,:,:,:,3) = rhoM*uM*vM - B1M*B2M
    Fx(:,:,:,:,4) = TM*uM - KM*B1M
    Fx(:,:,:,:,5) = 0
    Fx(:,:,:,:,6) = uM*B2M - vM*B1M
    
    Fy(:,:,:,:,1) = uG(:,:,:,:,3)
    Fy(:,:,:,:,2) = rhoM*uM*vM - B1M*B2M
    Fy(:,:,:,:,3) = rhoM*vM**2 + SM - B2M**2
    Fy(:,:,:,:,4) = TM*vM - KM*B2M
    Fy(:,:,:,:,5) = vM*B1M - uM*B2M
    Fy(:,:,:,:,6) = 0
    
    do i = 1,Nx
        do j = 1,Ny
            do d = 2,dimPk
                do n = 1,NumEq
                    do i1 = 1,NumGLP
                        do j1 = 1,NumGLP
                            du(i,j,d,n) = du(i,j,d,n) + 0.25d0*weight(i1)*weight(j1)*(Fx(i,j,i1,j1,n)*phixG(i1,j1,d) + Fy(i,j,i1,j1,n)*phiyG(i1,j1,d))
                        end do
                    end do
                end do
            end do
        end do
    end do
    
    ! The x-Flux
    UR = 0
    UL = 0
    do i = 0,Nx
        do j = 1,Ny
            do d = 1,dimPk
                do n = 1,NumEq
                    UR(i,j,:,n) = UR(i,j,:,n) + uh(i,j,d,n)*phiG(NumGLP,:,d)
                    UL(i + 1,j,:,n) = UL(i + 1,j,:,n) + uh(i + 1,j,d,n)*phiG(1,:,d)
                end do
            end do
        end do
    end do
    
    rhoMR = uR(:,:,:,1)
    uMR = uR(:,:,:,2)/rhoMR
    vMR = uR(:,:,:,3)/rhoMR
    EMR = uR(:,:,:,4)
    B1MR = uR(:,:,:,5)
    B2MR = uR(:,:,:,6)
    
    pMR = gamma1*(EMR - 0.5d0*rhoMR*(uMR**2 + vMR**2) - 0.5d0*(B1MR**2 + B2MR**2))
    
    SMR = pMR + 0.5d0*(B1MR**2 + B2MR**2)
    TMR = EMR + SMR
    KMR = uMR*B1MR + vMR*B2MR
    
    FR(:,:,:,1) = UR(:,:,:,2)
    FR(:,:,:,2) = rhoMR*uMR**2 + SMR - B1MR**2
    FR(:,:,:,3) = rhoMR*uMR*vMR - B1MR*B2MR
    FR(:,:,:,4) = TMR*uMR - KMR*B1MR
    FR(:,:,:,5) = 0
    FR(:,:,:,6) = uMR*B2MR - vMR*B1MR
    
    rhoML = uL(:,:,:,1)
    uML = uL(:,:,:,2)/rhoML
    vML = uL(:,:,:,3)/rhoML
    EML = uL(:,:,:,4)
    B1ML = uL(:,:,:,5)
    B2ML = uL(:,:,:,6)
    
    pML = gamma1*(EML - 0.5d0*rhoML*(uML**2 + vML**2) - 0.5d0*(B1ML**2 + B2ML**2))
    
    SML = pML + 0.5d0*(B1ML**2 + B2ML**2)
    TML = EML + SML
    KML = uML*B1ML + vML*B2ML
    
    FL(:,:,:,1) = UL(:,:,:,2)
    FL(:,:,:,2) = rhoML*uML**2 + SML - B1ML**2
    FL(:,:,:,3) = rhoML*uML*vML - B1ML*B2ML
    FL(:,:,:,4) = TML*uML - KML*B1ML
    FL(:,:,:,5) = 0
    FL(:,:,:,6) = uML*B2ML - vML*B1ML
    
    ! calculate F hat
    do i = 0,Nx
        do j = 1,Ny
            do j1 = 1,NumGLP
                call eigenvalueMm(SRmax,SRmin,UR(i,j,j1,1),UR(i,j,j1,2),UR(i,j,j1,3),UR(i,j,j1,4),UR(i,j,j1,5),UR(i,j,j1,6),1,0)
                call eigenvalueMm(SLmax,SLmin,UL(i + 1,j,j1,1),UL(i + 1,j,j1,2),UL(i + 1,j,j1,3),UL(i + 1,j,j1,4),UL(i + 1,j,j1,5),UL(i + 1,j,j1,6),1,0)
                SR = max(SRmax,SLmax)
                SL = min(SRmin,SLmin)
                FR1 = FL(i + 1,j,j1,:)
                FL1 = FR(i,j,j1,:)
                UR1 = UL(i + 1,j,j1,:)
                UL1 = UR(i,j,j1,:)
                if (flux_type == 1) then
                    call LF_Flux
                else if (flux_type == 2) then
                    call HLL_Flux
                else if (flux_type == 3) then
                    direction = 1
                    call HLLC_Flux
                end if
                Fxhat(i,j,j1,:) = Fhat1
            end do
        end do
    end do
    
    do i = 1,Nx
        do j = 1,Ny
            do d = 1,dimPk
                do n = 1,NumEq
                    do j1 = 1,NumGLP
                        du(i,j,d,n) = du(i,j,d,n) - (0.5d0/hx)*weight(j1)*(Fxhat(i,j,j1,n)*phiG(NumGLP,j1,d) - Fxhat(i - 1,j,j1,n)*phiG(1,j1,d))
                    end do
                end do
            end do
        end do
    end do
                
    
    ! The y-Flux
    UU = 0
    UD = 0
    do i = 1,Nx
        do j = 0,Ny
            do d = 1,dimPk
                do n = 1,NumEq
                    UU(i,j,:,n) = UU(i,j,:,n) + uh(i,j,d,n)*phiG(:,NumGLP,d)
                    UD(i,j + 1,:,n) = UD(i,j + 1,:,n) + uh(i,j + 1,d,n)*phiG(:,1,d)
                end do
            end do
        end do
    end do
    
    rhoMU = UU(:,:,:,1)
    uMU = UU(:,:,:,2)/rhoMU
    vMU = UU(:,:,:,3)/rhoMU
    EMU = UU(:,:,:,4)
    B1MU = UU(:,:,:,5)
    B2MU = UU(:,:,:,6)
    
    pMU = gamma1*(EMU - 0.5d0*rhoMU*(uMU**2 + vMU**2) - 0.5d0*(B1MU**2 + B2MU**2))
    
    SMU = pMU + 0.5d0*(B1MU**2 + B2MU**2)
    TMU = EMU + SMU
    KMU = uMU*B1MU + vMU*B2MU
    
    FU(:,:,:,1) = UU(:,:,:,3)
    FU(:,:,:,2) = rhoMU*uMU*vMU - B1MU*B2MU
    FU(:,:,:,3) = rhoMU*vMU**2 + SMU - B2MU**2
    FU(:,:,:,4) = TMU*vMU - KMU*B2MU
    FU(:,:,:,5) = vMU*B1MU - uMU*B2MU
    FU(:,:,:,6) = 0
    
    rhoMD = UD(:,:,:,1)
    uMD = UD(:,:,:,2)/rhoMD
    vMD = UD(:,:,:,3)/rhoMD
    EMD = UD(:,:,:,4)
    B1MD = UD(:,:,:,5)
    B2MD = UD(:,:,:,6)
    
    pMD = gamma1*(EMD - 0.5d0*rhoMD*(uMD**2 + vMD**2) - 0.5d0*(B1MD**2 + B2MD**2))
    
    SMD = pMD + 0.5d0*(B1MD**2 + B2MD**2)
    TMD = EMD + SMD
    KMD = uMD*B1MD + vMD*B2MD
    
    FD(:,:,:,1) = UD(:,:,:,3)
    FD(:,:,:,2) = rhoMD*uMD*vMD - B1MD*B2MD
    FD(:,:,:,3) = rhoMD*vMD**2 + SMD - B2MD**2
    FD(:,:,:,4) = TMD*vMD - KMD*B2MD
    FD(:,:,:,5) = vMD*B1MD - uMD*B2MD
    FD(:,:,:,6) = 0
    
    ! calculate F hat
    do i = 1,Nx
        do j = 0,Ny
            do i1 = 1,NumGLP
                call eigenvalueMm(SRmax,SRmin,UU(i,j,i1,1),UU(i,j,i1,2),UU(i,j,i1,3),UU(i,j,i1,4),UU(i,j,i1,5),UU(i,j,i1,6),0,1)
                call eigenvalueMm(SLmax,SLmin,UD(i,j + 1,i1,1),UD(i,j + 1,i1,2),UD(i,j + 1,i1,3),UD(i,j + 1,i1,4),UD(i,j + 1,i1,5),UD(i,j + 1,i1,6),0,1)
                SR = max(SRmax,SLmax)
                SL = min(SRmin,SLmin)
                FR1 = FD(i,j + 1,i1,:)
                FL1 = FU(i,j,i1,:)
                UR1 = UD(i,j + 1,i1,:)
                UL1 = UU(i,j,i1,:)
                if (flux_type == 1) then
                    call LF_Flux
                else if (flux_type == 2) then
                    call HLL_Flux
                else if (flux_type == 3) then
                    direction = 2
                    call HLLC_Flux
                end if
                Fyhat(i,j,i1,:) = Fhat1
            end do
        end do
    end do
    
    do i = 1,Nx
        do j = 1,Ny
            do d = 1,dimPk
                do n = 1,NumEq
                    do i1 = 1,NumGLP
                        du(i,j,d,n) = du(i,j,d,n) - (0.5d0/hy)*weight(i1)*(Fyhat(i,j,i1,n)*phiG(i1,NumGLP,d) - Fyhat(i,j - 1,i1,n)*phiG(i1,1,d))
                    end do
                end do
            end do
        end do
    end do
    
    do d = 1,dimPk
        du(:,:,d,:) = du(:,:,d,:)/mm(d)
    end do
    
    end subroutine Lh