PROGRAM MAIN

    INTEGER , PARAMETER :: fp = selected_real_kind( 15 , 307 )

    REAL(fp) :: m_e , m_p , q_e , q_p
    REAL(fp) :: k , G , Fe , Fg , E
    
    REAL(fp) , DIMENSION(3) :: Electron_pos , Electron_vel , Electron_acc , Electron_acc2
    REAL(fp) , DIMENSION(3) :: Proton_pos , Proton_vel , Proton_acc , Proton_acc2

    REAL(fp) :: T , dT
    REAL(fp) :: dx , dy , dz , r

    INTEGER :: n , count

    m_e = 9.1093837015E-31_fp
    m_p = 1.67262192369E-27_fp

    q_e = -1.602176634E-19_fp
    q_p = 1.602176634E-19_fp

    k = 8.987551782E9_fp
    G = 6.6743E-11_fp

    Electron_pos(1) = 1_fp
    Electron_pos(2) = 0_fp
    Electron_pos(3) = 0_fp

    Electron_vel(1) = 0_fp
    Electron_vel(2) = 14_fp
    Electron_vel(3) = 0_fp

    Proton_pos(1) = 0_fp
    Proton_pos(2) = 0_fp
    Proton_pos(3) = 0_fp

    Proton_vel(1) = 0_fp
    Proton_vel(2) = 0_fp
    Proton_vel(3) = 0_fp

    T = 1E+0_fp
    dT = 1E-5_fp

    n = NINT( T/dT )

    OPEN( unit = 1 , file = "data.csv" )
    OPEN( unit = 2 , file = "data2.csv" )

    DO count = 1 , n

        WRITE( 1 , * ) Electron_pos , Proton_pos

        dx = Electron_pos(1) - Proton_pos(1)
        dy = Electron_pos(2) - Proton_pos(2)
        dz = Electron_pos(3) - Proton_pos(3)

        r = sqrt( dx**2 + dy**2 + dz**2 )

        Fe = k * q_e * q_p / r**2
        Fg = -G * m_e * m_p / r**2

        Electron_acc(1) = ( Fe + Fg ) * cos(atan2( sqrt( dy**2 + dz**2 ) , dx )) / m_e
        Electron_acc(2) = ( Fe + Fg ) * cos(atan2( sqrt( dx**2 + dz**2 ) , dy )) / m_e
        Electron_acc(3) = ( Fe + Fg ) * cos(atan2( sqrt( dx**2 + dy**2 ) , dz )) / m_e

        Proton_acc(1) = ( -Fe - Fg ) * cos(atan2( sqrt( dy**2 + dz**2 ) , dx )) / m_p
        Proton_acc(2) = ( -Fe - Fg ) * cos(atan2( sqrt( dx**2 + dz**2 ) , dy )) / m_p
        Proton_acc(3) = ( -Fe - Fg ) * cos(atan2( sqrt( dx**2 + dy**2 ) , dz )) / m_p

        E = 0.5 * m_e * dot_product( Electron_vel , Electron_vel ) + &
        & 0.5 * p_e * dot_product( Proton_vel , Proton_vel ) + &
        & k * q_e * q_p / r - G * m_e * m_p / r

        WRITE( 2 , * ) E

        Electron_pos = Electron_pos + Electron_vel * dT + 0.5 * Electron_acc * dT**2

        Proton_pos = Proton_pos + Proton_vel * dT + 0.5 * Proton_acc * dT**2
        
        Electron_acc2(1) = ( Fe + Fg ) * cos(atan2( sqrt( dy**2 + dz**2 ) , dx )) / m_e
        Electron_acc2(2) = ( Fe + Fg ) * cos(atan2( sqrt( dx**2 + dz**2 ) , dy )) / m_e
        Electron_acc2(3) = ( Fe + Fg ) * cos(atan2( sqrt( dx**2 + dy**2 ) , dz )) / m_e

        Proton_acc2(1) = ( -Fe - Fg ) * cos(atan2( sqrt( dy**2 + dz**2 ) , dx )) / m_p
        Proton_acc2(2) = ( -Fe - Fg ) * cos(atan2( sqrt( dx**2 + dz**2 ) , dy )) / m_p
        Proton_acc2(3) = ( -Fe - Fg ) * cos(atan2( sqrt( dx**2 + dy**2 ) , dz )) / m_p

        Electron_vel = Electron_vel + 0.5 * ( Electron_Acc + Electron_Acc ) * dT

        Proton_vel = Proton_vel + 0.5 * ( Proton_Acc + Proton_Acc ) * dT

    END DO

    STOP

END PROGRAM MAIN