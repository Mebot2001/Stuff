PROGRAM MAIN

    ! A THING will be an array of 8 elements.
        ! (1) Mass
        ! (2) Charge
        ! (3) X-position
        ! (4) Y-position
        ! (5) Z-position
        ! (6) X-velocity
        ! (7) Y-velocity
        ! (8) Z-velocity
    !

    REAL*16 Electron(8) , Proton(8) , dT , G , k , r , Fe , Fg , dx , dy , dz , T
    INTEGER*8 count , n

    G = 6.6743E-11
    k = 8.987551782E9

    Electron(1) = 9.1093837015E-31
    Electron(2) = -1.602176634E-19

    Proton(1) = 1.67262192369E-27
    Proton(2) = 1.602176634E-19

    Electron(3) = 1
    Electron(4) = 0
    Electron(5) = 0

    Electron(6) = 0
    Electron(7) = 14
    Electron(8) = 0

    Proton(6) = 0
    Proton(7) = 0
    Proton(8) = 0

    Proton(3) = 0
    Proton(4) = 0
    Proton(5) = 0

    T = 1E-4
    dT = 1E-9

    n = T/dT

    OPEN( unit = 1 , file = "data.csv" )
    OPEN( unit = 2 , file = "data2.csv" )

    DO count = 1 , n

        WRITE( 1 , * ) Electron(3) , Electron(4) , Electron(5) , Proton(3) , Proton(4) , Proton(5)

        dx = Electron(3) - Proton(3)
        dy = Electron(4) - Proton(4)
        dz = Electron(5) - Proton(5)
        
        r = sqrt( dx**2 + dy**2 + dz**2 )

        Fe = k * Electron(2) * Proton(2) / r**2
        Fg = -G * Electron(1) * Proton(1) / r**2

        WRITE( 2 , * ) Fe , Fg

        Electron(3) = Electron(3) + Electron(6) * dT + cos(atan2(sqrt(dy**2+dz**2),dx)) * ( Fe + Fg ) * dT**2 / Electron(1) / 2
        Electron(4) = Electron(4) + Electron(7) * dT + cos(atan2(sqrt(dx**2+dz**2),dy)) * ( Fe + Fg ) * dT**2 / Electron(1) / 2
        Electron(5) = Electron(5) + Electron(8) * dT + cos(atan2(sqrt(dx**2+dy**2),dz)) * ( Fe + Fg ) * dT**2 / Electron(1) / 2

        Electron(6) = Electron(6) + cos(atan2(sqrt(dy**2+dz**2),dx)) * ( Fe + Fg ) * dT / Electron(1)
        Electron(7) = Electron(7) + cos(atan2(sqrt(dx**2+dz**2),dy)) * ( Fe + Fg ) * dT / Electron(1)
        Electron(8) = Electron(8) + cos(atan2(sqrt(dx**2+dy**2),dz)) * ( Fe + Fg ) * dT / Electron(1)

        Proton(3) = Proton(3) + Proton(6) * dT - cos(atan2(sqrt(dy**2+dz**2),dx)) * ( Fe + Fg ) * dT**2 / Proton(1) / 2
        Proton(4) = Proton(3) + Proton(6) * dT - cos(atan2(sqrt(dx**2+dz**2),dy)) * ( Fe + Fg ) * dT**2 / Proton(1) / 2
        Proton(5) = Proton(3) + Proton(6) * dT - cos(atan2(sqrt(dx**2+dy**2),dz)) * ( Fe + Fg ) * dT**2 / Proton(1) / 2

        Proton(6) = Proton(6) - cos(atan2(sqrt(dy**2+dz**2),dx)) * ( Fe + Fg ) * dT / Proton(1)
        Proton(7) = Proton(7) - cos(atan2(sqrt(dx**2+dz**2),dy)) * ( Fe + Fg ) * dT / Proton(1)
        Proton(8) = Proton(8) - cos(atan2(sqrt(dx**2+dy**2),dz)) * ( Fe + Fg ) * dT / Proton(1)
        
    END DO

    STOP

END PROGRAM MAIN
