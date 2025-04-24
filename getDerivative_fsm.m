function rhs = getDerivative_fsm(delta_t, m, x1, x2, x3, x4, x5, x6, theta, phi, psi, T) 
    grav = 9.81;
    r1 = x2;
    r2 = -1/m*(cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi))*T;
    r3 = x4;
    r4 = -1/m*(cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi))*T;
    r5 = x6;
    r6 = grav - 1/m*cos(phi)*cos(theta)*T;
    rhs = [r1; r2; r3; r4; r5; r6];

