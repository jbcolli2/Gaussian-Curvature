from dolfin import *


def GC_Problems(prob, N):
    # cutoff = pow(N,2.0);
    cutoff = 0;
    xtol = 1e-18;

    #u(x,y) = x^4 + y^2
    if(prob == 1):
        x0 = -.57; y0 = -.57; x1 = .57; y1 = .57;
        exact = Expression('sqrt(1-pow(x[0],2) - pow(x[1],2))');
        gx = Expression('-x[0]*pow(1-pow(x[0],2) - pow(x[1],2), -1.0/2.0)');
        gy = Expression('-x[1]*pow(1-pow(x[0],2) - pow(x[1],2), -1.0/2.0)');
        K = 1;

        return (x0, y0, x1, y1, exact, gx, gy, K);




    elif(prob == 2):
        ep = .5;
        x0 = -.57; y0 = -.57; x1 = .57; y1 = .57;
        exact = Expression('eps*pow(x[0],2)/2 + eps*pow(x[1],2)/2', eps = ep);
        gx = Expression('eps*x[0]', eps=ep);
        gy = Expression('eps*x[1]', eps = ep);
        K = .345;
        f = Expression('-Kparam + pow(eps,2)*pow( pow(eps,2) * ( pow(x[0],2) + pow(x[1],2) ) + 1, -2.0)', eps=ep, Kparam = K)

        return (x0, y0, x1, y1, exact,f, gx, gy, K);


    elif(prob == 3):
        ep = 1;
        x0 = 0; y0 = 0; x1 = 1; y1 = 1;
        exact = Expression('-eps*( cos(x[0]*pi)+cos(x[1]*pi) )/(pi*pi)', eps = ep);
        gx = Expression('(eps/pi)*(sin(pi*x[0])', eps=ep);
        gy = Expression('(eps/pi)*(sin(pi*x[1])', eps = ep);
        K = 1;
        f = Expression(('-eps * (-eps * cos(pi * x[0]) * pi * pi - \
            eps * cos(pi * x[1]) * pi * pi) + eps * eps * cos(pi * x[0]) *\
             cos(pi * x[1]) * pow(0.1e1 + eps * eps * pow(pi, -0.2e1) * pow(sin(pi * x[0]), 0.2e1)\
              + eps * eps * pow(pi, -0.2e1) * pow(sin(pi * x[1]), 0.2e1), -0.2e1) - Kparam', eps = ep, Kparam = K)

        return (x0, y0, x1, y1, exact,f, gx, gy, K);




    elif(prob == 4):
        ep = 1;
        x0 = 0; y0 = 0; x1 = 1; y1 = 1;
        exact = Expression('-eps*( cos(x[0]*pi)+cos(x[1]*pi) )/(pi*pi)', eps = ep);
        gx = Expression('(eps/pi)*(sin(pi*x[0])', eps=ep);
        gy = Expression('(eps/pi)*(sin(pi*x[1])', eps = ep);
        K = 1;
        f = Expression(('-eps * (-eps * cos(pi * x[0]) * pi * pi - \
            eps * cos(pi * x[1]) * pi * pi) + eps * eps * cos(pi * x[0]) *\
             cos(pi * x[1]) * pow(0.1e1 + eps * eps * pow(pi, -0.2e1) * pow(sin(pi * x[0]), 0.2e1)\
              + eps * eps * pow(pi, -0.2e1) * pow(sin(pi * x[1]), 0.2e1), -0.2e1) - Kparam', eps = ep, Kparam = K)

        return (x0, y0, x1, y1, exact,f, gx, gy, K);


  
