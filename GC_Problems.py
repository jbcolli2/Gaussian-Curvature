from dolfin import *


def GC_Problems(prob):
    # cutoff = pow(N,2.0);
    cutoff = 0;
    xtol = 1e-18;

    #u(x,y) = x^4 + y^2
    if(prob == 1):
        x0 = -.57; y0 = -.57; x1 = .57; y1 = .57;
        exact = Expression('sqrt(1-pow(x[0],2) - pow(x[1],2))');
        gx = Expression('-x[0]/exac', exac = exact);
        gy = Expression('-x[1]/exac', exac = exact);
        K = 1;

        return (x0, y0, x1, y1, exact, gx, gy, K);



    elif(prob == 2):
        x0 = -.57; y0 = -.57; x1 = .57; y1 = .57;
        exact = Expression('-sqrt(1-pow(x[0],2) - pow(x[1],2))');
        gx = Expression('x[0]/exac', exac = exact);
        gy = Expression('x[1]/exac', exac = exact);
        K = 1.0;

        return (x0, y0, x1, y1, exact, gx, gy, K);




    elif(prob == 3):
        x0 = -.57; y0 = -.57; x1 = .57; y1 = .57;
        exact = Expression('sqrt(1-pow(x[0],2))');
        gx = Expression('-x[0]/exac', exac = exact);
        gy = Expression('0.0', exac = exact);
        K = 0.0;

        return (x0, y0, x1, y1, exact, gx, gy, K);



    elif(prob == 4):
        x0 = 0; y0 = 0; x1 = 1; y1 = 1;
        exact = Expression('sin(pi*x[0])*sin(1.5*pi*x[1])');
        gx = Expression('pi*cos(pi*x[0])*sin(1.5*pi*x[1])', exac = exact);
        gy = Expression('1.5*pi*sin(pi*x[0])*cos(1.5*pi*x[1])', exac = exact);
        K = 0.5;

        return (x0, y0, x1, y1, exact, gx, gy, K);

    elif(prob == 5):
        x0 = 0; y0 = 0; x1 = 1; y1 = 1;
        exact = Expression('0.0');
        gx = exact
        gy = exact
        K = 0.01;

        return (x0, y0, x1, y1, exact, gx, gy, K);

    elif(prob == 6):
        x0 = 0; y0 = 0; x1 = 1; y1 = 1;
        exact = Expression('x[0]*x[0] + x[1]*x[1]');
        gx = Expression('2.0*x[0]', ex = exact)
        gy = Expression('2.0*x[1]', ex = exact)
        K = 0.01;

        return (x0, y0, x1, y1, exact, gx, gy, K);


