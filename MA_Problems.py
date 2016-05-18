from dolfin import *


def Problems(prob, N):
    # cutoff = pow(N,2.0);
    cutoff = 0;
    xtol = 1e-18;

    #u(x,y) = x^4 + y^2
    if(prob == 1):
        x0 = 0; y0 = 0; x1 = 1; y1 = 1;
        exact = Expression('pow(x[0],4.0) + pow(x[1],2.0)');
        f = Expression('24.0*pow(x[0],2.0)');
        gx = Expression('4.0*pow(x[0],3.0)');
        gy = Expression('2.0*x[1]');

        return (x0, y0, x1, y1, exact, f, gx, gy);

    #u(x,y) = exp(.5*(x^2+y^2))
    elif(prob == 2):
        x0 = 0; y0 = 0; x1 = 1; y1 = 1;
        exact = Expression('exp( 0.5*(pow(x[0],2.0) + pow(x[1],2.0)) )');
        f = Expression('exp( (pow(x[0],2.0) + pow(x[1],2.0)) )* (pow(x[0],2.0) + pow(x[1],2.0)+1.0)');
        gx = Expression('x[0]*exp( 0.5*(pow(x[0],2.0) + pow(x[1],2.0)) )');
        gy = Expression('x[1]*exp( 0.5*(pow(x[0],2.0) + pow(x[1],2.0)) )');

        return (x0, y0, x1, y1, exact, f, gx, gy);



    #u(x,y) = (1/3)(4x^2 + 4y^2)^(3/4)
    # Full domain, function cutoff
    elif(prob == 5):
        x0 = 0; y0 = 0; x1 = 1; y1 = 1;
        exact = Expression('(1.0/3.0)*pow(4*pow(x[0],2.0) + 4*pow(x[1],2.0),(3.0/4.0))');
        class Sing_f(Expression):
            def eval(self, value, x):
                temp = pow(pow(x[0],2.0) + pow(x[1],2.0), (-1.0/2.0));
                if(abs(x[0]) +abs(x[1]) < xtol):
                    value[0] = cutoff
                else:
                    value[0] = temp

        class Sing_gx(Expression):
            def eval(self, value, x):
                temp = 2*x[0]*pow(4*pow(x[0],2.0) + 4*pow(x[1],2.0), (-1.0/4.0));
                if(abs(x[0]) +abs(x[1]) < xtol):
                    value[0] = cutoff
                else:
                    value[0] = temp

        class Sing_gy(Expression):
            def eval(self, value, x):
                temp = 2*x[1]*pow(4*pow(x[0],2.0) + 4*pow(x[1],2.0), (-1.0/4.0));
                if(abs(x[0]) +abs(x[1]) < xtol):
                    value[0] = cutoff
                else:
                    value[0] = temp

        f = Sing_f();
        gx = Sing_gx();
        gy = Sing_gy();

        return (x0, y0, x1, y1, exact, f, gx, gy);


    #u(x,y) = -sqrt(2 - x^2 - y^2)
    # Full domain, function cutoff
    elif(prob == 6):
        x0 = 0; y0 = 0; x1 = 1; y1 = 1;
        exact = Expression('-pow(2.0 - pow(x[0],2.0) - pow(x[1],2.0), 1.0/2.0)');
        class Sing_f1(Expression):
            def eval(self, value, x):
                temp = 2.0*pow(2.0 - pow(x[0],2.0) - pow(x[1],2.0), -2.0);
                if(abs(x[0] - 1) < xtol and abs(x[1] - 1) < xtol):
                    value[0] = cutoff
                else:
                    value[0] = temp
        f = Sing_f1()
        class Sing_gx1(Expression):
            def eval(self, value, x):
                temp = x[0]*pow(2.0 - pow(x[0],2.0) - pow(x[1],2.0), -1.0/2.0);
                if(abs(x[0] - 1) < xtol and abs(x[1] - 1) < xtol):
                    value[0] = cutoff
                else:
                    value[0] = temp
        gx = Sing_gx1()
        class Sing_gy1(Expression):
            def eval(self, value, x):
                temp = x[1]*pow(2.0 - pow(x[0],2.0) - pow(x[1],2.0), -1.0/2.0);
                if(abs(x[0] - 1) < xtol and abs(x[1] - 1) < xtol):
                    value[0] = cutoff
                else:
                    value[0] = temp
        gy = Sing_gy1()

        return (x0, y0, x1, y1, exact, f, gx, gy);



    #u(x,y) = abs(x)
    # Full domain, full function
    elif(prob == 7):
        x0 = -1; y0 = -1; x1 = 1; y1 = 1;
        class Sing_u2(Expression):
            def eval(self, value, x):
                if(x[0] < 0.0):
                    value[0] = -x[0]
                else:
                    value[0] = x[0]
        exact = Sing_u2()
        f = Expression('0.0')
        class Sing_gx2(Expression):
            def eval(self, value, x):
                if(x[0] < 0.0):
                    value[0] = -1.0
                else:
                    value[0] = 1.0
        gx = Sing_gx2()
        gy = Expression('0.0')

        return (x0, y0, x1, y1, exact, f, gx, gy);



    #u(x,y) = x/x^2 piecewise function
    elif(prob == 8):
        x0 = -1; y0 = -1; x1 = 1; y1 = 1;
        class Sing_u2(Expression):
            def eval(self, value, x):
                if(x[0] < 0.0):
                    value[0] = -x[0]
                else:
                    value[0] = x[0]**2
        exact = Sing_u2()
        f = Expression('0.0')
        class Sing_gx2(Expression):
            def eval(self, value, x):
                if(x[0] < 0.0):
                    value[0] = -1.0
                else:
                    value[0] = 2.0*x[0]
        gx = Sing_gx2()
        gy = Expression('0.0')

        return (x0, y0, x1, y1, exact, f, gx, gy);



    #u(x,y) = sqrt(x^2 + y^2)
    # numerical Dirac delta function
    elif(prob == 9):
        x0 = -1; y0 = -1; x1 = 1; y1 = 1;
        exact = Expression('sqrt(pow(x[0],2.0) + pow(x[1],2.0))')
        class Sing_f3(Expression):
            def eval(self, value, x):
                if(abs(x[0]) < cutoff and abs(x[1]) < cutoff):
                    value[0] = 4*cutoff;
                else:
                    value[0] = 0.0;
        f = Sing_f3()
        class Sing_gx3(Expression):
            def eval(self, value, x):
                if(abs(x[0]) < xtol and abs(x[1]) < xtol):
                    value[0] = cutoff;
                else:
                    value[0] = x[0]/(sqrt(x[0]**2 + x[1]**2))

        class Sing_gy3(Expression):
            def eval(self, value, x):
                if(abs(x[0]) < xtol and abs(x[1]) < xtol):
                    value[0] = cutoff;
                else:
                    value[0] = x[1]/(sqrt(x[0]**2 + x[1]**2))

        gx = Sing_gx3()
        gy = Sing_gy3()

        return (x0, y0, x1, y1, exact, f, gx, gy);

