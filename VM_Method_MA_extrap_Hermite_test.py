from dolfin import *
from MA_Problems import *
import numpy as np
import sys;
from VM_Utilities import *
from VM_Solver import *


set_log_level(50)

# Values of N for the mesh
params = np.array([4, 8, 16, 32]);
# params = np.array([16]);

L = len(params);
e = np.zeros([L, 1]);
ratio = np.zeros([L, 1]);

p = 2;

ep = np.array([1e-1]);

for ii in range(L):
    N = params[ii];
    print(N);



    # Define Problem
    # 1. u(x,y) = x^4 + y^2
    # 2. u(x,y) = exp(.5*(x^2+y^2))
    # 3. u(x,y) = (1/3)(4x^2 + 4y^2)^(3/4)
    # #       Full domain, no function cutoff
    # 4. u(x,y) = (1/3)(4x^2 + 4y^2)^(3/4)
    # #       cutoff domain, no function cutoff
    # 5. u(x,y) = (1/3)(4x^2 + 4y^2)^(3/4)
    # #       Full domain, function cutoff
    # 6. u(x,y) = -sqrt(2 - x^2 - y^2)
    # #       Full domain, function cutoff
    # 7. u(x,y) = abs(x)
    # 8. u(x,y) = x/x^2 piecewise function
    # 9. u(x,y) = sqrt(x^2 + y^2)
    # #       numerical Dirac delta function
    prob = 1;
    (x0, y0, x1, y1, exact, f, gx, gy) = Problems(prob, N);

    exact = Expression('eps*(x[0]*x[0] + x[1]*x[1])*0.5', eps = ep[-1])
    dexact = Expression('(x[0]*x[0] + x[1]*x[1])*0.5')
    gx = Expression('eps*x[0]', eps = ep[-1])
    gy = Expression('eps*x[1]', eps = ep[-1])
    f = Expression('eps*eps', eps = ep[-1]);

    # Create mesh and define function space
    mesh = RectangleMesh(Point(x0,y0),Point(x1,y1),N,N)
    V = FunctionSpace(mesh, 'Lagrange', p)
    MixedV = MixedFunctionSpace([V,V,V,V]);

    finemesh = RectangleMesh(Point(x0,y0),Point(x1,y1),N,N)
    fineV = FunctionSpace(finemesh, 'Lagrange', p)
    fineMixedV = MixedFunctionSpace([fineV,fineV,fineV,fineV]);



    SetupUtilities(mesh, x0, x1, y0, y1);



    ds = Create_dsMeasure()






    ##### Loop through epsilon values and solve ####################
    w = Function(MixedV);
    ep_err = [];
    sols = []; solsSxx0 = []; solsSxy0 = []; solsSyy0 = []; solsw = [];
    for epjj in ep:
        print('Epsilon = ', epjj)

        w = ForwardProblemWithBoundary(MixedV,ds, epjj, w, exact, f, gx, gy)

        # Print out the error at this value of epsilon
        (Sxx0,Sxy0,Syy0,u) = w.split(deepcopy=True);
        ep_err.append( np.sqrt( assemble(abs(exact-u)**2*dx) ) );
        print('Run finished at epsilon = ', epjj)
        print('L2 error = ', ep_err[-1])
        print ''

        sols.append(u); solsSxx0.append(Sxx0); solsSxy0.append(Sxy0); solsSyy0.append(Syy0);
        solsw.append(w);


    SetupUtilities(finemesh, x0, x1, y0, y1);
    ds = Create_dsMeasure();

    wfine = interpolate(w,fineMixedV)
    dudep_n = SolveDavidenko(fineMixedV, ds, w, ep[-1], exact, f, gx, gy);





    error = abs(dexact - dudep_n) ** 2 * dx
    e[ii] = np.sqrt(assemble(error))

    if (ii > 0):
        ratio[ii] = np.log(e[ii - 1] / e[ii]) / np.log(2)

# plot(abs(exact-u))
# interactive()
print("Error: ", e)
print("Ratio: ", ratio)
