from dolfin import *
from MA_Problems import *
import numpy as np
import sys;
from VM_Utilities import *
from VM_Solver import *


set_log_level(50)

# Values of N for the mesh
params = np.array([4, 8, 16, 32]);
params = np.array([16]);

L = len(params);
e = np.zeros([L, 1]);
ratio = np.zeros([L, 1]);

p = 2;

ep = np.array([1, 5e-3, 1e-5]);

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



    # Create mesh and define function space
    mesh = RectangleMesh(Point(x0,y0),Point(x1,y1),N,N)
    V = FunctionSpace(mesh, 'Lagrange', p)
    MixedV = MixedFunctionSpace([V,V,V,V]);



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



    un = sols[-1]; unm1 = sols[-1 -1]; unm2 = sols[-1 -2];
    Sxxn = solsSxx0[-1]; Sxxnm1 = solsSxx0[-1 -1]; Sxxnm2 = solsSxx0[-1 -2];
    Sxyn = solsSxy0[-1]; Sxynm1 = solsSxy0[-1 -1]; Sxynm2 = solsSxy0[-1 -2];
    Syyn = solsSyy0[-1]; Syynm1 = solsSyy0[-1 -1]; Syynm2 = solsSyy0[-1 -2];
    epn = ep[-1]; epnm1 = ep[-1 -1]; epnm2 = ep[-1 -2];
    wn = solsw[-1]; wnm1 = solsw[-1 - 1];


    dudep_n = SolveDavidenko(MixedV, ds, wn, epn, exact, f, gx, gy);
    dudep_nm1 = SolveDavidenko(MixedV, ds, wnm1, epnm1, exact, f, gx, gy);






    #Hermite cubic
    A11 = unm1;
    A12 = dudep_nm1;
    A22 = (un - unm1)/(epn - epnm1);
    A32 = dudep_n;
    A13 = (A22 - A12)/(epn - epnm1);
    A23 = (A32 - A22)/(epn - epnm1);
    A14 = (A23 - A13)/(epn - epnm1);

    epeval = 0;
    u = A11 + (epeval-epnm1)*A12 + ((epeval-epnm1)**2)*A13 + ((epeval-epnm1)**2)*(epeval-epn)*A14;






    error = abs(exact - u) ** 2 * dx
    e[ii] = np.sqrt(assemble(error))

    if (ii > 0):
        ratio[ii] = np.log(e[ii - 1] / e[ii]) / np.log(2)

# plot(abs(exact-u))
# interactive()
print("Error: ", e)
print("Ratio: ", ratio)
