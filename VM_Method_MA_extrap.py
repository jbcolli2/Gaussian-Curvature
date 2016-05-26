from dolfin import *
from MA_Problems import *
import numpy as np
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

ep = np.array([1, 1e-2, 1e-3 ,6e-4, 5e-4, 4e-4, 3e-4]);

for ii in range(L):
    N = params[ii];
    print(N);



    # Define Problem
    # 1. u(x,y) = x^4 + y^2
    # 2. u(x,y) = exp(.5*(x^2+y^2))
    # 5. u(x,y) = (1/3)(4x^2 + 4y^2)^(3/4)
    # #       Full domain, function cutoff
    # 6. u(x,y) = -sqrt(2 - x^2 - y^2)
    # #       Full domain, function cutoff
    # 7. u(x,y) = abs(x)
    # 8. u(x,y) = x/x^2 piecewise function
    # 9. u(x,y) = sqrt(x^2 + y^2)
    # #       numerical Dirac delta function
    prob = 7;
    (x0, y0, x1, y1, exact, f, gx, gy) = Problems(prob, N);



    # Create mesh and define function space
    mesh = RectangleMesh(Point(x0,y0),Point(x1,y1),N,N)
    V = FunctionSpace(mesh, 'Lagrange', p)
    MixedV = MixedFunctionSpace([V,V,V,V]);


    SetupUtilities(mesh, x0, x1, y0, y1);



    ds = Create_dsMeasure()






    ##### Loop through epsilon values and solve ####################
    w = Function(MixedV);
    sols = []; ep_err = []; solsw = [];
    for epjj in ep:
        print('Epsilon = ', epjj)

        w = ForwardProblem(MixedV,ds, epjj, w, exact, f, gx, gy)

        # Print out the error at this value of epsilon
        (Sxx,Sxy,Syy,u) = w.split(deepcopy=True);
        ep_err.append( np.sqrt( assemble(abs(exact-u)**2*dx) ) );
        print('Run finished at epsilon = ', epjj)
        print('L2 error = ', ep_err[-1])
        print ''

        # Store solution at this value of epsilon
        sols.append(u); solsw.append(w);




    # Quadratic Extrapolation in epsilon
    u = (-ep[-1 - 1]) * (-ep[-1]) / ((ep[-1 - 2] - ep[-1 - 1]) * (ep[-1 - 2] - ep[-1])) * sols[-1 - 2] + (-ep[
        -1 - 2]) * (-ep[-1]) / ((ep[-1 - 1] - ep[-1 - 2]) * (ep[-1 - 1] - ep[-1])) * sols[-1 - 1] + \
        (-ep[-1 - 1]) * (-ep[-1 - 2]) / ((ep[-1] - ep[-1 - 1]) * (ep[-1] - ep[-1 - 2])) * sols[-1]

    w = (-ep[-1 - 1]) * (-ep[-1]) / ((ep[-1 - 2] - ep[-1 - 1]) * (ep[-1 - 2] - ep[-1])) * solsw[-1 - 2] + (-ep[
        -1 - 2]) * (-ep[-1]) / ((ep[-1 - 1] - ep[-1 - 2]) * (ep[-1 - 1] - ep[-1])) * solsw[-1 - 1] + \
        (-ep[-1 - 1]) * (-ep[-1 - 2]) / ((ep[-1] - ep[-1 - 1]) * (ep[-1] - ep[-1 - 2])) * solsw[-1]

    # Cubic Extrapolation in epsilon
    # u = (-ep[-1-2])*(-ep[-1-1])*(-ep[-1])/((ep[-1-3]-ep[-1-1])*(ep[-1-3]-ep[-1])*(ep[-1-3]-ep[-1-2]))*sols[-1-3] + \
    #     (-ep[-1-3])*(-ep[-1-1])*(-ep[-1])/((ep[-1-2]-ep[-1-1])*(ep[-1-2]-ep[-1])*(ep[-1-2]-ep[-1-3]))*sols[-1-2] + \
    #     (-ep[-1-2])*(-ep[-1-3])*(-ep[-1])/((ep[-1-1]-ep[-1-3])*(ep[-1-1]-ep[-1])*(ep[-1-1]-ep[-1-2]))*sols[-1-1] + \
    #     (-ep[-1-2])*(-ep[-1-1])*(-ep[-1-3])/((ep[-1-0]-ep[-1-1])*(ep[-1-0]-ep[-1-3])*(ep[-1-0]-ep[-1-2]))*sols[-1-0]
    w = interpolate(w[0],MixedV)
    w = ForwardProblem(MixedV,ds, 1e-4, w[0], exact, f, gx, gy)
    (Sxx,Sxy,Syy,u) = w.split(deepcopy=True);

    error = abs(exact - u) ** 2 * dx
    e[ii] = np.sqrt(assemble(error))

    if (ii > 0):
        ratio[ii] = np.log(e[ii - 1] / e[ii]) / np.log(2)

# plot(abs(exact-u))
# interactive()
print("Error: ", e)
print("Ratio: ", ratio)
