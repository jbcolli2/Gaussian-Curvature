from dolfin import *
from GC_Problems import *
import numpy as np
from VM_Utilities import *
from VM_Solver import *

set_log_level(20)


#Values of N for the mesh
params = np.array([4, 8, 16,32]);
# params = np.array([32]);

L = len(params);
e = np.zeros([L,1]);
ratio = np.zeros([L,1]);

p = 2;

ep = np.array([0,1, 1e-1, 5e-2, 2e-2, 1e-2, 7e-3, 5e-3, 2e-3, 1e-3, 7e-4, 6e-4, 5.5e-4, 2e-4]);
ep = np.array([0]);

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
    (x0, y0, x1, y1, exact, gx, gy, K) = GC_Problems(prob, N);



    # Create mesh and define function space
    mesh = RectangleMesh(Point(x0,y0),Point(x1,y1),N,N)
    V = FunctionSpace(mesh, 'Lagrange', p)
    MixedV = MixedFunctionSpace([V,V,V,V]);


    SetupUtilities(mesh, x0, x1, y0, y1);



    ds = Create_dsMeasure()
    
    
    class MixedExact(Expression):
        def eval(self, values, x):
            values[0] = ( x[1]**2 - 1)/( pow( 1 - x[0]**2, 3.0/2.0) )
            values[1] = ( x[0]*x[1])/( pow( 1 - x[0]**2, 3.0/2.0) )
            values[2] = ( x[0]**2 - 1)/( pow( 1 - x[0]**2, 3.0/2.0) )
            values[3] = sqrt(1-pow(x[0],2) - x[1]**2)
        def value_shape(self):
            return (4,)


    u = Function(MixedV)
    ex = MixedExact();
    u.interpolate(ex)



    ##### Loop through epsilon values and solve ####################
    w = Function(MixedV);
    w = u;
    ep_err = [];

    for epjj in ep:
        print('Epsilon = ',epjj)

        w = ForwardProblem_GC(MixedV,K,ds, epjj, w, exact, gx, gy)

        (Sxx,Sxy,Syy,u) = w.split(deepcopy=True);

        ep_err.append(np.sqrt(assemble(abs(exact-u)**2*dx)));
        print('Run finished at epsilon = ', epjj)
        print('L2 error = ', ep_err[-1])
        print ''
  
  


  



    (Sxx,Sxy,Syy,u) = w.split(deepcopy=True);

    error = abs(exact-u)**2*dx
    e[ii] = np.sqrt(assemble(error))
    
    if(ii > 0):
        ratio[ii] = np.log(e[ii-1]/e[ii])/np.log(2)
 
# plot(abs(exact-u))
# interactive()
print("Error: ", e)
print("Ratio: ", ratio)
