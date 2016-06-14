from dolfin import *
from GC_Problems import *
import numpy as np
from VM_Utilities import *
from VM_Solver import *

set_log_level(20)


#Values of N for the mesh
params = np.array([4, 8,16,32]);
params = np.array([32]);

L = len(params);
e = np.zeros([L,1]);
ratio = np.zeros([L,1]);

p = 2;

ep = np.array([1, 5e-1, 1e-1, 5e-2, 1e-2,5e-3, 1e-3,5e-4,3e-4,2e-4]);
# ep = np.array([1]);
# ep = -ep;
# ep = np.array([1, 1e-1]);

for ii in range(L):
    N = params[ii];
    print(N);



    prob = 4;
    (x0, y0, x1, y1, exact, gx, gy, K) = GC_Problems(prob, N);
    K = 0.0


    # Create mesh and define function space
    mesh = RectangleMesh(Point(x0,y0),Point(x1,y1),N,N)
    V = FunctionSpace(mesh, 'Lagrange', p)
    MixedV = MixedFunctionSpace([V,V,V,V]);


    SetupUtilities(mesh, x0, x1, y0, y1);



    ds = Create_dsMeasure()
    
    
    # class MixedExact(Expression):
    #     def eval(self, values, x):
    #         values[0] = ( x[1]**2 - 1.1)/( pow( 1 - x[0]**2, 3.0/2.0) )
    #         values[1] = ( x[0]*x[1])/( pow( 1.2 - x[0]**2, 3.0/2.0) )
    #         values[2] = ( x[0]**2 - 1)/( pow( 1 - x[0]**2, 3.1/2.0) )
    #         values[3] = sqrt(.8-pow(x[0],2) - x[1]**2) + .2
    #     def value_shape(self):
    #         return (4,)


    # u = Function(MixedV)
    # ex = MixedExact();
    # u.interpolate(ex)



    ##### Loop through epsilon values and solve ####################
    w = Function(MixedV);
    # w = u;
    ep_err = [];

    for epjj in ep:
        print('Epsilon = ',epjj)

        w, prob, solv = ForwardProblem_GC(MixedV,K,ds, epjj, w, exact, gx, gy)

        (Sxx,Sxy,Syy,u) = w.split(deepcopy=True);

        ep_err.append(np.sqrt(assemble(abs(exact-u)**2*dx)));
        print('Run finished at epsilon = ', epjj)
        print('L2 error = ', ep_err[-1])
        print ''

        # PlotToFile(u, 'Epsilon = ' + epjj.__str__(), 'file')
  
  
        bcv = DirichletBC(MixedV.sub(3), exact, Dir_boundary)
        bcxx = DirichletBC(MixedV.sub(0), epjj, EW_boundary)
        bcyy = DirichletBC(MixedV.sub(2), epjj, NS_boundary)
        bc = [bcxx,bcyy,bcv]

        F = F_Form_GC(MixedV,K,ds,epjj,gx,gy);
        R = EvalResidual(F, bc, w)
  
    w, prob, solv = ForwardProblem_GC(MixedV,K+.01,ds, epjj, w, exact, gx, gy)
    w, prob, solv = ForwardProblem_GC(MixedV,K+.02,ds, epjj, w, exact, gx, gy)
    w, prob, solv = ForwardProblem_GC(MixedV,K+.05,ds, epjj, w, exact, gx, gy)
    w, prob, solv = ForwardProblem_GC(MixedV,K+.2,ds, epjj, w, exact, gx, gy)
    w, prob, solv = ForwardProblem_GC(MixedV,K+.3,ds, epjj, w, exact, gx, gy)
    w, prob, solv = ForwardProblem_GC(MixedV,K+.4,ds, epjj, w, exact, gx, gy)
    w, prob, solv = ForwardProblem_GC(MixedV,K+.5,ds, epjj, w, exact, gx, gy)
    w, prob, solv = ForwardProblem_GC(MixedV,K+.55,ds, epjj, w, exact, gx, gy)
    # w, prob, solv = ForwardProblem_GC(MixedV,K+.57,ds, epjj, w, exact, gx, gy)


    (Sxx,Sxy,Syy,u) = w.split(deepcopy=True);

    error = abs(exact-u)**2*dx
    e[ii] = np.sqrt(assemble(error))
    
    if(ii > 0):
        ratio[ii] = np.log(e[ii-1]/e[ii])/np.log(2)

    # plot(exact-u)
    # interactive()
 
# plot(abs(exact-u))
# interactive()
print("Error: ", e)
print("Ratio: ", ratio)
