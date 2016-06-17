from dolfin import *
from GC_Problems import *
import numpy as np
from VM_Utilities import *
from VM_Solver import *

set_log_level(20)


#Values of N for the mesh
params = np.array([4, 8, 16,32]);
params = np.array([16]);

L = len(params);
e = np.zeros([L,1]);
ratio = np.zeros([L,1]);

p = 2;

ep = np.array([ 1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 0]);
ep = np.array([1]);
ep = -ep;
# ep = np.array([1, 1e-1]);

for ii in range(L):
    N = params[ii];
    print(N);



    prob = 5;
    (x0, y0, x1, y1, exact, gx, gy, K) = GC_Problems(prob, N);
    K = 1.0;


    # Create mesh and define function space
    mesh = RectangleMesh(Point(x0,y0),Point(x1,y1),N,N)
    V = FunctionSpace(mesh, 'Lagrange', p)
    MixedV = MixedFunctionSpace([V,V,V,V]);


    SetupUtilities(mesh, x0, x1, y0, y1);



    ds = Create_dsMeasure()
    
    
    # class MixedExact(Expression):
    #     def eval(self, values, x):
    #         values[0] = ( x[1]**2 - 1.0)/( pow( 1 - x[0]**2 - x[1]**2, 3.0/2.0) )
    #         values[1] = ( -x[0]*x[1])/( pow( 1 - x[0]**2- x[1]**2, 3.0/2.0) )
    #         values[2] = ( x[0]**2 - 1.0)/( pow( 1 - x[0]**2- x[1]**2, 3.1/2.0) )
    #         values[3] = sqrt(1.0-pow(x[0],2) - x[1]**2)
    #     def value_shape(self):
    #         return (4,)


    # u = Function(MixedV)
    # ex = MixedExact();
    # u.interpolate(ex)



    ##### Loop through epsilon values and solve ####################
    w = Function(MixedV);
    # w = u;
    ep_err = [];

    
    myinitial = Function(MixedV)
    (Sxx,Sxy,Syy,u) = myinitial.split();
    Sxx.vector()[:] = 1.0;
    Syy.vector()[:] = 1.0;
    u.vector()[:] = 1.0;
    feninitial = Function(MixedV)
    for epjj in ep:
        print 'Epsilon = ',epjj
        #Define boundary conditions
        bcv = DirichletBC(MixedV.sub(3), exact, Dir_boundary)
        bcxx = DirichletBC(MixedV.sub(0), epjj, EW_boundary)
        bcyy = DirichletBC(MixedV.sub(2), epjj, NS_boundary)
        bc = [bcxx,bcyy,bcv]
        bcvh = DirichletBC(MixedV.sub(3), exact, Dir_boundary)
        bcxxh = DirichletBC(MixedV.sub(0), epjj, EW_boundary)
        bcyyh = DirichletBC(MixedV.sub(2), epjj, NS_boundary)
        bch = [bcxxh,bcyyh,bcvh]
        for bcii in bch:
            bcii.homogenize();


        F = F_Form_GC(MixedV, K, ds, epjj, gx, gy);

        # Solve the problem

        # myinitial.assign(NewtonIteration(MixedV, feninitial, F, bc, bch));

        print 'Initial Residual = ', EvalResidual(F, bc, feninitial).norm('l2');
        R = action(F,feninitial);
        DR = derivative(R, feninitial);
        bcv = DirichletBC(MixedV.sub(3), exact, Dir_boundary)
        bcxx = DirichletBC(MixedV.sub(0), epjj, EW_boundary)
        bcyy = DirichletBC(MixedV.sub(2), epjj, NS_boundary)
        bc = [bcxx,bcyy,bcv]
        problem = NonlinearVariationalProblem(R,feninitial,bc,DR);
        solver = NonlinearVariationalSolver(problem);
        solver.parameters['newton_solver']['maximum_iterations'] = 5;
        solver.parameters['newton_solver']['absolute_tolerance'] = 1e-10;
        solver.solve();

        F = action(F,myinitial)
        solve(F == 0, myinitial, bc, solver_parameters={"newton_solver":
                                        {"relative_tolerance": 1e-6}})



        print 'Computed solutions: Fenics = ', feninitial.vector().norm('l2'), ',  Mine = ', myinitial.vector().norm('l2');
        print ''


  
  

