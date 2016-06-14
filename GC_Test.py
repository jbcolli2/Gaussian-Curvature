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

ep = np.array([ 1]);
ep = -ep;
# ep = np.array([1, 1e-1]);

for ii in range(L):
    N = params[ii];
    print(N);



    prob = 1;
    (x0, y0, x1, y1, exact, gx, gy, K) = GC_Problems(prob, N);



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

    print('Epsilon = ',1)

    
    #Define boundary conditions
    bcv = DirichletBC(MixedV.sub(3), exact, Dir_boundary)
    bcxx = DirichletBC(MixedV.sub(0), 1, EW_boundary)
    bcyy = DirichletBC(MixedV.sub(2), 1, NS_boundary)
    bc = [bcxx,bcyy,bcv]

    #Define the weak residual F
    (Sxx, Sxy, Syy, u) = TrialFunction(MixedV)
    (muxx, muxy, muyy, v) = TestFunction(MixedV)

    F = inner(Sxx,muxx)*dx + 2*inner(Sxy,muxy)*dx + inner(Syy,muyy)*dx;
    F += inner(Dx(u,0), Dx(muxx,0))*dx + inner(Dx(u,0), Dx(muxy,1))*dx;
    F += inner(Dx(u,1), Dx(muxy,0))*dx + inner(Dx(u,1), Dx(muyy,1))*dx;

    F += ( inner(Dx(Sxx,0), Dx(v,0)) + inner(Dx(Sxy,0), Dx(v,1)))*dx;
    F += ( inner(Dx(Sxy,1), Dx(v,0)) + inner(Dx(Syy,1), Dx(v,1)))*dx;

    F += (((Sxx*Syy - Sxy*Sxy)-(1 + (Dx(u,0)**2 + Dx(u,1)**2))**(2)) * K)*v*dx;

    F -= (-gy*muxy*ds(1) + gx*muxy*ds(2) + gy*muxy*ds(3) - gx*muxy*ds(4));

    # Solve the problem
    w0 = Function(MixedV)
    initial = Function(MixedV)
    print 'Initial Residual = ', EvalResidual(F, bc, initial).norm('l2');
    R = action(F,initial);
    DR = derivative(R, initial);
    problem = NonlinearVariationalProblem(R,initial,bc,DR);
    solver = NonlinearVariationalSolver(problem);
    solver.parameters['newton_solver']['maximum_iterations'] = 4;
    solver.parameters['newton_solver']['absolute_tolerance'] = 1e-10;
    solver.solve();

    mysol = NewtonIteration(MixedV, w0, F, bc);

    print 'Computed solutions: Fenics = ', initial.vector().norm('l2'), ',  Mine = ', mysol.vector().norm('l2');


  
  

