from dolfin import *
from GC_Problems import *
import numpy as np
from VM_Utilities import *
from VM_Solver import *

set_log_level(20)


#Values of N for the mesh
params = np.array([4, 8, 16,32]);
# params = np.array([4]);

L = len(params);
e = np.zeros([L,1]);
ratio = np.zeros([L,1]);

p = 2;

ep = np.array([1, 1e-2, 1e-4, 0]);
ep = -ep;
# ep = np.array([1, 1e-1]);

for ii in range(L):
    N = params[ii];
    print(N);






    # Create mesh and define function space
    x0 = 0; x1 = 1; y0 = 0; y1 = 1;
    mesh = RectangleMesh(Point(x0,y0),Point(x1,y1),N,N)
    V = FunctionSpace(mesh, 'Lagrange', p)
    exact = Expression('sin(pi*x[0])*sin(pi*x[1])', domain=mesh)
    f = Expression('2*pi*pi*sin(pi*x[0])*sin(pi*x[1]) - pow( sin(pi*x[0])*sin(pi*x[1]), 2.0 )', domain=mesh)


    
    



    ##### Loop through epsilon values and solve ####################
    bc = DirichletBC(V, 0.0, Dir_boundary)
    u = TrialFunction(V);
    v = TestFunction(V);

    F = inner(grad(u), grad(v))*dx - u*u*v*dx - f*v*dx;

    u0 = interpolate(Expression(' 1.5*sin(pi*2*x[0])*sin(pi*x[1])'), V);

    u0 = NewtonIteration(V, u0, F, bc)

    
    # R = action(F,u0);
    # DR = derivative(R, u0);
    # problem = NonlinearVariationalProblem(R,u0,bc,DR);
    # solver = NonlinearVariationalSolver(problem);
    # # solver.parameters['newton_solver']['absolute_tolerance'] = 1e-9
    # prm = solver.parameters
    # solver.solve();


    # s = SystemAssembler(DR,R,bc)
    # b = Vector();
    # s.assemble(b);





    error = abs(exact-u0)**2*dx
    e[ii] = np.sqrt(assemble(error))
    
    if(ii > 0):
        ratio[ii] = np.log(e[ii-1]/e[ii])/np.log(2)

    # plot(exact-u)
    # interactive()
 
# plot(abs(exact-u))
# interactive()
print("Error: ", e)
print("Ratio: ", ratio)
