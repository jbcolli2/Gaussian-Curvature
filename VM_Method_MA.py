from dolfin import *
import numpy as np

set_log_level(16)


#Values of N for the mesh
params = np.array([4, 8,16,32]);
# params = np.array([8]);

L = len(params);
e = np.zeros([L,1]);
ratio = np.zeros([L,1]);

p = 2;

ep = np.array([1, 1e-1, 1e-2, 1e-3]);
# ep = np.logspace(0,-3,5)

for ii in range(L):
    N = params[ii];
    print(N);
    
    # Create mesh and define function space
    x0 = -1.0; y0 = -1.0; x1 = 1.0; y1 = 1.0;
    mesh = RectangleMesh(Point(x0,y0),Point(x1,y1),N,N)
    V = FunctionSpace(mesh, 'Lagrange', p)
    MixedV = MixedFunctionSpace([V,V,V,V]);



    # Boundaries for the W_xx, W_yy and V spaces
    def Wxx_boundary(x, on_boundary):
        return near(x[0],x0) or near(x[0],x1);
    def Wyy_boundary(x, on_boundary):
        return near(x[1],y0) or near(x[1],y1);
    def V_boundary(x, on_boundary):
        return near(x[1],y0) or near(x[1],y1) or near(x[0],x0) or near(x[0],x1);


    # Boundaries data for integrating <g,\mu>
    class Left(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[0], x0)

    class Right(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[0], x1)

    class Bottom(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[1], y0)

    class Top(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[1], y1)
    left = Left();
    right = Right();
    top = Top();
    bottom = Bottom();
    
    
    
    # Set facet functions and define boundary measures
    boundaries = FacetFunction("size_t", mesh)
    boundaries.set_all(0)
    left.mark(boundaries, 1)
    top.mark(boundaries, 2)
    right.mark(boundaries, 3)
    bottom.mark(boundaries, 4)
    
    ds = Measure("ds")[boundaries]
    
    
    
    
   ################ Problem Definition ############################
    cutoff = pow(N,2.0);
    xtol = 1e-7;

    # exact = Expression('pow(x[0],4.0) + pow(x[1],2.0)');
    # f = Expression('24.0*pow(x[0],2.0)');
    # gx = Expression('4.0*pow(x[0],3.0)');
    # gy = Expression('2.0*x[1]');


    # exact = Expression('exp( 0.5*(pow(x[0],2.0) + pow(x[1],2.0)) )');
    # f = Expression('exp( (pow(x[0],2.0) + pow(x[1],2.0)) )* (pow(x[0],2.0) + pow(x[1],2.0)+1.0)');
    # gx = Expression('x[0]*exp( 0.5*(pow(x[0],2.0) + pow(x[1],2.0)) )');
    # gy = Expression('x[1]*exp( 0.5*(pow(x[0],2.0) + pow(x[1],2.0)) )');

    # exact = Expression('(1.0/3.0)*pow(4*pow(x[0],2.0) + 4*pow(x[1],2.0),(3.0/4.0))');
    # f = Expression('pow(pow(x[0],2.0) + pow(x[1],2.0), (-1.0/2.0))');
    # gx = Expression('2*x[0]*pow(4*pow(x[0],2.0) + 4*pow(x[1],2.0), (-1.0/4.0))');
    # gy = Expression('2*x[1]*pow(4*pow(x[0],2.0) + 4*pow(x[1],2.0), (-1.0/4.0))');

    # exact = Expression('(1.0/5.0)*pow(4*pow(x[0],2.0) + 4*pow(x[1],2.0),(5.0/4.0))');
    # f = Expression('12.0*pow(pow(x[0],2.0) + pow(x[1],2.0), (1.0/2.0))');
    # gx = Expression('2*x[0]*pow(4*pow(x[0],2.0) + 4*pow(x[1],2.0), (1.0/4.0))');
    # gy = Expression('2*x[1]*pow(4*pow(x[0],2.0) + 4*pow(x[1],2.0), (1.0/4.0))');

    # exact = Expression('(1.0/3.0)*pow(4*pow(x[0],2.0) + 4*pow(x[1],2.0),(3.0/4.0))');
    # class Sing_f(Expression):
    #     def eval(self, value, x):
    #         temp = pow(pow(x[0],2.0) + pow(x[1],2.0), (-1.0/2.0));
    #         if(abs(x[0]) < xtol and abs(x[1]) < xtol):
    #             value[0] = cutoff
    #         else:
    #             value[0] = temp
    #
    # class Sing_gx(Expression):
    #     def eval(self, value, x):
    #         temp = 2*x[0]*pow(4*pow(x[0],2.0) + 4*pow(x[1],2.0), (-1.0/4.0));
    #         if(abs(x[0]) < xtol and abs(x[1]) < xtol):
    #             value[0] = cutoff
    #         else:
    #             value[0] = temp
    #
    # class Sing_gy(Expression):
    #     def eval(self, value, x):
    #         temp = 2*x[1]*pow(4*pow(x[0],2.0) + 4*pow(x[1],2.0), (-1.0/4.0));
    #         if(abs(x[0]) < xtol and abs(x[1]) < xtol):
    #             value[0] = cutoff
    #         else:
    #             value[0] = temp
    #
    # f = Sing_f();
    # gx = Sing_gx();
    # gy = Sing_gy();

    # exact = Expression('-pow(2.0 - pow(x[0],2.0) - pow(x[1],2.0), 1.0/2.0)');
    # class Sing_f1(Expression):
    #     def eval(self, value, x):
    #         temp = 2.0*pow(2.0 - pow(x[0],2.0) - pow(x[1],2.0), -2.0);
    #         if(abs(x[0] - 1) < xtol and abs(x[1] - 1) < xtol):
    #             value[0] = cutoff
    #         else:
    #             value[0] = temp
    # f = Sing_f1()
    # class Sing_gx1(Expression):
    #     def eval(self, value, x):
    #         temp = x[0]*pow(2.0 - pow(x[0],2.0) - pow(x[1],2.0), -1.0/2.0);
    #         if(abs(x[0] - 1) < xtol and abs(x[1] - 1) < xtol):
    #             value[0] = cutoff
    #         else:
    #             value[0] = temp
    # gx = Sing_gx1()
    # class Sing_gy1(Expression):
    #     def eval(self, value, x):
    #         temp = x[1]*pow(2.0 - pow(x[0],2.0) - pow(x[1],2.0), -1.0/2.0);
    #         if(abs(x[0] - 1) < xtol and abs(x[1] - 1) < xtol):
    #             value[0] = cutoff
    #         else:
    #             value[0] = temp
    # gy = Sing_gy1()

    # class Sing_u2(Expression):
    #     def eval(self, value, x):
    #         if(abs(x[0]) < xtol):
    #             value[0] = -x[0]
    #         else:
    #             value[0] = x[0]
    # exact = Sing_u2()
    # f = Expression('0.0')
    # class Sing_gx2(Expression):
    #     def eval(self, value, x):
    #         if(abs(x[0]) < xtol):
    #             value[0] = -1.0
    #         else:
    #             value[0] = 1.0
    # gx = Sing_gx2()
    # gy = Expression('0.0')

    class Sing_u2(Expression):
        def eval(self, value, x):
            if(abs(x[0]) < xtol):
                value[0] = sin(x[0])
            else:
                value[0] = x[0]**2
    exact = Sing_u2()
    f = Expression('0.0')
    class Sing_gx2(Expression):
        def eval(self, value, x):
            if(abs(x[0]) < xtol):
                value[0] = cos(x[0])
            else:
                value[0] = 2*x[0]
    gx = Sing_gx2()
    gy = Expression('0.0')

    # exact = Expression('sqrt(pow(x[0],2.0) + pow(x[1],2.0))')
    # class Sing_f3(Expression):
    #     def eval(self, value, x):
    #         if(abs(x[0]) < cutoff and abs(x[1]) < cutoff):
    #             value[0] = 0*4*cutoff;
    #         else:
    #             value[0] = 0.0;
    # f = Sing_f3()
    # class Sing_gx3(Expression):
    #     def eval(self, value, x):
    #         if(abs(x[0]) < xtol and abs(x[1]) < xtol):
    #             value[0] = cutoff;
    #         else:
    #             value[0] = x[0]/(sqrt(x[0]**2 + x[1]**2))
    #
    # class Sing_gy3(Expression):
    #     def eval(self, value, x):
    #         if(abs(x[0]) < xtol and abs(x[1]) < xtol):
    #             value[0] = cutoff;
    #         else:
    #             value[0] = x[1]/(sqrt(x[0]**2 + x[1]**2))
    #
    # gx = Sing_gx3()
    # gy = Sing_gy3()




    ##### Loop through epsilon values and solve ####################
    w = Function(MixedV);
    bcv = DirichletBC(MixedV.sub(3), exact, V_boundary)
    for epii in ep:
        print('Epsilon = ',epii)

        bcxx = DirichletBC(MixedV.sub(0), epii, Wxx_boundary)
        bcyy = DirichletBC(MixedV.sub(2), epii, Wyy_boundary)

        bc = [bcxx,bcyy,bcv]
    
        # Define variational problem
        (Sxx, Sxy, Syy, u) = TrialFunction(MixedV)
        (muxx, muxy, muyy, v) = TestFunction(MixedV)
    
        F = inner(Sxx,muxx)*dx + 2*inner(Sxy,muxy)*dx + inner(Syy,muyy)*dx;
        F += inner(Dx(u,0), Dx(muxx,0))*dx + inner(Dx(u,0), Dx(muxy,1))*dx;
        F += inner(Dx(u,1), Dx(muxy,0))*dx + inner(Dx(u,1), Dx(muyy,1))*dx;

        if(epii != 0):
            F += epii*( inner(Dx(Sxx,0), Dx(v,0)) + inner(Dx(Sxy,0), Dx(v,1)))*dx;
            F += epii*( inner(Dx(Sxy,1), Dx(v,0)) + inner(Dx(Syy,1), Dx(v,1)))*dx;

        # Determinant term/Nonlinear term
        F += inner(Sxx*Syy,v)*dx - inner(Sxy*Sxy,v)*dx;

        F -= (f*v*dx - gy*muxy*ds(1) + gx*muxy*ds(2) + gy*muxy*ds(3) - gx*muxy*ds(4));



        # Solve problem

        R = action(F,w);
        DR = derivative(R, w);
        problem = NonlinearVariationalProblem(R,w,bc,DR);
        solver = NonlinearVariationalSolver(problem);
        prm = solver.parameters
        # prm["nonlinear_solver"] = "snes"
        # prm['snes_solver']['report'] = False;
        # prm["newton_solver"]["absolute_tolerance"] = 1e-13;
        # prm["newton_solver"]["relative_tolerance"] = 1e-13;
        # prm["newton_solver"]["linear_solver"] = "gmres"
        # prm["newton_solver"]["krylov_solver"]["absolute_tolerance"] = 1E-9
        # prm["newton_solver"]["krylov_solver"]["relative_tolerance"] = 1E-7
        # prm["newton_solver"]["krylov_solver"]["maximum_iterations"] = 1000
        # prm["newton_solver"]["krylov_solver"]["monitor_convergence"] = False
        # prm["newton_solver"]["krylov_solver"]["nonzero_initial_guess"] = False
        # prm["newton_solver"]["krylov_solver"]["gmres"]["restart"] = 40
        # prm["newton_solver"]["preconditioner"] = "ilu" # default is "ilu"
        # prm["newton_solver"]["krylov_solver"]["preconditioner"]["structure"]\
        # = "same_nonzero_pattern"
        # prm["newton_solver"]["krylov_solver"]["preconditioner"]["ilu"]["fill_level"] =0
        solver.solve();

        (Sxx,Sxy,Syy,u) = w.split(deepcopy=True);

        error = abs(exact-u)**2*dx
        print('At epsilon = ', epii, ' error = ', np.sqrt(assemble(error)))
  
  
  
  



    (Sxx,Sxy,Syy,u) = w.split(deepcopy=True);

    error = abs(exact-u)**2*dx
    u0 = project(exact,V)
    grad_error = inner(nabla_grad(u0) - nabla_grad(u), nabla_grad(u0) - nabla_grad(u))*dx
    e[ii] = np.sqrt(assemble(error))
    
    if(ii > 0):
        ratio[ii] = np.log(e[ii-1]/e[ii])/np.log(2)
 
# plot(abs(exact-u))
# interactive()
print("Error: ", e)
print("Ratio: ", ratio)
