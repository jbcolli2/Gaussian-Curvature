from dolfin import *
import numpy as np

set_log_level(16)

# Values of N for the mesh
params = np.array([4, 8, 16, 32]);
# params = np.array([1e-2, 1e-4, 1e-6, 1e-8]);
# params = np.array([8]);

L = len(params);
e = np.zeros([L, 1]);
ratio = np.zeros([L, 1]);

p = 2;

ep = np.array([1, 1e-1, 1e-2, 1e-3]);

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
    prob = 5;
    (x0, y0, x1, y1, exact, f, gx, gy) = Problems(prob, N);



    # Create mesh and define function space
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

    # exact = Expression('-pow(2.0 - pow(x[0],2.0) - pow(x[1],2.0), 1.0/2.0)');
    cutoff = pow(N,2.0);
    xtol = 1e-7;
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

    # exact = Expression('sqrt(pow(x[0],2.0) + pow(x[1],2.0))')
    # class Sing_f3(Expression):
    #     def eval(self, value, x):
    #         if(abs(x[0]) < cutoff and abs(x[1]) < cutoff):
    #             value[0] = 4*cutoff;
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




    ##### Loop through epsilon values and solve ####################
    w = Function(MixedV);
    sols = [];
    bcv = DirichletBC(MixedV.sub(3), exact, V_boundary)
    for epii in ep:
        print('Epsilon = ', epii)

        bcxx = DirichletBC(MixedV.sub(0), epii, Wxx_boundary)
        bcyy = DirichletBC(MixedV.sub(2), epii, Wyy_boundary)

        bc = [bcxx, bcyy, bcv]

        # Define variational problem
        (Sxx, Sxy, Syy, u) = TrialFunction(MixedV)
        (muxx, muxy, muyy, v) = TestFunction(MixedV)

        F = inner(Sxx, muxx) * dx + 2 * inner(Sxy, muxy) * dx + inner(Syy, muyy) * dx;
        F += inner(Dx(u, 0), Dx(muxx, 0)) * dx + inner(Dx(u, 0), Dx(muxy, 1)) * dx;
        F += inner(Dx(u, 1), Dx(muxy, 0)) * dx + inner(Dx(u, 1), Dx(muyy, 1)) * dx;

        if(epii != 0):
            F += epii * (inner(Dx(Sxx, 0), Dx(v, 0)) + inner(Dx(Sxy, 0), Dx(v, 1))) * dx;
            F += epii * (inner(Dx(Sxy, 1), Dx(v, 0)) + inner(Dx(Syy, 1), Dx(v, 1))) * dx;

        F += inner(Sxx * Syy, v) * dx - inner(Sxy * Sxy, v) * dx;

        F -= (f * v * dx - gy * muxy * ds(1) + gx * muxy * ds(2) + gy * muxy * ds(3) - gx * muxy * ds(4));

        # Solve problem

        R = action(F, w);
        DR = derivative(R, w);
        problem = NonlinearVariationalProblem(R, w, bc, DR);
        solver = NonlinearVariationalSolver(problem);
        solver.solve();

        (Sxx,Sxy,Syy,u) = w.split(deepcopy=True);

        error = abs(exact-u)**2*dx
        print('At epsilon = ', epii, ' error = ', np.sqrt(assemble(error)))

        Sxxtemp, Sxytemp, Syytemp, utemp = w.split(deepcopy=True);
        sols.append(utemp);

    # Quadratic Extrapolation in epsilon

    u = (-ep[-1 - 1]) * (-ep[-1]) / ((ep[-1 - 2] - ep[-1 - 1]) * (ep[-1 - 2] - ep[-1])) * sols[-1 - 2] + (-ep[
        -1 - 2]) * (-ep[-1]) / ((ep[-1 - 1] - ep[-1 - 2]) * (ep[-1 - 1] - ep[-1])) * sols[-1 - 1] + \
        (-ep[-1 - 1]) * (-ep[-1 - 2]) / ((ep[-1] - ep[-1 - 1]) * (ep[-1] - ep[-1 - 2])) * sols[-1]

    # Cubic Extrapolation in epsilon

    # u = (-ep[-1-2])*(-ep[-1-1])*(-ep[-1])/((ep[-1-3]-ep[-1-1])*(ep[-1-3]-ep[-1])*(ep[-1-3]-ep[-1-2]))*sols[-1-3] + \
    #     (-ep[-1-3])*(-ep[-1-1])*(-ep[-1])/((ep[-1-2]-ep[-1-1])*(ep[-1-2]-ep[-1])*(ep[-1-2]-ep[-1-3]))*sols[-1-2] + \
    #     (-ep[-1-2])*(-ep[-1-3])*(-ep[-1])/((ep[-1-1]-ep[-1-3])*(ep[-1-1]-ep[-1])*(ep[-1-1]-ep[-1-2]))*sols[-1-1] + \
    #     (-ep[-1-2])*(-ep[-1-1])*(-ep[-1-3])/((ep[-1-0]-ep[-1-1])*(ep[-1-0]-ep[-1-3])*(ep[-1-0]-ep[-1-2]))*sols[-1-0]




    # (Sxx,Sxy,Syy,u) = w.split(deepcopy=True);
    #     Sxx = w.sub(0,deepcopy=true);
    #     Sxy = w.sub(1);
    #     Syy = w.sub(2);
    #     u = w.sub(3);
    error = abs(exact - u) ** 2 * dx
    u0 = project(exact, V)
    grad_error = inner(nabla_grad(u0) - nabla_grad(u), nabla_grad(u0) - nabla_grad(u)) * dx
    e[ii] = np.sqrt(assemble(error))

    if (ii > 0):
        ratio[ii] = np.log(e[ii - 1] / e[ii]) / np.log(2)

# plot(abs(exact-u))
# interactive()
print("Error: ", e)
print("Ratio: ", ratio)
