from dolfin import *
from MA_Problems import *
import numpy as np
import sys;

set_log_level(16)

# Values of N for the mesh
params = np.array([4, 8, 16, 32]);
params = np.array([16]);

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
    # 8. u(x,y) = x/x^2 piecewise function
    # 9. u(x,y) = sqrt(x^2 + y^2)
    # #       numerical Dirac delta function
    prob = 1;
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






    ##### Loop through epsilon values and solve ####################
    w = Function(MixedV);
    sols = []; solsSxx = []; solsSxy = []; solsSyy = []; solsw = [];
    bcv = DirichletBC(MixedV.sub(3), exact, V_boundary)
    for epii in ep:
        print('Epsilon = ', epii)

        bcxx = DirichletBC(MixedV.sub(0), 0.0, Wxx_boundary)
        bcyy = DirichletBC(MixedV.sub(2), 0.0, Wyy_boundary)

        bc = [bcxx,bcyy,bcv]

        # Define variational problem
        (Sxx, Sxy, Syy, u) = TrialFunction(MixedV)
        (muxx, muxy, muyy, v) = TestFunction(MixedV)

        F = inner(Sxx,muxx)*dx + 2*inner(Sxy,muxy)*dx + inner(Syy,muyy)*dx;

        F += inner(Dx(u,0), Dx(muxx,0))*dx + inner(Dx(u,0), Dx(muxy,1))*dx;
        F += inner(Dx(u,1), Dx(muxy,0))*dx + inner(Dx(u,1), Dx(muyy,1))*dx;

        if(epii != 0):
            F += epii*muxx*dx + epii*muyy*dx;

            F += epii*( inner(Dx(Sxx,0), Dx(v,0)) + inner(Dx(Sxy,0), Dx(v,1)))*dx;
            F += epii*( inner(Dx(Sxy,1), Dx(v,0)) + inner(Dx(Syy,1), Dx(v,1)))*dx;

            F += inner(epii*Sxx + epii*Syy,v)*dx + epii*epii*v*dx;

        # Determinant term/Nonlinear term
        F += inner(Sxx*Syy,v)*dx - inner(Sxy*Sxy,v)*dx;

        F -= (f*v*dx - gy*muxy*ds(1) + gx*muxy*ds(2) + gy*muxy*ds(3) - gx*muxy*ds(4));

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
        Sxxtemp += epii; Syytemp += epii;
        sols.append(utemp); solsSxx.append(Sxxtemp); solsSxy.append(Sxytemp); solsSyy.append(Syytemp);
        solsw.append(w);



    un = sols[-1]; unm1 = sols[-1 -1]; unm2 = sols[-1 -2];
    Sxxn = solsSxx[-1]; Sxxnm1 = solsSxx[-1 -1]; Sxxnm2 = solsSxx[-1 -2];
    Sxyn = solsSxy[-1]; Sxynm1 = solsSxy[-1 -1]; Sxynm2 = solsSxy[-1 -2];
    Syyn = solsSyy[-1]; Syynm1 = solsSyy[-1 -1]; Syynm2 = solsSyy[-1 -2];
    epn = ep[-1]; epnm1 = ep[-1 -1]; epnm2 = ep[-1 -2];
    wn = solsw[-1]; wnm1 = solsw[-1 - 1];

    dFdep_n = muxx*dx + muyy*dx + (inner(Dx(Sxxn, 0), Dx(v, 0)) + inner(Dx(Sxyn, 0), Dx(v, 1))) * dx + 0*muxy;
    dFdep_n += (inner(Dx(Sxyn, 1), Dx(v, 0)) + inner(Dx(Syyn, 1), Dx(v, 1))) * dx;
    dFdep_n += inner(Sxxn+Syyn+epn,v)*dx;

    dFdep_nm1 = muxx*dx + muyy*dx + (inner(Dx(Sxxnm1, 0), Dx(v, 0)) + inner(Dx(Sxynm1, 0), Dx(v, 1))) * dx + 0*muxy;
    dFdep_nm1 += (inner(Dx(Sxynm1, 1), Dx(v, 0)) + inner(Dx(Syynm1, 1), Dx(v, 1))) * dx;
    dFdep_nm1 += inner(Sxxnm1+Syynm1+epnm1,v)*dx;

    dwdep_n = Function(MixedV);
    dwdep_nm1 = Function(MixedV);
    R = action(F,dwdep_n)
    dFdu = derivative(R,dwdep_n)
    dFdu = replace(dFdu,{dwdep_n: wn})
    solve(dFdu == -dFdep_n, dwdep_n, bcs=bc);
    R = action(F,dwdep_nm1)
    dFdu = derivative(R,dwdep_nm1)
    dFdu = replace(dFdu,{dwdep_nm1: wnm1})
    solve(dFdu == -dFdep_nm1, dwdep_nm1, bcs=bc);

    dudep_n = dwdep_n.sub(3,deepcopy=True);
    dudep_nm1 = dwdep_nm1.sub(3,deepcopy=True);






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
