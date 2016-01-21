from dolfin import *
import numpy as np
import cmath
from IPython.core.debugger import Tracer;

set_log_level(16)


#Values of N for the mesh
params = np.array([4, 8,16,32,64]);
# params = np.array([1e-2, 1e-4, 1e-6, 1e-8]);
# params = np.array([8]);

L = len(params);
e = np.zeros([L,1]);
errorIm = np.zeros([L,1]);
ratio = np.zeros([L,1]);

p = 2;

ep = np.array([1, 1e-2, 1e-4, 1e-5, 1e-6]);



for ii in range(L):
    N = params[ii];
    print(N);


    # ep = np.logspace(np.log10(1),np.log10(params[ii]),5);
    
    # Create mesh and define function space
    mesh = UnitSquareMesh(N, N)
    V = FunctionSpace(mesh, 'Lagrange', p)
    MixedV = MixedFunctionSpace([V,V,V,V]);
    MixedVComplex = MixedFunctionSpace([V,V,V,V,V,V,V,V]);
    

    
    # Boundaries for the W_xx, W_yy and V spaces
    def Wxx_boundary(x, on_boundary):
        return near(x[0],0.0) or near(x[0],1.0);
    def Wyy_boundary(x, on_boundary):
        return near(x[1],0.0) or near(x[1],1.0);
    def V_boundary(x, on_boundary):
        return near(x[1],0.0) or near(x[1],1.0) or near(x[0],0.0) or near(x[0],1.0);
  
  
    # Boundaries data for integrating <g,\mu>
    class Left(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[0], 0.0)

    class Right(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[0], 1.0)
    
    class Bottom(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[1], 0.0)
    
    class Top(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[1], 1.0)
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

    exact = Expression('pow(x[0],4.0) + pow(x[1],2.0)',domain=mesh);
    f = Expression('24.0*pow(x[0],2.0)');
    gx = Expression('4.0*pow(x[0],3.0)');
    gy = Expression('2.0*x[1]');


    # exact = Expression('exp( 0.5*(pow(x[0],2.0) + pow(x[1],2.0)) )');
    # f = Expression('exp( (pow(x[0],2.0) + pow(x[1],2.0)) )* (pow(x[0],2.0) + pow(x[1],2.0)+1.0)');
    # gx = Expression('x[0]*exp( 0.5*(pow(x[0],2.0) + pow(x[1],2.0)) )');
    # gy = Expression('x[1]*exp( 0.5*(pow(x[0],2.0) + pow(x[1],2.0)) )');



    solRe = [Function(V),Function(V),Function(V)]
    solIm = [Function(V),Function(V),Function(V)]


    ##### Loop through epsilon values and solve for real solution ####################
    w = Function(MixedV);
    sols = [];
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

        F += epii*( inner(Dx(Sxx,0), Dx(v,0)) + inner(Dx(Sxy,0), Dx(v,1)))*dx;
        F += epii*( inner(Dx(Sxy,1), Dx(v,0)) + inner(Dx(Syy,1), Dx(v,1)))*dx;

        F += inner(Sxx*Syy,v)*dx - inner(Sxy*Sxy,v)*dx;

        F -= (f*v*dx - gy*muxy*ds(1) + gx*muxy*ds(2) + gy*muxy*ds(3) - gx*muxy*ds(4));



        # Solve problem

        R = action(F,w);
        DR = derivative(R, w);
        problem = NonlinearVariationalProblem(R,w,bc,DR);
        solver = NonlinearVariationalSolver(problem);
        solver.solve();


    solRe[0] = w.sub(3,deepcopy=True);





    ################### Solve for other two complex solutions ######################
    epComplex = np.array([ep[-1],cmath.rect(ep[-1],2*cmath.pi/3.0), cmath.rect(ep[-1],4*cmath.pi/3.0)]);
    wRe = w
    wIm = interpolate(Expression(('0.0', '0.0', '0.0', '0.0')),MixedV)
    w = Function(MixedVComplex);
    assignerkk = FunctionAssigner(MixedVComplex, [V,V,V,V,V,V,V,V])
    assignerkk.assign(w,[interpolate(wRe.sub(0),V),interpolate(wRe.sub(1),V),interpolate(wRe.sub(2),V),interpolate(wRe.sub(3),V),\
                         interpolate(wIm.sub(0),V),interpolate(wIm.sub(1),V),interpolate(wIm.sub(2),V),interpolate(wIm.sub(3),V),])

    bcvRe = DirichletBC(MixedVComplex.sub(3), exact, V_boundary)
    bcvIm = DirichletBC(MixedVComplex.sub(7), 0.0, V_boundary)

    for kk, epkk in enumerate(epComplex[1:len(epComplex)]):
        print('Epsilon = ',epkk)
        eplin = np.linspace(epComplex[kk],epComplex[kk+1],1)

        for jj, epjj in enumerate(eplin):

            bcxxRe = DirichletBC(MixedVComplex.sub(0), epjj.real, Wxx_boundary)
            bcyyRe = DirichletBC(MixedVComplex.sub(2), epjj.real, Wyy_boundary)
            bcxxIm = DirichletBC(MixedVComplex.sub(4), epjj.imag, Wxx_boundary)
            bcyyIm = DirichletBC(MixedVComplex.sub(6), epjj.imag, Wyy_boundary)

            bcRe = [bcxxRe,bcyyRe,bcvRe]
            bcIm = [bcxxIm, bcyyIm, bcvIm]
            bc = [bcxxRe,bcyyRe,bcvRe,bcxxIm,bcyyIm,bcvIm]

            # Define variational problem
            (SxxRe, SxyRe, SyyRe, uRe, SxxIm, SxyIm, SyyIm, uIm) = TrialFunction(MixedVComplex)
            (muxxRe, muxyRe, muyyRe, vRe, muxxIm, muxyIm, muyyIm, vIm) = TestFunction(MixedVComplex)

            # inner(Sxx,muxx)*dx
            F = (SxxRe + SxxIm)*muxxRe*dx + (SxxRe - SxxIm)*muxxIm*dx
            #2*inner(Sxy,muxy)*dx
            F += 2*(SxyRe + SxyIm)*muxyRe*dx + 2*(SxyRe - SxyIm)*muxyIm*dx
            # inner(Syy,muyy)*dx
            F += (SyyRe + SyyIm)*muyyRe*dx + (SyyRe - SyyIm)*muyyIm*dx

            #inner(Dx(u,0), Dx(muxx,0))*dx
            F += (Dx(uRe,0) + Dx(uIm,0))*Dx(muxxRe,0)*dx + (Dx(uRe,0) - Dx(uIm,0))*Dx(muxxIm,0)*dx
            # inner(Dx(u,0), Dx(muxy,1))*dx
            F += (Dx(uRe,0) + Dx(uIm,0))*Dx(muxyRe,1)*dx + (Dx(uRe,0) - Dx(uIm,0))*Dx(muxyIm,1)*dx
            #inner(Dx(u,1), Dx(muxy,0))*dx
            F += (Dx(uRe,1) + Dx(uIm,1))*Dx(muxyRe,0)*dx + (Dx(uRe,1) - Dx(uIm,1))*Dx(muxyIm,0)*dx
            #inner(Dx(u,1), Dx(muyy,1))*dx
            F += (Dx(uRe,1) + Dx(uIm,1))*Dx(muyyRe,1)*dx + (Dx(uRe,1) - Dx(uIm,1))*Dx(muyyIm,1)*dx


            # epsilon*( inner(Dx(Sxx,0), Dx(v,0))
            F += (epjj.real*(Dx(SxxRe,0) + Dx(SxxIm,0)) + epjj.imag*(Dx(SxxRe,0) - Dx(SxxIm,0))) * Dx(vRe,0)*dx
            F += (epjj.real*(Dx(SxxRe,0) - Dx(SxxIm,0)) - epjj.imag*(Dx(SxxRe,0) + Dx(SxxIm,0))) * Dx(vIm,0)*dx
            # inner(Dx(Sxy,0), Dx(v,1)))*dx
            F += (epjj.real*(Dx(SxyRe,0) + Dx(SxyIm,0)) + epjj.imag*(Dx(SxyRe,0) - Dx(SxyIm,0))) * Dx(vRe,1)*dx
            F += (epjj.real*(Dx(SxyRe,0) - Dx(SxyIm,0)) - epjj.imag*(Dx(SxyRe,0) + Dx(SxyIm,0))) * Dx(vIm,1)*dx
            # epsilon*( inner(Dx(Sxy,1), Dx(v,0))
            F += (epjj.real*(Dx(SxyRe,1) + Dx(SxyIm,1)) + epjj.imag*(Dx(SxyRe,1) - Dx(SxyIm,1))) * Dx(vRe,0)*dx
            F += (epjj.real*(Dx(SxyRe,1) - Dx(SxyIm,1)) - epjj.imag*(Dx(SxyRe,1) + Dx(SxyIm,1))) * Dx(vIm,0)*dx
            # inner(Dx(Syy,1), Dx(v,1)))*dx
            F += (epjj.real*(Dx(SyyRe,1) + Dx(SyyIm,1)) + epjj.imag*(Dx(SyyRe,1) - Dx(SyyIm,1))) * Dx(vRe,1)*dx
            F += (epjj.real*(Dx(SyyRe,1) - Dx(SyyIm,1)) - epjj.imag*(Dx(SyyRe,1) + Dx(SyyIm,1))) * Dx(vIm,1)*dx

            # inner(Sxx*Syy,v)*dx
            F += ( SxxRe*(SyyRe + SyyIm) + SxxIm*(SyyRe - SyyIm) )*vRe*dx
            F += ( SxxRe*(SyyRe - SyyIm) - SxxIm*(SyyRe + SyyIm) )*vIm*dx
            # - inner(Sxy*Sxy,v)*dx
            F -= ( SxyRe*(SxyRe + SxyIm) + SxyIm*(SxyRe - SxyIm) )*vRe*dx
            F += ( SxyRe*(SxyRe - SxyIm) - SxyIm*(SxyRe + SxyIm) )*vIm*dx

            # (f*v*dx - gy*muxy*ds(1) + gx*muxy*ds(2) + gy*muxy*ds(3) - gx*muxy*ds(4))
            F -= (f*vRe*dx - gy*muxyRe*ds(1) + gx*muxyRe*ds(2) + gy*muxyRe*ds(3) - gx*muxyRe*ds(4));
            F -= (f*vIm*dx - gy*muxyIm*ds(1) + gx*muxyIm*ds(2) + gy*muxyIm*ds(3) - gx*muxyIm*ds(4));



            # Solve problem

            R = action(F,w);
            DR = derivative(R, w);
            problem = NonlinearVariationalProblem(R,w,bc,DR);
            solver = NonlinearVariationalSolver(problem);
            solver.solve();

        solRe[kk+1] = w.sub(3,deepcopy=True)
        solIm[kk+1] = w.sub(7,deepcopy=True)













    # (Sxx,Sxy,Syy,u) = w.split(deepcopy=True);
#     Sxx = w.sub(0,deepcopy=true);
#     Sxy = w.sub(1);
#     Syy = w.sub(2);
#     u = w.sub(3);
    u = (solRe[0] + solRe[1] + solRe[2])
    error = abs(exact-u)**2*dx
    errorIm[ii] = np.sqrt(assemble(abs(u)**2*dx))
    u0 = project(exact,V)
    grad_error = inner(nabla_grad(u0) - nabla_grad(u), nabla_grad(u0) - nabla_grad(u))*dx
    e[ii] = np.sqrt(assemble(error)) - (4.0/3)
    
    if(ii > 0):
        ratio[ii] = np.log(e[ii-1]/e[ii])/np.log(2)
 
# plot(exact)
# interactive()
print("Error: ", e)
print("Ratio: ", ratio)
