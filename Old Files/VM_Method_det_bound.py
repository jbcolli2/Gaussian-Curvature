from dolfin import *
import numpy

params = numpy.array([2,4,8,16,32]);
params = numpy.array([2]);
L = len(params);
e = numpy.zeros([L,1]);
ratio = numpy.zeros([L,1])
for ii in range(L):
    N = params[ii]
    epsilon = 1;
    
    # Create mesh and define function space
    mesh = UnitSquareMesh(N, N)
    V = FunctionSpace(mesh, 'Lagrange', 1)
    MixedV = MixedFunctionSpace([V,V,V,V]);
    

    
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
    
    
    
    
    # Problem data
#     epsilon = 1.0;
#     exact = Expression('0.5*x[0]*(x[0]-1) + 0.5*x[1]*(x[1]-1)',domain=mesh);
#     f = Expression('1.0')
#     g = exact;
#     gx = Expression('x[0]-0.5');
#     gy = Expression('x[1]-0.5');
    
    epsilon = 1.0;
    exact = Expression('sin(pi*x[0])*sin(pi*x[1]) + .5*ep*(pow(x[0],2) + pow(x[1],2))',ep=epsilon, domain=mesh );
      
    f = Expression('-0.4e1 * ep* sin(pi * x[0]) * pow(pi, 0.4e1) * sin(pi * x[1]) + \
     pow(-sin(pi * x[0]) * 0.3141592654e1 * pi * sin(pi * x[1]) + ep, 0.2e1) - \
      pow(cos(pi * x[0]), 0.2e1) * pow(pi, 0.4e1) * pow(cos(pi * x[1]), 0.2e1)',ep=epsilon);
#     f = Expression('-4.0*pow(pi,4.0)*sin(pi*x[0])*sin(pi*x[1])',ep=epsilon);
   
    g = exact;
    gx = Expression('pi*cos(pi*x[0])*sin(pi*x[1]) + ep*x[0]', ep=epsilon);
    gy = Expression('pi*sin(pi*x[0])*cos(pi*x[1]) + ep*x[1]',ep=epsilon);


#     epsilon = 0.01;
#     exact = Expression('exp(pow(x[0],2.0) + pow(x[1],2.0))',ep=epsilon);
#     
#     f = Expression('4.0*exp(2.0*(pow(x[0],2.0) + pow(x[1],2.0)))* \
#         (2.0*(pow(x[0],2.0) + pow(x[1],2.0)) + 1.0)',ep=epsilon);
#       
#     gx = Expression('2.0*x[0]*exp(pow(x[0],2.0) + pow(x[1],2.0))', ep=epsilon);
#     gy = Expression('2.0*x[1]*exp(pow(x[0],2.0) + pow(x[1],2.0))',ep=epsilon);

    
    bcxx = DirichletBC(MixedV.sub(0), 0.0, Wxx_boundary)
    bcyy = DirichletBC(MixedV.sub(2), 0.0, Wyy_boundary)
    bcv = DirichletBC(MixedV.sub(3), 0.0, V_boundary)
    bc = [bcxx,bcyy,bcv]
    
    # Define variational problem
    (Sxx, Sxy, Syy, u) = TrialFunction(MixedV)
    (muxx, muxy, muyy, v) = TestFunction(MixedV)
    
    F = inner(Sxx,muxx)*dx + 2*inner(Sxy,muxy)*dx + inner(Syy,muyy)*dx;
    F += inner(Dx(u,0), Dx(muxx,0))*dx + inner(Dx(u,0), Dx(muxy,1))*dx;
    F += inner(Dx(u,1), Dx(muxy,0))*dx + inner(Dx(u,1), Dx(muyy,1))*dx;
    
    F += epsilon*( inner(Dx(Sxx,0), Dx(v,0)) + inner(Dx(Sxy,0), Dx(v,1)))*dx;
    F += epsilon*( inner(Dx(Sxy,1), Dx(v,0)) + inner(Dx(Syy,1), Dx(v,1)))*dx;
    
    F += inner((Sxx+epsilon)*(Syy+epsilon),v)*dx - inner(Sxy*Sxy,v)*dx;
    
    F -= (f*v*dx - gy*muxy*ds(1) + gx*muxy*ds(2) + gy*muxy*ds(3) - gx*muxy*ds(4));

    F -= -epsilon*muxx*dx - epsilon*muyy*dx;
    F -= -Dx(g,0)*Dx(muxx,0)*dx;
    F -= -Dx(g,0)*Dx(muxy,1)*dx - Dx(g,1)*Dx(muxy,0)*dx;
    F -= -Dx(g,1)*Dx(muyy,1)*dx;
    
    DF = (Sxx+epsilon)*(Syy+epsilon)*v*dx - Sxy*Sxy*v*dx;
    
    
    # Solve problem
    w = Function(MixedV);
#     solve(a==F,w,bc);
    R = action(F,w);
    DR = derivative(R, w);
    problem = NonlinearVariationalProblem(R,w,bc,DR);
    solver = NonlinearVariationalSolver(problem);
    solver.solve();
    
    (Sxx,Sxy,Syy,u) = w.split(deepcopy=True);
    Sxx += epsilon;
    Syy += epsilon;
    u += g;
#     Sxx = w.sub(0,deepcopy=true);
#     Sxy = w.sub(1);
#     Syy = w.sub(2);
#     u = w.sub(3);
    error = abs(exact-u)**2*dx
    e[ii] = numpy.sqrt(assemble(error))
    
    if(ii > 0):
        ratio[ii] = numpy.log(e[ii-1]/e[ii])/numpy.log(2)
 
# plot(Sxx)
# interactive()
print("Error: ", e)
print("Ratio: ", ratio)
w = interpolate(Expression(('0.0','0.0','0.0','0.0')),MixedV)
R = action(DF,w);
Raa = assemble(R).array();
