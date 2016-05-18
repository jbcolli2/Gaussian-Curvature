from dolfin import *
import numpy

L = 1;
e = numpy.zeros([L,1]);
ratio = numpy.zeros([L,1])
for ii in range(L):
    N = 2
    
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
#     exact = Expression('0.5*x[0]*(x[0]-1) + 0.5*x[1]*(x[1]-1)');
#     f = Expression('1.0')
#     gx = Expression('x[0]-0.5');
#     gy = Expression('x[1]-0.5');
    
#     epsilon = 0.25;
#     exact = Expression('sin(pi*x[0])*sin(pi*x[1]) + .5*ep*(pow(x[0],2) + pow(x[1],2))',ep=epsilon);
#     
#     f = Expression('-0.4e1 * ep* sin(pi * x[0]) * pow(pi, 0.4e1) * sin(pi * x[1]) + \
#      pow(-sin(pi * x[0]) * 0.3141592654e1 * pi * sin(pi * x[1]) + ep, 0.2e1) - \
#       pow(cos(pi * x[0]), 0.2e1) * pow(pi, 0.4e1) * pow(cos(pi * x[1]), 0.2e1)',ep=epsilon);
#       
#     gx = Expression('pi*cos(pi*x[0])*sin(pi*x[1]) + ep*x[0]', ep=epsilon);
#     gy = Expression('pi*sin(pi*x[0])*cos(pi*x[1]) + ep*x[1]',ep=epsilon);



    exact = Expression('pow(x[0],4.0) + pow(x[1],2.0)');
    
    f = Expression('1.0');
      
    gx = Expression('4.0*pow(x[0],3.0)');
    gy = Expression('2.0*x[1]');

    
  
  
  
    # Solve actual problem with initial condition w
    epsilon = .035
    bcxx = DirichletBC(V, epsilon, Wxx_boundary)
    bcyy = DirichletBC(V, epsilon, Wyy_boundary)
    bcv = DirichletBC(V, exact, V_boundary)
    bc = [bcxx,bcyy,bcv]
    
    # Define variational problem
#     (Sxx, Sxy, Syy, u) = TrialFunction(MixedV)
#     (muxx, muxy, muyy, v) = TestFunction(MixedV)
#     
#     F = inner(Sxx,muxx)*dx + 2*inner(Sxy,muxy)*dx + inner(Syy,muyy)*dx;
#     F += inner(Dx(u,0), Dx(muxx,0))*dx + inner(Dx(u,0), Dx(muxy,1))*dx;
#     F += inner(Dx(u,1), Dx(muxy,0))*dx + inner(Dx(u,1), Dx(muyy,1))*dx;
#     
#     F += epsilon*( inner(Dx(Sxx,0), Dx(v,0)) + inner(Dx(Sxy,0), Dx(v,1)))*dx;
#     F += epsilon*( inner(Dx(Sxy,1), Dx(v,0)) + inner(Dx(Syy,1), Dx(v,1)))*dx;
#     
#     F += inner(Sxx*Syy,v)*dx - inner(Sxy*Sxy,v)*dx;
#     
#     F -= (f*v*dx - gy*muxy*ds(1) + gx*muxy*ds(2) + gy*muxy*ds(3) - gx*muxy*ds(4));
    f = Expression('exp(x[0])');
    bcxx = DirichletBC(V, f, Wxx_boundary)
    bcyy = DirichletBC(V, epsilon, Wyy_boundary)
    bcv = DirichletBC(V, exact, V_boundary)

    Sxx = TrialFunction(V);
    v = TestFunction(V);
    axx = Sxx*v*dx;
    Axx = assemble(axx);
    
    L = f*v*dx;
    F = assemble(L)
    
    Axx = Axx.array();
    F = F.array();
    u = Function(V);
    solve(axx==L,u,bcxx)
    U = u.vector().array();
    
    c = mesh.coordinates();
    for ii in range(mesh.num_vertices()):
        print('u(%8g,%8g) = %g' % (c[ii][0], c[ii][1], u(c[ii][0],c[ii][1])))
    
    
    # Solve problem
# plot(Sxx)
# interactive()
# print("Error: ", e)
# print("Ratio: ", ratio)
