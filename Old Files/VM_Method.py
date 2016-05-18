from dolfin import *
import numpy

L = 5;
e = numpy.zeros([L,1]);
ratio = numpy.zeros([L,1])
for ii in range(L):
    N = pow(2,ii+1)
    epsilon = 1;
    
    # Create mesh and define function space
    mesh = UnitSquareMesh(N, N)
    V = FunctionSpace(mesh, 'Lagrange', 1)
    MixedV = MixedFunctionSpace([V,V]);
    

    
    # Boundaries for the W_xx, W_yy and V spaces
    def Wxx_boundary(x, on_boundary):
        return near(x[0],0.0) or near(x[0],1.0);
    def Wyy_boundary(x, on_boundary):
        return near(x[1],0.0) or near(x[1],1.0);
    def V_boundary(x, on_boundary):
        return near(x[1],0.0) or near(x[1],1.0) or near(x[0],0.0) or near(x[0],1.0);
#     bc = DirichletBC(MixedV.sub(2), 1.0, Wyy_boundary)
#     bc = DirichletBC(V, 0.0, Wyy_boundary)
    bcxx = DirichletBC(MixedV.sub(0), 1.0, Wxx_boundary)

  
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
    epsilon = 1.0;
    exact = Expression('.5*x[0]*(x[0]-1)')
    f = Expression('0.0')
    gx = Expression('x[0]-0.5');
    gy = Expression('x[1]-0.5');
    bcv = DirichletBC(MixedV.sub(1), exact, V_boundary)
    bc = [bcxx,bcv]
    
    # Define variational problem
    (Sxx, u) = TrialFunction(MixedV)
    (v, w) = TestFunction(MixedV)
    
    a = inner(Dx(Sxx,0),Dx(v,0))*dx + inner(Dx(u,0), Dx(w,0))*dx + inner(Sxx,w)*dx;
    L = f*v*dx + inner(gx,w)*ds(1) - inner(gx,w)*ds(3);
#     a = inner(Sxx,muxx)*dx + 2*inner(Sxy,muxy)*dx + inner(Syy,muyy)*dx;
#     a += inner(Dx(u,0), Dx(muxx,0))*dx + inner(Dx(u,0), Dx(muxy,1))*dx;
#     a += inner(Dx(u,1), Dx(muxy,0))*dx + inner(Dx(u,1), Dx(muyy,1))*dx;
#      
#     a += epsilon*( inner(Dx(Sxx,0), Dx(v,0)) + inner(Dx(Sxy,0), Dx(v,1)))*dx;
#     a += epsilon*( inner(Dx(Sxy,1), Dx(v,0)) + inner(Dx(Syy,1), Dx(v,1)))*dx;
#      
#      
#     L = f*v*dx - gy*muxy*ds(1) + gx*muxy*ds(2) + gy*muxy*ds(3) - gx*muxy*ds(4);


    
    
    # Solve problem
    sol = Function(MixedV);
    solve(a == L, sol, bc)
    (Sxx, u) = sol.split(deepcopy=True);
    error = abs(exact - u)**2*dx
    e[ii] = numpy.sqrt(assemble(error))
    
    if(ii > 0):
        ratio[ii] = numpy.log(e[ii-1]/e[ii])/numpy.log(2)
 
# plot(Syy)
# interactive()
print("Error: ", e)
print("Ratio: ", ratio)
