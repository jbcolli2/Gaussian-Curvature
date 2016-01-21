

from dolfin import *
import numpy

params = numpy.array([2,4,8,16,32]);
# params = numpy.array([4]);
L = len(params);
e = numpy.zeros([L,1]);
ratio = numpy.zeros([L,1])
for ii in range(L):
    N = params[ii];
    
    # Create mesh and define function space
    mesh = UnitSquareMesh(N, N)
    V = FunctionSpace(mesh, 'Lagrange', 2)
    MixedV = MixedFunctionSpace([V,V]);
    

    
    # Boundaries for the W_xx, W_yy and V spaces
    def V_boundary(x, on_boundary):
        return near(x[1],0.0) or near(x[1],1.0) or near(x[0],0.0) or near(x[0],1.0);
    def Wxx_boundary(x, on_boundary):
        return near(x[0],0.0) or near(x[0],1.0);

    
    
    
    
    
    
    
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



    exactRe = Expression('sin(pi*x[0])*sin(pi*x[1])')
    exactIm = Expression('x[0]*(x[0]-1.0)*x[1]*(x[1]-1.0)')
    fRe = Expression('2.0*pow(pi,2.0)*sin(pi*x[0])*sin(pi*x[1])');
    fIm = Expression('-2*x[1]*(x[1]-1.0) - 2*x[0]*(x[0]-1.0)');
    bc = DirichletBC(MixedV, Expression(('0.0','0.0')), V_boundary)
    # g = interpolate(exact,V);
    
    
    uRe,uIm = TrialFunction(MixedV);
    vRe, vIm = TestFunction(MixedV);
    
    # a = inner(grad(u),grad(v))*dx
    a = (Dx(uRe,0)+Dx(uIm,0))*Dx(vRe,0)*dx + (Dx(uRe,0)-Dx(uIm,0))*Dx(vIm,0)*dx
    a += (Dx(uRe,1)+Dx(uIm,1))*Dx(vRe,1)*dx + (Dx(uRe,1)-Dx(uIm,1))*Dx(vIm,1)*dx
    L = (fRe+fIm)*vRe*dx + (fRe-fIm)*vIm*dx
  
    u = Function(MixedV)
    solve(a == L, u, bc)
    uRe,uIm = u.split(deepcopy=True);
    error = abs(exactRe-uRe)**2*dx + abs(exactIm-uIm)**2*dx
    e[ii] = numpy.sqrt(assemble(error))
    
    if(ii > 0):
        ratio[ii] = numpy.log(e[ii-1]/e[ii])/numpy.log(2)
    
#     K = assemble(a);
#     F = assemble(L);
#     bc.apply(K,F)
#     K = K.array();
#     F = F.array();
#
#     cInds = bc.get_boundary_values().keys();
#     fInds = range(K.shape[1]);
#     for x in cInds:
#         fInds.remove(x);
#
#     Kf = numpy.zeros([len(fInds),len(fInds)]);
#     Ff = numpy.zeros(len(fInds));
#     for ii in range(len(fInds)):
#         Ff[ii] = F[fInds[ii]];
#         for jj in range(len(fInds)):
#             Kf[ii,jj] = K[fInds[ii],fInds[jj]];
# #     K = K.array();
# #     F = F.array();
#
#     onefunc = interpolate(Expression('1.0'),V)
#     R = action(a,onefunc)
#     R = assemble(R);
#     bc.apply(R)
    # Solve problem
# plot(Sxx)
# interactive()
print("Error: ", e)
print("Ratio: ", ratio)
