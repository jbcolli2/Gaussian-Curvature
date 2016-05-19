from dolfin import *
from MA_Problems import *
import numpy as np
import cmath
import pdb
from IPython.core.debugger import Tracer;
from VM_Utilities import *
from VM_Solver import *


set_log_level(50)


#Values of N for the mesh
params = np.array([4, 8,16,32]);
# params = np.array([16]);

L = len(params);
e = np.zeros([L,1]);
errorIm = np.zeros([L,1]);
ratio = np.zeros([L,1]);

p = 2;

ep = np.array([1, 1e-2, 5e-3]);
loopsteps = 5


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
    prob = 7;
    (x0, y0, x1, y1, exact, f, gx, gy) = Problems(prob, N);



    # Create mesh and define function space
    mesh = RectangleMesh(Point(x0,y0),Point(x1,y1),N,N)
    V = FunctionSpace(mesh, 'Lagrange', p)
    MixedV = MixedFunctionSpace([V,V,V,V]);
    MixedVComplex = MixedFunctionSpace([V,V,V,V,V,V,V,V]);
    

    
    SetupUtilities(mesh, x0, x1, y0, y1);



    ds = Create_dsMeasure()
    
    
    
    




    solRe = [Function(V)]
    solIm = [Function(V)]


    ##### Loop through epsilon values and solve for real solution ####################
    w = Function(MixedV);
    sols = []; ep_err = [];
    for epjj in ep:
        print('Running epsilon = ',epjj)

        w = ForwardProblem(MixedV,ds, epjj, w, exact, f, gx, gy)


        (Sxx,Sxy,Syy,u) = w.split(deepcopy=True);

        ep_err.append( np.sqrt( assemble(abs(exact-u)**2*dx) ) );
        print('Run finished at epsilon = ', epjj)
        print('L2 error = ', ep_err[-1])



    solRe[0] = u;





    ################### Solve for other two complex solutions ######################

    fracs = np.linspace(0,2,loopsteps+1);
    epComplex = np.zeros(len(fracs)-1,dtype=complex)
    for ll,fracll in enumerate(fracs[0:len(fracs)-1]):
        epComplex[ll] = cmath.rect(ep[-1],fracll*cmath.pi);

    wRe = w
    wIm = interpolate(Expression(('0.0', '0.0', '0.0', '0.0')),MixedV)
    w = Function(MixedVComplex);
    assignerkk = FunctionAssigner(MixedVComplex, [V,V,V,V,V,V,V,V])
    assignerkk.assign(w,[interpolate(wRe.sub(0),V),interpolate(wRe.sub(1),V),interpolate(wRe.sub(2),V),interpolate(wRe.sub(3),V),\
                         interpolate(wIm.sub(0),V),interpolate(wIm.sub(1),V),interpolate(wIm.sub(2),V),interpolate(wIm.sub(3),V),])

    for kk, epkk in enumerate(epComplex[1:len(epComplex)]):
        solRe.append(Function(V))
        solIm.append(Function(V))
        print('Loop Epsilon = ',epkk)

        substeps = 2;
        eplinRe = np.linspace(epComplex[kk].real, epComplex[kk+1].real,substeps);
        eplinIm = np.linspace(epComplex[kk].imag, epComplex[kk+1].imag,substeps);
        eplin = np.zeros(substeps, dtype=complex);
        for ll in range(0,substeps):
            eplin[ll] = complex(eplinRe[ll], eplinIm[ll]);



        for jj, epjj in enumerate(eplin[1:len(eplin)]):
            print('Loop Line Epsilon = ',epjj)

            w = ForwardComplexProblem(MixedVComplex,ds, epjj, w, exact, f, gx, gy)


        solRe[kk+1] = w.sub(3,deepcopy=True)
        solIm[kk+1] = w.sub(7,deepcopy=True)













    # (Sxx,Sxy,Syy,u) = w.split(deepcopy=True);
#     Sxx = w.sub(0,deepcopy=true);
#     Sxy = w.sub(1);
#     Syy = w.sub(2);
#     u = w.sub(3);
    u = (solRe[0])
    for ll in range(1,len(epComplex)):
        u = u + solRe[ll]
    u = u/len(epComplex);

    error = abs(exact-u)**2*dx
    errorIm[ii] = np.sqrt(assemble(abs(u)**2*dx))
    u0 = project(exact,V)
    grad_error = inner(nabla_grad(u0) - nabla_grad(u), nabla_grad(u0) - nabla_grad(u))*dx
    e[ii] = np.sqrt(assemble(error))
    
    if(ii > 0):
        ratio[ii] = np.log(e[ii-1]/e[ii])/np.log(2)
 
# plot(exact)
# interactive()
print("Error: ", e)
print("Ratio: ", ratio)
