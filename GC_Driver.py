from dolfin import *
from GC_Problems import *
import numpy as np
from VM_Utilities import *
from VM_Solver import *

set_log_level(20)


#Values of N for the mesh
params = np.array([4, 8,16,32]);
params = np.array([32]);

L = len(params);
e = np.zeros([L,1]);
ratio = np.zeros([L,1]);

p = 2;

ep = np.array([1, 1e-1, 1e-2, 1e-3, 1e-4]);
Karr = np.array([0, 0.05, 0.15, 0.25, .32, .45, .55, .65, .75, .85, .95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6])
# ep = np.array([1]);
Karr = np.array([2, 2.05, 2.1, 2.515, 2.52, 2.54, 2.55, 2.56, 2.57, 2.58, 2.59])
ep = -ep;
# ep = np.array([1, 1e-1]);

prob = 1;
(x0, y0, x1, y1, exact, gx, gy, K) = GC_Problems(prob);
K = Karr[0]

for ii in range(L):
    N = params[ii];
    print(N);

    # Create mesh and define function space
    mesh = RectangleMesh(Point(x0,y0),Point(x1,y1),N,N)
    V = FunctionSpace(mesh, 'Lagrange', p)
    MixedV = MixedFunctionSpace([V,V,V,V]);


    SetupUtilities(mesh, x0, x1, y0, y1);



    ds = Create_dsMeasure()
    

    saved = raw_input('Would you like to load saved function? (y or n)')
    
    if(saved == 'n'):
        ##### Loop through epsilon values and solve ####################
        w = Function(MixedV);
        ep_err = [];

        for epjj in ep:
            print('Epsilon = ',epjj)

            w = ForwardProblem_GC(MixedV,K,ds, epjj, w, exact, gx, gy)


            #Compute error for this epsilon
            (Sxx,Sxy,Syy,u) = w.split(deepcopy=True);

            ep_err.append(np.sqrt(assemble(abs(exact-u)**2*dx)));
            print('Run finished at epsilon = ', epjj)
            print('L2 error = ', ep_err[-1])
            print ''

            # PlotToFile(u, 'Epsilon = ' + epjj.__str__(), 'file')
      
      
        # Now start changing the curvature
        for Kii in Karr[1:]:
            print 'Running curvature K = ', Kii
            w = ForwardProblem_GC(MixedV,Kii,ds, epjj, w, exact, gx, gy)
            print ''


        (Sxx,Sxy,Syy,u) = w.split(deepcopy=True);

        error = abs(exact-u)**2*dx
        e[ii] = np.sqrt(assemble(error))
        
        if(ii > 0):
            ratio[ii] = np.log(e[ii-1]/e[ii])/np.log(2)

        # Save solution to file
        file_Sxx = File('Sxx.xml');
        file_Sxy = File('Sxy.xml');
        file_Syy = File('Syy.xml');
        file_u = File('u.xml');  

        file_Sxx << Sxx; file_Sxy << Sxy; file_Syy << Syy; file_u << u;

    else:
        Sxx = Function(V,'Sxx.xml')
        Sxy = Function(V,'Sxy.xml')
        Syy = Function(V,'Syy.xml')
        u = Function(V,'u.xml')
        w = Function(MixedV);

        func_assign = FunctionAssigner([V,V,V,V],MixedV);
        func_assign.assign([Sxx,Sxy,Syy,u],w);

        Karr = np.array([.25,.35,.45,.55,.56])
        epjj = ep[-1];

        for Kii in Karr[0:]:
            print 'Running curvature K = ', Kii
            w = ForwardProblem_GC(MixedV,Kii,ds, epjj, w, exact, gx, gy)
            print ''

    # plot(exact-u)
    # interactive()
 
# plot(abs(exact-u))
# interactive()
print("Error: ", e)
print("Ratio: ", ratio)
