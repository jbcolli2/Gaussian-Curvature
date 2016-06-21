from dolfin import *
from GC_Problems import *
import numpy as np
from VM_Utilities import *
from VM_Solver import *
from scipy.interpolate import InterpolatedUnivariateSpline
import scipy.io

set_log_level(20)


#Values of N for the mesh
params = np.array([4, 8,16,32]);
params = np.array([8]);

L = len(params);
e = np.zeros([L,1]);
ratio = np.zeros([L,1]);

p = 2;

ep = np.array([1, 1e-1, 1e-2, 1e-3, 5e-4]);
Karr = np.array([.01])
# Karr = np.array([0]);
# ep = np.array([1]);
# ep = -ep;
# ep = np.array([1, 1e-1]);

prob = 5;
(x0, y0, x1, y1, exact, gx, gy, K) = GC_Problems(prob);



saved = raw_input('Would you like to load saved function? (y or n)')


for ii in range(L):
    N = params[ii];
    print(N);

    # Create mesh and define function space
    mesh = RectangleMesh(Point(x0,y0),Point(x1,y1),N,N)
    V = FunctionSpace(mesh, 'Lagrange', p)
    MixedV = MixedFunctionSpace([V,V,V,V]);


    SetupUtilities(mesh, x0, x1, y0, y1);



    ds = Create_dsMeasure()
    

    
    
    if(saved == 'n' or saved == ''):
        K = Karr[0]
        ##### Loop through epsilon values and solve ####################
        w = Function(MixedV);
        ep_err = [];

        for epjj in ep:
            print('Epsilon = ',epjj)

            bc = GetBC(MixedV, exact, epjj);

            F = F_Form_GC(MixedV, K, ds, epjj, gx, gy);
            w,temp = NewtonIteration(MixedV, w, F, bc, True);


            #Compute error for this epsilon
            (Sxx,Sxy,Syy,u) = w.split(deepcopy=True);

            ep_err.append(np.sqrt(assemble(abs(exact-u)**2*dx)));
            print('Run finished at epsilon = ', epjj)
            print('L2 error = ', ep_err[-1])
            print ''

            # PlotToFile(u, 'Epsilon = ' + epjj.__str__(), 'file')
      
        A = EvalJacobian(F,bc,w);
        s = np.linalg.svd(A.array(), compute_uv=False)
        s_val = np.min(s);
      
        # Now start changing the curvature
        sing = [];
        Kmax = []; Kout = []
        for jj,Kii in enumerate(Karr[1:]):
            print 'Running curvature K = ', Kii
            F = F_Form_GC(MixedV, Kii, ds, epjj, gx, gy);
            w, s = NewtonIteration(MixedV, w, F, bc);
            (Sxx,Sxy,Syy,u) = split(w);
            if(Kii > 0.35):
                Kout.append(Kii);
                if(s < 0):
                    sing.append(sing[-1]);
                else:
                    sing.append(s);

                if(len(Kout) > 2):
                    f = InterpolatedUnivariateSpline(sing,Kout[0:], k = 2)
                    Kmax.append(f(0));
                    print 'Kmax = ', f(0);

            print 'Norm of solution = ', assemble(u*u*dx);
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

        scipy.io.savemat('sval.mat', dict(x=sing, y = Karr))
    else:
        folder = 'g=0_N=8';
        Kfile = 'K=0_73';
        print 'Loading functions...'
        Sxx = Function(V); Sxx_file = File('Data Files/'+folder+'/Sxx'+'_'+folder+'_'+Kfile + '.xml'); Sxx_file >> Sxx;
        Sxy = Function(V); Sxy_file = File('Data Files/'+folder+'/Sxy'+'_'+folder+'_'+Kfile+ '.xml'); Sxy_file >> Sxy;
        Syy = Function(V); Syy_file = File('Data Files/'+folder+'/Syy'+'_'+folder+'_'+Kfile+ '.xml'); Syy_file >> Syy;
        u = Function(V); u_file = File('Data Files/'+folder+'/u'+'_'+folder+'_'+Kfile+ '.xml'); u_file >> u;
        w = Function(MixedV);

        func_assign = FunctionAssigner(MixedV,[V,V,V,V]);
        func_assign.assign(w,[Sxx,Sxy,Syy,u]);


        Karr = np.linspace(0.73,0.78, 200);
        epjj = ep[-1];
        bc = GetBC(MixedV, exact, epjj);
        sing = [];
        Kmax = []; Kout = []

        for jj,Kii in enumerate(Karr):
            print 'Running curvature K = ', Kii
            F = F_Form_GC(MixedV, Kii, ds, epjj, gx, gy);
            w, s = NewtonIteration(MixedV, w, F, bc);
            (Sxx,Sxy,Syy,u) = split(w);
            if(Kii > 0.35):
                Kout.append(Kii);
                if(s < 0):
                    sing.append(sing[-1]);
                else:
                    sing.append(-s);

                if(len(Kout) > 8):
                    f = InterpolatedUnivariateSpline(sing[-1 -9:],Kout[-1-9:], k = 3)
                    Kmax.append(f(0));
                    print 'Kmax = ', f(0);
                elif(len(Kout) > 3):
                    f = InterpolatedUnivariateSpline(sing,Kout, k = 3)
                    Kmax.append(f(0));
                    print 'Kmax = ', f(0);
                else:
                    Kmax.append(0);

            print 'Norm of solution = ', assemble(u*u*dx);
            print ''

        if(saved == 's'):
            # Save solution to file
            Sxx, Sxy, Syy, u = w.split(deepcopy=True);
            Kfileout = 'K=endgame'
            file_Sxx = File('Sxx'+'_'+folder+'_'+Kfileout+'.xml');
            file_Sxy = File('Sxy.xml');
            file_Syy = File('Syy.xml');
            file_u = File('u.xml');  

            file_Sxx << Sxx; file_Sxy << Sxy; file_Syy << Syy; file_u << u;

        scipy.io.savemat('sval.mat', dict(sing= sing, K = Kout, Kmax = Kmax))

        error = abs(exact-u)**2*dx
        e[ii] = np.sqrt(assemble(error))


    # plot(exact-u)
    # interactive()
 


# plot(abs(exact-u))
# interactive()
print("Error: ", e)
print("Ratio: ", ratio)
