from dolfin import *
from GC_Problems import *
import numpy as np
from VM_Utilities import *
from VM_Solver import *
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy import interpolate
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



# Various test problems I have predefined in file GC_Problems
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
    

    
    #If initial input is n or just hit enter?

    # Fix curvature K at a small value
    # Track epsilon from 1 to as small a value as possible by solving PDE at each value of epsilon
    #    and using solution from previous value of eps as initial iterate to Newton's method
    #    for next value of epsilon.

    # Values for epsilon are stored in ep.
    if(saved == 'n' or saved == ''):
        K = Karr[0]
        ##### Loop through epsilon values and solve ####################
        w = Function(MixedV);
        ep_err = [];

        for epjj in ep:
            print('Epsilon = ',epjj)

            bc = GetBC(MixedV, exact, epjj);

            F = F_Form_GC(MixedV, K, ds, epjj, gx, gy);
            w,temp = NewtonIteration(MixedV, w, F, bc, False);


            #Compute error for this epsilon
            (Sxx,Sxy,Syy,u) = w.split(deepcopy=True);

            ep_err.append(np.sqrt(assemble(abs(exact-u)**2*dx)));
            print('Run finished at epsilon = ', epjj)
            print('L2 error = ', ep_err[-1])
            print ''

            # PlotToFile(u, 'Epsilon = ' + epjj.__str__(), 'file')
        
        # Final solution for single value of K and smallest eps is stored in w.
      









        # Now start changing the curvature.
        #   ***** If only looking at a single value for curvature, this loop is not used!!!! *****
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
        if(saved == 'n'):
            file_Sxx = File('Sxx.xml');
            file_Sxy = File('Sxy.xml');
            file_Syy = File('Syy.xml');
            file_u = File('u.xml');  

            file_Sxx << Sxx; file_Sxy << Sxy; file_Syy << Syy; file_u << u;

        # scipy.io.savemat('sval.mat', dict(x=sing, y = Karr))



























    # Here is where we track the curvature parameter, using the same general method as with epsilon
    else:
        # Load solutions from file

        folder = 'g=0_N=8_2';
        Kfile = 'K=0_74';
        print 'Loading functions...'
        
        Sxx = Function(V); Sxx_file = File('Data Files/'+folder+'/Sxx'+'_'+folder+'_'+Kfile + '.xml'); Sxx_file >> Sxx;
        Sxy = Function(V); Sxy_file = File('Data Files/'+folder+'/Sxy'+'_'+folder+'_'+Kfile+ '.xml'); Sxy_file >> Sxy;
        Syy = Function(V); Syy_file = File('Data Files/'+folder+'/Syy'+'_'+folder+'_'+Kfile+ '.xml'); Syy_file >> Syy;
        u = Function(V); u_file = File('Data Files/'+folder+'/u'+'_'+folder+'_'+Kfile+ '.xml'); u_file >> u;
        w = Function(MixedV);

        func_assign = FunctionAssigner(MixedV,[V,V,V,V]);
        func_assign.assign(w,[Sxx,Sxy,Syy,u]);


        





        # Solve PDE for various values of the curvature

        # Values of K to track with
        Karr = np.linspace(0.74,0.78, 150);
        # Smallest value of epsilon
        epjj = ep[-1];


        # Initialize things
        bc = GetBC(MixedV, exact, epjj);
        sing = [];
        Kmax = []; Kout = []



        #Loop through values of K
        for jj,Kii in enumerate(Karr):
            # Solve PDE
            print 'Running curvature K = ', Kii
            F = F_Form_GC(MixedV, Kii, ds, epjj, gx, gy);
            # Note: w is solution of PDE, s is smallest singular value of Jacobian at final Newton iteration
            # Thought: 
            w, s = NewtonIteration(MixedV, w, F, bc);
            (Sxx,Sxy,Syy,u) = split(w);


            # This conditional is needed because s(K) is not a function until K > 0.35
            if(Kii > 0.35):
                Kout.append(Kii);
                if(s < 0):
                    sing.append(sing[-1]+1e-5);
                else:
                    sing.append(-s);

                if(len(Kout) > 4):
                    f = interpolate.PchipInterpolator(sing[-1-3:],Kout[-1-3:],True)
                    # z = np.polyfit(sing[-1-20:],Kout[-1-20:],3);
                    # f = np.poly1d(z);
                    Kmax.append(f(0));
                    print 'Kmax = ', f(0);
                # elif(len(Kout) > 3):
                #     f = InterpolatedUnivariateSpline(sing,Kout, k = 3)
                #     Kmax.append(f(0));
                #     print 'Kmax = ', f(0);
                else:
                    Kmax.append(0);

            print 'Norm of solution = ', assemble(u*u*dx);
            print ''

        if(saved == 's'):
            # Save solution to file
            Sxx, Sxy, Syy, u = w.split(deepcopy=True);
            Kfileout = 'K=0_74'
            file_Sxx = File('Sxx'+'_'+folder+'_'+Kfileout+'.xml');
            file_Sxy = File('Sxy'+'_'+folder+'_'+Kfileout+'.xml');
            file_Syy = File('Syy'+'_'+folder+'_'+Kfileout+'.xml');
            file_u = File('u'+'_'+folder+'_'+Kfileout+'.xml');  

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
