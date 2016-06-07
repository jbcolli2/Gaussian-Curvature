from dolfin import *
import numpy as np




x0 = 0; x1 = 1; y0 = 0; y1 = 1; mesh = ();
def SetupUtilities(pmesh, px0, px1, py0, py1):
    global x0; global x1; global y0; global y1; global mesh;
    mesh = pmesh;
    x0 = px0;
    x1 = px1;
    y0 = py0;
    y1 = py1;


############## Boundary functions ##########
def EW_boundary(x, on_boundary):
    return near(x[0],x0) or near(x[0],x1);
def NS_boundary(x, on_boundary):
    return near(x[1],y0) or near(x[1],y1);
def Dir_boundary(x, on_boundary):
    return near(x[1],y0) or near(x[1],y1) or near(x[0],x0) or near(x[0],x1);


# Create measures for each side of cube
def Create_dsMeasure():
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

    return ds;









def PlotToFile(sol, title, filename):
    viz = plot(sol, title, axes = True)
    viz.write_png(filename)







#################################
  #
  # Solver
  #
#####################
def SetParameters(prm):
    prm["nonlinear_solver"] = "snes"
    prm['snes_solver']['report'] = False;
    prm["newton_solver"]["absolute_tolerance"] = 1e-13;
    prm["newton_solver"]["relative_tolerance"] = 1e-13;
    prm["newton_solver"]["linear_solver"] = "gmres"
    prm["newton_solver"]["krylov_solver"]["absolute_tolerance"] = 1E-9
    prm["newton_solver"]["krylov_solver"]["relative_tolerance"] = 1E-7
    prm["newton_solver"]["krylov_solver"]["maximum_iterations"] = 1000
    prm["newton_solver"]["krylov_solver"]["monitor_convergence"] = False
    prm["newton_solver"]["krylov_solver"]["nonzero_initial_guess"] = False
    prm["newton_solver"]["krylov_solver"]["gmres"]["restart"] = 40
    prm["newton_solver"]["preconditioner"] = "ilu" # default is "ilu"
    prm["newton_solver"]["krylov_solver"]["preconditioner"]["structure"]\
    = "same_nonzero_pattern"
    prm["newton_solver"]["krylov_solver"]["preconditioner"]["ilu"]["fill_level"] =0








def NewtonIteration(V, w0, form, bc):
    max_iters = 30;
    abstol = 1e-9;
    reltol = 1e-6;

    it_count = 0;

    wk = Function(V);
    wk.assign(w0);
    wkp1 = Function(V);
    wkp1.vector()[:] += 100.0;
    abserr = np.linalg.norm( (wk.vector().array() - wkp1.vector().array()), np.Inf );
    relerr = abserr/np.linalg.norm(wkp1.vector().array());
    while((abserr > abstol or relerr > reltol) and it_count < max_iters):
        it_count += 1;
        dw = TrialFunction(V);

        Fwk = action(form, wk);
        Feval_k = np.linalg.norm( assemble(Fwk).array(), ord=np.Inf )
        DF = derivative(Fwk, wk, dw);   # derivative of F evaluated at wk.  dw is trial function to be solved for

        A,b = assemble_system(DF, (-1.0)*Fwk, bc);

        dw = Function(V);
        solve(A, dw.vector(), b);

        lam = 1.0;
        print 'Old Function value = ', Feval_k

        for ii in range(1):
            wkp1.assign(wk);
            wkp1.vector().axpy(lam, dw.vector());
            Fwkp1 = action(form, wkp1);
            Feval_kp1 = np.linalg.norm( assemble(Fwkp1).array(), ord=np.Inf )
            print 'Function value = ', Feval_kp1

            lam = lam/2;


        abserr = np.linalg.norm( (wk.vector().array() - wkp1.vector().array()), np.Inf );
        relerr = abserr/np.linalg.norm(wkp1.vector().array(), np.Inf);

        

        print 'Iteration Number: ', it_count
        print 'Absolute Error = ', abserr
        print 'Relative Error = ', relerr
        print ''

        wk.assign(wkp1);


    return wk;




#################################
  #
  # Forms
  #
#####################

def F_Form_MA(MixedV, ds, ep, f, gx, gy):
    (Sxx, Sxy, Syy, u) = TrialFunction(MixedV)
    (muxx, muxy, muyy, v) = TestFunction(MixedV)

    F = inner(Sxx,muxx)*dx + 2*inner(Sxy,muxy)*dx + inner(Syy,muyy)*dx;
    F += inner(Dx(u,0), Dx(muxx,0))*dx + inner(Dx(u,0), Dx(muxy,1))*dx;
    F += inner(Dx(u,1), Dx(muxy,0))*dx + inner(Dx(u,1), Dx(muyy,1))*dx;

    if(ep != 0):
        F += ep*( inner(Dx(Sxx,0), Dx(v,0)) + inner(Dx(Sxy,0), Dx(v,1)))*dx;
        F += ep*( inner(Dx(Sxy,1), Dx(v,0)) + inner(Dx(Syy,1), Dx(v,1)))*dx;

    # Determinant term/Nonlinear term
    F += inner(Sxx*Syy,v)*dx - inner(Sxy*Sxy,v)*dx;

    F -= (f*v*dx - gy*muxy*ds(1) + gx*muxy*ds(2) + gy*muxy*ds(3) - gx*muxy*ds(4));

    return F;



def F_Form_GC(MixedV, K, ds, ep, gx, gy):
    (Sxx, Sxy, Syy, u) = TrialFunction(MixedV)
    (muxx, muxy, muyy, v) = TestFunction(MixedV)

    F = inner(Sxx,muxx)*dx + 2*inner(Sxy,muxy)*dx + inner(Syy,muyy)*dx;
    F += inner(Dx(u,0), Dx(muxx,0))*dx + inner(Dx(u,0), Dx(muxy,1))*dx;
    F += inner(Dx(u,1), Dx(muxy,0))*dx + inner(Dx(u,1), Dx(muyy,1))*dx;

    if(ep != 0):
        F += ep*( inner(Dx(Sxx,0), Dx(v,0)) + inner(Dx(Sxy,0), Dx(v,1)))*dx;
        F += ep*( inner(Dx(Sxy,1), Dx(v,0)) + inner(Dx(Syy,1), Dx(v,1)))*dx;

    # Determinant term/Nonlinear term
    F += (((Sxx*Syy - Sxy*Sxy)*(1 + (Dx(u,0)**2 + Dx(u,1)**2))**(-2)) - K)*v*dx;


    F -= (-gy*muxy*ds(1) + gx*muxy*ds(2) + gy*muxy*ds(3) - gx*muxy*ds(4));

    return F;






def F_FormWithBoundary_MA(MixedV, ds, ep, f, gx, gy):
    (Sxx, Sxy, Syy, u) = TrialFunction(MixedV)
    (muxx, muxy, muyy, v) = TestFunction(MixedV)

    F = inner(Sxx,muxx)*dx + 2*inner(Sxy,muxy)*dx + inner(Syy,muyy)*dx;
    F += inner(Dx(u,0), Dx(muxx,0))*dx + inner(Dx(u,0), Dx(muxy,1))*dx;
    F += inner(Dx(u,1), Dx(muxy,0))*dx + inner(Dx(u,1), Dx(muyy,1))*dx;

    if(ep != 0):
        F += ep*( inner(Dx(Sxx,0), Dx(v,0)) + inner(Dx(Sxy,0), Dx(v,1)))*dx;
        F += ep*( inner(Dx(Sxy,1), Dx(v,0)) + inner(Dx(Syy,1), Dx(v,1)))*dx;
        F += ep*(muxx + muyy)*dx;
        F += ep*(Sxx + Syy)*v*dx + ep*ep*v*dx;

    # Determinant term/Nonlinear term
    F += inner(Sxx*Syy,v)*dx - inner(Sxy*Sxy,v)*dx;

    F -= (f*v*dx - gy*muxy*ds(1) + gx*muxy*ds(2) + gy*muxy*ds(3) - gx*muxy*ds(4));

    return F;







def F_ComplexForm_MA(MixedVComplex, ds, ep, f, gx, gy):
    # Define variational problem
    (SxxRe, SxyRe, SyyRe, uRe, SxxIm, SxyIm, SyyIm, uIm) = TrialFunction(MixedVComplex)
    (muxxRe, muxyRe, muyyRe, vRe, muxxIm, muxyIm, muyyIm, vIm) = TestFunction(MixedVComplex)

    epRe = ep.real;
    epIm = ep.imag;
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


    if(epIm != 0 and epRe != 0):
        # epsilon*( inner(Dx(Sxx,0), Dx(v,0))
        F += (epRe*(Dx(SxxRe,0) + Dx(SxxIm,0)) + epIm*(Dx(SxxRe,0) - Dx(SxxIm,0))) * Dx(vRe,0)*dx
        F += (epRe*(Dx(SxxRe,0) - Dx(SxxIm,0)) - epIm*(Dx(SxxRe,0) + Dx(SxxIm,0))) * Dx(vIm,0)*dx
        # inner(Dx(Sxy,0), Dx(v,1)))*dx
        F += (epRe*(Dx(SxyRe,0) + Dx(SxyIm,0)) + epIm*(Dx(SxyRe,0) - Dx(SxyIm,0))) * Dx(vRe,1)*dx
        F += (epRe*(Dx(SxyRe,0) - Dx(SxyIm,0)) - epIm*(Dx(SxyRe,0) + Dx(SxyIm,0))) * Dx(vIm,1)*dx
        # epsilon*( inner(Dx(Sxy,1), Dx(v,0))
        F += (epRe*(Dx(SxyRe,1) + Dx(SxyIm,1)) + epIm*(Dx(SxyRe,1) - Dx(SxyIm,1))) * Dx(vRe,0)*dx
        F += (epRe*(Dx(SxyRe,1) - Dx(SxyIm,1)) - epIm*(Dx(SxyRe,1) + Dx(SxyIm,1))) * Dx(vIm,0)*dx
        # inner(Dx(Syy,1), Dx(v,1)))*dx
        F += (epRe*(Dx(SyyRe,1) + Dx(SyyIm,1)) + epIm*(Dx(SyyRe,1) - Dx(SyyIm,1))) * Dx(vRe,1)*dx
        F += (epRe*(Dx(SyyRe,1) - Dx(SyyIm,1)) - epIm*(Dx(SyyRe,1) + Dx(SyyIm,1))) * Dx(vIm,1)*dx
    elif(epRe != 0):
        # epsilon*( inner(Dx(Sxx,0), Dx(v,0))
        F += (epRe*(Dx(SxxRe,0) + Dx(SxxIm,0))) * Dx(vRe,0)*dx
        F += (epRe*(Dx(SxxRe,0) - Dx(SxxIm,0))) * Dx(vIm,0)*dx
        # inner(Dx(Sxy,0), Dx(v,1)))*dx
        F += (epRe*(Dx(SxyRe,0) + Dx(SxyIm,0))) * Dx(vRe,1)*dx
        F += (epRe*(Dx(SxyRe,0) - Dx(SxyIm,0))) * Dx(vIm,1)*dx
        # epsilon*( inner(Dx(Sxy,1), Dx(v,0))
        F += (epRe*(Dx(SxyRe,1) + Dx(SxyIm,1))) * Dx(vRe,0)*dx
        F += (epRe*(Dx(SxyRe,1) - Dx(SxyIm,1))) * Dx(vIm,0)*dx
        # inner(Dx(Syy,1), Dx(v,1)))*dx
        F += (epRe*(Dx(SyyRe,1) + Dx(SyyIm,1))) * Dx(vRe,1)*dx
        F += (epRe*(Dx(SyyRe,1) - Dx(SyyIm,1))) * Dx(vIm,1)*dx
    elif(epIm != 0):
        # epsilon*( inner(Dx(Sxx,0), Dx(v,0))
        F += (epIm*(Dx(SxxRe,0) - Dx(SxxIm,0))) * Dx(vRe,0)*dx
        F += (epIm*(Dx(SxxRe,0) + Dx(SxxIm,0))) * Dx(vIm,0)*dx
        # inner(Dx(Sxy,0), Dx(v,1)))*dx
        F += (epIm*(Dx(SxyRe,0) - Dx(SxyIm,0))) * Dx(vRe,1)*dx
        F += (epIm*(Dx(SxyRe,0) + Dx(SxyIm,0))) * Dx(vIm,1)*dx
        # epsilon*( inner(Dx(Sxy,1), Dx(v,0))
        F += (epIm*(Dx(SxyRe,1) - Dx(SxyIm,1))) * Dx(vRe,0)*dx
        F += (epIm*(Dx(SxyRe,1) + Dx(SxyIm,1))) * Dx(vIm,0)*dx
        # inner(Dx(Syy,1), Dx(v,1)))*dx
        F += (epIm*(Dx(SyyRe,1) - Dx(SyyIm,1))) * Dx(vRe,1)*dx
        F += (epIm*(Dx(SyyRe,1) + Dx(SyyIm,1))) * Dx(vIm,1)*dx


    # inner(Sxx*Syy,v)*dx
    F += ( SxxRe*(SyyRe + SyyIm) + SxxIm*(SyyRe - SyyIm) )*vRe*dx
    F += ( SxxRe*(SyyRe - SyyIm) - SxxIm*(SyyRe + SyyIm) )*vIm*dx
    # - inner(Sxy*Sxy,v)*dx
    F -= ( SxyRe*(SxyRe + SxyIm) + SxyIm*(SxyRe - SxyIm) )*vRe*dx
    F -= ( SxyRe*(SxyRe - SxyIm) - SxyIm*(SxyRe + SxyIm) )*vIm*dx

    # (f*v*dx - gy*muxy*ds(1) + gx*muxy*ds(2) + gy*muxy*ds(3) - gx*muxy*ds(4))
    F -= (f*vRe*dx - gy*muxyRe*ds(1) + gx*muxyRe*ds(2) + gy*muxyRe*ds(3) - gx*muxyRe*ds(4));
    F -= (f*vIm*dx - gy*muxyIm*ds(1) + gx*muxyIm*ds(2) + gy*muxyIm*ds(3) - gx*muxyIm*ds(4));

    return F;





def dFdEps_Form(MixedV, ds, ep, f, gx, gy):
    (Sxx, Sxy, Syy, u) = TrialFunction(MixedV)
    (muxx, muxy, muyy, v) = TestFunction(MixedV)

    F = (inner(Dx(Sxx,0), Dx(v,0)) + inner(Dx(Sxy,0), Dx(v,1)))*dx;
    F += (inner(Dx(Sxy,1), Dx(v,0)) + inner(Dx(Syy,1), Dx(v,1)))*dx;

    F += (muxx + muyy)*dx + (Sxx+Syy + 2*ep)*v*dx;


    return F;


