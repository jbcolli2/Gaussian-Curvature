from dolfin import *
from VM_Utilities import *



###################################
#
#   Solver
#
############################

def ForwardProblem_MA(MixedV,ds, ep, initial, exact, f, gx, gy):
    bcv = DirichletBC(MixedV.sub(3), exact, Dir_boundary)
    bcxx = DirichletBC(MixedV.sub(0), ep, EW_boundary)
    bcyy = DirichletBC(MixedV.sub(2), ep, NS_boundary)
    bc = [bcxx,bcyy,bcv]

    # Define variational problem
    F = F_Form_MA(MixedV, ds, ep, f, gx, gy);

    # Solve problem
    R = action(F,initial);
    DR = derivative(R, initial);
    problem = NonlinearVariationalProblem(R,initial,bc,DR);
    solver = NonlinearVariationalSolver(problem);
    solver.solve();

    return initial;




def ForwardProblem_GC(MixedV,K,ds, ep, initial, exact, gx, gy):
    bcv = DirichletBC(MixedV.sub(3), exact, Dir_boundary)
    bcxx = DirichletBC(MixedV.sub(0), ep, EW_boundary)
    bcyy = DirichletBC(MixedV.sub(2), ep, NS_boundary)
    bc = [bcxx,bcyy,bcv]

    # Define variational problem
    F = F_Form_GC(MixedV, K, ds, ep, gx, gy);

    # initial = NewtonIteration(MixedV, initial, F, bc);
    #Solve problem
    R = action(F,initial);
    DR = derivative(R, initial);
    problem = NonlinearVariationalProblem(R,initial,bc,DR);
    solver = NonlinearVariationalSolver(problem);
    # solver.parameters['newton_solver']['absolute_tolerance'] = 1e-9
    prm = solver.parameters
    solver.solve();


    Fw = action(F,initial);
    DR = derivative(Fw,initial);
    A,b = assemble_system(DR,Fw,bc);
    print "My computed residual = ", np.linalg.norm(b, 2)

    return initial, problem;





def ForwardProblemWithBoundary(MixedV,ds, ep, initial, exact, f, gx, gy):
    bcv = DirichletBC(MixedV.sub(3), exact, Dir_boundary)
    bcxx = DirichletBC(MixedV.sub(0), 0.0, EW_boundary)
    bcyy = DirichletBC(MixedV.sub(2), 0.0, NS_boundary)
    bc = [bcxx,bcyy,bcv]

    # Define variational problem
    F = F_FormWithBoundary(MixedV, ds, ep, f, gx, gy);

    # Solve problem
    R = action(F,initial);
    DR = derivative(R, initial);
    problem = NonlinearVariationalProblem(R,initial,bc,DR);
    solver = NonlinearVariationalSolver(problem);
    solver.solve();

    return initial;




def ForwardComplexProblem(MixedVComplex,ds, ep, initial, exact, f, gx, gy):
    bcvRe = DirichletBC(MixedVComplex.sub(3), exact, Dir_boundary)
    bcvIm = DirichletBC(MixedVComplex.sub(7), 0.0, Dir_boundary)
    bcxxRe = DirichletBC(MixedVComplex.sub(0), ep.real, EW_boundary)
    bcyyRe = DirichletBC(MixedVComplex.sub(2), ep.real, NS_boundary)
    bcxxIm = DirichletBC(MixedVComplex.sub(4), ep.imag, EW_boundary)
    bcyyIm = DirichletBC(MixedVComplex.sub(6), ep.imag, NS_boundary)

    bcRe = [bcxxRe,bcyyRe,bcvRe]
    bcIm = [bcxxIm, bcyyIm, bcvIm]
    bc = [bcxxRe,bcyyRe,bcvRe,bcxxIm,bcyyIm,bcvIm]

    # Define variational problem
    F = F_ComplexForm(MixedVComplex, ds, ep, f, gx, gy);

    # Solve problem
    R = action(F,initial);
    DR = derivative(R, initial);
    problem = NonlinearVariationalProblem(R,initial,bc,DR);
    solver = NonlinearVariationalSolver(problem);
    solver.solve();

    return initial;





def SolveDavidenko(MixedV, ds, w, ep, exact, f, gx, gy):
    bcv = DirichletBC(MixedV.sub(3), exact, Dir_boundary)
    bcxx = DirichletBC(MixedV.sub(0), 0.0, EW_boundary)
    bcyy = DirichletBC(MixedV.sub(2), 0.0, NS_boundary)
    bc = [bcxx,bcyy,bcv]

    dFdEps = action(dFdEps_Form(MixedV, ds, ep, f, gx, gy), w);
    F_Eps = F_FormWithBoundary(MixedV, ds, ep, f, gx, gy);

    dWdEps = Function(MixedV);
    R = action(F_Eps,dWdEps)
    dFdu = derivative(R,dWdEps)
    dFdu = replace(dFdu,{dWdEps: w})
    solve(dFdu == -dFdEps, dWdEps, bcs=bc);

    dudEps = dWdEps.sub(3,deepcopy=True);

    return dudEps;
