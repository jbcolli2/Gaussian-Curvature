from dolfin import *
from VM_Utilities import *



###################################
#
#   Solver
#
############################

def ForwardProblem(MixedV,ds, ep, initial, exact, f, gx, gy):
    bcv = DirichletBC(MixedV.sub(3), exact, Dir_boundary)
    bcxx = DirichletBC(MixedV.sub(0), ep, EW_boundary)
    bcyy = DirichletBC(MixedV.sub(2), ep, NS_boundary)
    bc = [bcxx,bcyy,bcv]

    # Define variational problem
    F = F_Form(MixedV, ds, ep, f, gx, gy);

    # Solve problem
    R = action(F,initial);
    DR = derivative(R, initial);
    problem = NonlinearVariationalProblem(R,initial,bc,DR);
    solver = NonlinearVariationalSolver(problem);
    solver.solve();

    return initial;



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
