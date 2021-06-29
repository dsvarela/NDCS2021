H = [1 -1; -1 2]; 
f = [-2; -6];
A = [1 1; -1 2; 2 1];
b = [2; 2; 3];
Aeq = [1 1];
beq = 0;

[x,fval,exitflag,output,lambda] = quadprog(H,f,A,b, Aeq, beq);


cvx_begin
    variable uopt(2)
    minimize (quad_form(uopt,H/2) + f'*uopt)
    subject to
        A*uopt <= b
        Aeq*uopt ==  beq; 
    cvx_solver('gurobi')
    cvx_solver_settings('NumericFocus',2)
cvx_end