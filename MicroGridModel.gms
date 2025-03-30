$onText
========================================================================
For more details please refer to
Lecture 5 : Course "Optimisation in Modern Electricity Markets"
------------------------------------------------------------------------
Model type: LP
Short Description: Capacity Expansion Model: Benders decomposition
Normally Distributed Demand Scenarios
------------------------------------------------------------------------
Developed by
Dr. Hannes Hobbie
Postdoctoral Researcher at the Chair of Energy Economics, TU Dresden
Extended by
Philipp Lahmeyer
========================================================================
$offText

Scalar
cc_multi  multiplier of construction costs /1.2/
without_sh   whether to use shared plant(0) or not(1) /0/
;


Set
t        timesteps /t1*t24/
i        technologies /i1*i3/

;
*Micro grid addition
*********************************
Set
s         base scenarios /s1*s5/
ss(s)     relevant scenarios /s1,s2,s3/
z         combined scenarios /z1*z27/
nodes     nodes /mg1, mg2, mg3, pp1/
mg(nodes) microgrids /mg1*mg3/
sh(nodes) shared power plants /pp1/
i_mg(i)   technologies in microgrids /i1, i2/
i_sh(i)   shared technologies  /i3/;

Parameter
d_s(t,s)     Base demand scenario
pb_s(ss)     Scenario probabilities /s1 0.3, s2 0.55, s3 0.15/
*pb_s(ss)     Scenario probabilities /s4 0.3, s2 0.55, s5 0.15/
*pb_s(ss)     Scenario probabilities /s1 0.2, s2 0.5, s4 0.25/
*pb_s(ss)     Scenario probabilities /s3 0.2, s2 0.6, s5 0.2/
mg_size(mg)  size multiplier per microgri /mg1 1.0, mg2 1.2, mg3 2/
d_z(z, mg)   actual scenarios as combination from base scenarios
;


$onEmbeddedCode Connect:
- CSVReader:
    file: scenarios.csv
    name: d_z
    indexColumns: 1
    valueColumns: "2:lastCol"
    
- GAMSWriter:
    symbols: all
$offEmbeddedCode


*********************************

Set
iter     iterations / iter1 * iter1000 /
nu(iter) current iteration
k(iter)  all iterations until current iteration
;

Parameters
oc(i)       operating cost /i1 15.5, i2 42, i3 58.4/
cc(i)       capacity cost /i1 500, i2 384, i3 301/
pc(i)       provision costs for unused capacity /i1 0, i2 0, i3 5/
voll        Value of lost load EUR per MWh /400/

d(mg, t)

pb
pb_mg(mg)

ls_s(t,z,mg)

zm_l(iter)          level value master problem objective
zs_l(z,iter)        level value subproblem objective
c_l(iter,nodes,i)         level value capacity variables
lambda_l(z,iter,nodes, i)  dual values of capacity fixation
alpha_l(z,iter)     level value alpha variables
lb(iter)            lower bound
ub(iter)            upper bound
conv(iter)          convergence level
;

cc(i) = cc(i) * cc_multi;

execute_load  'demand.gdx' d_s ;

Scalar
epsilon / 0.00001 /
EndOfWhile / 0 /
iteration
;

Variable
ZM          objective master problem
ZS          objective subproblem
ALPHA(z)    Benders cuts
;

Positive variables
G(nodes, i, t)
G_mg(nodes,i,t)
G_sh(nodes, mg, i, t) Power to cover microgrid demand
LS(mg, t)
C(nodes, i)  installed capacity MW
;


Equations
master_problem
sub_problem
benders_cut
alpha_start
benders_fixation
energy_balance_mg
gen_constraint
cap_constrained_mg
cap_constrained_sh
init_constraint
gen_sum
;


* MASTER PROBLEM

master_problem..         ZM =e= sum(i_mg, sum(mg, C(mg, i_mg))*cc(i_mg)) + sum(i_sh, sum(sh, C(sh, i_sh) * cc(i_sh))) + sum(z, ALPHA(z));
cap_constrained_mg(mg).. sum(i_sh, C(mg,i_sh)) =e= 0;
cap_constrained_sh..     sum((sh,i_mg), C(sh,i_mg)) + sum((sh, i_sh), C(sh, i_sh)) * without_sh =e= 0;
init_constraint(mg)..    sum(i_mg, C(mg, i_mg)) =g= smax((t,ss), d_s(t,ss));
benders_cut(k,z)..       ALPHA(z) - sum((nodes, i),lambda_l(z,k,nodes,i)*(C(nodes, i)-c_l(k,nodes,i))) =g= zs_l(z,k);
alpha_start(z)..         ALPHA(z) =g= 0 ;


* SUBPROBLEM
sub_problem..                    ZS =e= pb *(sum((t,i), sum(nodes, G(nodes,i,t))*oc(i))
                                        +sum(i,sum((sh,t), C(sh,i) - G(sh, i, t)) * pc(i))
                                        +sum((mg, t), voll*LS(mg,t))) ;
benders_fixation(nodes, i)..     C(nodes, i) =e= sum(nu, c_l(nu,nodes,i)) ;
energy_balance_mg(mg, t)..       sum(i, G_mg(mg,i,t) + sum(sh, G_sh(sh, mg, i, t))) + LS(mg,t) =e= d(mg, t)  ;
gen_constraint(nodes,i,t)..      G(nodes,i,t) =l= C(nodes, i) ;
gen_sum(nodes, i, t)..           G(nodes,i,t) =e= sum(mg, G_sh(nodes, mg, i, t)) + G_mg(nodes,i,t);


* build and solve initialisation
nu('iter1') = YES ;
iteration = 1 ;
model Initialise / master_problem, alpha_start, init_constraint / ;
solve Initialise minimizing ZM using lp ;

* fix capacities
c_l(nu,mg,i) = C.l(mg,i) ;

* record objective value of initialisation
zm_l(nu) = ZM.l ;

* calculate lower bound
lb(nu) = zm_l(nu) ;

* record alpha
alpha_l(z,nu) = ALPHA.l(z) ;


*build master and subproblem model
model Masterproblem / master_problem, alpha_start, benders_cut, cap_constrained_mg, cap_constrained_sh / ;
model Subproblem / sub_problem, benders_fixation, energy_balance_mg, gen_constraint, gen_sum/ ;


While(EndOfWhile = 0,
    Loop(z,
* Set demand scenario and probabilities
        Loop(mg,
            d(mg,t) = sum(ss$(ord(ss)=d_z(z, mg)), d_s(t,ss))*mg_size(mg);
            pb_mg(mg) = sum(ss$(ord(ss)=d_z(z, mg)), pb_s(ss));
        );
   
        pb = prod(ss, pb_s(ss));

    solve Subproblem minimizing ZS using lp ;

* record dual value of capacity fixation
    lambda_l(z,nu,nodes, i) = benders_fixation.m(nodes,i) ;

* record objective value of subproblem
    zs_l(z,nu) = ZS.l ;

*record variables
    ls_s(t,z,mg) = LS.l(mg,t);
    ) ;
* calculate upper bound
    ub(nu) = zm_l(nu) + sum(z, zs_l(z,nu)) - sum(z, alpha_l(z,nu)) ;

* convergence check
    conv(nu) = abs(ub(nu) - lb(nu)) / abs(ub(nu)) ;
    If (sum(nu, conv(nu)) < epsilon , EndOfWhile = 1 ;

* update information
        Else iteration = iteration + 1 ;
        nu(iter) = NO ;
        nu(iter) $ (ord(iter) = iteration ) = YES ;
        k(iter) $ (ord(iter) = iteration - 1) = YES ;

        solve Masterproblem minimizing ZM using lp ;

*fix capacity values
        c_l(nu,nodes,i) = C.l(nodes,i) ;

* calculate objective value of master problem
        zm_l(nu) = ZM.l ;

* record alpha
       alpha_l(z,nu) = ALPHA.l(z) ;

* calculate lower bound
       lb(nu) = zm_l(nu) ;
   );
   If (iteration = 100, EndOfWhile = 1 ) ;
 );

display C.l, ls_s, conv, c_l;
*should be i1 58.557,    i2 21.710,    i3 15.085

execute_unload 'results.gdx';
