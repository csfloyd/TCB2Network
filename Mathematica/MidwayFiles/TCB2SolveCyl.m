(* ::Package:: *)

ClearAll["Global`*"];

args = $CommandLine[[3;;]];
Print[args];
outputDir = args[[3]];
var = ToExpression[args[[2]]];
 
protocol = "p"; (* "o", "p" *)


(* chemical parameters *)
kAct =     5*10^1; (* mM^-1 s^-1 activation rate of calcium *)
kInact =   5*10^0; (* s^-1 inactivation rate of calcium *)
kTrap =    3*10^2; (* mM^-1 s^-1 trapping rate of calcium by DMNP *)
kRel =     7*10^2; (* s^-1 release rate of calcium by DMNP when light is on*)
kIBind =   1*10^-1; (* mM^-1 s^-1 binding rate of inactivated TCB2 *)
kIUnbind = 10^0; (* s^-1 unbinding rate of inactivated TCB2 *)
kABind =   1*10^-1; (* mM^-1 s^-1 binding rate of activated TCB2 *)
kAUnbind = 0; (* s^-1 unbinding rate of activated TCB2 *)
\[Beta]  = 0.5; (* degradation factor *)

CSat = 5; (* mM, saturation concentration of bound TCB2 *)

DT = 1; (* \[Mu]m^2/s, diffusion constant for TCB2 *)
DD = 25; (* \[Mu]m^2/s, diffusion constant for DMNP *)
DC = 500 ;(* \[Mu]m^2/s, diffusion constant for calcium *) 

(* initial values *)
CDI0 =  2.5; (* mM *)
CDA0 =  0; (* mM *)
CBI0 =  0.1; (* mM *)
CBA0 =  0; (* mM *)
CC0 =   0; (* mM *)
CDst0 = 20; (* mM *)
CD0 =   20; (* mM *)





(* domain values *)
rMin = 1;
rMax = 500; (* \[Mu]m *)
tMax = 1000; (* s *)
pts =  250; 

If[protocol == "p",
	lent = 1;
	cyct = 30;
	,
	lent = 150;
	cyct = 1000;
];

(* light activation values *)
width = 5; (* \[Mu]m, width of gaussian beam *)
r0 =    var; (* \[Mu]m, radius of circle *)
len =   lent; (* s *)
cyc =   cyct + len;
del =   0.25; (* s *)

(* mechanical values *)
\[Lambda]0 = 3; (* \[Mu]m^2/s, 1st Lame modulus slope *)
\[Mu]0 = 3; (* \[Mu]m^2/s, 2nd Lame modulus slope *)
gMin = 0.4; (* shrinking factor of activated TCB2 *)


(* import functions *)
Get["/project/svaikunt/csfloyd/TCB2/src/TCB2FunctionsCyl.m"]

Print["here"];

ret=AbsoluteTiming[NDSolve[{pdes,bcs,ics},
{CDI[r,t],CDA[r,t],CBI[r,t],CBA[r,t],CC[r,t],CD[r,t],CDst[r,t],u[r,t]},
{r,rMin,rMax},{t,0,tMax},
Method->{"MethodOfLines","TemporalVariable"->t,
"SpatialDiscretization"->{"TensorProductGrid","MinPoints"->120,"MaxPoints"->pts,"DifferenceOrder"->6},
Method->{"IDA","ImplicitSolver"->{"GMRES"}}},
PrecisionGoal->8,
AccuracyGoal->8,
MaxSteps->1000000]];

Print[ret[[1]]];
sol = ret[[2]];


CDISol[r_,t_]:=Evaluate[CDI[r,t]/.sol];
CDASol[r_,t_]:=Evaluate[CDA[r,t]/.sol];
CBISol[r_,t_]:=Evaluate[CBI[r,t]/.sol];
CBASol[r_,t_]:=Evaluate[CBA[r,t]/.sol];
CCSol[r_,t_]:=Evaluate[CC[r,t]/.sol];
CDSol[r_,t_]:=Evaluate[CD[r,t]/.sol];
CDstSol[r_,t_]:=Evaluate[CDst[r,t]/.sol];
uSol[r_,t_]:=Evaluate[u[r,t]/.sol];

interpList = {
        InterpFunction[CDISol,rMin,rMax,tMax],
	InterpFunction[CDASol,rMin,rMax,tMax],
	InterpFunction[CBISol,rMin,rMax,tMax],
	InterpFunction[CBASol,rMin,rMax,tMax],
	InterpFunction[CCSol,rMin,rMax,tMax],
	InterpFunction[CDSol,rMin,rMax,tMax],
	InterpFunction[CDstSol,rMin,rMax,tMax],
	InterpFunction[uSol,rMin,rMax,tMax]
	};
	
Export[outputDir <> "/sol.mx", interpList];
 

 

