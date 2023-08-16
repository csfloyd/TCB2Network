(* ::Package:: *)

ClearAll["Global`*"];

args = $CommandLine[[3;;]];
Print[args];
outputDir = args[[3]];
var = ToExpression[args[[2]]];

 

(* chemical parameters *)
kAct =     2*10^2; (* mM^-1 s^-1 activation rate of calcium *)
kInact =   10^2; (* s^-1 inactivation rate of calcium *)
kTrap =    3*10^2; (* mM^-1 s^-1 trapping rate of calcium by DMNP *)
kRel =     7*10^2; (* s^-1 release rate of calcium by DMNP when light is on*)
kIBind =   10^0; (* mM^-1 s^-1 binding rate of inactivated TCB2 *)
kIUnbind = 10^(var[[1]]); (* s^-1 unbinding rate of inactivated TCB2 *)
kABind =   10^0; (* mM^-1 s^-1 binding rate of activated TCB2 *)
kAUnbind = 0; (* s^-1 unbinding rate of activated TCB2 *)
\[Beta]        = 1; (* degradation factor *)

CSat = 5; (* mM, saturation concentration of bound TCB2 *)

DT = 10^(var[[2]]); (* \[Mu]m^2/s, diffusion constant for TCB2 *)
DD = 50; (* \[Mu]m^2/s, diffusion constant for DMNP *)
DC = 500;(* \[Mu]m^2/s, diffusion constant for calcium *) 


(* initial values *)
CDI0 =  2.5; (* mM *)
CDA0 =  0; (* mM *)
CBI0 =  0.1; (* mM *)
CBA0 =  0; (* mM *)
CC0 =   0; (* mM *)
CDst0 = 20; (* mM *)
CD0 =   20; (* mM *)

(* domain values *)
xMax = 300; (* \[Mu]m *)
yMax = 300; (* \[Mu]m *)
tMax = 150; (* s *)
pts =  100; 

(* light activation values *)
width = 5; (* \[Mu]m, width of gaussian beam *)
r0 =    75; (* \[Mu]m, radius of circle *)
cyc =   20;
len =   5; (* s *)
del =   1; (* s *)

(* mechanical values *)
\[Lambda]0 = 3;(* \[Mu]m^2/s, 1st Lame modulus slope *)
\[Mu]0 = 3; (* \[Mu]m^2/s, 2nd Lame modulus slope *)
gMin = 0.4; (* shrinking factor of activated TCB2 *)


(* import functions *)
Get["/project/svaikunt/csfloyd/TCB2/src/TCB2Functions.m"]

ret = AbsoluteTiming[
    NDSolve[{pdes, bcs, ics},
   {CDI[x, y, t], CDA[x, y, t], CBI[x, y, t], CBA[x, y, t], 
    CC[x, y, t], CD[x, y, t], CDst[x, y, t], ux[x, y, t], uy[x, y, t]},
   {x, 0, xMax},
   {y, 0, yMax},
   {t, 0, tMax},
   Method -> {"MethodOfLines", "TemporalVariable" -> t, 
     "SpatialDiscretization" -> {"TensorProductGrid", 
	    "MinPoints" -> 60, "MaxPoints" -> pts, "DifferenceOrder" -> 6},
     Method -> {"IDA", "ImplicitSolver" -> {"GMRES"}}},
	    PrecisionGoal -> 6,
	    AccuracyGoal -> 8,
	    MaxSteps->1000000]
		     ];

Print[ret[[1]]];
sol = ret[[2]];


CDISol[x_,y_,t_]:=Evaluate[CDI[x,y,t]/.sol];
CDASol[x_,y_,t_]:=Evaluate[CDA[x,y,t]/.sol];
CBISol[x_,y_,t_]:=Evaluate[CBI[x,y,t]/.sol];
CBASol[x_,y_,t_]:=Evaluate[CBA[x,y,t]/.sol];
CCSol[x_,y_,t_]:=Evaluate[CC[x,y,t]/.sol];
CDSol[x_,y_,t_]:=Evaluate[CD[x,y,t]/.sol];
CDstSol[x_,y_,t_]:=Evaluate[CDst[x,y,t]/.sol];
uxSol[x_,y_,t_]:=Evaluate[ux[x,y,t]/.sol];
uySol[x_,y_,t_]:=Evaluate[uy[x,y,t]/.sol];

interpList = {
	InterpFunction[CDISol,xMax,yMax,tMax],
	InterpFunction[CDASol,xMax,yMax,tMax],
	InterpFunction[CBISol,xMax,yMax,tMax],
	InterpFunction[CBASol,xMax,yMax,tMax],
	InterpFunction[CCSol,xMax,yMax,tMax],
	InterpFunction[CDSol,xMax,yMax,tMax],
	InterpFunction[CDstSol,xMax,yMax,tMax],
	InterpFunction[uxSol,xMax,yMax,tMax],
	InterpFunction[uySol,xMax,yMax,tMax]
	};
	
Export[outputDir <> "/sol.mx", interpList];
 

 

