(* ::Package:: *)

ClearAll["Global`*"];

args = $CommandLine[[3;;]];
Print[args];
outputDir = args[[3]];
var = ToExpression[args[[2]]];
 
protocol = "o"; (* "o", "p" *)

kinFac = 1;
(* chemical parameters *)
kAct =     5*10^1; (* mM^-1 s^-1 activation rate of calcium *)
kInact =   5*10^0; (* s^-1 inactivation rate of calcium *)
kTrap =    3*10^2; (* mM^-1 s^-1 trapping rate of calcium by DMNP *)
kRel =     7*10^2; (* s^-1 release rate of calcium by DMNP when light is on*)
kETrap =   2.7*10^3;(* mM^-1 s^-1 trapping rate of calcium by EGTA *)
kERel =    0.5; (* s^-1 release rate of calcium by EGTA *)
kIBind =   kinFac * 1*10^-1; (* mM^-1 s^-1 binding rate of inactivated TCB2 *)
kIUnbind = kinFac * 1^1; (* s^-1 unbinding rate of inactivated TCB2 *)
kABind =   kinFac * 1*10^-1; (* mM^-1 s^-1 binding rate of activated TCB2 *)
kAUnbind = 0; (* s^-1 unbinding rate of activated TCB2 *)
\[Beta]  = 0.9; (* degradation factor *)

CSat = 5; (* mM, saturation concentration of bound TCB2 *)
phiFac = CSat * 3;

DT = 1; (* \[Mu]m^2/s, diffusion constant for TCB2 *)
DD = 25; (* \[Mu]m^2/s, diffusion constant for DMNP *)
DE = 180; (* \[Mu]m^2/s, diffusion constant for EGTA *)
DC = 500 ;(* \[Mu]m^2/s, diffusion constant for calcium *) 

DFunc[CBtot_]:=1;
(*DFunc[CBtot_]:=((1-CBtot/phiFac)/(1+CBtot/phiFac))^2;*)
(*DFunc[CBtot_]:=1-CBtot/(2 * CSat);*)


(* initial values *)
CDI0 =  0*2.5; (* mM *)
CDA0 =  0; (* mM *)
CBI0 =  0.0000001; (* mM *)
CBA0 =  0; (* mM *)
CC0 =   0; (* mM *)
CDst0 = 0*20; (* mM *)
CD0 =   20; (* mM *)
CE0 =   0; (* mM *)
CEst0 = var; (* mM *)

(* domain values *)
rMin = 1;
rMax = 500; (* \[Mu]m *)
tMax = 100; (* s *)
pts =  250; 

If[protocol == "p",
	lent = 1;
	cyct = 30;
	,
	lent = 5;
	cyct = 1000;
];

(* light activation values *)
width = 5; (* \[Mu]m, width of gaussian beam *)
r0 =    75; (* \[Mu]m, radius of circle *)
len =   lent; (* s *)
cyc =   cyct + len;
del =   0.25; (* s *)

(* mechanical values *)
\[Lambda]0 = 3; (* \[Mu]m^2/s, 1st Lame modulus slope *)
\[Mu]0 = 3; (* \[Mu]m^2/s, 2nd Lame modulus slope *)
gMin = 0.4; (* shrinking factor of activated TCB2 *)
advectionBool = False;

(* import functions *)
Get["/project/svaikunt/csfloyd/TCB2/src/TCB2FunctionsCyl.m"]

Print["here"];

ret=AbsoluteTiming[NDSolve[{pdes,bcs,ics},
{CDI[r,t],CDA[r,t],CBI[r,t],CBA[r,t],CC[r,t],CD[r,t],CDst[r,t],CE[r,t],CEst[r,t],u[r,t]},
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
CESol[r_,t_]:=Evaluate[CE[r,t]/.sol];
CEstSol[r_,t_]:=Evaluate[CEst[r,t]/.sol];
uSol[r_,t_]:=Evaluate[u[r,t]/.sol];

interpList = {
        InterpFunction[CDISol,rMin,rMax,tMax],
	InterpFunction[CDASol,rMin,rMax,tMax],
	InterpFunction[CBISol,rMin,rMax,tMax],
	InterpFunction[CBASol,rMin,rMax,tMax],
	InterpFunction[CCSol,rMin,rMax,tMax],
	InterpFunction[CDSol,rMin,rMax,tMax],
	InterpFunction[CDstSol,rMin,rMax,tMax],
	InterpFunction[CESol,rMin,rMax,tMax],
	InterpFunction[CEstSol,rMin,rMax,tMax],
	InterpFunction[uSol,rMin,rMax,tMax]
	};
	
Export[outputDir <> "/sol.mx", interpList];
 

 

