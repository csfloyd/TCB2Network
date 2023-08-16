(* ::Package:: *)

ClearAll["Global`*"];

args = $CommandLine[[3;;]];
Print[args];
outputDir = args[[3]];
var = ToExpression[args[[2]]];
 
CB0=1.0; (*mM*)

(*domain values*)
rMin=0.1; (*\[Mu]m*)
rMax=30; (*\[Mu]m*)
tMax=250; (*s*)
pts=300;


release = 2*tMax; (* turn off pulsing *)
flat = 0.2; (* homogeneity after illuminated region *)
aOff = 0.5; (* degree of accumulation *)
n = 2; (* sharpness of accumulation *)
v = 1.0;

(*light activation values*)
width=1; (*\[Mu]m,width of gaussian beam*)
r0=15; (*\[Mu]m,radius of circle*)
cyc=1;
len=0.25; (*s*)
del=0.025; (*s*)

(*mechanical values*)
\[Lambda]0=3;(*\[Mu]m^2/s,1st Lame modulus slope*)
\[Mu]0=3; (*\[Mu]m^2/s,2nd Lame modulus slope*)
gMin=var; (*shrinking factor of activated TCB2*)



(* import functions *)
Get["/project/svaikunt/csfloyd/TCB2/src/TCB2FunctionsND.m"]

retC=AbsoluteTiming[NDSolve[{pdesC,bcs,ics},{u[r,t]},{r,rMin,rMax},{t,0,tMax},
Method->{"MethodOfLines","TemporalVariable"->t,
"SpatialDiscretization"->{"TensorProductGrid","MinPoints"->400,"MaxPoints"->600,
"DifferenceOrder"->6},
Method->{"IDA","ImplicitSolver"->{"GMRES"}}},
PrecisionGoal->6,
AccuracyGoal->12,
MaxSteps->1000000]];

Print[retC[[1]]];
solC = retC[[2]];

retL=AbsoluteTiming[NDSolve[{pdesL,bcs,ics},{u[r,t]},{r,rMin,rMax},{t,0,tMax},
Method->{"MethodOfLines","TemporalVariable"->t,
"SpatialDiscretization"->{"TensorProductGrid","MinPoints"->400,"MaxPoints"->600,
"DifferenceOrder"->6},
Method->{"IDA","ImplicitSolver"->{"GMRES"}}},
PrecisionGoal->6,
AccuracyGoal->12,
MaxSteps->1000000]];

Print[retL[[1]]];
solL = retL[[2]];


uSolC[r_,t_]:=Evaluate[u[r,t]/.solC];
uSolL[r_,t_]:=Evaluate[u[r,t]/.solL];

interpList = {
	InterpFunction[uSolC,rMin,rMax,tMax],
	InterpFunction[uSolL,rMin,rMax,tMax]
	};
	
Export[outputDir <> "/sol.mx", interpList];
 

 

