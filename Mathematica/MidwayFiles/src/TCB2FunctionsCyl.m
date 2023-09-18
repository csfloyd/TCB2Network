(* ::Package:: *)

BS[CB_]:=CSat^2/4-(CB-CSat/2)^2;(*concentration of available binding sites*)
gaussian[r_,r0_,s_]:=Exp[-(1/2) ((r-r0)/s)^2];
spaceFunc[r_,r0_,s_]:=1/2*(1-Tanh[(r-r0)/s]);


sigmoid[x_]:=1/(1+Exp[-x]);
smoothBump[x_,a_,xa_,xb_]:=sigmoid[(x-xa)/a]-sigmoid[(x-xb)/a]
waveFunc[t_,cyc_,len_,del_]:=Sum[smoothBump[t,del,cyc*i,cyc*i+len],{i,0,Floor[tMax/cyc]}];
trange = Range[0,tMax,0.01];
wPoints = waveFunc[trange,cyc,len,del];
data = Transpose[{trange,wPoints}];
iFun = Interpolation[data, InterpolationOrder -> 1];
gamFunc[r_,t_]:=spaceFunc[r,r0,width]*iFun[t];
(*external light field*)

pA[CBA_,CBI_]:=CBA/(CBA+CBI);
g[CBA_,CBI_]:=-(1-gMin)*pA[CBA,CBI];



RDI[r_,t_,CDI_,CDA_,CBI_,CBA_,CC_,CD_,CDst_,CE_,CEst_]:=-kAct*CDI*CC+kInact*CDA-kIBind*CDI*BS[CBA+CBI]+kIUnbind*CBI;
RDA[r_,t_,CDI_,CDA_,CBI_,CBA_,CC_,CD_,CDst_,CE_,CEst_]:=kAct*CDI*CC-kInact*CDA-kABind*CDA*BS[CBA+CBI]+kAUnbind*CBA;

RBI[r_,t_,CDI_,CDA_,CBI_,CBA_,CC_,CD_,CDst_,CE_,CEst_]:=-kAct*CBI*CC+kInact*CBA+kIBind*CDI*BS[CBA+CBI]-kIUnbind*CBI;
RBA[r_,t_,CDI_,CDA_,CBI_,CBA_,CC_,CD_,CDst_,CE_,CEst_]:=kAct*CBI*CC-kInact*CBA+kABind*CDA*BS[CBA+CBI]-kAUnbind*CBA;

RC[r_,t_,CDI_,CDA_,CBI_,CBA_,CC_,CD_,CDst_,CE_,CEst_]:=-kAct*(CDI+CBI)*CC+kInact*(CDA+CBA)-kTrap*CDst*CC+kRel*gamFunc[r,t]*CD-kETrap*CEst*CC+kERel*CE;
RD[r_,t_,CDI_,CDA_,CBI_,CBA_,CC_,CD_,CDst_,CE_,CEst_]:=kTrap*CDst*CC-kRel*gamFunc[r,t]*CD;
RDst[r_,t_,CDI_,CDA_,CBI_,CBA_,CC_,CD_,CDst_,CE_,CEst_]:=-kTrap*CDst*CC+kRel*gamFunc[r,t]*CD*\[Beta];
RE[r_,t_,CDI_,CDA_,CBI_,CBA_,CC_,CD_,CDst_,CE_,CEst_]:=-kETrap*CEst*CC+kERel*CE;
REst[r_,t_,CDI_,CDA_,CBI_,CBA_,CC_,CD_,CDst_,CE_,CEst_]:=-kETrap*CEst*CC-kERel*CE;


gradpA[CBA_,CBI_]:=(gMin-1)/(CBA+CBI)*((1-pA[CBA,CBI])*D[CBA,r]-pA[CBA,CBI]*D[CBI,r]);
force[r_,u_,CBA_,CBI_]:=2 \[Mu]0*(D[CBA+CBI,r]*(D[u,r]-1/2 g[CBA,CBI])+(CBA+CBI)*(D[u,{r,2}]+1/r D[u,r]-u/r^2-1/2 gradpA[CBA,CBI]))+\[Lambda]0*(D[CBA+CBI,r]*(u/r+D[u,r]-g[CBA,CBI])+(CBA+CBI)*(D[u,{r,2}]+1/r D[u,r]-u/r^2-gradpA[CBA,CBI]));


If[advectionBool,
pdes={
D[CDI[r,t],t]==RDI[r,t,CDI[r,t],CDA[r,t],CBI[r,t],CBA[r,t],CC[r,t],CD[r,t],CDst[r,t],CE[r,t],CEst[r,t]]+(DT*DFunc[CBI[r,t]+CBA[r,t]]*D[r*D[CDI[r,t],r],r]/r)+DT*D[DFunc[CBI[r,t]+CBA[r,t]],r]*D[CDI[r,t],r],
D[CDA[r,t],t]==RDA[r,t,CDI[r,t],CDA[r,t],CBI[r,t],CBA[r,t],CC[r,t],CD[r,t],CDst[r,t],CE[r,t],CEst[r,t]]+(DT*DFunc[CBI[r,t]+CBA[r,t]]*D[r*D[CDA[r,t],r],r]/r)+DT*D[DFunc[CBI[r,t]+CBA[r,t]],r]*D[CDA[r,t],r],
 
(*D[CBI[r,t],t]==RBI[r,t,CDI[r,t],CDA[r,t],CBI[r,t],CBA[r,t],CC[r,t],CD[r,t],CDst[r,t],CE[r,t],CEst[r,t]]-D[force[r,u[r,t],CBA[r,t],CBI[r,t]] * CBI[r,t],r] -(force[r,u[r,t],CBA[r,t],CBI[r,t]] * CBI[r,t])/r,
D[CBA[r,t],t]==RBA[r,t,CDI[r,t],CDA[r,t],CBI[r,t],CBA[r,t],CC[r,t],CD[r,t],CDst[r,t],CE[r,t],CEst[r,t]]-D[force[r,u[r,t],CBA[r,t],CBI[r,t]] * CBA[r,t],r] -(force[r,u[r,t],CBA[r,t],CBI[r,t]] * CBA[r,t])/r,*)
D[CBI[r,t],t]==RBI[r,t,CDI[r,t],CDA[r,t],CBI[r,t],CBA[r,t],CC[r,t],CD[r,t],CDst[r,t],CE[r,t],CEst[r,t]]-force[r,u[r,t],CBA[r,t],CBI[r,t]]*(D[CBI[r,t],r] +CBI[r,t]/r),
D[CBA[r,t],t]==RBA[r,t,CDI[r,t],CDA[r,t],CBI[r,t],CBA[r,t],CC[r,t],CD[r,t],CDst[r,t],CE[r,t],CEst[r,t]]-force[r,u[r,t],CBA[r,t],CBI[r,t]]*(D[CBA[r,t],r] +CBA[r,t]/r),

D[CC[r,t],t]==RC[r,t,CDI[r,t],CDA[r,t],CBI[r,t],CBA[r,t],CC[r,t],CD[r,t],CDst[r,t],CE[r,t],CEst[r,t]]+(DC*DFunc[CBI[r,t]+CBA[r,t]]*D[r*D[CC[r,t],r],r]/r)+DC*D[DFunc[CBI[r,t]+CBA[r,t]],r]*D[CC[r,t],r],
D[CD[r,t],t]==RD[r,t,CDI[r,t],CDA[r,t],CBI[r,t],CBA[r,t],CC[r,t],CD[r,t],CDst[r,t],CE[r,t],CEst[r,t]]+(DD*DFunc[CBI[r,t]+CBA[r,t]]*D[r*D[CD[r,t],r],r]/r)+DD*D[DFunc[CBI[r,t]+CBA[r,t]],r]*D[CD[r,t],r],
D[CDst[r,t],t]==RDst[r,t,CDI[r,t],CDA[r,t],CBI[r,t],CBA[r,t],CC[r,t],CD[r,t],CDst[r,t],CE[r,t],CEst[r,t]]+(DD*DFunc[CBI[r,t]+CBA[r,t]]*D[r*D[CDst[r,t],r],r]/r)+DD*D[DFunc[CBI[r,t]+CBA[r,t]],r]*D[CDst[r,t],r],
D[CE[r,t],t]==RE[r,t,CDI[r,t],CDA[r,t],CBI[r,t],CBA[r,t],CC[r,t],CD[r,t],CDst[r,t],CE[r,t],CEst[r,t]]+(DE*DFunc[CBI[r,t]+CBA[r,t]]*D[r*D[CE[r,t],r],r]/r)+DE*D[DFunc[CBI[r,t]+CBA[r,t]],r]*D[CE[r,t],r],
D[CEst[r,t],t]==REst[r,t,CDI[r,t],CDA[r,t],CBI[r,t],CBA[r,t],CC[r,t],CD[r,t],CDst[r,t],CE[r,t],CEst[r,t]]+(DE*DFunc[CBI[r,t]+CBA[r,t]]*D[r*D[CEst[r,t],r],r]/r)+DE*D[DFunc[CBI[r,t]+CBA[r,t]],r]*D[CEst[r,t],r],

D[u[r,t],t] == force[r,u[r,t],CBA[r,t],CBI[r,t]]
};
, (* advection off *)
pdes={
D[CDI[r,t],t]==RDI[r,t,CDI[r,t],CDA[r,t],CBI[r,t],CBA[r,t],CC[r,t],CD[r,t],CDst[r,t],CE[r,t],CEst[r,t]]+DT*DFunc[CBI[r,t]+CBA[r,t]]*D[r*D[CDI[r,t],r],r]/r+DT*D[DFunc[CBI[r,t]+CBA[r,t]],r]*D[CDI[r,t],r],
D[CDA[r,t],t]==RDA[r,t,CDI[r,t],CDA[r,t],CBI[r,t],CBA[r,t],CC[r,t],CD[r,t],CDst[r,t],CE[r,t],CEst[r,t]]+DT*DFunc[CBI[r,t]+CBA[r,t]]*D[r*D[CDA[r,t],r],r]/r+DT*D[DFunc[CBI[r,t]+CBA[r,t]],r]*D[CDA[r,t],r],

D[CBI[r,t],t] == RBI[r,t,CDI[r,t],CDA[r,t],CBI[r,t],CBA[r,t],CC[r,t],CD[r,t],CDst[r,t],CE[r,t],CEst[r,t]] ,
D[CBA[r,t],t] == RBA[r,t,CDI[r,t],CDA[r,t],CBI[r,t],CBA[r,t],CC[r,t],CD[r,t],CDst[r,t],CE[r,t],CEst[r,t]],

D[CC[r,t],t]==RC[r,t,CDI[r,t],CDA[r,t],CBI[r,t],CBA[r,t],CC[r,t],CD[r,t],CDst[r,t],CE[r,t],CEst[r,t]]+DC*DFunc[CBI[r,t]+CBA[r,t]]*D[r*D[CC[r,t],r],r]/r+DC*D[DFunc[CBI[r,t]+CBA[r,t]],r]*D[CC[r,t],r],
D[CD[r,t],t]==RD[r,t,CDI[r,t],CDA[r,t],CBI[r,t],CBA[r,t],CC[r,t],CD[r,t],CDst[r,t],CE[r,t],CEst[r,t]]+DD*DFunc[CBI[r,t]+CBA[r,t]]*D[r*D[CD[r,t],r],r]/r+DD*D[DFunc[CBI[r,t]+CBA[r,t]],r]*D[CD[r,t],r],
D[CDst[r,t],t]==RDst[r,t,CDI[r,t],CDA[r,t],CBI[r,t],CBA[r,t],CC[r,t],CD[r,t],CDst[r,t],CE[r,t],CEst[r,t]]+DD*DFunc[CBI[r,t]+CBA[r,t]]*D[r*D[CDst[r,t],r],r]/r+DD*D[DFunc[CBI[r,t]+CBA[r,t]],r]*D[CDst[r,t],r],
D[CE[r,t],t]==RE[r,t,CDI[r,t],CDA[r,t],CBI[r,t],CBA[r,t],CC[r,t],CD[r,t],CDst[r,t],CE[r,t],CEst[r,t]]+(DE*DFunc[CBI[r,t]+CBA[r,t]]*D[r*D[CE[r,t],r],r]/r)+DE*D[DFunc[CBI[r,t]+CBA[r,t]],r]*D[CE[r,t],r],
D[CEst[r,t],t]==REst[r,t,CDI[r,t],CDA[r,t],CBI[r,t],CBA[r,t],CC[r,t],CD[r,t],CDst[r,t],CE[r,t],CEst[r,t]]+(DE*DFunc[CBI[r,t]+CBA[r,t]]*D[r*D[CEst[r,t],r],r]/r)+DE*D[DFunc[CBI[r,t]+CBA[r,t]],r]*D[CEst[r,t],r],

D[u[r,t],t] == force[r,u[r,t],CBA[r,t],CBI[r,t]]
};
];


ics={
CDI[r,0]==CDI0,
CDA[r,0]==CDA0,

CBI[r,0]==CBI0,
CBA[r,0]==CBA0,

CC[r,0]==CC0,
CD[r,0]==CD0,
CDst[r,0]==CDst0,
CE[r,0]==CE0,
CEst[r,0]==CEst0,

u[r,0] == 0
};


If[advectionBool,
bcs={

Derivative[1,0][CDI][rMin,t]==0,
Derivative[1,0][CDI][rMax,t]==0,

Derivative[1,0][CDA][rMin,t]==0,
Derivative[1,0][CDA][rMax,t]==0,

(*Derivative[1,0][CBA][rMin,t]==0,
Derivative[1,0][CBA][rMax,t]==0,

Derivative[1,0][CBI][rMin,t]==0,
Derivative[1,0][CBI][rMax,t]==0,*)

Derivative[1,0][CC][rMin,t]==0,
Derivative[1,0][CC][rMax,t]==0,

Derivative[1,0][CD][rMin,t]==0,
Derivative[1,0][CD][rMax,t]==0,

Derivative[1,0][CDst][rMin,t]==0,
Derivative[1,0][CDst][rMax,t]==0,

Derivative[1,0][CE][rMin,t]==0,
Derivative[1,0][CE][rMax,t]==0,

Derivative[1,0][CEst][rMin,t]==0,
Derivative[1,0][CEst][rMax,t]==0,

u[rMin,t] == 0,
u[rMax, t] == 0,

Derivative[1,0][u][rMin,t] == 0,
Derivative[1,0][u][rMax,t] == 0};
,(* advection off*)
bcs ={
Derivative[1,0][CDI][rMin,t]==0,
Derivative[1,0][CDI][rMax,t]==0,

Derivative[1,0][CDA][rMin,t]==0,
Derivative[1,0][CDA][rMax,t]==0,

Derivative[1,0][CC][rMin,t]==0,
Derivative[1,0][CC][rMax,t]==0,

Derivative[1,0][CD][rMin,t]==0,
Derivative[1,0][CD][rMax,t]==0,

Derivative[1,0][CDst][rMin,t]==0,
Derivative[1,0][CDst][rMax,t]==0,

Derivative[1,0][CE][rMin,t]==0,
Derivative[1,0][CE][rMax,t]==0,

Derivative[1,0][CEst][rMin,t]==0,
Derivative[1,0][CEst][rMax,t]==0,

u[rMin,t] == 0,
u[rMax, t] == 0,

Derivative[1,0][u][rMin,t] == 0,
Derivative[1,0][u][rMax,t] == 0};
];




pdesCO={
D[CDI[r,t],t]==RDI[r,t,CDI[r,t],CDA[r,t],CBI[r,t],CBA[r,t],CC[r,t],CD[r,t],CDst[r,t]]+DT*D[r*(CDI[r,t]),{r,2}]/r,
D[CDA[r,t],t]==RDA[r,t,CDI[r,t],CDA[r,t],CBI[r,t],CBA[r,t],CC[r,t],CD[r,t],CDst[r,t]]+DTA*D[r*(CDA[r,t]),{r,2}]/r,

     
D[CBI[r,t],t] == RBI[r,t,CDI[r,t],CDA[r,t],CBI[r,t],CBA[r,t],CC[r,t],CD[r,t],CDst[r,t]],
D[CBA[r,t],t] == RBA[r,t,CDI[r,t],CDA[r,t],CBI[r,t],CBA[r,t],CC[r,t],CD[r,t],CDst[r,t]],

D[CC[r,t],t]==RC[r,t,CDI[r,t],CDA[r,t],CBI[r,t],CBA[r,t],CC[r,t],CD[r,t],CDst[r,t]]+DC*D[r*CC[r,t],{r,2}]/r,
D[CD[r,t],t]==RD[r,t,CDI[r,t],CDA[r,t],CBI[r,t],CBA[r,t],CC[r,t],CD[r,t],CDst[r,t]]+DD*D[r*CD[r,t],{r,2}]/r,
D[CDst[r,t],t]==RDst[r,t,CDI[r,t],CDA[r,t],CBI[r,t],CBA[r,t],CC[r,t],CD[r,t],CDst[r,t]]+DD*D[r*CDst[r,t],{r,2}]/r
};


icsCO={
CDI[r,0]==CDI0,
CDA[r,0]==CDA0,

CBI[r,0]==CBI0,
CBA[r,0]==CBA0,

CC[r,0]==CC0,
CD[r,0]==CD0,
CDst[r,0]==CDst0
};


bcsCO={

Derivative[1,0][CDI][rMin,t]==0,
Derivative[1,0][CDI][rMax,t]==0,

Derivative[1,0][CDA][rMin,t]==0,
Derivative[1,0][CDA][rMax,t]==0,

Derivative[1,0][CC][rMin,t]==0,
Derivative[1,0][CC][rMax,t]==0,

Derivative[1,0][CD][rMin,t]==0,
Derivative[1,0][CD][rMax,t]==0,

Derivative[1,0][CDst][rMin,t]==0,
Derivative[1,0][CDst][rMax,t]==0

};


InterpFunction[funcN_,rMinN_,rMaxN_,tMaxN_]:=Module[{func = funcN, rMin = rMinN, rMax = rMaxN, tMax = tMaxN},
	intDat=Table[
	{{r,t},Evaluate[func[r,t][[1]]]}
	,{r,rMin,rMax,1.},{t,0.,tMax,1.}];
	dn=Flatten[intDat,1];
	fN=Interpolation[dn,InterpolationOrder->2]
]
