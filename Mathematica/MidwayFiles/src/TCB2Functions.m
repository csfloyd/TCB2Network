(* ::Package:: *)

BS[CB_]:=CSat^2/4-(CB-CSat/2)^2;(*concentration of available binding sites*)
gaussian[r_,r0_,s_]:=Exp[-(1/2) ((r-r0)/s)^2];
spaceFunc[r_,r0_,s_]:=1/2*(1-Tanh[(r-r0)/s]);

rFunc[x_,y_,x0_,y0_]:=Sqrt[(x-x0)^2+(y-y0)^2];

sigmoid[x_]:=1/(1+Exp[-x]);
smoothBump[x_,a_,xa_,xb_]:=sigmoid[(x-xa)/a]-sigmoid[(x-xb)/a]
waveFunc[t_,cyc_,len_,del_]:=Sum[smoothBump[t,del,cyc*i,cyc*i+len],{i,0,Floor[tMax/cyc]}];
trange = Range[0,tMax,0.01];
wPoints = waveFunc[trange,cyc,len,del];
data = Transpose[{trange,wPoints}];
iFun = Interpolation[data, InterpolationOrder -> 1];
gamFunc[x_,y_,t_]:=spaceFunc[rFunc[x,y,xMax/2,yMax/2],r0,width]*iFun[t];
(*external light field*)

pA[CBA_,CBI_]:=CBA/(CBA+CBI);
g[CBA_,CBI_]:=-(1-gMin)*pA[CBA,CBI];



RDI[x_,y_,t_,CDI_,CDA_,CBI_,CBA_,CC_,CD_,CDst_]:=-kAct*CDI*CC+kInact*CDA-kIBind*CDI*BS[CBA+CBI]+kIUnbind*CBI;
RDA[x_,y_,t_,CDI_,CDA_,CBI_,CBA_,CC_,CD_,CDst_]:=kAct*CDI*CC-kInact*CDA-kABind*CDA*BS[CBA+CBI]+kAUnbind*CBA;

RBI[x_,y_,t_,CDI_,CDA_,CBI_,CBA_,CC_,CD_,CDst_]:=-kAct*CBI*CC+kInact*CBA+kIBind*CDI*BS[CBA+CBI]-kIUnbind*CBI;
RBA[x_,y_,t_,CDI_,CDA_,CBI_,CBA_,CC_,CD_,CDst_]:=kAct*CBI*CC-kInact*CBA+kABind*CDA*BS[CBA+CBI]-kAUnbind*CBA;

RC[x_,y_,t_,CDI_,CDA_,CBI_,CBA_,CC_,CD_,CDst_]:=-kAct*(CDI+CBI)*CC+kInact*(CDA+CBA)-kTrap*CDst*CC+kRel*gamFunc[x,y,t]*CD;
RD[x_,y_,t_,CDI_,CDA_,CBI_,CBA_,CC_,CD_,CDst_]:=kTrap*CDst*CC-kRel*gamFunc[x,y,t]*CD;
RDst[x_,y_,t_,CDI_,CDA_,CBI_,CBA_,CC_,CD_,CDst_]:=-kTrap*CDst*CC+kRel*gamFunc[x,y,t]*CD*\[Beta];


gradpA[ind_,CBA_,CBI_]:=(gMin-1)/(CBA+CBI)*((1-pA[CBA,CBI])*D[CBA,ind]-pA[CBA,CBI]*D[CBI,ind]);
strainxx[ux_,uy_]:=D[ux,x];
strainxy[ux_,uy_]:=1/2 (D[ux,y]+D[uy,x]);
strainyx[ux_,uy_]:=strainxy[ux,uy];
strainyy[ux_,uy_]:=D[uy,y];
divu[ux_,uy_]:=D[ux,x]+D[uy,y];
gradiTr[ind_,ux_,uy_]:=D[divu[ux,uy],ind];

lapui[ui_]:=D[ui,{x,2}]+D[ui,{y,2}]
forcex[x_,y_,ux_,uy_,CBA_,CBI_]:=2 \[Mu]0*(((D[CBA+CBI,x]*strainxx[ux,uy]+D[CBA+CBI,y]*strainxy[ux,uy])-1/2 g[CBA,CBI]*D[CBA+CBI,x])+(CBA+CBI)*(1/2 (gradiTr[x,ux,uy]+lapui[ux])-1/2 gradpA[x,CBA,CBI]))+\[Lambda]0*(D[CBA+CBI,x]*(divu[ux,uy]-g[CBA,CBI])+(CBA+CBI)*(gradiTr[x,ux,uy]-gradpA[x,CBA,CBI]));

forcey[x_,y_,ux_,uy_,CBA_,CBI_]:=2 \[Mu]0*(((D[CBA+CBI,x]*strainyx[ux,uy]+D[CBA+CBI,y]*strainyy[ux,uy])-1/2 g[CBA,CBI]*D[CBA+CBI,y])+(CBA+CBI)*(1/2 (gradiTr[y,ux,uy]+lapui[uy])-1/2 gradpA[y,CBA,CBI]))+\[Lambda]0*(D[CBA+CBI,y]*(divu[ux,uy]-g[CBA,CBI])+(CBA+CBI)*(gradiTr[y,ux,uy]-gradpA[y,CBA,CBI]));


pdes={
D[CDI[x,y,t],t]==RDI[x,y,t,CDI[x,y,t],CDA[x,y,t],CBI[x,y,t],CBA[x,y,t],CC[x,y,t],CD[x,y,t],CDst[x,y,t]]+DT*Laplacian[CDI[x,y,t],{x,y}],
D[CDA[x,y,t],t]==RDA[x,y,t,CDI[x,y,t],CDA[x,y,t],CBI[x,y,t],CBA[x,y,t],CC[x,y,t],CD[x,y,t],CDst[x,y,t]]+DT*Laplacian[CDA[x,y,t],{x,y}],
     
D[CBI[x, y, t], t] == 
 RBI[x, y, t, CDI[x, y, t], CDA[x, y, t], CBI[x, y, t], CBA[x, y, t], CC[x, y, t], CD[x, y, t], CDst[x, y, t]] ,
D[CBA[x, y, t], t] == 
 RBA[x, y, t, CDI[x, y, t], CDA[x, y, t], CBI[x, y, t], CBA[x, y, t],  CC[x, y, t], CD[x, y, t], CDst[x, y, t]],

D[CC[x,y,t],t]==RC[x,y,t,CDI[x,y,t],CDA[x,y,t],CBI[x,y,t],CBA[x,y,t],CC[x,y,t],CD[x,y,t],CDst[x,y,t]]+DC*Laplacian[CC[x,y,t],{x,y}],
D[CD[x,y,t],t]==RD[x,y,t,CDI[x,y,t],CDA[x,y,t],CBI[x,y,t],CBA[x,y,t],CC[x,y,t],CD[x,y,t],CDst[x,y,t]]+DD*Laplacian[CD[x,y,t],{x,y}],
D[CDst[x,y,t],t]==RDst[x,y,t,CDI[x,y,t],CDA[x,y,t],CBI[x,y,t],CBA[x,y,t],CC[x,y,t],CD[x,y,t],CDst[x,y,t]]+DD*Laplacian[CDst[x,y,t],{x,y}],

D[ux[x, y, t], t] == forcex[x, y, ux[x, y, t], uy[x, y, t], CBA[x, y, t], CBI[x, y, t]],
D[uy[x, y, t], t] == forcey[x, y, ux[x, y, t], uy[x, y, t], CBA[x, y, t], CBI[x, y, t]]
};


ics={
CDI[x,y,0]==CDI0,
CDA[x,y,0]==CDA0,

CBI[x,y,0]==CBI0,
CBA[x,y,0]==CBA0,

CC[x,y,0]==CC0,
CD[x,y,0]==CD0,
CDst[x,y,0]==CDst0,

ux[x,y,0] == 0,
uy[x,y,0] == 0
};


bcs={

Derivative[1,0,0][CDI][0,y,t]==0,
Derivative[1,0,0][CDI][xMax,y,t]==0,
Derivative[0,1,0][CDI][x,0,t]==0,
Derivative[0,1,0][CDI][x,yMax,t]==0,

Derivative[1,0,0][CDA][0,y,t]==0,
Derivative[1,0,0][CDA][xMax,y,t]==0,
Derivative[0,1,0][CDA][x,0,t]==0,
Derivative[0,1,0][CDA][x,yMax,t]==0,

Derivative[1,0,0][CC][0,y,t]==0,
Derivative[1,0,0][CC][xMax,y,t]==0,
Derivative[0,1,0][CC][x,0,t]==0,
Derivative[0,1,0][CC][x,yMax,t]==0,

Derivative[1,0,0][CD][0,y,t]==0,
Derivative[1,0,0][CD][xMax,y,t]==0,
Derivative[0,1,0][CD][x,0,t]==0,
Derivative[0,1,0][CD][x,yMax,t]==0,

Derivative[1,0,0][CDst][0,y,t]==0,
Derivative[1,0,0][CDst][xMax,y,t]==0,
Derivative[0,1,0][CDst][x,0,t]==0,
Derivative[0,1,0][CDst][x,yMax,t]==0,

ux[0, y, t] == 0,
ux[xMax, y, t] == 0,
ux[x, 0, t] == 0,
ux[x, yMax, t] == 0,

Derivative[1, 0, 0][ux][0, y, t] == 0,
Derivative[1, 0, 0][ux][xMax, y, t] == 0,
Derivative[0, 1, 0][ux][x, 0, t] == 0,
Derivative[0, 1, 0][ux][x, yMax, t] == 0,

uy[0, y, t] == 0,
uy[xMax, y, t] == 0,
uy[x, 0, t] == 0,
uy[x, yMax, t] == 0,

Derivative[1, 0, 0][uy][0, y, t] == 0,
Derivative[1, 0, 0][uy][xMax, y, t] == 0,
Derivative[0, 1, 0][uy][x, 0, t] == 0,
Derivative[0, 1, 0][uy][x, yMax, t] == 0
};


pdesCO={
D[CDI[x,y,t],t]==RDI[x,y,t,CDI[x,y,t],CDA[x,y,t],CBI[x,y,t],CBA[x,y,t],CC[x,y,t],CD[x,y,t],CDst[x,y,t]]+DT*Laplacian[CDI[x,y,t],{x,y}],
D[CDA[x,y,t],t]==RDA[x,y,t,CDI[x,y,t],CDA[x,y,t],CBI[x,y,t],CBA[x,y,t],CC[x,y,t],CD[x,y,t],CDst[x,y,t]]+DT*Laplacian[CDA[x,y,t],{x,y}],

D[CBI[x,y,t],t]==RBI[x,y,t,CDI[x,y,t],CDA[x,y,t],CBI[x,y,t],CBA[x,y,t],CC[x,y,t],CD[x,y,t],CDst[x,y,t]],
D[CBA[x,y,t],t]==RBA[x,y,t,CDI[x,y,t],CDA[x,y,t],CBI[x,y,t],CBA[x,y,t],CC[x,y,t],CD[x,y,t],CDst[x,y,t]],

D[CC[x,y,t],t]==RC[x,y,t,CDI[x,y,t],CDA[x,y,t],CBI[x,y,t],CBA[x,y,t],CC[x,y,t],CD[x,y,t],CDst[x,y,t]]+DC*Laplacian[CC[x,y,t],{x,y}],
D[CD[x,y,t],t]==RD[x,y,t,CDI[x,y,t],CDA[x,y,t],CBI[x,y,t],CBA[x,y,t],CC[x,y,t],CD[x,y,t],CDst[x,y,t]]+DD*Laplacian[CD[x,y,t],{x,y}],
D[CDst[x,y,t],t]==RDst[x,y,t,CDI[x,y,t],CDA[x,y,t],CBI[x,y,t],CBA[x,y,t],CC[x,y,t],CD[x,y,t],CDst[x,y,t]]+DD*Laplacian[CDst[x,y,t],{x,y}]
};


icsCO={
CDI[x,y,0]==CDI0,
CDA[x,y,0]==CDA0,

CBI[x,y,0]==CBI0,
CBA[x,y,0]==CBA0,

CC[x,y,0]==CC0,
CD[x,y,0]==CD0,
CDst[x,y,0]==CDst0
};


bcsCO={

Derivative[1,0,0][CDI][0,y,t]==0,
Derivative[1,0,0][CDI][xMax,y,t]==0,
Derivative[0,1,0][CDI][x,0,t]==0,
Derivative[0,1,0][CDI][x,yMax,t]==0,

Derivative[1,0,0][CDA][0,y,t]==0,
Derivative[1,0,0][CDA][xMax,y,t]==0,
Derivative[0,1,0][CDA][x,0,t]==0,
Derivative[0,1,0][CDA][x,yMax,t]==0,

Derivative[1,0,0][CC][0,y,t]==0,
Derivative[1,0,0][CC][xMax,y,t]==0,
Derivative[0,1,0][CC][x,0,t]==0,
Derivative[0,1,0][CC][x,yMax,t]==0,

Derivative[1,0,0][CD][0,y,t]==0,
Derivative[1,0,0][CD][xMax,y,t]==0,
Derivative[0,1,0][CD][x,0,t]==0,
Derivative[0,1,0][CD][x,yMax,t]==0,

Derivative[1,0,0][CDst][0,y,t]==0,
Derivative[1,0,0][CDst][xMax,y,t]==0,
Derivative[0,1,0][CDst][x,0,t]==0,
Derivative[0,1,0][CDst][x,yMax,t]==0
};

InterpFunction[funcN_,xMaxN_,yMaxN_,tMaxN_]:=Module[{func = funcN, xMax = xMaxN, yMax = yMaxN, tMax = tMaxN},
	intDat=Table[
	{{x,y,t},Evaluate[func[x,y,t][[1]]]}
	,{x,0.,xMax,2.},{y,0.,xMax,2.},{t,0.,tMax,1.}];
	dn=Flatten[Flatten[intDat,1],1];
	fN=Interpolation[dn,InterpolationOrder->2]
]

