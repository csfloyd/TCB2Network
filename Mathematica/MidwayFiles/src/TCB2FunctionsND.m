(* ::Package:: *)

gaussian[r_,r0_,s_]:=Exp[-(1/2) ((r-r0)/s)^2];
spaceFunc[r_,r0_,s_]:=1/2*(1-Tanh[(r-r0)/s]);


sigmoid[x_]:=1/(1+Exp[-x]);
smoothBump[x_,a_,xa_,xb_]:=sigmoid[(x-xa)/a]-sigmoid[(x-xb)/a]
waveFunc[t_,cyc_,len_,del_]:=Sum[smoothBump[t,del,cyc*i,cyc*i+len],{i,0,Floor[tMax/cyc]}] * smoothBump[t,5*del,0,release];
trange = Range[0,tMax,0.05];
wPoints = waveFunc[trange,cyc,len,del];
data = Transpose[{trange,wPoints}];
iFun = Interpolation[data, InterpolationOrder -> 2];
gamFunc[r_,t_]:=spaceFunc[r,r0,width]*iFun[t];
g[r_,t_]:=-(1-gMin)*gamFunc[r,t];
(*g[r_,t_]:=gMin*gamFunc[r,t];*)
gradr[r_,t_]:=Evaluate[D[g[rp,t],rp]/.rp->r];

f[r_,t_,flat_,aOff_,n_,v_]:=(v*iFun[t]+(1-v))*(1-flat)*spaceFunc[r,r0,width]*(aOff+(1-aOff)*(r/r0)^n) +flat;
CB[r_,t_]:=f[r,t,flat,aOff,n,v];
(*CB[r_,t_]:=CB0;*)


forceC[r_,t_,u_]:=2 \[Mu]0*(D[CB[r,t],r]*(D[u,r]-1/2 g[r,t])+CB[r,t]*(D[u,{r,2}]+1/r D[u,r]-u/r^2-1/2 gradr[r,t]))+\[Lambda]0*(D[CB[r,t],r]*(u/r+D[u,r]-g[r,t])+CB[r,t]*(D[u,{r,2}]+1/r D[u,r]-u/r^2-gradr[r,t]));

forceL[r_,t_,u_]:=(\[Lambda]0 + 2 \[Mu]0)*(D[CB[r,t],r]*(D[u,r]- g[r,t])+CB[r,t]*(D[u,{r,2}]-gradr[r,t]));



pdesC={
D[u[r,t],t] == forceC[r,t,u[r,t]]
};

pdesL={
D[u[r,t],t] == forceL[r,t,u[r,t]]
};


ics={
u[r,0] == 0
};


bcs={

u[rMin,t] == 0,
u[rMax, t] == 0,

Derivative[1,0][u][rMin,t] == 0,
Derivative[1,0][u][rMax,t] == 0

};


InterpFunction[funcN_,rMinN_,rMaxN_,tMaxN_]:=Module[{func = funcN, rMin = rMinN, rMax = rMaxN, tMax = tMaxN},
	intDat=Table[
	{{r,t},Evaluate[func[r,t][[1]]]}
	,{r,rMin,rMax,0.5},{t,0.,tMax,0.005}];
	dn=Flatten[intDat,1];
	fN=Interpolation[dn,InterpolationOrder->2]
]
