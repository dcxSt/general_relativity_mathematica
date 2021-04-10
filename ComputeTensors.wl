(* Author: Stephen Fay
	dcxstephen@live.co.uk or
	stephen.fay@mail.mcgill.ca	*)

BeginPackage["ComputeTensors`"]

(*Nonzero::usage="true if input is not zero" *) (* dont think i need this , delete it at some point *)
ComputeChristoffel::usage="takes [g,x] compute the christoffel symbols a geometry"
ComputeRicciTensor::usage="takes [g,x] computes the ricci tensor"
ComputeRiemann::usage="takes [g,x] computes and returns riemann tensor, this is identical to Maloney's implementation";
ComputeRicciScalar::usage="computes the ricci scalar, takes [g,x]"
ComputeEinsteinTensor::usage="ComputeEinsteinTensor[g,x] where g is an nxn matrix and x is a n-vector will return the nxn einstein tensor with lowered indices"
DisplayMetric::usage="takes [metric, coords] displays nonzeor metric elements"
DisplayChristoffelSymbols::usage="takes [christoffel,coords] display nonzero christoffel symbols"
DisplayRiemann::usage="takes [riemannTensor,x] displays nonzero terms of riemann tensor"
DisplayRicciTensor::usage="takes [ricciTensor,x] displays nonzero terms"
DisplayEinstein::usage="takes [einsteinTensor,x] displays nonzero term"
GeodesicEquations::usage="takes [christoffel, x] displays general parameter geodesic equations"
GeodesicEquationsAffine::usage="takes [christoffel, x] displays affine parameter geodesic equations"
$schwarzschildCoords::usage="spherical coordinates are used"
$schwarzschildMetric::usage="the scwarzschild metric in spherical coordinates"
$antiDeSitterCoords::usage="spherical coordinates"
$antiDeSitterMetric::usage="anti de-Sitter space metric"
$frwFlatSphericalCoords::usage="spherical coords"
$frwFlatSphericalMetric::usage="metric of flat frw space, in spherical coordinates"
$frwFlatCartesianCoords::usage="cartesian coordinates"
$frwFlatCartesianMetric::usage="metric of flat frw space in cartesian coordinates"
$frwReducedCircumferenceSphericalCoords::usage="spherical coordinates"
$frwReducedCircumferenceSphericalMetric::usage="positive or negative curvature frw reduced circumference coordinates, does not cover whole space"
$frwHypersphericalPositiveCurvatureCoords::usage="Spherical coords"
$frwHypersphericalPositiveCurvatureMetric::usage="hyperspherical coordinates for positive curvature (k>0) metric"
$frwHypersphericalNegativeCurvatureCoords::usage="spherical coords"
$frwHypersphericalNegativeCurvatureMetric::usage="hyperspherical coordinates for negative curvature (k<0) metric"
$wormholeCoords::usage="spherical coordinates"
$wormholeMetric::usage="the metric of a simple wormhole connecting two assymtotically minkowski spaces, assuming the function f(r) -> r^2 for as r -> infty"


(* defines some metrics and their coordinate systems *)
antiDeSitterCoords=Block[{t,r,\[Theta],\[Phi]},{t,r,\[Theta],\[Phi]}];
antiDeSitterMetric=Block[{r,l,\[Theta]},DiagonalMatrix[{-(1+r^2/l^2) , (1+r^2/l^2)^(-1) , r^2 , Sin[\[Theta]]^2 * r^2}]];

schwarzschildCoords=Block[{t,r,\[Theta],\[Phi]},{t,r,\[Theta],\[Phi]}];
schwarzschildMetric=Block[{r,\[Theta],R},DiagonalMatrix[{-(1-R/r),(1-R/r)^-1,r^2,r^2*Sin[\[Theta]^2]}]];


frwFlatSphericalCoords=Block[{t,r,\[Theta],\[Phi]},{t,r,\[Theta],\[Phi]}];
frwFlatSphericalMetric=Block[{r,a,t,\[Theta]},{-1,a[t]^2,r^2*a[t]^2,r^2*a[t]^2*Sin[\[Theta]]^2}];
(*special case k=0 of reducedcircumference*)

frwFlatCartesianCoords=Block[{t,x,y,z},{t,x,y,z}];
frwFlatCartesianMetric=Block[{t,a},{-1,a[t]^2,a[t]^2,a[t]^2}];

frwReducedCircumferenceSphericalCoords=Block[{t,r,\[Theta],\[Phi]},{t,r,\[Theta],\[Phi]}];
frwReducedCircumferenceSphericalMetric=Block[{a,r,k,\[Theta]},{-1,a[t]^2,a[t]^2/(1-k*r^2),a[t]^2*r^2,a[t]^2*r^2*Sin[\[Theta]]^2}];

frwHypersphericalPositiveCurvatureCoords=Block[{t,r,\[Theta],\[Phi]},{t,r,\[Theta],\[Phi]}];
frwHypersphericalPositiveCurvatureMetric=Block[{r,a,t,\[Theta],k},DiagonalMatrix[{-1,a[t]^2,a[t]^2*Sin[r*Sqrt[k]]^2/k,a[t]^2*Sin[r*Sqrt[k]]^2/k*Sin[\[Theta]]^2}]];

frwHypersphericalNegativeCurvatureCoords=Block[{t,r,\[Theta],\[Phi]},{t,r,\[Theta],\[Phi]}];
frwHypersphericalNegativeCurvatureMetric=Block[{r,a,t,k,\[Theta]},DiagonalMatrix[{-1,a[t]^2,a[t]^2*Sinh[r*Sqrt[k]]^2/Abs[k],a[t]^2*Sinh[r*Sqrt[k]]^2/Abs[k]*Sin[\[Theta]]^2}]];

wormholeCoords=Block[{t,r,\[Theta],\[Phi]},{t,r,\[Theta],\[Phi]}];
wormholeMetric=Block[{r,\[Theta],f0},DiagonalMatrix[{-1,1,f0[r],Sin[\[Theta]]^2*f0[r]}]];

Begin["`Private`"]

(* helper *)
Nonzero[num_]:=Not[StringMatchQ[ToString[num],"0"] || StringMatchQ[ToString[num],"0."]];

(* compute stuff *)
ComputeChristoffel[metric_,x_]:=
	Block[{Dim, InvMetric, Christoffel,
		alpha, beta, gamma},
		Dim = Length[x];
		InvMetric = Simplify[Inverse[metric]];
		Christoffel = 
			Table[ D[metric[[gamma,alpha]],x[[beta]]]
				+ D[metric[[beta,gamma]],x[[alpha]]]
				- D[metric[[alpha,beta]],x[[gamma]]],
				{gamma,Dim},{alpha,Dim},{beta,Dim}];
			(* the lower index part of christoffel symbol*)
		Christoffel = (1/2) InvMetric . Christoffel;
		Simplify[Christoffel]]

(* Maloney's implementation *)
ComputeRiemann[metric_,x_]:=
	Block[{Dim,Chr,Riemann, sigma,rho,alpha,beta,gamma},
		Dim=Length[x];
		Chr = ComputeChristoffel[metric,x];
		(* Return the Riemann Tensor *)
		Simplify[Table[
			D[Chr[[sigma,alpha,gamma]],x[[beta]]]
			- D[Chr[[sigma,beta,gamma]],x[[alpha]]]
			+ Sum[
				Chr[[rho,alpha,gamma]]*Chr[[sigma,rho,beta]]
				- Chr[[rho,beta,gamma]]*Chr[[sigma,rho,alpha]],
			{rho,Dim}],{alpha,Dim},{beta,Dim},{gamma,Dim},{sigma,Dim}]]]

(* contract indices 2 and 4 *)
ComputeRicciTensor[metric_,x_]:=
	Block[{Dim,Riemann,Ricci,i,j,k},
		Dim=Length[x];
		Riemann=ComputeRiemann[metric,x];
		(* Return Ricci Tensor *)
		Simplify[Table[Sum[Riemann[[i,k,j,k]],{k,Dim}],
			{i,Dim},{j,Dim}]]]

(* Contracts the ricci tensor with the inverse metric *)
ComputeRicciScalar[metric_,x_]:=
	Block[{Dim,RicciTensor,alpha,beta},
		Dim=Length[x];
		RicciTensor=ComputeRicciTensor[metric,x];
		Simplify[Sum[Inverse[metric][[alpha,beta]]*RicciTensor[[alpha,beta]],
			{alpha,Dim},{beta,Dim}]]]

ComputeEinsteinTensor[metric_,x_]:=
	Block[ {Dim, Ricci, RicciScalar, alpha, beta},
		Dim = Length[x];
		Ricci = ComputeRicciTensor[metric,x];
        RicciScalar = Sum[ Inverse[metric][[alpha,beta]]
                                 Ricci[[alpha,beta]],
                                 {alpha,Dim}, {beta,Dim}];
        (* Return Einstein tensor: *)
        Simplify[Ricci - (1/2) RicciScalar metric] ]

(* display stuff *)
DisplayMetric[metric_,x_] := 
	Block[{i,j,dim},
		dim = Length[x];
		For[i=1, i<=dim, i++,
			For[j=1, j<=dim, j++,
				If[Nonzero[metric[[i,j]]],
					Print[StringForm[
						"\!\(\*SubscriptBox[\(g\), \(`1`\\\ `2`\)]\)=`3`", x[[i]], x[[j]], metric[[i,j]]]]]]]]

DisplayChristoffelSymbols[christoffel_, x_] := 
  Block[{alpha, beta, gamma, dim},
   dim = Length[x];
   For[alpha = 1, alpha <= dim, alpha++,
    For[beta = 1, beta <= dim, beta++,
     For[gamma = 1, gamma <= dim, gamma++,
      If[Nonzero[christoffel[[alpha,beta,gamma]]],
       Print[
        StringForm[
         "\!\(\*SubscriptBox[SuperscriptBox[\(\[CapitalGamma]\), \
\(`1`\)], \(`2`\\\ `3`\)]\)=`4`", x[[alpha]], x[[beta]], 
         x[[gamma]], christoffel[[alpha,beta,gamma]]]]]]]]];

DisplayRiemann[riemann_,x_] := 
	Block[{alpha,beta,gamma,sigma,dim},
	dim=Length[x];
	   For[alpha = 1, alpha <= dim, alpha++,
		For[beta = 1, beta <= dim, beta++,
		 For[gamma = 1, gamma <= dim, gamma++,
		  For[sigma = 1, sigma <= dim, sigma++, 
		   If[Nonzero[riemann[[alpha,beta,gamma,sigma]]],
			Print[
			 StringForm[
          "\!\(\*SuperscriptBox[SubscriptBox[\(R\), \(`1`\\\ `2`\\\ \
`3`\)], \(`4`\)]\)=`5`", x[[alpha]], x[[beta]], x[[gamma]], 
          x[[sigma]], riemann[[alpha,beta,gamma,sigma]] ]]]]]]]];

DisplayRicciTensor[ricci_,x_] := 
	Block[{Dim,i,j}, 
		Dim=Length[x];
		For[i = 1, i <= Dim, i++,
			For[j = 1, j <= Dim, j++,
			 If[Nonzero[ricci[[i,j]]], 
			  Print[StringForm[
				"\!\(\*SubscriptBox[\(R\), \(`1`\\\ `2`\)]\)=`3`", 
				x[[i]], x[[j]], ricci[[i,j]]]]]]]];

DisplayEinstein[einstein_,x_] := 
	Block[{alpha,beta,Dim}, 
	Dim=Length[x];
	For[alpha = 1, alpha <= Dim, alpha++,
		For[beta = 1, beta <= Dim, beta++,
		 If[Nonzero[einstein[[alpha,beta]]],
		  Print[
		   StringForm["\!\(\*SubscriptBox[\(G\), \(`1`\\\ `2`\)]\)=`3`", 
			x[[alpha]], x[[beta]], einstein[[alpha,beta]]]]]]]];

(* takes christoffel symbols and metric, returns geodesic equations*)
(* x should have componants written as functions of s, e.g. x[t[s],y[s]], s is the worldline parameter *)
GeodesicEquations[Chr_,x_] := 
	Block[{alpha,beta,gamma,Dim},
		Dim=Length[x];
		Print[D[x[[2]],s]];
		Print[x]
		Simplify[
			Table[
				\[Lambda]*D[x[[gamma]],s] == 
				D[D[x[[gamma]],s],s] + 
				Sum[Chr[[gamma,alpha,beta]]*D[x[[alpha]],s]*D[x[[beta]],s],
					{alpha,Dim},{beta,Dim}],{gamma,Dim}]]]

GeodesicEquationsAffine[Chr_,x_] := 
	Block[{alpha,beta,gamma,Dim},
		Dim=Length[x];
		Simplify[
			Table[
				0 == D[D[x[[gamma]],s],s] + 
				Sum[Chr[[gamma,alpha,beta]]*D[x[[alpha]],s]*D[x[[beta]],s],
					{alpha,Dim},{beta,Dim}],{gamma,Dim}]]]

End[]

EndPackage[];
