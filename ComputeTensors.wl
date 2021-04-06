(* Author: Stephen Fay
	dcxstephen@live.co.uk or
	stephen.fay@mail.mcgill.ca	*)

BeginPackage["ComputeTensors`"]

(*Nonzero::usage="true if input is not zero" *) (* dont think i need this , delete it at some point *)
ComputeChristoffel::usage="takes [g,x] compute the christoffel symbols a geometry"
ComputeRicciTensor::usage="takes [g,x] computes the ricci tensor"
ComputeRicciTensorAlt::usage="takes [g,x] computes ricci tensor, contracts indices 1 and 3 instead of 2 and 4"
ComputeRiemann::usage="takes [g,x] computes and returns riemann tensor, this is identical to Maloney's implementation";
ComputeRiemannAlt::usage="takes [g,x] computes riemann tensor, alternative to Maloney's implementation, they disagree on the schwarzschild metric"
ComputeRicciScalar::usage="computes the ricci scalar, takes [g,x]"
ComputeEinsteinTensor::usage="ComputeEinsteinTensor[g,x] where g is an nxn matrix and x is a n-vector will return the nxn einstein tensor with lowered indices"
ComputeEinsteinTensorAlt::usage="alternative computation of einstein tensor"
DisplayMetric::usage="takes [metric, coords] displays nonzeor metric elements"
DisplayChristoffelSymbols::usage="takes [christoffel,coords] display nonzero christoffel symbols"
DisplayRiemann::usage="takes [riemannTensor,x] displays nonzero terms of riemann tensor"
DisplayRicciTensor::usage="takes [ricciTensor,x] displays nonzero terms"
DisplayEinstein::usage="takes [einsteinTensor,x] displays nonzero term"
GeodesicEquations::usage="takes [christoffel, x] displays general parameter geodesic equations"
GeodesicEquationsAffine::usage="takes [christoffel, x] displays affine parameter geodesic equations"
$schwarzschildCoords::usage="polar coordinates are used"
$schwarzschildMetric::usage="the scwarzschild metric in polar coordinates"
$antiDeSitterCoords::usage="polar coordinates"
$antiDeSitterMetric::usage="anti de-Sitter space metric"

(* defines some metrics and their coordinate systems *)
antiDeSitterCoords=Block[{t,r,\[Theta],\[Phi]},{t,r,\[Theta],\[Phi]}];
antiDeSitterMetric=Block[{r,l,\[Theta]},DiagonalMatrix[{-(1+r^2/l^2) , (1+r^2/l^2)^(-1) , r^2 , Sin[\[Theta]]^2 * r^2}]];

schwarzschildCoords=Block[{t,r,\[Theta],\[Phi]},{t,r,\[Theta],\[Phi]}];
schwarzschildMetric=Block[{r,\[Theta],R},DiagonalMatrix[{-(1-R/r),(1-R/r)^-1,r^2,r^2*Sin[\[Theta]^2]}]];

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


(* Alternative implementation, this yeilds different results for the schwarzschild metric than Maloney's implementation *)
ComputeRiemannAlt[metric_,x_]:=
	Block[{Dim,Christoffel,Riemann,
		sigma, mu, nu, alpha, beta, gamma},
		Dim = Length[x];
		Christoffel = ComputeChristoffel[metric,x];
		Riemann = Table[ D[Christoffel[[sigma,alpha,nu]],x[[mu]]]
                    + Sum[Christoffel[[gamma,alpha,nu]] Christoffel[[sigma,gamma,mu]],{gamma,Dim} ],{sigma,Dim}, {alpha,Dim}, {mu,Dim}, {nu,Dim} ];
           	(* antisymmetrize Riemann tensor: *)
        Riemann = Table[ Riemann[[sigma,alpha,mu,nu]] - Riemann[[sigma,alpha,nu,mu]],{sigma,Dim}, {alpha,Dim},{mu,Dim}, {nu,Dim} ];
		Riemann];

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

(* in this implementation we contract 1st index with 3rd, in maloney's we contract the 2nd with the 4th *)
ComputeRicciTensorAlt[metric_,x_]:=
	Block[{Dim,Riemann,Ricci,i,j,k},
		Dim=Length[x];
		Riemann=ComputeRiemannAlt[metric,x];
		(* Return Ricci Tensor *)
		Simplify[Table[Sum[Riemann[[i,j,i,k]],{i,Dim}],
			{j,Dim},{k,Dim}]]]

(* Maloney's implementation, contract indices 2 and 4 *)
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

ComputeEinsteinTensorAlt[metric_,x_]:=
	Block[ {Dim,Ricci, RicciScalar, i, j},
		Dim = Length[x];
		Ricci = ComputeRicciTensorAlt[metric,x];
		RicciScalar = Sum[ Inverse[metric][[i,j]], {i,Dim},{j,Dim}];
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
`3`\)], \(`4`\)]\)=`5`", coords[[alpha]], coords[[beta]], coords[[gamma]], 
          coords[[sigma]], riemann[[alpha,beta,gamma,sigma]] ]]]]]]]];

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
	Block[{i,j,Dim}, 
	Dim=Length[x];
	For[i = 1, i <= Dim, i++,
		For[j = 1, j <= Dim, j++,
		 If[Nonzero[einstein[[i,j]]],
		  Print[
		   StringForm["\!\(\*SubscriptBox[\(G\), \(`1`\\\ `2`\)]\)=`3`", 
			coords[[i]], coords[[j]], einstein[[i,j]]]]]]]];

(* takes christoffel symbols and metric, returns geodesic equations*)
(* x should have componants written as functions of s, e.g. x[t[s],x[s]], s is the worldline parameter *)
GeodesicEquations[Chr_,x_] := 
	Block[{alpha,beta,gamma,Dim},
		Dim=Length[x];
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
