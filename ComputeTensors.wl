(* Author: Stephen Fay
	dcxstephen@live.co.uk or
	stephen.fay@mail.mcgill.ca	*)

BeginPacakge["ComputeTensors`"]

Nonzero::usage="true if input is not zero"
ComputeRiemann::usage="takes [g,x] computes and returns riemann tensor, this is identical to Maloney's implementation";
ComputeRiemannAlt::usage="takes [g,x] computes riemann tensor, alternative to Maloney's implementation, they disagree on the schwarzschild metric"
ComputeRicciScalar::usage="computes the ricci scalar, takes [g,x]"
ComputeEinsteinTensor::usage="ComputeEinsteinTensor[g,x] where g is an nxn matrix and x is a n-vector will return the nxn einstein tensor with lowered indices"
ComputeChristoffel::usage="compute the christoffel symbols a geometry"
ComputeRicciTensor::usage="takes [g,x] computes the ricci tensor"
ComputeRicciTensorAlt::usage="takes [g,x] computes ricci tensor, contracts indices 1 and 3 instead of 2 and 4"
DisplayChristoffelSymbols::usage="takes [christoffel,coords] display nonzero christoffel symbols"
DisplayRiemann::usage="takes [riemannTensor,x] displays nonzero terms of riemann tensor"
DisplayRicciTensor::usage="takes [ricciTensor,x] displays nonzero terms"
DisplayEinstein::usage="takes [einsteinTensor,x] displays nonzero term"
$schwarzschildCoords::usage="polar coordinates are used"
$schwarzschildMetric::usage="the scwarzschild metric in polar coordinates"



(* defines some metrics and their coordinate systems *)
schwarzschildCoords={t,r,\[Theta],\[Phi]};
schwarzschildMetric=DiagonalMatrix[{-(1-R)/r,(1-R/r)^-1,r^2,r^2*Sin[\[Theta]^2]}];

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
		Christoffel = Simplify[Christoffel];
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
        Riemann = Table[ Riemann[[sigma,alpha,mu,nu]] - Riemann[[sigma,alpha,nu,mu]],{sigma,Dim}, {alpha,Dim},{mu,Dim}, {nu,Dim} ]];

(* Maloney's implementation *)
ComputeRiemann[metric_,x_]:=
	Block[{Dim,Chr,Riemann,
		l,i,k,j,m},
		Dim=Length[x];
		Chr = ComputeChristoffel[metric,x];
		(* Return the Riemann Tensor *)
		Simplify[Table[
			D[Chr[[l,i,k]],x[[j]]]
			- D[Chr[[l,j,k]],x[[i]]]
			+ Sum[Chr[[m,i,k]]*Chr[[l,m,j]]
				- Chr[[m,j,k]]*Chr[[l,m,i]],
			{m,Dim}],{i,Dim},{j,Dim},{k,Dim},{l,Dim}]]]

(* in this implementation we contract 1st index with 3rd, in maloney's we contract the 2nd with the 4th *)
ComputeRicciTensorAlt[metric_,x_]:=
	Block[{Dim,Riemann,Ricci,i,j,k},
		Dim=Length[x];
		Riemann=ComputeRiemann[metric,x];
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
	Block[{Dim,RicciTensor,i,j},
		Dim=Length[x];
		RicciTensor=ComputeRicciTensor[metric,x];
		Simplify[Sum[Inverse[metric][[i,j]]*RicciTensor[[i,j]],
			{i,Dim},{j,Dim}]]]

ComputeEinsteinTensor[metric_,x_]:=
	Block[ {Dim, Ricci, RicciScalar,
		i, j},
		Dim = Length[x];
		Ricci = ComputeRicciTensor[metric,x];
        RicciScalar = Sum[ Inverse[metric][[i,j]]
                                 Ricci[[i,j]],
                                 {i,Dim}, {j,Dim}];
        (* Return Einstein tensor: *)
        Simplify[Ricci - (1/2) RicciScalar metric] ]


(* display stuff *)
DisplayChristoffelSymbols[christoffel_, x_] := 
  Block[{i, j, k, dim},
   dim = Length[x];
   For[i = 1, i <= dim, i++,
    For[j = 1, j <= dim, j++,
     For[k = 1, k <= dim, k++,
      If[Nonzero[christoffel[[i]][[j]][[k]]],
       Print[
        StringForm[
         "\!\(\*SubscriptBox[SuperscriptBox[\(\[CapitalGamma]\), \
\(`1`\)], \(`2`\\\ `3`\)]\)=`4`", x[[i]], x[[j]], 
         x[[k]], christoffel[[i]][[j]][[k]]]]]]]]];

DisplayRiemann[riemann_,x_] := 
	Block[{i,j,k,l,dim},
	dim=Length[x];
	   For[i = 1, i <= dim, i++,
		For[j = 1, j <= dim, j++,
		 For[k = 1, k <= dim, k++,
		  For[l = 1, l <= dim, l++, 
		   If[Nonzero[riemann[[i]][[j]][[k]][[l]]],
			Print[
			 StringForm[
          "\!\(\*SuperscriptBox[SubscriptBox[\(R\), \(`1`\\\ `2`\\\ \
`3`\)], \(`4`\)]\)=`5`", coords[[i]], coords[[j]], coords[[k]], 
          coords[[l]], riemann[[i]][[j]][[k]][[l]] ]]]]]]]];

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
	For[i = 1, i <= dim, i++,
		For[j = 1, j <= dim, j++,
		 If[Nonzero[einstein[[i,j]]],
		  Print[
		   StringForm["\!\(\*SubscriptBox[\(G\), \(`1`\\\ `2`\)]\)=`3`", 
			coords[[i]], coords[[j]], einstein[[i,j]]]]]]]];

End[]

EndPackage[];
