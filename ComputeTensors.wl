(* Author: Stephen Fay
	dcxstephen@live.co.uk or
	stephen.fay@mail.mcgill.ca	*)

BeginPacakge["ComputeTensors`"]

ComputeEinsteinTensor::usage="ComputeEinsteinTensor[g,x] where g is an nxn matrix and x is a n-vector will return the nxn einstein tensor with lowered indices"

Begin["`Private`"]

ComputeEinsteinTensor[metric_,x_]:=
  Block[ {Dim,InvMetric, Christoffel, Riemann,
          Ricci, RicciScalar,
          sigma, mu, nu, alpha, beta, gamma},
          Dim = Length[x];
          InvMetric = Simplify[Inverse[metric]];
          Christoffel =
            Table[ D[metric[[gamma,alpha]],x[[beta]]]
                 + D[metric[[beta,gamma]],x[[alpha]]]
           		    - D[metric[[alpha,beta]],x[[gamma]]],
           	   	 {gamma,Dim}, {alpha,Dim}, {beta,Dim} ];
           	 (* The lower index part of Christoffel symbols *)
          Christoffel = Simplify[Christoffel];
          Christoffel = (1/2) InvMetric . Christoffel;
             (* The Christoffel symbols *)
          Christoffel = Simplify[Christoffel];
          Riemann = 
             Table[ D[Christoffel[[sigma,alpha,nu]],x[[mu]]]
                    + Sum[Christoffel[[gamma,alpha,nu]]
                            Christoffel[[sigma,gamma,mu]],
                          {gamma,Dim} ],
                    {sigma,Dim}, {alpha,Dim}, {mu,Dim}, {nu,Dim} ];
           	(* antisymmetrize Riemann tensor: *)
          Riemann = Table[ Riemann[[sigma,alpha,mu,nu]]
                         - Riemann[[sigma,alpha,nu,mu]],
                           {sigma,Dim}, {alpha,Dim},
                           {mu,Dim}, {nu,Dim} ];
          Ricci = Table[ Sum[Riemann[[sigma,alpha,sigma,beta]],
                             {sigma,Dim}],
                         {alpha,Dim}, {beta,Dim} ];
          RicciScalar = Sum[ InvMetric[[alpha,beta]]
                                 Ricci[[alpha,beta]],
                                 {alpha,Dim}, {beta,Dim} ];
          (* Return Einstein tensor: *)
          Ricci - (1/2) RicciScalar metric ]


End[]


EndPackage[];
