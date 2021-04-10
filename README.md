# General Relativity Package

A lightweight mathematica package for computing : Christoffel symbols, Riemann Tensor, Ricci tensor, Ricci scalar, Einstein Tensor given an arbitrary metric and coordinate system in any number of dimensions. 

In the ExampleUsage2 notebook you will find code to **compute the geodesic equations of an arbitrary metric** for both and affine and non-affine parameter. Curently this feature is not included in package, but I'll add it as soon as I figure it out (I'm a Mathematica noob).

**Todo** Einstein Equations, raising and lowering indices of tensors for General relativity.


### Quickstart 

1. Clone the repository or download ComputeTensors.wl 
2. Place mathematica notebook in the same folder. 
3. Run 
	```
	SetDirectory[NotebookDirectory[]];
	<< ComputeTensors.wl
	```

5. Define coordinates and a metric, you can do this yourself
	```
	coords = {t,x,y,z};
	metric = DiagonalMatrix[{-1,a[t]^2,a[t]^2,a[t]^2}];
	```

	Or you can import them from a preset
	```
	coords = antiDeSitterCoords;
	metric = antiDeSitterMetric;
	```

6. Compute and display the christoffel symbols and the einstein tensor.
	```
	christoffel = ComputeChristoffelSymbols[metric,coords];
	DisplayChristoffelSymbols[christoffel,coords];
	einsteinTensor = ComputeEinsteinTensor[metric,coords];
	DisplayEinstein[einteinTensor,coords];
	```


