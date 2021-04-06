# General Relativity Package

Lightweight mathematica package for computing : Christoffel symbols, Riemann Tensor, Ricci tensor, Ricci scalar, Einstein Tensor, Einstein Equations, raising and lowering indices of tensors for General relativity in any number of dimensions.

### Quickstart 

1. Clone the repository or download ComputeTensors.wl 
2. Place mathematica notebook in the same folder. 
3. Run 

	SetDirectory[NotebookDirectory[]];

4. Run 

	<< ComputeTensors.wl

5. Define a coordinates and a metric

	coords={t,x,y,z};
	metric=DiagonalMatrix[{-1,a[t]^2,a[t]^2,a[t]^2}];

6. Display the einstein tensor.

	einsteinTensor=Simplify[ComputeEinsteinTensor[metric,coords]]

### Functions


