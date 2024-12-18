
##Principle


The ISO 25178-2:2012(E) states that

*“The observed area is calculated as a function of scale by a series of
virtual tiling exercises covering the measured surface in a patchwork fashion. The areas of the tiles, or patches,
represent the areal scales of observation. The tiling exercises are repeated with tiles of progressively smaller areas to
determine the observed areas to determine the observed areas as a function of the areal scales of observation.*

Then, the function log(relative area)=f(log(element area)) can be determined.
The area-scale fractal analysis complexity parameter *Asfc* is defined as a thousand times minus the line slope:.

<p style="text-align:center;"><img src="../media/Asfc.png" alt="25178-2:2012(E) Asfc" width="500px"/></p>

##Determining the steepest part of the curve in a robust way
It can be reasonably stated that most of relative area plots are S-shaped, even if the trail can be more linear
than curved. Instead of arbitrarily defining a “central” region where to calculate the maximum of the slope, the whole
curve is fitted with a family *Tn*, of monotonic functions.
\[
   Tn(x) = y_0 +a \cdot \tanh \left( \dfrac{x -x_0}{b} \right)^n
\]

In most cases \(n=1\) leads to an accurate fit, especially in the linear part of the curve, where the studied surface exhibits
self-similarity properties.

The proposed model is parsimonious — only four parameters have to be determined, with a
least square procedure for instance — and then doesn’t suffer from overfitting. Some surfaces have better fit results with
\(n=2\) but it doesn’t change much the *Asfc* value. The location \(x_s\) of the steepest part is analytically determined by
canceling the second derivative of *Tn*. A simple calculus yields:

\[
   x_s = x_0 +b \cdot \text{atanh} \left( \sqrt{\dfrac{x -x_0}{b}} \right)
\]

when \(n=1\), it reduces to \(x_s = x_0\).

To address the computation efficiency question — induced by the series of tiling exercises — the surface itself
isn’t tiled as explained in ISO 25178. The different scales are those of the grids that are used to discretize the surface.

The grids are chosen regular with lateral steps from *(hx, hy)* for the finest grid to *(Hx, Hy)* for the largest grid. The finest
grid is the one of the original surface, and the coarsest grid is an 8 × 8 grid.

The original grid is the only one for which the surface points match the grid points, for the other grids the surface heights are obtained by interpolation.

The intermediate grids are defined so that the element size has a geometric progression, hence its location is evenly spaced on a logarithm axis.

Three surface interpolation methods are tested: linear (linear Finite Elements), interpolant cubic splines, and Hermite
(Hermite Finite Elements)

<p style="text-align:center;"><img src="../media/Asfc_interp.png" alt="Asfc interpolations" width="500px"/></p>

The chosen interpolation methods are: left, linear FE - middle, Hermite FE -  right, cubic spline. “dof” stands for “degree of
freedom”. In the upper part, one can see that from left to right, the smoothness increases. Below, calling *u* the unknown, it can be
seen that increasing smoothness requires to take into account more derivatives. As the spline method guarantees curvature continuity,
the unknown and its derivatives must be computed globally.

The surface heights are computed on the different grids with one of the three interpolation procedures that are
proposed – linear (finite element style), Hermite (finite element style) and cubic spline. 128 grids are used to draw the
relative area curve.
Provided the interpolating function *f(x,y)*, the area of a surface element defined on \([x_i,x_i+hx]\times[y_j,y_j+hy]\) is obtained
by integration of

\[
   \sqrt{ 1 + \left(\dfrac{\partial f}{\partial x}\right)^2 + \left(\dfrac{\partial f}{\partial y}\right)^2 }
\]

Such an expression being barely easy to analytically integrate, the Gauss integrating method is used with two points in each direction,
which ensures the exact integration of a third degree polynomial on both directions.








