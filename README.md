# Heat-Equation-2D-Douglas-Rachford
Numerical solution for 2D heat equation in a rectangular domain using Douglas‐Rachford splitting method

$$\frac{\partial u}{\partial t} = \Delta u + f(x, y), \quad (x, y) \in \Omega$$

with initial condition

$$u|_{t = 0} = \varphi(x, y)$$

and Dirichlet boundary conditions

$$u = \mu(x, y, t), \quad (x, y) \in \partial \Omega$$

## Structure
1. **main.cpp** - main file
2. **Problems.cpp** - implementation of functions $f$, $\varphi$ and $\mu$
3. **Results.txt** - results of a numerical solution with the parameters used for this solution

## Compile
```
g++ main.cpp Problems.cpp -o main.exe
```
