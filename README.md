# BoundaryScattering
Codebase for the manuscript: _Thermal Resistance at Twist Boundaries and Heterointerfaces_
## Scattering Scripts
#### ArrayScattering Class
General initialization of an interfacial dislocation array. Includes methods to calculate grain boundary energy and the scattering rate from dislocation strain.
#### AMMTransport Class
Defines an generic acoustic mismatch at an interface. Includes methods to describe the phonon transmissivity and thermal boundary conductance limited by an interfacial acoustic mismatch. 
##### HetAMMTransport Class
Child class of AMMTransport for the acoustic mismatch at a heterointerface, where the acoustic mismatch stems from a change in the stiffness matrix and density of the material.
### Specific Grain Boundary Types
__TiltScattering__ : Instantiates ArrayScattering object with tilt character. Defines strain fields for 1D array of edge dislocations. Instantiates AMMTransport object with acoustic mismatch coming from a rotation of the Christoffel matrix with the axis of rotation in the direction of dislocation spacing.  
__TwistScattering__ : Instantiates ArrayScattering object with twist character. Defines strain fields for 2D grid of screw dislocations. Instantiates AMMTransport object with acoustic mismatch coming from a rotation of the Christoffel matrix with the axis of rotation normal to the boundary plane.  
__HetIntScattering__ : Instantiates ArrayScattering object with semicoherent heterointerface character. Defines strain fields for 2D grid of edge-type misfit dislocations. Instantiates HetAMMTransport object with acoustic mismatch coming from a step-function change in the stiffness matrix and density at the interface.  

### Additional Support Scripts
__ThermalTransport__ : Defines methods for calculating transport properties including the spectral relaxation time, phonon transmissivity, thermal boundary conductance, and lattice thermal conductivity (for a specified grain size)  
__ScatteringPlots__ : Defines plotting functions for phonon diffraction plots, convergence plots, and spectral relaxation time and transport properties  
__helper__ : helper functions for generic vector operations and interfacing with pymatgen  
