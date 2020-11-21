# BoundaryScattering
Codebase for the manuscript: _Thermal Resistance at Twist Boundaries and Heterointerfaces_
## Scattering Scripts
### ArrayScattering Class
General initialization of an interfacial dislocation array. Includes methods to calculate grain boundary energy and the scattering rate from dislocation strain.
### AMMTranposrt Class
Defines an generic acoustic mismatch at an interface. Includes methods to describe the phonon transmissivity and thermal boundary conductance limited by an interfacial acoustic mismatch.
### HetAMMTransport Class
Child class of AMMTransport for the acoustic mismatch at a heteorinterface, where the acoustic mismatch stems from a change in the stiffness matrix and density of the material.
