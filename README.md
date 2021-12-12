# Fiber Orientation Tools

By Charles L. Tucker III, Department of Mechanical Science and Engineering, University of Illinois at Urbana-Champaign

This toolkit contains Matlab functions for modeling flow-induced fiber orientation and mechanical properties for discontinuous-fiber composite materials.  It was developed as a supplement to the book **_Fundamentals of Fiber Orientation: Description, Measurement and Prediction_** by C. L. Tucker, Hanser, Munich, 2022.  Errata for the book will also be posted here.

This software is distributed under an MIT license; see the LICENSE file for details.  It can be cited as:
- C. L. Tucker, Fiber Orientation Tools, https://https://github.com/charlestucker3/Fiber-Orientation-Tools, v1.0.0, 2021.

Included are functions to calculate:
- Flow-induced fiber orientation using a variety of tensor-based models: Jeffrey, Folgar-Tucker or anisotropic rotary diffusion (Phelps-Tucker 5-constant, iARD, pARD, MRD and Wang/2-constant models), and slow kinetics (SRF, RSC and RPR).
- A range of closure approximations for the fourth-order orientation tensor: linear, quadratic, hybrid, orthotropic, natural and IBOF.
- The full orientation distribution function under flow, either planar or 3-D orientation (Jeffery, Folgar-Tucker and ARD).
- Evolution of the fiber length distribution using the Phelps-Tucker model.
- Linear elastic stiffness via mean-field models (Halpin-Tsai, dilute Eshelby, Mori-Tanaka and Lielens/double inclusion) for either unidirectional or distributed fiber orientation.
- Thermal stress and thermal expansion tensors for any stiffness model.
- The stiffness of a composite plate whose orientation varies across its thickness, using classical laminated plate theory.  

A number of smaller utility functions are also included.  

To use the tools, click the green **Code** button, and select **Download ZIP**.  Unzip the file and follow the installation suggestions in the _Documentation_ folder.  
