# Fiber Orientation Tools

By Charles L. Tucker III, Department of Mechanical Science and Engineering, University of Illinois at Urbana-Champaign

This toolkit contains Matlab functions for modeling flow-induced fiber orientation and mechanical properties for discontinuous-fiber composite materials.  It was developed as a supplement to the book **_Fundamentals of Fiber Orientation: Description, Measurement and Prediction_** by C. L. Tucker, Hanser, Munich, 2022.  Errata for the book is available at https://github.com/charlestucker3/Fundamentals-of-Fiber-Orientation-errata.

Also included are functions referenced in the paper "Planar Fiber Orientation: Jeffery, Non-Orthotropic Closures and Reconstructing Distribution Functions," submitted to the _Journal of Non-Newtonian Fluid Mechanics_, June, 2022.  

This software is distributed under an MIT license; see the LICENSE file for details.  It can be cited as:
- C. L. Tucker, Fiber Orientation Tools, https://github.com/charlestucker3/Fiber-Orientation-Tools, v1.0.0, 2021.

Included are functions to calculate:
- Flow-induced fiber orientation using a variety of tensor-based models: Jeffrey, Folgar-Tucker or anisotropic rotary diffusion (Phelps-Tucker 5-constant, iARD, pARD, MRD and Wang/2-constant models), and slow kinetics (SRF, RSC and RPR).
- A range of closure approximations for the fourth-order orientation tensor: linear, quadratic, hybrid, orthotropic, natural and IBOF in 3-D; Bingham, natural, ellipse-radius, orthotropic D and non-orthotropic F for planar orientation.
- The full orientation distribution function, transient or steady state, for either planar or 3-D orientation (Jeffery, Folgar-Tucker and ARD).
- Evolution of the fiber length distribution using the Phelps-Tucker model.
- Linear elastic stiffness via mean-field models (Halpin-Tsai, dilute Eshelby, Mori-Tanaka and Lielens/double inclusion) for either unidirectional or distributed fiber orientation.
- Thermal stress and thermal expansion tensors for any stiffness model.
- The stiffness of a composite plate whose orientation varies across its thickness, using classical laminated plate theory. 
- Reconstruction of planar orientation distributions using the Jeffery, Bingham, ellipse radius (ER), and fourth-order maximum entropy functions.   

A number of smaller utility functions are also included.  Example scripts show how the functions are used in many common applications, and reproduce some of the figures and calculations from the book.   

To use the tools, click the green **Code** button, and select **Download ZIP**.  Unzip the file and follow the installation suggestions in the _Documentation_ folder.  

Python users may want to explore https://github.com/nilsmeyerkit/fiberoripy and several packages by https://github.com/JulianKarlBauer.
