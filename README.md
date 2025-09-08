# HPReAn
High-Pressure Resolvent Analysis
HPReAn is designed as a flexible tool for performing resolvent analyses under both ideal-gas and real-fluid thermodynamic conditions. It includes built-in wrappers for the CoolProp and RefProp libraries. The code is thoroughly commented for readability, and the repository contains a complete description, user instructions, and a setup guide.
The resolvent operator can be interrogated via a triplet-sweep (in temporal frequency and wavenumbers) or along a user-specified parameter space. While the code is currently configured for Poiseuille flow, it can be easily adapted to other wall-bounded flow configurations (e.g., Couette flow) by modifying the initial and boundary conditions.
The weighted operator formulation follows the definition presented in this paper, although the code is structured to allow easy customization if needed.
\texttt{HPReAn} depends on two additional packages that require installation beforehand: the High-Pressure Compressible Flow Solver (HPCFS), available at https://github.com/marc-bernades/HPCFS, and the High-Pressure Linear Stability Analysis (HPLSA), found at https://github.com/marc-bernades/HPLSA, which provide supporting functionality and thermodynamic models required by HPReAn.

- Main_Resolvent_HP is the main file to run the analytical resolvent analysis for a Poiseuille flow under high-pressure (or ideal-gas) thermodynamics. The baseflow can be (1) data-driven or (2) laminar/iterative and the triplet (wavenumbers and phase speed) can be swept or fixed at single configuration.
- Main_Postprocess and Main_Resolvent_HP_postprocess are the main scripts that call the functions and utilities to interpret and analyse the results (amplification maps, profiles, baseflows, scale motions, spectrums, etc.
