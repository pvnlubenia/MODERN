==========================================================================================

  MODERN: Mass action Ordinary Differential Equation system to chemical Reaction Network

==========================================================================================

Matlab was used to develop the function used here.
The function ode_to_crn returns the chemical reaction network corresponding to a system of ordinary differential equations (ODEs) with mass action kinetics. If the system does not satisfy the Hars-Toth criterion, a message appears stating what needs to be added in a flux. Furthermore, the output variable 'ode' allows the user to view the complete system of ODEs with all the species and fluxes listed in the 'species' and 'flux' fields, respectively.



=================================
How to fill out 'ode' structure
=================================

'ode' is the input for the function ode_to_crn. It is a structure, representing the system of ODEs, with the following fields:

   - id: name of the model
   - species: a list of all species in the system; this is left blank since incorporated into the function is a step which compiles all species used in the ODEs
   - flux: a list of all fluxes in the system; this is left blank since incorporated into the function is a step which compiles all fluxes used in the ODEs
   - reaction: a list of all reactions/fluxes in the system, each with the following subfields:
        - id: has the following further subfields:
             - var: a string representing the flux number of the kinetics
             - eq: a string representing the kinetic equation
        - species: a list of strings representing the species in the flux
        - order: a list of numbers representing the kinetic order of each species in the flux in the left to right direction (listed  in the same order of the species)
   - equation: a list of all ODEs in the system, each with the following subfields:
        - id: has the following further subfields:
             - var: a string representing the dependent variable (state variable) of the ODE
             - eq: a string representing the ODE flux balance equation
        - flux: a list of fluxes representing the fluxes in the ODE
        - coeff: a list of numbers representing the coefficient of each flux in the ODE in the left to right direction (listed  in the same order of the fluxes)

To fill out the 'ode' structure, write a string for 'ode.id': this is just to put a name to the set of ODEs. To add fluxes and equations to the network, use the functions addFlux and addEquation, respectively, where both have the output 'ode'.

   addFlux
      - OUTPUT: Returns a structure called 'ode' with added field 'reaction' with subfields 'id' (with further subfields 'var' and 'flux'), 'species', and 'kinetic' (see README.txt for details). The output variable 'ode' allows the user to view the ordinary differential equations with the added reaction.
      - INPUTS
           - ode: a structure, representing the ordinary differential equations
           - var: variable representing the flux (string)
           - flux: a visual representation of the flux (string)
           - species: list of species involved in the reaction (cell within a cell)
           - kinetic: list of kinetic orders (exponents) of the species (cell within a cell)

   addEquation
      - OUTPUT: Returns a structure called 'ode' with added field 'equation' with subfields 'id' (with further subfields 'var' and 'flux'), 'flux', and 'coeff' (see README.txt for details). The output variable 'ode' allows the user to view the ordinary differential equations with the added equation.
      - INPUTS
           - ode: a structure, representing the ordinary differential equations
           - var: variable representing the state variable (string)
           - eq: a visual representation of the flux balance expression (string)
           - flux: list of fluxes involved in the ODE (cell within a cell)
           - coeff: list of coefficients (stoichiometry) of the fluxes (cell within a cell)

   * Make sure the functions addFlux and addEquation are in the same folder/path being used as the current working directory.

Note that it is assumed that the ODE system is a MASS ACTION SYSTEM.



========
Examples
========

4 examples are included in this folder, all of which came from [1].



===================
Contact Information
===================

For questions, comments, and suggestions, feel free to contact me at pvnlubenia@yahoo.co.uk.


- Patrick Lubenia (18 July 2022)



=========
Reference
=========

   [1] Chellaboina V, Bhat S, Haddad W, Bernstein D (2009) Modeling and analysis of mass-action kinetics. IEEE Control Syst 29(4):60-78. https://doi.org/10.1109/MCS.2009.932926