Outline
==============

Rabacus is a `Python <http://www.python.org>`_  package for performing 
analytic radiative transfer calculations in simple geometries relevant to 
cosmology and astrophysics.   It also contains tools to calculate cosmological
quantities such as the power spectrum and mass function.  
The core of the package is a series of single zone solvers
that can be used to find equilibrium ionization fractions and/or temperatures 
as a function of hydrogen and helium density.  
These single zone solvers are used in the geometric solvers to compute 
solutions in a number of simple geometries.   
In this section we will outline some of its unique features and review some
of the relevant basic physics. 
In the following sections we provide usage examples. 


Features
-----------------

  * **Units**: All physical quantites in Rabacus have a magnitude and
    a unit through use of the `Quantities
    <https://github.com/python-quantities/python-quantities>`_ python package.
  * **Helium**: All solvers in Rabacus take accout of both hydrogen and helium.
  * **Accelerated**: The core functionality of Rabacus is available in Python 
    wrapped Fortran 90 which is thousands of times faster than a pure Python
    implementation. 
  * **Recombination Radiation**: Rabacus can solve for the transfer of diffuse
    recombination radiation. 





Case A Fraction
-----------------------

Atomic rates play a central role in Rabacus.  When creating a chemistry 
or cooling rates object, the user must specify if case A or case B 
recombination rates should be used for each recombing species 
{``HII``, ``HeII``, ``HeIII``}.

.. note:: 

   The recombination rate to all atomic levels is called the case A rate.  
   Recombinations to the ground state can cause ionizations and so a case B 
   rate is defined that negelects these recombinations, 

   .. math::
      \alpha_{\rm B} = \alpha_{\rm A} - \alpha_{\rm 1}

To allow for a continuous transition between case A and case B, we define a 
variable called the case A fraction, :math:`f_{\rm A}`, for each recombining 
ionic species {``HII``, ``HeII``, ``HeIII``}
such that the returned recombination rate is, 

.. math:: 

   \log \alpha = \log \alpha_{\rm A} f_{\rm A} + 
   \log \alpha_{\rm B} (1-f_{\rm A})

