Fabber models for ASL-MRI
=========================

These models use the Fabber_
Bayesian model fitting framework [1]_ to implement a number of models
for Arterial Spin Labelling MRI (ASL-MRI).

.. note::
    If you have ASL data that you are looking to process you should
    start with the BASIL_ toolset or the OXASL_ both of 
    which use ``FABBER_ASL`` internally.

Getting FABBER_ASL
------------------

The ASL models are included in FSL_. We
stongly recommend version 6.0.1 or later.

If you need an updated version of the model which has not yet been released to
FSL, you will either need to 
`build from source <https://fabber-core.readthedocs.io/en/latest/building.html#building-new-or-updated-model-libraries>`_ 
using an existing FSL 6.0.1 or later installation, or download 
the pre-built `Fabber bundle <https://fabber-core.readthedocs.io/en/latest/getting.html#standalone-fabber-distribution>`_ 
which contains the latest ASL release alongside other models in a standalone package.

Models included
---------------

The resting state ASL model
~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is the most common ASL model. In it simplest form it implements the
well-known Buxton model, however it can incorporate additional features
such as an arterial component, exchange and dispersion models and partial
volume correction. It is the main model used in the OXASL_ and BASIL_
pipelines.

This model is selected using ``--model=aslrest``.
                
The multiphase model
~~~~~~~~~~~~~~~~~~~~

This model is designed for processing multiphase ASL data. It performs the equivalent
of label-control subtraction, resulting in a data set which is suitable for processing
using the ``aslrest`` model. It is used by the multiphase plugin for the OXASL_ pipeline.

This model is selected using ``--model=asl_multiphase``.

The QUASAR model
~~~~~~~~~~~~~~~~

This model is intended for processing data from the QUASAR ASL sequence.

This model is selected using ``--model=quasar``.

The Turbo-QUASAR model
~~~~~~~~~~~~~~~~~~~~~~

This model is intended for processing data from the Turbo-QUASAR ASL sequence.

This model is selected using ``--model=turboquasar``.

The Saturation-Recovery model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This model is designed for the saturation recovery curve calibration method.

This model is selected using ``--model=satrecov``.

The multi-TE model
~~~~~~~~~~~~~~~~~~

This model is for ASL data captured at a series of different TE values. It is used
by the multi-TE plugin for the OXASL_ pipeline.

This model is selected using ``--model=asl_multite``.

The Buxton model
~~~~~~~~~~~~~~~~

This model implements only the basic Buxton kinetic model. It has been superceded
by the more generic ``aslrest`` model and is kept only for historical compatibility.

The 2-compartment model
~~~~~~~~~~~~~~~~~~~~~~~

This model has been superceded by the exchange options in the ``aslrest`` model.

Examples
--------

References
----------

.. [1] *Chappell, M.A., Groves, A.R., Woolrich, M.W., "Variational Bayesian
   inference for a non-linear forward model", IEEE Trans. Sig. Proc., 2009,
   57(1), 223–236.*

.. _Fabber: https://fabber-core.readthedocs.io/

.. _FSL: https://fsl.fmrib.ox.ac.uk/fsl/

.. _BASIL: https://asl-docs.readthedocs.io/

.. _OXASL: https://oxasl.readthedocs.io/
