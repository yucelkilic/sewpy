Source Extractor Wrapper for Python
===================================

``sewpy`` is a lightweight Python wrapper for
`SExtractor <https://www.astromatic.net/software/sextractor>`_,
allowing you to run it as if it were a native Python module.

Example
-------

.. code-block:: python

    import sewpy

    sew = sewpy.SEW(
        params=["X_IMAGE", "Y_IMAGE", "FLUX_RADIUS(3)", "FLAGS"],
        config={"DETECT_MINAREA": 10, "PHOT_FLUXFRAC": "0.3, 0.5, 0.8"}
    )

    out = sew("myimage.fits")
    print(out["table"])  # Returns an astropy.table.Table

This simple interface hides the complexity of temporary files and system calls,
while still giving you full access to SExtractor’s powerful features.

Motivation
----------

Why `yet another <https://pypi.org/project/pysex/>`_
`SExtractor wrapper <https://gitorious.org/pysextractor>`_?

Because ``sewpy`` was designed to be:

* Fully compatible with **Python ≥ 3.8**
* Based on the `Astropy <https://www.astropy.org>`_ ecosystem (especially ``astropy.table``)
* Using the standard ``logging`` module (captures both stdout & stderr from SExtractor)
* Managing temporary input/output files automatically via ``tempfile``
  *(but lets you keep them if you want)*
* Providing convenient support for SExtractor’s ``ASSOC`` feature
  *(e.g. appending new measurements to an existing catalog)*

Demonstrations
--------------

The **examples/** directory contains ready-to-run demos that can be executed
without installing the package. They provide a quick overview of what ``sewpy`` can do.

Installation
------------

Modern installation (PEP 517 build system):

.. code-block:: bash

    # Clone the repository
    git clone https://github.com/yucelkilic/sewpy.git
    cd sewpy

    # Install in editable mode using pip (recommended)
    pip install -e .

Or simply from PyPI (when available):

.. code-block:: bash

    pip install sewpy

Dependencies
------------

* Python ≥ 3.8
* `Astropy <https://www.astropy.org>`_ ≥ 7.1
* NumPy ≥ 2.2
* SExtractor (installed and available in your system PATH)

Development and Documentation
-----------------------------

To learn more about how to install and use ``sewpy``, visit the documentation:
*(The original ReadTheDocs link may be outdated; an updated version is available on GitHub.)*

* GitHub repository: https://github.com/yucelkilic/sewpy
* Documentation: https://sewpy.readthedocs.io (legacy version)

License
-------

Licensed under the GNU General Public License v3 (GPLv3).
See the LICENSE file for details.