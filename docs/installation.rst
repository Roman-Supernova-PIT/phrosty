.. highlight:: shell

============
Installation
============

Prerequisites
-------------

Currently, phrosty is designed inside a container built from the the Roman Supernova PIT environment defined at https://github.com/Roman-Supernova-PIT/environment .

You can pull the necessary container from one of the following two sources:

* ``registry.nersc.gov/m4385/rknop/roman-snpit-env:cuda-dev``
* ``docker.io/rknop/roman-snpit-env:cuda-dev``

Because phrosty (and other libraries it depends on, such as snappl) is under heavy development, it's possible that the latest container will not work properly with phrosty at any given moment.

.. _running-snpit-container:
Running the SNPIT container
---------------------------

TODO.


Stable release
--------------

(NOTE: this is not currently supported, we don't have up-to-date pip releases.)

To install phrosty, run this command in your terminal:

.. code-block:: console

    $ pip install phrosty

This is the preferred method to install phrosty, as it will always install the most recent stable release.

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/


From sources
------------

The sources for phrosty can be downloaded from the `Github repo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git@github.com:Roman-Supernova-PIT/phrosty.git

Or download the `tarball`_:

.. code-block:: console

    $ curl -OJL https://github.com/Roman-Supernova-PIT/phrosty/tarball/master

Once you have a copy of the source, you can install it with:

.. code-block:: console

    $ pip install .

If you would like to do an editable install:

.. code-block:: console

    $ pip install -e .
    $ pip install -e .[docs]  # install document build packages during install


.. _Github repo: https://github.com/Roman-Supernova-PIT/phrosty
.. _tarball: https://github.com/Roman-Supernova-PIT/phrosty/tarball/master
