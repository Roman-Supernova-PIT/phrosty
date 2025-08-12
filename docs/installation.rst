.. highlight:: shell

============
Installation
============

.. _phrosty-installation-prerequisites:
Prerequisites
-------------

Currently, phrosty is designed inside a container built from the the `Roman Supernova PIT environment <https://github.com/Roman-Supernova-PIT/environment>`_.

You can pull the necessary container image from one of the following two sources:

* ``registry.nersc.gov/m4385/rknop/roman-snpit-env:cuda-dev``
* ``docker.io/rknop/roman-snpit-env:cuda-dev``

Because phrosty (and other libraries it depends on, such as snappl) is under heavy development, it's possible that the latest container will not work properly with phrosty at any given moment.

Stable release
--------------

phrosty is under heavy development and does not currently have stable releases.  Eventually it will be released on pypi either under the name ``phrosty`` or ``roman-snpit-phrosty`` (depending on what's available).

.. _install-from-sources:
From sources
------------

Currently, the only way to install phrosty is to download it from the `github repo <https://github.com/Roman-Supernova-PIT/phrosty>`_.  Clone it with::

    git clone https://github.com/Roman-Supernova-PIT/phrosty.git

(you can also clone it via the ``git@`` code link if you know what you're doing.)

.. _running-snpit-container:
Running the SNPIT container
---------------------------

To use phrosty inside the container, you will need to run it with ``docker`` or ``podman``, and bind-mount the directory where you've cloned phrosty.  Assuming you are running from the directory where you ran the ``git clone`` command in the previous step, that command will look something like::

  docker run --mount type=bind,source=$PWD/phrosty,target=/phrosty \
    -it docker.io/rknop/roman-snpit-env:cuda-dev \
    /bin/bash

In pratctice, however, you will almost certainly need to bind-mount other directories so that phrosty will have data to work with.  For a couple of examples, see :ref:`running-tests` and :ref:`example-usage`.




.. _Github repo: https://github.com/Roman-Supernova-PIT/phrosty
.. _tarball: https://github.com/Roman-Supernova-PIT/phrosty/tarball/master
