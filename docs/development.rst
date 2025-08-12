.. highlight:: shell

===========
Development
===========

If you're one of the phrosty developers, or are otherwise interested in contributing, here is some useful information.


Running Tests
-------------

Running all of the tests requires an NVIDIA GPU with enough memory.  We are able to run them on 40GB NVIDIA GPUs, a GPU with only 12GB is not enough.  (TODO: figure out the actual cutoff.)

To run the tests, make sure to run the SNPIT container as described in :ref:`running-snpit-container`.  Inside the container, cd into the `/phrosty/phrosty/tests` directory and run::

  SNPIT_CONFIG=phrosty_test_config.yaml pytest -v

