internal_autodockvina_contestant
================================

.. image:: https://travis-ci.org/cookiecutter/cookiecutter-pycustomdock.png
   :target: https://travis-ci.org/cookiecutter/cookiecutter-pycustomdock
   :alt: Latest Travis CI build status

Containerized in-house AutoDock Vina contestant for `D3R CELPP <https://drugdesigndata.org/about/celpp>`_ competition. 

Usage
-----

.. code-block:: bash

   cd /tmp
   # get challenge data
   challdir="1-get_challenge_data/"
   mkdir -p $challdir
   singularity run getchallengedata.py --unpackdir $challdir -f ~/ftp.config

   # protein prep
   protdir="2-protein_prep/"
   mkdir $protdir
   singularity run internal_autodockvina_contestant_protein_prep.py --challengedata $challdir --prepdir $protdir
   
   # ligand prep
   ligdir="3-ligand_prep/"
   mkdir $ligdir
   singularity run internal_autodockvina_contestant_ligand_prep.py --challengedata $challdir --prepdir $ligdir

   #dock
   dockdir="4-docking/"
   mkdir $dockdir
   singularity run internal_autodockvina_contestant_dock.py --protsciprepdir $protdir --ligsciprepdir $ligdir --outdir $dockdir

   #upload results
   packdir="5-pack_docking_results"
   mkdir $packdir
   singularity run packdockingresults.py --dockdir $dockdir --packdir $packdir --challengedata $challdir -f ~/ftp.config


Building the container
----------------------

Build Requirements
^^^^^^^^^^^^^^^^^^

* Vagrant https://www.vagrantup.com/

* Virtual Box https://www.virtualbox.org/

* Binary of 64-bit Linux distribution of Chimera (tested with version `1.13 <https://www.cgl.ucsf.edu/chimera/cgi-bin/secure/chimera-get.py?file=linux_x86_64/chimera-1.13-linux_x86_64.bin>`_) https://www.cgl.ucsf.edu/chimera

The following commands spin up a `Virtual Box <https://www.virtualbox.org>`_ virtual machine via `Vagrant <https://www.vagrantup.com>`_ with `Singularity <https://www.sylabs.io>`_ installed. A `Makefile <https://www.gnu.org/software/make/manual/make.html>`_ is then used to create the `Singularity <https://www.sylabs.io>`_ Container runnable on any machine that can run `Singularity <https://www.sylabs.io>`_. 


.. code-block:: bash

   git clone https://github.com/drugdata/internal_autodockvina_contestant.git
   cd internal_autodockvina_contestant
   #
   # Be sure to download 64-bit Linux version of Chimera and put binary
   # in source tree directory
   #
   vagrant up
   vagrant ssh
   cd /vagrant
   make singularity

Compatibility
-------------

License
-------

See LICENSE.txt_

Authors
-------

`internal_autodockvina_contestant` was written by `Jeff Wagner <j5wagner@ucsd.edu>`_.
