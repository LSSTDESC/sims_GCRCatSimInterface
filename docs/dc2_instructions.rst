==========================
InstanceCatalog generation
==========================

Before generating InstanceCatalogs
----------------------------------

- Before generating InstanceCatalogs, go into the `data/` directory
  of this repo and run the `get_sed_mags.py` script.  That will generate
  the `data/CatSimMagGrid.txt` file, which contains the data necessary
  to fit protoDC2 galaxies to SEDs from the CatSim library.

- Before generating InstanceCatalogs, you need to run the script
  ``bin.src/create_agn_db.py``.  This will create a ``sqlite`` database
  of AGN variability parameters associated with the protoDC2 galaxies.
  This database will need to be created in a space where all of the
  InstanceCatalog generation processes can connect to it.  You can
  specify the location of the database using the ``--out_dir`` and
  ``--out_file`` command line options.  You can use the ``--mbh_cut``
  to limit the mass of black holes that are treated as AGN and
  the ``--m_i_cut`` to put a limit on the observed i-band magnitude of
  AGNs that will be simulating.  Running the script with::

      python bin.src/create_agn_db.py --mbh_cut 7.0 --m_i_cut 30.0

  Results in a 240MB database file.
