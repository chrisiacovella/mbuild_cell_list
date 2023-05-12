Installation
===============

pip install
-----------
Check out the source from the github repository:

    $ git clone https://github.com/chrisiacovella/mbuild_cell_list.git

In the top level of hoomdxml_reader directory, use pip to install:

    $ pip install -e .

This requires mbuild and numpy to be installed. To create an environment named mbuild_cell_list
with this necessary packages, run the following from the top level of the hoomdxml_reader directory.

    $ conda env create -f environment.yml
    
Building the documentation
--------------------------

This packag uses `sphinx <https://www.sphinx-doc.org/en/master/index.html>`_ to build its documentation. To build the docs locally, run the following while in the ``docs`` directory::

    $ pip install -r requirements.yaml
    $ make html

To view documentation, open ``$INSTALL_PATH/hoomdxml_reader/docs/_build/html/index.html`` in your web browser.

Usage
===========

Basic usage
------------
The cell list can be populated in one of two ways.  (1) For a given Compound inserted into the cell list, it will be inserted into a cell based on the Compound center-of-mass. (2) The individual Particles that make up a Compound are inserted into the cell.

Let us first start by construct a box of ethane molecules.

.. code:: ipython3
    import mbuild as mb
    import numpy as np
    import mbuild_cell_list as mbcl
    from mbuild.lib.molecules.ethane import Ethane

    system_box = mb.Box([2,2,2])
    ethane_box = mb.fill_box(Ethane(), n_compounds=100, overlap=0.22,
                       edge=0.15, sidemax=2, box=system_box, seed=12345)




To insert Compounds based on their center of mass, we use the 'insert_compound_position' function after initializing an empty cell list.

.. code:: ipython3
    
    # examine the structure of the ethane_box Compound
    ethane_box.print_hierarchy()
    '''
    Compound, 800 particles, 700 bonds, 100 children
    └── [Ethane x 100], 8 particles, 7 bonds, 2 children
        └── [CH3 x 2], 4 particles, 3 bonds, 4 children
            ├── [C x 1], 1 particles, 4 bonds, 0 children
            └── [H x 3], 1 particles, 1 bonds, 0 children
    '''
    # initialize the cell list with 4 cells in each direction
    ethane_cell_list = mbcl.CellList(box=system_box, n_cells=[4, 4, 4])

    # looping over ethane_box.children will give us the ethane Compounds
    for child in ethane_box.children:
        ethane_cell_list.insert_compound_position(child)

With the cell list initialized and Compounds added, we can examine the contents:

.. code:: ipython3
    # print the total number of cells
    print(f'# cells: {ethane_cell_list.n_cells_total}')

    # print number of cells in each direction
    print(f'# cells in each direction: {ethane_cell_list.n_cells}')

    # print size of cells in each direction
    print(f'length of cells in each direction: {ethane_cell_list.cell_sizes}')

    # print periodicity
    print(f'periodicity: {ethane_cell_list.periodicity}')

*output*:

.. code:: ipython3
    # cells: 64
    # cells in each direction: [4 4 4]
    length of cells in each direction: [0.5 0.5 0.5]
    periodicity: [ True  True  True]

Examine the contents of the cells (note the 'id' will be different each time executed):

.. code:: ipython3

    # examine the contents of the cell list
    c = 0
    print(f'# members in cell {c}: {len(ethane_cell_list.members(c))}')
    print(f'members in cell {c}: {ethane_cell_list.members(c)}')

    # examine the contents of the 27 neighboring cells
    print(f'# members in neighboring cell {c}: {len(ethane_cell_list.neighbor_members(c))}')

*output*:

.. code:: ipython3
    # members in cell 0: 1
    members in cell 0: [<Ethane 8 particles, 7 bonds, non-periodic, id: 5691885888>]
    # members in neighboring cell 0: 33

We can additionally return the minimum image shifting along with each neighbor. This allows us to quicky calculate the distance between compound
.. code:: ipython3
    neigh_im = ethane_cell_list.neighbor_members_and_min_image_shift(c)

    print(f'First neighbor Compound: {neigh_im[0][0]}')
    print(f'The minimum image shifting for the compound: {neigh_im[0][1]}')

    # we can use this to calculate the distance between two compounds
    compounds_in_cell = ethane_cell_list.members(c)
    distance = np.linalg.norm(compounds_in_cell[0].pos-neigh_im[0][0].pos-neigh_im[0][1]*ethane_cell_list.box.lengths)

    print(f'Distance between first compound and first neighbor: {distance.round(4)} nm')
*output*:

.. code:: ipython3
    First neighbor Compound: <Ethane 8 particles, 7 bonds, non-periodic, id: 5684693552>
    The minimum image shifting for the compound: [-1  0  0]
    Distance between first compound and first neighbor: 0.6805 nm

To insert the individual particles that make up a Compound into the cell list, we can  use the 'insert_compound_particles' function. Note, unlike above, we do not need to traverse the hierarchy to individually insert Particles, this is automatically done for the provided Compound (since Particles is a well-defined level in the tree).

.. code:: ipython3

    ethane_particles_cell_list = mbcl.CellList(box=system_box, n_cells=[4, 4, 4])
    ethane_particles_cell_list.insert_compound_particles(ethane_box)


Note, we cannot use both the insert_compound_position and insert_compound_particles function with the same cell list instance. Howeer, since Particles are simply just Compounds, we can use the insert_compound_position function to insert Particles alongside center-of-mass representations of a Compound, but the routine will not automatically traverse the hierarchy.

We can examine the contents the same as above
.. code:: ipython3

    c = 0
    print(f'# members in cell {c}: {len(ethane_particles_cell_list.members(c))}')
    print(f'members in cell {c}: {ethane_particles_cell_list.members(c)}')
    
*output*:

.. code:: ipython3

    # members in cell 0: 8
    members in cell 0: [<C pos=([0.2574 0.3175 0.2521]), 4 bonds, id: 5691886320>, <H pos=([0.2418 0.3457 0.15  ]), 1 bonds, id: 5691886464>, <H pos=([0.302  0.3989 0.3052]), 1 bonds, id: 5691886608>, <H pos=([0.1634 0.2928 0.297 ]), 1 bonds, id: 5691886752>, <C pos=([0.298  0.19   0.2107]), 4 bonds, id: 5691887040>, <H pos=([0.3135 0.1618 0.3127]), 1 bonds, id: 5691887184>, <H pos=([0.3769 0.1503 0.1502]), 1 bonds, id: 5691887328>, <H pos=([0.2041 0.1512 0.1768]), 1 bonds, id: 5691887472>]


