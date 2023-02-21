"""
Unit and regression test for the mbuild_cell_list package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import mbuild_cell_list as mbcl
import mbuild as mb
import numpy as np

def test_mbuild_cell_list_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "mbuild_cell_list" in sys.modules
    

def test_init_cell_list_basic():
    """Examine the minimal size cell list to ensure behavior is as expected"""
    box = mb.Box([3,3,3])
    cell_list = mbcl.CellList(box=box, n_cells=[3,3,3], periodicity=[True,True,True], box_min=[0,0,0])

    # check that variables are correctly set
    assert (cell_list.n_cells == np.array([3,3,3])).all()
    assert cell_list.n_cells_total == 27
    assert cell_list.box == box
    assert (cell_list.periodicity == np.array([True,True,True])).all()
    
    # examine each cell
    for c, cell in enumerate(cell_list.cells):
        assert len(cell_list.members(c)) == 0
        assert len(cell_list.neighbor_members(c)) == 0
        assert len(cell.neighbor_cells) == 26
    
    # check that the center of each cell is correctly determined.
    assert (cell_list.cells[0].pos == np.array([0.5,0.5,0.5])).all()
    assert (cell_list.cells[1].pos == np.array([1.5,0.5,0.5])).all()
    assert (cell_list.cells[2].pos == np.array([2.5,0.5,0.5])).all()
    assert (cell_list.cells[3].pos == np.array([0.5,1.5,0.5])).all()
    assert (cell_list.cells[4].pos == np.array([1.5,1.5,0.5])).all()
    assert (cell_list.cells[5].pos == np.array([2.5,1.5,0.5])).all()
    assert (cell_list.cells[6].pos == np.array([0.5,2.5,0.5])).all()
    assert (cell_list.cells[7].pos == np.array([1.5,2.5,0.5])).all()
    assert (cell_list.cells[8].pos == np.array([2.5,2.5,0.5])).all()
    assert (cell_list.cells[9].pos == np.array([0.5,0.5,1.5])).all()
    assert (cell_list.cells[10].pos == np.array([1.5,0.5,1.5])).all()
    assert (cell_list.cells[11].pos == np.array([2.5,0.5,1.5])).all()
    assert (cell_list.cells[12].pos == np.array([0.5,1.5,1.5])).all()
    assert (cell_list.cells[13].pos == np.array([1.5,1.5,1.5])).all()
    assert (cell_list.cells[14].pos == np.array([2.5,1.5,1.5])).all()
    assert (cell_list.cells[15].pos == np.array([0.5,2.5,1.5])).all()
    assert (cell_list.cells[16].pos == np.array([1.5,2.5,1.5])).all()
    assert (cell_list.cells[17].pos == np.array([2.5,2.5,1.5])).all()
    assert (cell_list.cells[18].pos == np.array([0.5,0.5,2.5])).all()
    assert (cell_list.cells[19].pos == np.array([1.5,0.5,2.5])).all()
    assert (cell_list.cells[20].pos == np.array([2.5,0.5,2.5])).all()
    assert (cell_list.cells[21].pos == np.array([0.5,1.5,2.5])).all()
    assert (cell_list.cells[22].pos == np.array([1.5,1.5,2.5])).all()
    assert (cell_list.cells[23].pos == np.array([2.5,1.5,2.5])).all()
    assert (cell_list.cells[24].pos == np.array([0.5,2.5,2.5])).all()
    assert (cell_list.cells[25].pos == np.array([1.5,2.5,2.5])).all()
    assert (cell_list.cells[26].pos == np.array([2.5,2.5,2.5])).all()

    """ Check that changing the number of cells produces the correct behavior"""
    cell_list = mbcl.CellList(box=box, n_cells=[5,3,7], periodicity=[True,True,True], box_min=[0,0,0])

    assert (cell_list.n_cells == np.array([5,3,7])).all()
    assert cell_list.n_cells_total == 105
    
    for c, cell in enumerate(cell_list.cells):
        assert len(cell_list.members(c)) == 0
        assert len(cell_list.neighbor_members(c)) == 0
        assert len(cell.neighbor_cells) == 26
        
def test_init_cell_list_periodicity_FFF():

    box = mb.Box([3,3,3])

    """ Check to see if changing the periodicity yields the correct number of neighboring cells for each cell."""
    cell_list = mbcl.CellList(box=box, n_cells=[3,3,3], periodicity=[False,False,False], box_min=[0,0,0])
    assert (cell_list.periodicity == np.array([False,False,False])).all()
    
    num_neighbors = [7, 26, 7, 7, 26, 7, 7, 26, 7, 7, 26, 7, 7, 26, 7, 7, 26, 7, 7, 26, 7, 7, 26, 7, 7, 26, 7 ]
    for c, cell in enumerate(cell_list.cells):
        assert len(cell_list.cells[c].neighbor_cells) == num_neighbors[c]

def test_init_cell_list_periodicity_FFT():

    box = mb.Box([3,3,3])
    cell_list = mbcl.CellList(box=box, n_cells=[3,3,3], periodicity=[False,False,True], box_min=[0,0,0])
    assert (cell_list.periodicity == np.array([False,False,True])).all()
    
    num_neighbors = [11, 26, 11, 11, 26, 11, 11, 26, 11, 11, 26, 11, 11, 26, 11, 11, 26, 11, 11, 26, 11, 11, 26, 11, 11, 26, 11 ]
    for c, cell in enumerate(cell_list.cells):
        assert len(cell_list.cells[c].neighbor_cells) == num_neighbors[c]

def test_init_cell_list_periodicity_FTT():

    box = mb.Box([3,3,3])
    cell_list = mbcl.CellList(box=box, n_cells=[3,3,3], periodicity=[False,True,True], box_min=[0,0,0])
    assert (cell_list.periodicity == np.array([False,True,True])).all()

    num_neighbors = [17, 26, 17, 17, 26, 17, 17, 26, 17, 17, 26, 17, 17, 26, 17, 17, 26, 17, 17, 26, 17, 17, 26, 17, 17, 26, 17 ]
    for c, cell in enumerate(cell_list.cells):
        assert len(cell_list.cells[c].neighbor_cells) == num_neighbors[c]


def test_init_cell_list_too_small():
    box = mb.Box([3,3,3])

    """ check that the code fails if we have any dimension less than 3."""
    with pytest.raises(Exception):
        cell_list = mbcl.CellList(box=box, n_cells=[2,3,3], periodicity=[True,True,True], box_min=[0,0,0])
    with pytest.raises(Exception):
        cell_list = mbcl.CellList(box=box, n_cells=[3,2,3], periodicity=[True,True,True], box_min=[0,0,0])
    with pytest.raises(Exception):
        cell_list = mbcl.CellList(box=box, n_cells=[3,3,2], periodicity=[True,True,True], box_min=[0,0,0])
    with pytest.raises(Exception):
        cell_list = mbcl.CellList(box=box, n_cells=[3,2,2], periodicity=[True,True,True], box_min=[0,0,0])
    with pytest.raises(Exception):
        cell_list = mbcl.CellList(box=box, n_cells=[2,3,2], periodicity=[True,True,True], box_min=[0,0,0])
    with pytest.raises(Exception):
        cell_list = mbcl.CellList(box=box, n_cells=[2,2,3], periodicity=[True,True,True], box_min=[0,0,0])
    with pytest.raises(Exception):
        cell_list = mbcl.CellList(box=box, n_cells=[2,2,2], periodicity=[True,True,True], box_min=[0,0,0])


 

    
def test_init_cell_list_check_shifting():
    box = mb.Box([3,3,3])

    """ Ensure that changing the box_min results in the correct shifting of each cell"""
    cell_list = mbcl.CellList(box=box, n_cells=[3,3,3], periodicity=[True,True,True], box_min=[1,0,0])
    assert (cell_list.cells[0].pos == np.array([1.5, 0.5, 0.5])).all()
    
    cell_list = mbcl.CellList(box=box, n_cells=[3,3,3], periodicity=[True,True,True], box_min=[1,1,0])
    assert (cell_list.cells[0].pos == np.array([1.5, 1.5, 0.5])).all()
    
    cell_list = mbcl.CellList(box=box, n_cells=[3,3,3], periodicity=[True,True,True], box_min=[1,1,1])
    assert (cell_list.cells[0].pos == np.array([1.5, 1.5, 1.5])).all()
    
    """ Check that the size of the cells are correct when changing the number of cells"""
    cell_list = mbcl.CellList(box=box, n_cells=[3,3,3], periodicity=[True,True,True], box_min=[0,0,0])
    assert (cell_list.cell_sizes == np.array([1.0,1.0,1.0])).all()

    cell_list = mbcl.CellList(box=box, n_cells=[6,6,6], periodicity=[True,True,True], box_min=[0,0,0])
    assert (cell_list.cell_sizes == np.array([0.5,0.5,0.5])).all()

    cell_list = mbcl.CellList(box=box, n_cells=[8,10,16], periodicity=[True,True,True], box_min=[0,0,0])
    assert (cell_list.cell_sizes == np.array([0.375, 0.3, 0.1875])).all()


def test_cell_containing():
    box = mb.Box([3,3,3])
    cell_list = mbcl.CellList(box=box, n_cells=[3,3,3], periodicity=[True,True,True], box_min=[0,0,0])

    #  when call the cell_containing function with the center of the cell, the function should return the same cell
    for c, cell in enumerate(cell_list.cells):
        assert cell_list.cell_containing(cell.pos) == c

    # we should raise an error if out of bounds
    with pytest.raises(Exception):
        cell_list.cell_containing([3.1,1.0,1.0])

    with pytest.raises(Exception):
        cell_list.cell_containing([1.0,3.1,1.0])

    with pytest.raises(Exception):
        cell_list.cell_containing([1.0,1.0,3.1])

    with pytest.raises(Exception):
        cell_list.cell_containing([-1.0,1.0,1.0])

    with pytest.raises(Exception):
        cell_list.cell_containing([1.0,-1.0,1.0])

    with pytest.raises(Exception):
        cell_list.cell_containing([1.0,1.0,-1.0])
        

    # same test as above, but let us check that shifting the box minimum also produces the correct cell
    cell_list = mbcl.CellList(box=box, n_cells=[3,3,3], periodicity=[True,True,True], box_min=[0.2,1.2,0.8])

    for c, cell in enumerate(cell_list.cells):
        assert cell_list.cell_containing(cell.pos) == c
    


        
def test_insert_compounds_by_com_basic():
    box = mb.Box([3,3,3])
 
    argon = mb.Compound(name='Ar', element='Ar', charge=0)
    
    
    cell_list = mbcl.CellList(box=box, n_cells=[3,3,3], periodicity=[True,True,True], box_min=[0,0,0])

  
    # add argon atoms at the center of each cell
    for c, cell in enumerate(cell_list.cells):
        temp = mb.clone(argon)
        temp.translate_to(cell.pos)
        cell_list.insert_compound_position(temp)
        assert len(cell_list.members(c)) == 1
    
    for c, cell in enumerate(cell_list.cells):
        assert len(cell_list.neighbor_members(c)) == 26
      
    # check to see we raise an exception if we are outside the bounds:
    cell_list = mbcl.CellList(box=box, n_cells=[3,3,3], periodicity=[True,True,True], box_min=[0,0,0])
    temp = mb.clone(argon)
    temp.translate_to([3.1,0.5,0.5])
    with pytest.raises(Exception):
        cell_list.insert_compound_position(temp)
        
    temp.translate_to([-3.1,0.5,0.5])
    with pytest.raises(Exception):
        cell_list.insert_compound_position(temp)
        
    cell_list = mbcl.CellList(box=box, n_cells=[3,3,3], periodicity=[True,True,True], box_min=[1,0,0])
    temp = mb.clone(argon)
    temp.translate_to([4.1,0.5,0.5])
    with pytest.raises(Exception):
        cell_list.insert_compound_position(temp)
        
    temp.translate_to([0.5,0.5,0.5])
    with pytest.raises(Exception):
        cell_list.insert_compound_position(temp)

def test_insert_compounds_by_com():
    box = mb.Box([3,3,3])
 
    argon = mb.Compound(name='Ar', element='Ar', charge=0)
    
    
    # add more particles to each cell. they will overlap but that is fine for this test
    # we will also add them to a compound called system for later checking
    system = mb.Compound()

    cell_list = mbcl.CellList(box=box, n_cells=[3,3,3], periodicity=[True,True,True], box_min=[0,0,0])

    for c, cell in enumerate(cell_list.cells):
        for j in range(0,3):
            temp = mb.clone(argon)
            temp.translate_to(cell.pos)
            system.add(temp)
            cell_list.insert_compound_position(temp)
        assert len(cell_list.members(c)) == 3

    for c, cell in enumerate(cell_list.cells):
        assert len(cell_list.neighbor_members(c)) == 78

    # if we try to add a compound to the cell list using underlying particles
    # we want the code to throw an exception if we used insert_compound_position previously
    with pytest.raises(Exception):
        cell_list.insert_compound_particles(system)
       
    # empty cells and check to ensure they were emptied
    cell_list.empty_cells()
    
    for c, cell in enumerate(cell_list.cells):
        assert len(cell_list.neighbor_members(c)) == 0
        assert len(cell_list.members(c)) == 0
    
    # if the cells are emptied, we can now insert compounds by the underlying particle

    cell_list.insert_compound_particles(system)
    for c, cell in enumerate(cell_list.cells):
        assert len(cell_list.members(c)) == 3
        assert len(cell_list.neighbor_members(c)) == 78

def check_cell_list_exceptions():
    box = mb.Box([3,3,3])
 
    argon = mb.Compound(name='Ar', element='Ar', charge=0)

    cell_list = mbcl.CellList(box=box, n_cells=[3,3,3], periodicity=[True,True,True], box_min=[0,0,0])
    for c, cell in enumerate(cell_list.cells):
        for j in range(0,3):
            temp = mb.clone(argon)
            temp.translate_to(cell.pos)
            cell_list.insert_compound_position(temp)

    with pytest.raises(Exception):
        cell_list._check_cell(cell_list.n_cells_total)
        
    assert cell_list._check_cell(0) == True
    
    with pytest.raises(Exception):
        cell_list.members(cell_list.n_cells_total)

    with pytest.raises(Exception):
        cell_list.neighbor_members(cell_list.n_cells_total)


def test_insert_compound_by_underlying_particles_basic():
    box = mb.Box([3,3,3])
 
    argon = mb.Compound(name='Ar', element='Ar', charge=0)
    
    
    # add more particles to each cell. they will overlap but that is fine for this test
    # we will also add them to a compound called system for later checking
    system = mb.Compound()
    coumpound_list = []
    cell_list = mbcl.CellList(box=box, n_cells=[3,3,3], periodicity=[True,True,True], box_min=[0,0,0])


    for c, cell in enumerate(cell_list.cells):
        for j in range(0,3):
            temp = mb.clone(argon)
            temp.translate_to(cell.pos)
            system.add(temp)
            coumpound_list.append(temp)

  
    cell_list.insert_compound_particles(system)
    for c, cell in enumerate(cell_list.cells):
        assert len(cell_list.members(c)) == 3
        assert len(cell_list.neighbor_members(c)) == 78

    # this should fail if we mix ways of initializing
    with pytest.raises(Exception):
        cell_list.insert_compound_position(argon)
     
       
    # empty cells and check to ensure they were emptied
    cell_list.empty_cells()
    
    for c, cell in enumerate(cell_list.cells):
        assert len(cell_list.neighbor_members(c)) == 0
        assert len(cell_list.members(c)) == 0
    
    # if the cells are emptied, we can now insert compounds by the center of mass
    for compound in coumpound_list:
          cell_list.insert_compound_position(compound)

    for c, cell in enumerate(cell_list.cells):
        assert len(cell_list.members(c)) == 3
        assert len(cell_list.neighbor_members(c)) == 78
    
    
