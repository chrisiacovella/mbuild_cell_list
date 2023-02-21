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
    

def test_init_cell_list():
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
    assert cell_list.cells[0].pos.all() == np.array([0.5,0.5,0.5]).all()
    assert cell_list.cells[1].pos.all() == np.array([1.5,0.5,0.5]).all()
    assert cell_list.cells[2].pos.all() == np.array([2.5,0.5,0.5]).all()
    assert cell_list.cells[3].pos.all() == np.array([0.5,1.5,0.5]).all()
    assert cell_list.cells[4].pos.all() == np.array([1.5,1.5,0.5]).all()
    assert cell_list.cells[5].pos.all() == np.array([2.5,1.5,0.5]).all()
    assert cell_list.cells[6].pos.all() == np.array([0.5,2.5,0.5]).all()
    assert cell_list.cells[7].pos.all() == np.array([1.5,2.5,0.5]).all()
    assert cell_list.cells[8].pos.all() == np.array([2.5,2.5,0.5]).all()
    assert cell_list.cells[9].pos.all() == np.array([0.5,0.5,1.5]).all()
    assert cell_list.cells[10].pos.all() == np.array([1.5,0.5,1.5]).all()
    assert cell_list.cells[11].pos.all() == np.array([2.5,0.5,1.5]).all()
    assert cell_list.cells[12].pos.all() == np.array([0.5,1.5,1.5]).all()
    assert cell_list.cells[13].pos.all() == np.array([1.5,1.5,1.5]).all()
    assert cell_list.cells[14].pos.all() == np.array([2.5,1.5,1.5]).all()
    assert cell_list.cells[15].pos.all() == np.array([0.5,2.5,1.5]).all()
    assert cell_list.cells[16].pos.all() == np.array([1.5,2.5,1.5]).all()
    assert cell_list.cells[17].pos.all() == np.array([2.5,2.5,1.5]).all()
    assert cell_list.cells[18].pos.all() == np.array([0.5,0.5,2.5]).all()
    assert cell_list.cells[19].pos.all() == np.array([1.5,0.5,2.5]).all()
    assert cell_list.cells[20].pos.all() == np.array([2.5,0.5,2.5]).all()
    assert cell_list.cells[21].pos.all() == np.array([0.5,1.5,2.5]).all()
    assert cell_list.cells[22].pos.all() == np.array([1.5,1.5,2.5]).all()
    assert cell_list.cells[23].pos.all() == np.array([2.5,1.5,2.5]).all()
    assert cell_list.cells[24].pos.all() == np.array([0.5,2.5,2.5]).all()
    assert cell_list.cells[25].pos.all() == np.array([1.5,2.5,2.5]).all()
    assert cell_list.cells[26].pos.all() == np.array([2.5,2.5,2.5]).all()

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

    """ Check to see if changing the periodicity yields the correct number of neighboring cells for each cell."""
    cell_list = mbcl.CellList(box=box, n_cells=[3,3,3], periodicity=[False,False,False], box_min=[0,0,0])
    assert (cell_list.periodicity == np.array([False,False,False])).all()

    assert len(cell_list.cells[0].neighbor_cells) == 7
    assert len(cell_list.cells[1].neighbor_cells) == 26
    assert len(cell_list.cells[2].neighbor_cells) == 7
    assert len(cell_list.cells[3].neighbor_cells) == 7
    assert len(cell_list.cells[4].neighbor_cells) == 26
    assert len(cell_list.cells[5].neighbor_cells) == 7
    assert len(cell_list.cells[6].neighbor_cells) == 7
    assert len(cell_list.cells[7].neighbor_cells) == 26
    assert len(cell_list.cells[8].neighbor_cells) == 7
    assert len(cell_list.cells[9].neighbor_cells) == 7
    assert len(cell_list.cells[10].neighbor_cells) == 26
    assert len(cell_list.cells[11].neighbor_cells) == 7
    assert len(cell_list.cells[12].neighbor_cells) == 7
    assert len(cell_list.cells[13].neighbor_cells) == 26
    assert len(cell_list.cells[14].neighbor_cells) == 7
    assert len(cell_list.cells[15].neighbor_cells) == 7
    assert len(cell_list.cells[16].neighbor_cells) == 26
    assert len(cell_list.cells[17].neighbor_cells) == 7
    assert len(cell_list.cells[18].neighbor_cells) == 7
    assert len(cell_list.cells[19].neighbor_cells) == 26
    assert len(cell_list.cells[20].neighbor_cells) == 7
    assert len(cell_list.cells[21].neighbor_cells) == 7
    assert len(cell_list.cells[22].neighbor_cells) == 26
    assert len(cell_list.cells[23].neighbor_cells) == 7
    assert len(cell_list.cells[24].neighbor_cells) == 7
    assert len(cell_list.cells[25].neighbor_cells) == 26
    assert len(cell_list.cells[26].neighbor_cells) == 7
    
    cell_list = mbcl.CellList(box=box, n_cells=[3,3,3], periodicity=[False,False,True], box_min=[0,0,0])
    assert (cell_list.periodicity == np.array([False,False,True])).all()

    assert len(cell_list.cells[0].neighbor_cells) == 11
    assert len(cell_list.cells[1].neighbor_cells) == 26
    assert len(cell_list.cells[2].neighbor_cells) == 11
    assert len(cell_list.cells[3].neighbor_cells) == 11
    assert len(cell_list.cells[4].neighbor_cells) == 26
    assert len(cell_list.cells[5].neighbor_cells) == 11
    assert len(cell_list.cells[6].neighbor_cells) == 11
    assert len(cell_list.cells[7].neighbor_cells) == 26
    assert len(cell_list.cells[8].neighbor_cells) == 11
    assert len(cell_list.cells[9].neighbor_cells) == 11
    assert len(cell_list.cells[10].neighbor_cells) == 26
    assert len(cell_list.cells[11].neighbor_cells) == 11
    assert len(cell_list.cells[12].neighbor_cells) == 11
    assert len(cell_list.cells[13].neighbor_cells) == 26
    assert len(cell_list.cells[14].neighbor_cells) == 11
    assert len(cell_list.cells[15].neighbor_cells) == 11
    assert len(cell_list.cells[16].neighbor_cells) == 26
    assert len(cell_list.cells[17].neighbor_cells) == 11
    assert len(cell_list.cells[18].neighbor_cells) == 11
    assert len(cell_list.cells[19].neighbor_cells) == 26
    assert len(cell_list.cells[20].neighbor_cells) == 11
    assert len(cell_list.cells[21].neighbor_cells) == 11
    assert len(cell_list.cells[22].neighbor_cells) == 26
    assert len(cell_list.cells[23].neighbor_cells) == 11
    assert len(cell_list.cells[24].neighbor_cells) == 11
    assert len(cell_list.cells[25].neighbor_cells) == 26
    assert len(cell_list.cells[26].neighbor_cells) == 11
    
    cell_list = mbcl.CellList(box=box, n_cells=[3,3,3], periodicity=[False,True,True], box_min=[0,0,0])
    assert (cell_list.periodicity == np.array([False,True,True])).all()

    assert len(cell_list.cells[0].neighbor_cells) == 17
    assert len(cell_list.cells[1].neighbor_cells) == 26
    assert len(cell_list.cells[2].neighbor_cells) == 17
    assert len(cell_list.cells[3].neighbor_cells) == 17
    assert len(cell_list.cells[4].neighbor_cells) == 26
    assert len(cell_list.cells[5].neighbor_cells) == 17
    assert len(cell_list.cells[6].neighbor_cells) == 17
    assert len(cell_list.cells[7].neighbor_cells) == 26
    assert len(cell_list.cells[8].neighbor_cells) == 17
    assert len(cell_list.cells[9].neighbor_cells) == 17
    assert len(cell_list.cells[10].neighbor_cells) == 26
    assert len(cell_list.cells[11].neighbor_cells) == 17
    assert len(cell_list.cells[12].neighbor_cells) == 17
    assert len(cell_list.cells[13].neighbor_cells) == 26
    assert len(cell_list.cells[14].neighbor_cells) == 17
    assert len(cell_list.cells[15].neighbor_cells) == 17
    assert len(cell_list.cells[16].neighbor_cells) == 26
    assert len(cell_list.cells[17].neighbor_cells) == 17
    assert len(cell_list.cells[18].neighbor_cells) == 17
    assert len(cell_list.cells[19].neighbor_cells) == 26
    assert len(cell_list.cells[20].neighbor_cells) == 17
    assert len(cell_list.cells[21].neighbor_cells) == 17
    assert len(cell_list.cells[22].neighbor_cells) == 26
    assert len(cell_list.cells[23].neighbor_cells) == 17
    assert len(cell_list.cells[24].neighbor_cells) == 17
    assert len(cell_list.cells[25].neighbor_cells) == 26
    assert len(cell_list.cells[26].neighbor_cells) == 17
    
    """ Check that changing the number of cells produces the correct behavior"""
    cell_list = mbcl.CellList(box=box, n_cells=[5,3,7], periodicity=[True,True,True], box_min=[0,0,0])

    assert (cell_list.n_cells == np.array([5,3,7])).all()
    assert cell_list.n_cells_total == 105
    
    for c, cell in enumerate(cell_list.cells):
        assert len(cell_list.members(c)) == 0
        assert len(cell_list.neighbor_members(c)) == 0
        assert len(cell.neighbor_cells) == 26

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

def test_insert_compounds():
    box = mb.Box([3,3,3])
    argon = mb.Compound(name='Ar', element='Ar', charge=0)
    
    system = mb.Compound()
    cell_list = mbcl.CellList(box=box, n_cells=[3,3,3], periodicity=[True,True,True], box_min=[0,0,0])

    #  when call the cell_containing function with the center of the cell, the function should return the same cell
    for c, cell in enumerate(cell_list.cells):
        assert cell_list.cell_containing(cell.pos) == c

    # same as above, but let us check that shifting the box minimum also produces the correct cell
    cell_list = mbcl.CellList(box=box, n_cells=[3,3,3], periodicity=[True,True,True], box_min=[0.2,1.2,0.8])

    for c, cell in enumerate(cell_list.cells):
        assert cell_list.cell_containing(cell.pos) == c

    cell_list = mbcl.CellList(box=box, n_cells=[3,3,3], periodicity=[True,True,True], box_min=[0,0,0])

    # add argon atoms at the center of each cell
    for c, cell in enumerate(cell_list.cells):
        temp = mb.clone(argon)
        temp.translate_to(cell.pos)
        system.add(temp)
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
        
    # add more particles to each cell. they will overlap but that is fine for this text
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

    # this should fail if we mix ways of initializing
    with pytest.raises(Exception):
        cell_list.insert_compound_position(argon)
        
    with pytest.raises(Exception):
        cell_list._check_cell(cell_list.n_cells_total+10)
    assert cell_list._check_cell(0) == True

