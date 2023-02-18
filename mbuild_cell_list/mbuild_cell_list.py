"""Simple cell list that is compatible with mBuild Compounds."""


__all__ = ["CellList"]

import mbuild as mb
import numpy as np


class Cell():
    """
    A generic container to hold the relevant
    data for each cell in the cell list.
    
    """
    def __init__(self):
        self._neighbor_cells = []
        self._members = []
        self._neighbor_members = []
        self._pos = np.array([0.0,0.0,0.0])
        
        
    @property
    def members(self):
        """Returns a list of all members in a cell."""
        return self._members
    
    @property
    def pos(self):
        """Returns the center of the cell as a numpy array."""
        return self._pos

    @property
    def neighbor_cells(self):
        """Returns a list of all cells that are neigbors of the current cell."""
        return self._neighbor_cells
        
    @property
    def neighbor_members(self):
        """Returns a list of all members of the cell and neighboring cells."""
        return self._neighbor_members


class CellList():
    """Cell list compatible with mbuild Compounds.
    The cell list can be constructed based on either the center of mass of a Compound
    or based on the particles contained in a Compound.

    Parameters
    ----------
    xyz :  np.ndarray, shape=3(,), dtype=float
    
    Returns
    ------
    c : int
        The cell containing the point
    """
    def __init__(self, box, n_cells=[3,3,3], periodicity=[True,True,True], box_min=[0,0,0]):
        self.box = box
        self._n_cells = np.array(n_cells,dtype=int)
        self._n_cells_total = np.prod(self._n_cells)
        

        if (self._n_cells < [3,3,3] ).any():
            raise Exception(f'The CellList must have at least 3 cells in each dimension, found: {n_cells}')
        
        self._box_min = box_min
        
        self._cell_sizes = np.array(self.box.lengths)/self._n_cells
        
        self.cells = []
        self._periodicity = periodicity
        
        self._init_full()
        self._from_particles = False
        self._from_com = False

    def _init_full(self):
        #initialize empty cells and calculate the center of each
        for i in range(0, self._n_cells_total):
            cell_temp = Cell()
            
            cell_temp._pos[0] = (i%self._n_cells[0])*self._cell_sizes[0]+self._box_min[0]+self._cell_sizes[0]/2.0
            cell_temp._pos[1] = (int(i/self._n_cells[0])%self._n_cells[1])*self._cell_sizes[1]+self._box_min[1]+self._cell_sizes[1]/2.0
            cell_temp._pos[2] = (int(i/(self._n_cells[0]*self._n_cells[1]))%self._n_cells[2])*self._cell_sizes[2]+self._box_min[2]+self._cell_sizes[2]/2.0

            self.cells.append(cell_temp)

        #define which cells are neighboring
        for k in range(0,self._n_cells[2]):
            for j in range(0,self._n_cells[1]):
                for i in range(0,self._n_cells[0]):
                    c = i + j*self._n_cells[0] + k*self._n_cells[0]*self._n_cells[1]
                    
                    start = []
                    end = []
                    
                    for kk in range(0,3):
                        if not self._periodicity[kk]:
                            if i == 0:
                                start.append(0)
                            else:
                                start.append(-1)
                            if i == self._n_cells[kk]-1:
                                end.append(0)
                            else:
                                end.append(1)
                        else:
                            start.append(-1)
                            end.append(1)
                    for z in range(start[2], end[2]+1):
                        for y in range(start[1], end[1]+1):
                            for x in range(start[0], end[0]+1):
                                cn = (i+x+self._n_cells[0])%self._n_cells[0] + \
                                 ((j+y+self._n_cells[1])%self._n_cells[1])*self._n_cells[0]+\
                                 ((k+z+self._n_cells[2])%self._n_cells[2])*self._n_cells[0]*self._n_cells[1]
                                if cn != c:
                                    self.cells[c]._neighbor_cells.append(cn)
    
    def cell_containing(self, xyz):
        """Return the cell that contains a given point in 3d space.

        Parameters
        ----------
        xyz :  np.ndarray, shape=3(,), dtype=float
        
        Returns
        ------
        c : int
            The cell containing the point
        """
        vals = np.array((np.array(xyz)-self._box_min)/self._cell_sizes, dtype=int)
        c = vals[0]+vals[1]*self._n_cells[0]+vals[2]*self._n_cells[0]*self._n_cells[1]
        
        return c
        
    
    def insert_compound_particles(self, compound):
        """This will look at the lowest level of the hierarchy of an mbuild Compound
        (i.e., the particles) and insert them  into the cell list.

        Parameters
        ----------
        compound :  mb.Compound
            An mbuild Compound whose particles will be inserted into the cell list.
        
        Returns
        ------
        """
        if self._from_particles == False and self._from_com == False:
            self._from_particles = True
        elif self._from_com == True:
            raise Exception('Cell list should be consistent in use of Compound center of mass or underlying particle positions, not mixing them.')

        
        if isinstance(compound, mb.Compound):
            for particle in compound.particles():
                c = self.cell_containing(particle.pos)
                if c < self._n_cells_total:
                    self.cells[c]._members.append(particle)
                    for neigh in self.cells[c]._neighbor_cells:
                        self.cells[neigh]._neighbor_members.append(particle)
                else:
                    msg = f'The particle is outside bounds of the box.\
                            \nposition: {particle.pos}\n\
                            box: {self.box.lengths}\n\
                            min box dimensions: {self.minbox}'
                    raise Exception(msg)

                    
    def insert_compound_position(self, compound):
        """This will insert an mbuild Compound into the cell list based upon the
        center-of-mass of the Compound (i.e., compound.pos).

        Parameters
        ----------
        compound :  mb.Compound
            An mbuild Compound that will be inserted into the cell list.
        
        Returns
        ------
        """
        if self._from_particles == False and self._from_com == False:
            self._from_com = True
        elif self._from_particles== True:
            raise Exception('Cell list should be consistent in use of Compound center of mass or underlying particle positions, not mixing them.')

        if isinstance(compound, mb.Compound):
            c = self.cell_containing(compound.pos)
            if c < self._n_cells_total:
                self.cells[c]._members.append(compound)
                for neigh in self.cells[c].neighbor_cells:
                    self.cells[neigh]._neighbor_members.append(compound)
            else:
                msg = f'The compound is outside bounds of the box.\
                        \nposition: {compound.pos}\n\
                        box: {self.box.lengths}\n\
                        min box dimensions: {self.minbox}'
                raise Exception(msg)
                
    def empty_cells(self):
        """Remove all members from the cell list.

        Parameters
        ----------

        Returns
        ------
        """
        for cell in self.cells:
            cell._members = []
            cell._neighbor_members = []
        #since it is empty we
        self._from_particles = False
        self._from_com = False
    def members(self, c):
        """Returns all members of a given cell .

        Parameters
        ----------
        c : int
            The cell of interest.
        Returns
        ------
        members : list, dtype=mb.Compound
            A list of all compounds that are within the cell.
        """
        if c < self._n_cells_total:
            return self.cells[c].members
        else:
            raise Exception(f'The cell requested {c} is out of bounds, total number of cells: {self._n_cells_total}')

    def neighbor_members(self, c):
        """Returns members of all neighboring cells.

        Parameters
        ----------
        c : int
            The cell of interest.
        Returns
        ------
        members : list, dtype=mb.Compound
            A list of all compounds that are within the cell.
        """
        if c < self._n_cells_total:
            return self.cells[c].neighbor_members
        else:
            raise Exception(f'The cell requested {c} is out of bounds, total number of cells: {self._n_cells_total}')

    @property
    def n_cells(self):
        """Returns a numpy array of the number of cells in each direction.
        Returns
        ------
        n_cells : np.array, dtype=int
            A numpy array of the number of cells in each x,y, and z direction.
        """
        return self._n_cells

    @property
    def n_cells_total(self):
        """Returns the total number of cells in each direction.
        Returns
        ------
        n_cells_total : int
            The total number of cells in the cell list.
        """
        return self._n_cells_total
        
    @property
    def periodicity(self):
        """Returns the periodicity in each direction.
        Returns
        ------
        periodicity : list, dtype=bool
            The total periodicity in each direction.
        """
        return self._periodicity

    @property
    def cell_sizes(self):
        """Returns a numpy array of the size of cells in each direction.
        Returns
        ------
        cell_sizes : np.array, dtype=float
            A numpy array of the size of the cells of cells in each x,y, and z direction.
        """
        return self._cell_sizes
