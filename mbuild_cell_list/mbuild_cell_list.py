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
        self._neighbor_cells_shift = {}
        self._ghost_cells = []
        
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
        
    @property
    def neighbor_cells_shift(self):
        """Returns a dictionary that defines how to shift the contents of a neighboring cell
        that exists across a periodic boundary, relative to this cell. The key of the dictionary
        corresponds to the numerical index of the neighboring cell."""
        return self._neighbor_cells_shift


class CellList():
    """Cell list compatible with mbuild Compounds.
    The cell list can be constructed based on either the center of mass of a Compound
    or based on the position of the particles contained within a Compound.
    """
    def __init__(self, box, n_cells=[3,3,3], periodicity=[True,True,True], box_min=[0.0,0.0,0.0]):
        """Initialize the cell list.
        Note this will initialize the full cell list where each cell has 26 neighbors when fully periodic.

        Parameters
        ----------
        box : list, length=3, dtype=float or mb.Box
            Either an mBuild Box or list of length=3 representing box lengths
        n_cells : list, length=3, dtype=int, default=[3,3,3]
            Number of cells in x,y,z dimensions, must be greater than 3
        periodicity, list, length=3, type=bool, default=[True,True,True]
            Periodicity in each box dimensions
        box_min, list, length=3, dtype=float, default=[0.0,0.0,0.0]
            Minimium position of the box.
        Returns
        ------
        """
        if isinstance(box, mb.Box):
            self._box = box
        else:
            assert len(box) == 3
            self._box = mb.Box(box)

        self._n_cells = np.array(n_cells,dtype=int)
        self._n_cells_total = np.prod(self._n_cells)
        

        if (self._n_cells < np.array([3,3,3]) ).any():
            raise Exception(f'The CellList must have at least 3 cells in each dimension, found: {n_cells}')
        
        self._box_min = np.array(box_min)
        
        self._cell_sizes = np.array(self._box.lengths)/self._n_cells
        
        self.cells = []
        self._periodicity = np.array(periodicity)
        
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
                    
                    
                    for kk, tmp in enumerate([i,j,k]):
                        
                        if not self._periodicity[kk]:
                            if tmp == 0:
                                start.append(0)
                            else:
                                start.append(-1)
                            if tmp == self._n_cells[kk]-1:
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
                                    
        for c, cell in enumerate(self.cells):
            for neigh in cell._neighbor_cells:
                dist = (self.cells[c].pos-self.cells[neigh].pos)/self._box.lengths
                flag = []
                for i in range(0,3):
                    flag.append(self._anint(dist[i]))
                cell.neighbor_cells_shift[neigh] = flag
    
    def cell_containing(self, xyz):
        """Return the cell that contains a given point in 3d space.

        Parameters
        ----------
        xyz :  np.ndarray, shape=(3), dtype=float
        
        Returns
        ------
        c : int
            The cell containing the point
        """
        if (np.array(xyz) <self._box_min).any():
            raise Exception('Particle outside bounds of the box.')
        if ((np.array(xyz) - self._box_min) > np.array(self._box.lengths)).any():
            raise Exception('Particle outside bounds of the box.')
            
        vals = np.array((np.array(xyz)-self._box_min)/self._cell_sizes, dtype=int)
        c = vals[0]+vals[1]*self._n_cells[0]+vals[2]*self._n_cells[0]*self._n_cells[1]
        
        return c
        
    def _shift(self, x):
        if x >= 1:
            return int(x)
        elif x < 0:
            return int(x-1)
        else:
            return 0
            
    def _anint(self, x):
        # function to construct periodic images
        if x >= 0.5:
            return 1
        elif x < -0.5:
            return -1
        else:
            return 0
            
    def _wrap_position(self, xyz):
        deltas = (np.array(xyz)-self._box_min)/np.array(self._box.lengths)
        xyz_shifted = np.array(xyz)
        for i, delta in enumerate(deltas):
            if self._periodicity[i]:
                xyz_shifted[i] = xyz[i] - self._shift(delta)*self._box.lengths[i]

        return xyz_shifted
        
    def _check_cell(self, c):
        # The cell_containing function  explicitly checks if a particle is within the
        # box bounds, and will raise an exception if we are outside the box.
        # We can actually have a particle that is out of bounds and still satisfy c < n_cell_total in some cases.
        # Checking to ensure that c < n_cells_total only ensure we don't have a seg-fault
        # This is likley unnecessary since cell_containing should catch an out of bounds particle before we check.
        if c < self._n_cells_total:
            return True
        else:
            raise Exception(f'Cell {c} is outside the bounds of the cell list.\n n_cell_total: {self._n_cells_total}')
            
    def insert_compound_particles(self, compound, wrap_pbc=False):
        """This will look at the lowest level of the hierarchy of an mbuild Compound
        (i.e., the particles) and insert them  into the cell list.

        Parameters
        ----------
        compound :  mb.Compound
            An mbuild Compound whose particles will be inserted into the cell list.
        wrap_pbc : bool, default=False
            If True, particle positions outside of the box bounds will be wrapped to the other side based on defined periodicity.
        Returns
        ------
        """
        if self._from_particles == False and self._from_com == False:
            self._from_particles = True
        elif self._from_com == True:
            raise Exception('Cell list should be consistent in use of Compound center of mass or underlying particle positions, not mixing them.')

        
        if isinstance(compound, mb.Compound):
            for particle in compound.particles():
                if wrap_pbc:
                    pos_shifted = self._wrap_position(particle.pos)
                    c = self.cell_containing(pos_shifted)
                else:
                    c = self.cell_containing(particle.pos)
                if self._check_cell(c):
                    self.cells[c]._members.append(particle)
                    for neigh in self.cells[c]._neighbor_cells:
                        self.cells[neigh]._neighbor_members.append((particle,c))
                   

                    
    def insert_compound_position(self, compound, wrap_pbc=False):
        """This will insert an mbuild Compound into the cell list based upon the
        center-of-mass of the Compound (i.e., compound.pos).

        Parameters
        ----------
        compound :  mb.Compound
            An mbuild Compound that will be inserted into the cell list.
        wrap_pbc : bool, default=False
            If True, particle positions outside of the box bounds will be wrapped to the other side based on defined periodicity.
        Returns
        ------
        """
        if self._from_particles == False and self._from_com == False:
            self._from_com = True
        elif self._from_particles== True:
            raise Exception('Cell list should be consistent in use of Compound center of mass or underlying particle positions, not mixing them.')

        if isinstance(compound, mb.Compound):
            if wrap_pbc:
                pos_shifted = self._wrap_position(compound.pos)
                c = self.cell_containing(pos_shifted)
            else:
                c = self.cell_containing(compound.pos)
            if self._check_cell(c):
                self.cells[c]._members.append(compound)
                for neigh in self.cells[c].neighbor_cells:
                    self.cells[neigh]._neighbor_members.append((compound, c))
                
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
        """Returns all members of a given cell.

        Parameters
        ----------
        c : int
            The cell of interest.
            
        Returns
        ------
        members : list, dtype=mb.Compound
            A list of all compounds that are within the cell.
        """
        if self._check_cell(c):
            return self.cells[c].members
 
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
        if self._check_cell(c):
            return [m[0] for m in self.cells[c].neighbor_members]
    
    def neighbor_members_and_min_image_shift(self, c):
        """Returns a list that contains members of all neighboring cells
        and how to shift those members to create a minimum image reconstruction
        relative to the cell of interest.

        Parameters
        ----------
        c : int
            The cell of interest.
            
        Returns
        ------
        (members, shift) : list, dtype=[mb.Compound, np.array shape=(3)]
            A list of all compounds that are within the cell.
        """
        if self._check_cell(c):
            tmp_list = []
            for neigh in self.cells[c].neighbor_members:
                shift = self.cells[c].neighbor_cells_shift[neigh[1]]
                tmp = [neigh[0], np.array(shift)]
                tmp_list.append(tmp)
            return tmp_list
    

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
        periodicity : np.array, dtype=bool
            The periodicity in each direction.
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
        
    @property
    def box(self):
        """Returns the box information used to initialize the cell list.
        Returns
        ------
        box : mb.Box
            An mbuild Box.
        """
        return self._box
