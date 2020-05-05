'''
MD MSD analysis
'''

import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as et

class MSD(object):

    def __init__(self, FileName, check, filecheck, initpos):
        self.FileName = FileName
        self.check = check   # Controls sample rate
        self.filecheck = filecheck   # Used to ensure initial position from 1st file used
        self.init_posarray = initpos

        # - Test Vals - #
        # self.L = np.matrix([[5, 0, 0], [3, 4, 0], [0, 0, 4]])  # basis trans cart -> slant
        # self.posarray = np.array([[[0.25, 0.3, 0.15], [0.8, 0.6, 0.85]]]) #[1 iter][2 particles][3d pos]
        # self.posxy = self.posarray[:,:,:-1]
        # self.init_posarray = np.array([[0.2, 0.25, 0.1], [0.75, 0.6, 0.9]])
        # self.init_posxy = self.init_posarray[:,:-1]

    '''
    Function to pars the lattice translation matrix (convert Cartesian to 'Slant' basis)

    Sets self lattice and inverse lattice matrices
    '''

    def parsL(self):

        tree = et.parse(self.FileName)
        root = tree.getroot()

        L = []

        for v in root.find('structure/crystal/varray[@name="basis"]'):
            readin = v.text
            basevec = [float(n) for n in readin.split()]
            L.append(basevec)

        self.L = np.asmatrix(L)
        self.Linverse = np.linalg.inv(self.L)

    '''
    Function to pars positions of all atoms at the beginning.

    Sets self initial positions 2d numpy array of [atom][coordinates(slant basis)]
    and aslo sets xy coord array for closest image checker in plane
    '''

    def parsInitPos(self):

        if self.filecheck == 1:   # First file needs initial position read in, stored and returned to use in next files
            tree = et.parse(self.FileName)
            root = tree.getroot()

            init_posarray = []

            for v in root.find('structure[@name="initialpos"]/varray[@name="positions"]'):
                readin = v.text
                atompos = [float(n) for n in readin.split()]
                init_posarray.append(atompos)

            self.init_posarray = np.array(init_posarray)
            self.init_posxy = self.init_posarray[:,:-1] # xy coord for image checker

            return self.init_posarray

        else:   # If not 1st file, use 1st file passed inital positions
            self.init_posxy = self.init_posarray[:,:-1]
            return False



    '''
    Function to pars positions of all atoms at the each iteration.

    Sets self positions 3d numpy array of[iteration][atom][coordinates(slant basis)]
    Also sets xy coord array for image checker in plane
    '''

    def parsPos(self):

        tree = et.parse(self.FileName)
        root = tree.getroot()

        positions = []

        for varray in root.findall('calculation/structure/varray'):
            # - each step
            posset = []
            for v in varray.findall('v'):
                # - positions in step
                readin = v.text
                coord = [float(n) for n in readin.split()]
                posset.append(coord)
            positions.append(posset)

        self.posarray = np.array(positions)
        self.posxy = self.posarray[:,:,:-1] # for xy plane image checker

    '''
    Create in-plane (xy) images to iterte to find closest in the plane

    Returns array of images
    '''

    def create_images(self, pos):
        images = np.zeros((9,3))
        pos = np.append(pos,0)
        images[0] = pos
        # right
        im = pos
        im[0] = im[0] + 1
        images[1] = im
        # diagonal upper right
        im[1] = im[1] + 1
        images[2] = im
        # upper
        im[0] = im[0] - 1
        images[3] = im
        # diagonal upper left
        im[0] = im[0] - 1
        images[4] = im
        # left
        im[1] = im[1] - 1
        images[5] = im
        # diagonal lower left
        im[1] = im[1] - 1
        images[6] = im
        # lower
        im[0] = im[0] + 1
        images[7] = im
        # diagonal lower right
        im[0] = im[0] + 1
        images[8] = im

        return images

    '''
    Determines the closest in plane particles using brute force image creation and checking
    the displacement squared. Only uses the xy components as z is orthogonal and can be done
    by modulus method.

    Returns the displacement squared of the xy coordinates (r^2 = x^2 + y^2; returning r^2)
    having used minimum image convention
    '''

    def closest_inplane(self, pos1, pos2):
        L_trans = np.transpose(self.L)

        pos1 = np.append(pos1,0) # ensure correct dimensions
        pos1 = np.ndarray.transpose(pos1) # transpose
        pos1 = L_trans * pos1[:, None] # cartesian
        pos1 = np.ndarray.transpose(pos1) # transpose back

        images = self.create_images(pos2) # gets image array
        images = np.ndarray.transpose(images) # transpose
        images = L_trans * images # cartesian
        images = np.ndarray.transpose(images) # transpose back

        min_displsq_xy = 10000
        # High value to ensure 1st calculated is always lower

        for i in range(len(images)):
            r12 = images[i] - pos1 # displacement of two atoms
            displsq_xy = (np.linalg.norm(r12))**2.0 # displ squared
            if displsq_xy < min_displsq_xy: # if lower than current min
                min_displsq_xy = displsq_xy # set to new min

        return min_displsq_xy

    '''
    Finds z component squared in MSD, by modulus method using
    orthogonality

    Returns z displacement squared
    '''

    def closest_z(self, step, i):
        closest_z_displ = np.mod((self.init_posarray[i][2]-self.posarray[step][i][2]+1.5), 1.0) - 0.5
        # Modulus method, 1.5, 1.0 and 0.5 are relative values as the
        # positions are stored by vasp between 0 and 1, with self.L being
        # the basis set for this.
        closest_z_displ_sq = (closest_z_displ)**2.0
        closest_z_displ_sq = ((self.L.item(2,2))**2.0) * closest_z_displ_sq
        # Converts to Euclidian distance

        return closest_z_displ_sq

    '''
    Finds the MSD of the the atoms, with same atoms of each type.

    NOTE: Currently the atom indices are hard coded; TO BE UPDATED

    Returns 2d array of the normalised MSDs at each iteration for each atom type
    '''

    def find_MSD(self, atoms):
        MSDs_norm = []   # array intiialisations
        MSDs_Mgnorm = []
        MSDs_Onorm = []
        MSDs_Hnorm = []
        # over all steps
        for step in range(len(self.posarray)):
            MSDstep_Mg = 0   # each step array initialisations
            MSDstep_O = 0
            MSDstep_H = 0
            #for step % certain number -> not every step finding MSD
            if ((step%self.check) == 0):
                # over all Mg atoms
                for i in range(0,atoms[0]):
                    closest_xy_displ_sq = self.closest_inplane(self.init_posxy[i], self.posxy[step][i]) # xy closest displ squared
                    closest_z_displ_sq = self.closest_z(step, i) # z closest displ squared
                    MSDstep_Mg = MSDstep_Mg + ((closest_xy_displ_sq + closest_z_displ_sq)/32.0) # total MSD for step counter
                # Over all O atoms
                for i in range(atoms[0],atoms[1]):
                    closest_xy_displ_sq = self.closest_inplane(self.init_posxy[i], self.posxy[step][i])
                    closest_z_displ_sq = self.closest_z(step, i) # z closest displ squared
                    MSDstep_O = MSDstep_O + ((closest_xy_displ_sq + closest_z_displ_sq)/64.0)
                # Over all H atoms
                for i in range(atoms[1],atoms[2]):
                    closest_xy_displ_sq = self.closest_inplane(self.init_posxy[i], self.posxy[step][i])
                    closest_z_displ_sq = self.closest_z(step, i) # z closest displ squared
                    MSDstep_H = MSDstep_H + ((closest_xy_displ_sq + closest_z_displ_sq)/64.0)
                # Appending results to relevant arrays
                MSDs_Mgnorm.append(MSDstep_Mg)
                MSDs_Onorm.append(MSDstep_O)
                MSDs_Hnorm.append(MSDstep_H)
        # Appending results to single array for returning.
        MSDs_norm.append(MSDs_Mgnorm)
        MSDs_norm.append(MSDs_Onorm)
        MSDs_norm.append(MSDs_Hnorm)

        return MSDs_norm

'''
MAIN FUNC
'''

def main():

    check = 10   # set to alter number of MSD checks; 100 means 100 MSD calcs in 10000 iterations
    FileName1 = 'vasprun-1.xml'   # file to check
    FileName2 = 'vasprun-2.xml'
    FileName3 = 'vasprun-3.xml'
    atoms = [32, 96, 160]   # [1-32 Mg, 33-96 O, 97-160 H]

    initpos = [0,0,0]
    # Initialise
    MSDs1 = MSD(FileName1, check, 1, initpos)
    MSDs1.parsL()
    initpos = MSDs1.parsInitPos()
    MSDs1.parsPos()

    MSDs2 = MSD(FileName2, check, 2, initpos)
    MSDs2.parsL()
    unpos = MSDs2.parsInitPos()
    MSDs2.parsPos()

    MSDs3 = MSD(FileName3, check, 3, initpos)
    MSDs3.parsL()
    unpos = MSDs3.parsInitPos()
    MSDs3.parsPos()

    # MSD calculation
    MSDs1_norm = MSDs1.find_MSD(atoms)
    MSDs2_norm = MSDs2.find_MSD(atoms)
    MSDs3_norm = MSDs3.find_MSD(atoms)

    MSDs_norm = np.hstack((MSDs1_norm, MSDs2_norm, MSDs3_norm))

    # Storing data in files to allow analysis without regathering.
    fileMg = open('p3m1_msd_Mg.txt', 'a')
    for i in range(len(MSDs_norm[0])):
        fileMg.write(str(MSDs_norm[0][i])+'\n')
    fileMg.close()

    fileO = open('p3m1_msd_O.txt', 'a')
    for i in range(len(MSDs_norm[1])):
        fileO.write(str(MSDs_norm[1][i])+'\n')
    fileO.close()

    fileH = open('p3m1_msd_H.txt', 'a')
    for i in range(len(MSDs_norm[2])):
        fileH.write(str(MSDs_norm[2][i])+'\n')
    fileH.close()

main()
