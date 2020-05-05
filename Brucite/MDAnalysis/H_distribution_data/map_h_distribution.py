import xml.etree.ElementTree as ET
import numpy as np

# User settings
###############
# Supercell size
NSUPX = 4
NSUPY = 4
NSUPZ = 2

# Grid size (from OUTCAR)
NGX=60
NGY=60
NGZ=42

# Number of atoms in supercell
NATOMS=160
# Types of atoms
ATTYPES = "Mg O H"
# Numbers of atom types
ATNUMS = "32 64 64"
# Index range of atoms to track
ATBEG=97
ATEND=160

# Equilibration period (in timesteps)
NEQUI=1000

# Number of runs, producing vasprun_X.xml files
NRUNS=1

# End of user input
###################


# Import vasprun_1.xml information
tree = ET.parse("vasprun_1.xml")
root = tree.getroot()

# Get unit cell and initial positions
latt_vec = []
pos_init = []
for entry in root.findall('structure'):
    if entry.get('name')=="initialpos":
        crystal_xml = entry.find('crystal')
        for c_entry in crystal_xml.findall('varray'):
            if c_entry.get('name')=="basis":
                for dir in range(3):
                    latt_vec.append(c_entry[dir].text)
        for p_entry in entry.findall('varray'):
            if p_entry.get('name')=="positions":
                for at in range(len(p_entry)):
                    pos_init.append(p_entry[at].text)

# Set up counting grid
counts = np.zeros((NGZ,NGY,NGX))
counts_prim = np.zeros((NGZ/NSUPZ,NGY/NSUPY,NGX/NSUPX))

# Add counts of atomic positions
nsteps = 0
for run in range(1,NRUNS+1):
    if run>1:
        tree = ET.parse("vasprun_"+str(run)+".xml")
        root = tree.getroot()
    for step in root.findall('calculation'):
        nsteps = nsteps + 1
        if nsteps >= NEQUI:
            positions = step.find('structure').find('varray')
            for atom in range(ATBEG-1,ATEND):
                atom_pos = positions[atom].text.split()
                ix = int(float(atom_pos[0])*NGX)%NGX
                iy = int(float(atom_pos[1])*NGY)%NGY
                iz = int(float(atom_pos[2])*NGZ)%NGZ
                counts[iz,iy,ix] = counts[iz,iy,ix] + 1

                px = ix%(NGX/NSUPX)
                py = iy%(NGY/NSUPY)
                pz = iz%(NGZ/NSUPZ)
                counts_prim[pz,py,px] = counts_prim[pz,py,px] + 1

# renormalise counts
counts = counts/nsteps
counts_prim = counts_prim/nsteps

counts_flattened = counts.flatten()
prim_flattened = counts_prim.flatten()

print "Number of geometry steps found: ",nsteps
print "Grid sites occupied throughout run: ",len(np.nonzero(counts_flattened)[0])

# write data to file
probcar = open("PROBCAR.vasp","w")
probcar.write("Probability distribution\n")
probcar.write("1.00\n")
probcar.write(latt_vec[0]+"\n")
probcar.write(latt_vec[1]+"\n")
probcar.write(latt_vec[2]+"\n")
probcar.write(ATTYPES+"\n")
probcar.write(ATNUMS+"\n")
probcar.write("Direct\n")
for at in pos_init:
    probcar.write(at+"\n")
probcar.write("\n")
probcar.write(str(NGX)+" "+str(NGY)+" "+str(NGZ)+"\n")
for i, pr in enumerate(counts_flattened):
    probcar.write(str(pr)+" ")
    if i%5==4: probcar.write("\n")
probcar.close()

# project back onto primitive unit cell
primcar = open("PROBPRIM.vasp","w")
primcar.write(str(NGX/NSUPX)+" "+str(NGY/NSUPY)+" "+str(NGZ/NSUPZ)+"\n")
for i, pr in enumerate(prim_flattened):
    primcar.write(str(pr)+" ")
    if i%5==4: primcar.write("\n")
primcar.close()
