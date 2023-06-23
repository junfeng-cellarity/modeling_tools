__author__ = 'jfeng1'
from openeye.oechem import *
import math,tempfile,shutil
from glide_utilities import MyTmpFile
import os
import subprocess

SCHRODINGER = os.getenv("SCHRODINGER")

def get_ligand_dimension(mol):
    cx = 0.0
    cy = 0.0
    cz = 0.0
    n = 0
    minx = None
    maxx = None
    miny = None
    maxy = None
    minz = None
    maxz = None
    for atm in mol.GetAtoms():
        # if not atm.IsHydrogen():
        x,y,z =  mol.GetCoords(atm)

        if minx is None or minx > x:
            minx = x
        if maxx is None or maxx < x:
            maxx = x

        if miny is None or miny > y:
            miny = y
        if maxy is None or maxy < y:
            maxy = y

        if minz is None or minz > z:
            minz = z
        if maxz is None or maxz < z:
            maxz = z

        cx += x
        cy += y
        cz += z
        n += 1

    xdim = (maxx-minx)/2
    ydim = (maxy-miny)/2
    zdim = (maxz-minz)/2

    max_dim = int(math.ceil(max(xdim,ydim,zdim)))

    cx/=float(n)
    cy/=float(n)
    cz/=float(n)

    return cx,cy,cz,max_dim


class GlideGrid:
    def __init__(self, ligand, gridfile, recep_file):
        # Each rigidly docked ligand or flexibly docked conformation has an associated length, L,
        # which can be defined as twice the distance from the ligand center to its farthest atom.
        # The required relationship between L and the lengths E and B of the outer and inner boxes for successful
        # placement of the ligand center anywhere within the inner box is:
        # E>=B+L

        self.cx,self.cy,self.cz, self.B = get_ligand_dimension(ligand)
        self.L = 2 * self.B
        self.E = self.B + self.L
        self.gridfile = gridfile
        self.recep_file = recep_file
        self.dict = {}

    def setKeyword(self,keyword, value):
        self.dict[keyword] = value

    def writeInput(self):
        dir = tempfile.mkdtemp(prefix="glide_grid")
        shutil.copy(self.recep_file,dir)
        os.chdir(dir)
        grid_input = "glide_grid.in"
        inputfile = open(grid_input,"w")
        inputfile.write("GRID_CENTER\t%5.3f,%5.3f,%5.3f\n"%(self.cx,self.cy,self.cz))
        inputfile.write("GRIDFILE\t%s\n"%self.gridfile)
        inputfile.write("INNERBOX\t%d,%d,%d\n"%(self.B,self.B,self.B))
        inputfile.write("OUTERBOX\t%d,%d,%d\n"%(self.E,self.E,self.E))
        inputfile.write("RECEP_FILE\t%s\n"%self.recep_file)
        for keyword in self.dict:
            inputfile.write("%s\t%s\n"%(keyword,self.dict[keyword]))
        inputfile.close()
        glide_command = os.path.join(SCHRODINGER,"glide")
        p = subprocess.Popen([glide_command,"-WAIT",grid_input],stdout=subprocess.PIPE)
        print p.communicate()[0]
        return

if __name__ == "__main__":
    ifs = oemolistream()
    ifs.open("ligand.sdf")
    mol = OEGraphMol()
    OEReadMolecule(ifs,mol)
    ifs.close()
    g = GlideGrid(mol,"grid.zip","receptor.maegz")
    g.writeInput()

