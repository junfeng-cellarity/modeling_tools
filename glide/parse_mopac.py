#!/usr/bin/env python
from enum import Enum
import parse
class Mode(Enum):
    GENERAL = 1
    EIGENVECTORS = 2
    ATOMIC_ORBITAL = 3
    SAR=4

class MOPAC:
    def __init__(self,filename):
        if filename is not None:
            lines = open(filename,"r").readlines()
            mode = Mode.GENERAL
            self.eigen_values = []
            orbital_values = []
            vector_line_count = 0
            self.atom_symbols = {}
            self.sar_dict = {}
            orbital_dict = {} #key is atom_id: {1:{'S':[],'PX':[],'PY':[]},2:{} ...}
            self.atom_orbitals_dict = {}
            self.eigen_vectors = {}
            orbital_list = []
            for line_no,line in enumerate(lines):
                buffer = line.strip()
                if not len(buffer)==0:
                    if buffer.startswith("*"):
                        continue
                    if set(buffer)=={'-'}:
                        continue
                    if buffer.startswith("E_homo:"):
                        p = parse.parse("E_homo:{homo} eV",buffer)
                        self.homo = p['homo'].strip()
                    if buffer.startswith(("E_lumo:")):
                        p = parse.parse("E_lumo:{lumo} eV",buffer)
                        self.lumo = p['lumo'].strip()
                    if buffer.startswith("NO. OF FILLED LEVELS"):
                        p = parse.parse("NO. OF FILLED LEVELS    =         {num_level}",buffer)
                        self.num_level = int(p['num_level'])
                    # if buffer.startswith('MULLIKEN POPULATION ANALYSIS'):
                    #     mode = Mode.ATOMIC_ORBITAL
                    #     continue
                    if buffer.startswith("EIGENVECTORS"):
                        mode = Mode.EIGENVECTORS
                        orbital_list.clear()
                        continue
                    if buffer.startswith("Atom       FE        FN        SE        SN        SR       ASP"):
                        mode = Mode.SAR
                        continue
                    if buffer.startswith("Totals") and mode == Mode.SAR:
                        mode = Mode.GENERAL
                        continue
                    if mode == mode.SAR:
                        #FE: Electrophilic Frontier Electron Density
                        #FN: Nucleophilic Frontier Electron Density
                        #SE: Electrophilic Superdelocalizability
                        #SN: Electrophilic Superdelocalizability
                        #SR: Radical Superdelocalizability
                        #ASP:Atom Self Polarizability
                        p = parse.parse("{atom_symbol}   {atom_id}     {FE}     {FN}    {SE}    {SN}    {SR}    {ASP}",buffer)
                        if p:
                            atom_id = p['atom_id']
                            fe = p['FE']
                            fn = p['FN']
                            se = p['SE']
                            sn = p['SN']
                            sr = p['SR']
                            asp = p['ASP']
                            self.sar_dict[atom_id] = p.named
                            continue

                    if mode == Mode.EIGENVECTORS:
                        print(vector_line_count,buffer)
                        if buffer.startswith("Root No."):
                            for a in buffer.split()[2:]:
                                orbital_list.append(int(a))
                            vector_line_count = 1
                            continue
                        if vector_line_count == 1:
                            vector_line_count += 1
                            #skip this line
                            continue
                        if vector_line_count == 2:
                            #read eigenvalues
                            for eg in buffer.split():
                                self.eigen_values.append(float(eg))
                            vector_line_count += 1
                            continue
                        if vector_line_count >=3:
                            if buffer.startswith("ATOMIC PROPERTIES"):
                                mode = Mode.GENERAL
                                num_orbitals = len(self.eigen_values)
                                print(self.eigen_values)
                                print(self.eigen_vectors)
                            else:
                                print(buffer,self.eigen_values)
                                args = buffer.split()
                                print(args)
                                orbital_name = args[0]
                                atom_symbol =  args[1]
                                atom_id = args[2]
                                if atom_id not in self.eigen_vectors:
                                    self.eigen_vectors[atom_id] = {}
                                if orbital_name not in self.eigen_vectors[atom_id]:
                                    self.eigen_vectors[atom_id][orbital_name] = []
                                for idx,arg in enumerate(args[3:]):
                                    self.eigen_vectors[atom_id][orbital_name].append(float(arg))
                    if mode == Mode.ATOMIC_ORBITAL:
                        p = []
                        p.append(parse.parse("{orbital_1} {atom_1} {atom_id_1}", " ".join(buffer.split())))
                        p.append(parse.parse("{orbital_1} {atom_1} {atom_id_1} {orbital_2} {atom_2} {atom_id_2}", " ".join(buffer.split())))
                        p.append(parse.parse("{orbital_1} {atom_1} {atom_id_1} {orbital_2} {atom_2} {atom_id_2} {orbital_3} {atom_3} {atom_id_3}", " ".join(buffer.split())))
                        p.append(parse.parse("{orbital_1} {atom_1} {atom_id_1} {orbital_2} {atom_2} {atom_id_2} {orbital_3} {atom_3} {atom_id_3} {orbital_4} {atom_4} {atom_id_4}", " ".join(buffer.split())))
                        p.append(parse.parse("{orbital_1} {atom_1} {atom_id_1} {orbital_2} {atom_2} {atom_id_2} {orbital_3} {atom_3} {atom_id_3} {orbital_4} {atom_4} {atom_id_4} {orbital_5} {atom_5}", " ".join(buffer.split())))
                        p.append(parse.parse("{orbital_1} {atom_1} {atom_id_1} {orbital_2} {atom_2} {atom_id_2} {orbital_3} {atom_3} {atom_id_3} {orbital_4} {atom_4} {atom_id_4} {orbital_5} {atom_5} {atom_id_5} {orbital_6} {atom_6} {atom_id_6}", " ".join(buffer.split())))
                        if (p[0] or p[1] or p[2] or p[3] or p[4] or p[5]) and set(lines[line_no+1].strip())=={'-'}:
                            orbitals = []
                            atom_ids = []
                            atom_symbols = []
                            num_args = int(len(buffer.split())/3)
                            for i in range(1,num_args+1):
                                orbitals.append(p[num_args-1]['orbital_%d'%i])
                                atom_symbols.append(p[num_args-1]['atom_%d'%i])
                                atom_ids.append(p[num_args-1]['atom_id_%d'%i])

                            # print(orbital_dict)
                        else:
                            args = buffer.split()
                            if len(args) > 3:
                                orbital_name = args[0]
                                atom_symbol = args[1]
                                atom_id = args[2]
                                if atom_id not in orbital_dict:
                                    orbital_dict[atom_id] = []
                                if orbital_name not in orbital_dict[atom_id]:
                                    orbital_dict[atom_id].append(orbital_name)
                                if atom_id not in self.atom_symbols:
                                    self.atom_symbols[atom_id] = atom_symbol
                                for idx,arg in enumerate(args[3:]):
                                    # print(buffer)
                                    # print(idx,orbitals,atom_ids)
                                    key1 = "%s_%s-%s_%s"%(orbital_name,atom_id,orbitals[idx],atom_ids[idx])
                                    key2 = "%s_%s-%s_%s"%(orbitals[idx],atom_ids[idx],orbital_name,atom_id)
                                    self.atom_orbitals_dict[key1] = float(arg)
                                    if orbital_name == orbitals[idx] and atom_id == atom_ids[idx]:
                                        orbital_values.append(arg)
                                        if len(orbital_values) == num_orbitals:
                                            mode = Mode.GENERAL
                                    else:
                                        self.atom_orbitals_dict[key2] = float(arg)
                        continue

            # self.homo = eigen_values[int(num_orbitals/2)]
            # self.lumo = eigen_values[int(num_orbitals/2)+1]
            self.orbital_dict = orbital_dict
            self.num_orbitals = num_orbitals
            self.atom_ids = self.orbital_dict.keys()

    def calculate_superdelocalizability(self):
        mol_occupied_orbitals = self.eigen_values[0:self.num_level]
        mol_unoccupied_orbitals = self.eigen_values[self.num_level:]
        self.SEs = {}
        num_atoms = len(self.atom_ids)
        atom_ids = list(self.atom_symbols.keys())
        print(self.orbital_dict)

        for atom_id in self.atom_ids:
            self.SEs[atom_id] = 0.0
            orbital_dict = self.eigen_vectors[atom_id]
            for idx,eigen_value in enumerate(mol_occupied_orbitals):
                for orbital in orbital_dict.keys():
                    atomic_orbitals = orbital_dict[orbital]
                    self.SEs[atom_id] += atomic_orbitals[idx]*atomic_orbitals[idx]/eigen_value*2
                    print(atomic_orbitals[idx],eigen_value)
        print(self.SEs)
        #     atomic_orbitals_1 = self.orbital_dict[atom_id_1]
        #     sum_orbitals = []
        #     for ao_1 in atomic_orbitals_1:
        #         total = []
        #         for atom_id_2 in self.atom_ids:
        #             atomic_orbitals_2 = self.orbital_dict[atom_id_2]
        #             for ao_2 in atomic_orbitals_2:
        #                 key = "%s_%s-%s_%s" % (ao_1, atom_id_1, ao_2, atom_id_2)
        #                 total.append(self.atom_orbitals_dict[key])
        #         print(total,len(total))
        #         sum_orbitals.append(sum(total))
        #     #print(sum_orbitals)
        #
        #     for moo in mol_occupied_orbitals:
        #         for aoo in sum_orbitals:
        #             self.SEs[atom_id_1] += (2*aoo*aoo/moo)
        # print(self.SEs)
        # print(self.atom_symbols)
        # self.SNs = {}

if __name__ == "__main__":
    mopac = MOPAC("/home/jfeng/mopac.out")
    mopac.calculate_superdelocalizability()


