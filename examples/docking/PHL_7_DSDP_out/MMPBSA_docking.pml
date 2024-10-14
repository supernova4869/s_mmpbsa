cmd.load(r"../PHL_7.pdbqt", "Protein")
cmd.save(r"MMPBSA_docking_PHL_7.pdb", selection="(Protein)", state=1)
cmd.load(r"../DSDP_out.pdbqt", "Ligand")
cmd.h_add("all")
cmd.save(r"LIG.mol2", selection="(Ligand)", state=1)
cmd.save(r"MMPBSA_docking_DSDP_out.pdb", selection="(Ligand)", state=0)
quit
