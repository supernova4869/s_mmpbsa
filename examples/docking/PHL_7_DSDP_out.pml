cmd.load(r"PHL_7.pdbqt", "Protein")
cmd.load(r"DSDP_out.pdbqt", "Ligand")
cmd.h_add("all")
cmd.save(r"\\?\F:\Programming\Rust\Projects\s_mmpbsa\examples\docking\LIG.mol2", selection="(Ligand)", state=1)
cmd.save(r"\\?\F:\Programming\Rust\Projects\s_mmpbsa\examples\docking\MMPBSA_PHL_7_DSDP_out.pdb", selection="(all)", state=0)
quit
