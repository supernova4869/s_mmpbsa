cmd.load("F:\Programming\Rust\Projects\s_mmpbsa\examples\docking\case1\PHL_7.pdbqt", "Protein")
cmd.load("F:\Programming\Rust\Projects\s_mmpbsa\examples\docking\case1\DSDP_out.pdbqt", "Ligand")
cmd.h_add("all")
cmd.save(r"F:\Programming\Rust\Projects\s_mmpbsa\examples\docking\case1\PHL_7.pdb", selection="(all)", state=0)
quit
