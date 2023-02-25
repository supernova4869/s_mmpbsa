#ifndef H_TABIPB_STRUCT_H
#define H_TABIPB_STRUCT_H

struct TABIPBInput {

    int mesh_flag_;
    double mesh_density_;
    double mesh_probe_radius_;
    
    double phys_temp_;
    double phys_eps_solute_;
    double phys_eps_solvent_;
    double phys_bulk_strength_;
    
    int tree_degree_;
    int tree_max_per_leaf_;
    double tree_theta_;
    
    int precondition_;
    int nonpolar_;
    
    int output_data_;

};

struct TABIPBOutput {

    double solvation_energy_;
    double coulombic_energy_;
    double free_energy_;

};

#endif
