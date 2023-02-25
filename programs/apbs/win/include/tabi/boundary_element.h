#ifndef H_TABIPB_TREECODE_STRUCT_H
#define H_TABIPB_TREECODE_STRUCT_H

#include "timer.h"
#include "particles.h"
#include "clusters.h"
#include "interaction_list.h"

struct Timers_BoundaryElement;
struct Timers;

class BoundaryElement
{
private:
    class Particles& particles_;
    class Clusters& clusters_;
    const class Tree& tree_;
    const class InteractionList& interaction_list_;
    const class Molecule& molecule_;
    const struct Params& params_;
    struct Timers_BoundaryElement& timers_;
    
    std::vector<double> potential_;
    
    long int num_iter_;
    double residual_;
    
    double solvation_energy_;
    double free_energy_;
    double coulombic_energy_;
    
    double pot_min_;
    double pot_max_;
    double pot_normal_min_;
    double pot_normal_max_;
    
    int gmres_(long int n, const double* b, double* x, long int restrt,
               double* work, long int ldw, double *h, long int ldh,
               long int& iter, double& residual);
    
    void matrix_vector(double alpha, const double* __restrict potential_old,
                       double beta,        double* __restrict potential_new);
                       
    void precondition_diagonal(double* z, double* r);
    void precondition_block(double* z, double* r);
    
    void particle_particle_interact(double* __restrict potential,
                              const double* __restrict potential_old,
            std::array<std::size_t, 2> target_node_particle_idxs,
            std::array<std::size_t, 2> source_node_particle_idxs);
    
    void particle_cluster_interact(double* __restrict potential,
            std::array<std::size_t, 2> target_node_particle_idxs, std::size_t source_node_idx);
                                   
    void cluster_particle_interact(double* __restrict potential,
            std::size_t target_node_idx, std::array<std::size_t, 2> source_node_particle_idxs);
            
    void cluster_cluster_interact(double* __restrict potential,
            std::size_t target_node_idx, std::size_t source_node_idx);
    
public:
    BoundaryElement(class Particles& particles, class Clusters& clusters,
             const class Tree& tree, const class InteractionList& interaction_list,
             const class Molecule& molecule, const struct Params& params,
             struct Timers_BoundaryElement& timers);
    ~BoundaryElement() = default;
    
    void run_GMRES();
    void finalize();
    
    friend std::array<double, 3> Output(const BoundaryElement&, const Timers&);
};


struct Timers_BoundaryElement
{
    Timer ctor;
    Timer run_GMRES;
    Timer finalize;

    Timer matrix_vector;
    Timer precondition;
    Timer particle_particle_interact;
    Timer particle_cluster_interact;
    Timer cluster_particle_interact;
    Timer cluster_cluster_interact;

    void print() const;
    std::string get_durations() const;
    std::string get_headers() const;

    Timers_BoundaryElement() = default;
    ~Timers_BoundaryElement() = default;
};

#endif /* H_TABIPB_TREECODE_STRUCT_H */
