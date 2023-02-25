#ifndef H_TABIPB_TREE_STRUCT_H
#define H_TABIPB_TREE_STRUCT_H

#include <array>
#include <cstddef>

#include "timer.h"
#include "params.h"
#include "particles.h"

class InteractionList;
struct Timers_Tree;

class Tree
{
private:
    class Particles& particles_;
    const struct Params& params_;
    struct Timers_Tree& timers_;

    std::size_t num_nodes_;
    std::size_t num_leaves_;
    
    std::size_t min_leaf_size_;
    std::size_t max_leaf_size_;
    std::size_t max_depth_;
    
    std::vector<std::size_t> node_num_particles_;
    std::vector<std::size_t> node_particles_begin_;
    std::vector<std::size_t> node_particles_end_;
    
    std::vector<std::size_t> leaves_;
    
    std::vector<double> node_x_min_;
    std::vector<double> node_y_min_;
    std::vector<double> node_z_min_;
    
    std::vector<double> node_x_max_;
    std::vector<double> node_y_max_;
    std::vector<double> node_z_max_;
    
    std::vector<double> node_x_mid_;
    std::vector<double> node_y_mid_;
    std::vector<double> node_z_mid_;
    
    std::vector<double> node_radius_;
    std::vector<std::size_t> node_num_children_;
    std::vector<std::size_t> node_children_idx_;
    std::vector<std::size_t> node_parent_idx_;
    std::vector<std::size_t> node_level_;
    
    void construct(std::size_t, std::size_t, std::size_t, std::size_t);
    
public:
    Tree(class Particles&, const struct Params&, struct Timers_Tree&);
    ~Tree() = default;
    
    std::size_t num_nodes() const { return num_nodes_; };
    const std::array<double, 12> node_particle_bounds(std::size_t node_idx) const;
    const std::array<std::size_t, 2> node_particle_idxs(std::size_t node_idx) const;
    const std::vector<std::size_t>& leaves() const { return leaves_; }
    
    friend class InteractionList;
};


struct Timers_Tree
{
    Timer ctor;
    
    void print() const;
    std::string get_durations() const;
    std::string get_headers() const;

    Timers_Tree() = default;
    ~Timers_Tree() = default;
};

#endif /* H_STRUCT_TREE_H */
