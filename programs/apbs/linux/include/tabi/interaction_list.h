#ifndef H_TABIPB_INTERACTION_LIST_STRUCT_H
#define H_TABIPB_INTERACTION_LIST_STRUCT_H

#include <cstddef>

#include "timer.h"
#include "tree.h"
#include "params.h"

struct Timers_InteractionList;

class InteractionList
{
private:
    const class Tree& tree_;
    const struct Params& params_;
    struct Timers_InteractionList& timers_;
    
    int size_check_;
    
    std::vector<std::vector<std::size_t>> particle_particle_;
    std::vector<std::vector<std::size_t>> particle_cluster_;
    std::vector<std::vector<std::size_t>> cluster_particle_;
    std::vector<std::vector<std::size_t>> cluster_cluster_;
    
    void build_BLTC_lists(std::size_t batch_idx, std::size_t node_idx);
    void build_BLDTT_lists(std::size_t target_node_idx, std::size_t source_node_idx);
    
public:
    InteractionList(const class Tree&, const struct Params&, struct Timers_InteractionList&);
    ~InteractionList() = default;
    
    const std::vector<std::size_t>& particle_particle(std::size_t idx) const { return particle_particle_[idx]; }
    const std::vector<std::size_t>& particle_cluster (std::size_t idx) const { return particle_cluster_ [idx]; }
    const std::vector<std::size_t>& cluster_particle (std::size_t idx) const { return cluster_particle_ [idx]; }
    const std::vector<std::size_t>& cluster_cluster  (std::size_t idx) const { return cluster_cluster_  [idx]; }
};


struct Timers_InteractionList
{
    Timer ctor;
    
    void print() const;
    std::string get_durations() const;
    std::string get_headers() const;

    Timers_InteractionList() = default;
    ~Timers_InteractionList() = default;
};

#endif /* H_TABIPB_INTERACTION_LIST_STRUCT_H */
