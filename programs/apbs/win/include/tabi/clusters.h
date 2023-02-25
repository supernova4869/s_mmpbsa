#ifndef H_TABIPB_CLUSTERS_STRUCT_H
#define H_TABIPB_CLUSTERS_STRUCT_H

#include <cstddef>

#include "timer.h"
#include "particles.h"
#include "tree.h"
#include "params.h"

struct Timers_Clusters;

class Clusters
{
private:
    const class Particles& particles_;
    const class Tree& tree_;
    const struct Params& params_;
    struct Timers_Clusters& timers_;

    int num_interp_pts_per_node_;
    int num_charges_per_node_;
    
    std::size_t num_interp_pts_;
    std::size_t num_charges_;

    std::vector<double> interp_x_;
    std::vector<double> interp_y_;
    std::vector<double> interp_z_;
    
    std::vector<double> interp_charge_;
    std::vector<double> interp_charge_dx_;
    std::vector<double> interp_charge_dy_;
    std::vector<double> interp_charge_dz_;
    
    std::vector<double> interp_potential_;
    std::vector<double> interp_potential_dx_;
    std::vector<double> interp_potential_dy_;
    std::vector<double> interp_potential_dz_;
    
    
public:
    Clusters(const class Particles&, const class Tree&, const struct Params&, struct Timers_Clusters&);
    ~Clusters() = default;
    
    void upward_pass();
    void downward_pass(double* potential);
    
    void clear_charges();
    void clear_potentials();
    
    std::size_t num_interp_pts_per_node() const { return num_interp_pts_per_node_; };
    std::size_t num_charges_per_node()    const { return num_charges_per_node_; };
    const std::array<std::size_t, 2> cluster_interp_pts_idxs(std::size_t node_idx) const;
    const std::array<std::size_t, 2> cluster_charges_idxs(std::size_t node_idx) const;
    
    const double* interp_x_ptr() const { return interp_x_.data(); };
    const double* interp_y_ptr() const { return interp_y_.data(); };
    const double* interp_z_ptr() const { return interp_z_.data(); };
    
    const double* interp_charge_ptr()    const { return interp_charge_.data(); };
    const double* interp_charge_dx_ptr() const { return interp_charge_dx_.data(); };
    const double* interp_charge_dy_ptr() const { return interp_charge_dy_.data(); };
    const double* interp_charge_dz_ptr() const { return interp_charge_dz_.data(); };
    
    double* interp_potential_ptr()    { return interp_potential_.data(); };
    double* interp_potential_dx_ptr() { return interp_potential_dx_.data(); };
    double* interp_potential_dy_ptr() { return interp_potential_dy_.data(); };
    double* interp_potential_dz_ptr() { return interp_potential_dz_.data(); };
    
    void compute_all_interp_pts();
    void copyin_to_device() const;
    void delete_from_device() const;
    
};


struct Timers_Clusters
{
    Timer ctor;
    Timer upward_pass;
    Timer downward_pass;
    Timer clear_charges;
    Timer clear_potentials;
    Timer compute_all_interp_pts;
    Timer copyin_to_device;
    Timer delete_from_device;
    
    void print() const;
    std::string get_durations() const;
    std::string get_headers() const;

    Timers_Clusters() = default;
    ~Timers_Clusters() = default;
};

#endif /* H_TABIPB_CLUSTERS_STRUCT_H */
