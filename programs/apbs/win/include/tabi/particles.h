#ifndef H_TABIPB_PARTICLES_STRUCT_H
#define H_TABIPB_PARTICLES_STRUCT_H

#include <vector>
#include <cstdlib>

#include "timer.h"
#include "molecule.h"
#include "params.h"

struct Timers_Particles;

class Particles
{
private:
    const class Molecule& molecule_;
    const struct Params& params_;
    struct Timers_Particles& timers_;
    
    std::size_t num_;
    std::size_t num_faces_;
    double surface_area_;
    
    std::vector<std::size_t> face_x_;
    std::vector<std::size_t> face_y_;
    std::vector<std::size_t> face_z_;
    
    std::vector<double> x_;
    std::vector<double> y_;
    std::vector<double> z_;
    
    std::vector<double> nx_;
    std::vector<double> ny_;
    std::vector<double> nz_;
    
    std::vector<double> area_;
    std::vector<double> source_term_;
    
    std::vector<double> target_charge_;
    std::vector<double> target_charge_dx_;
    std::vector<double> target_charge_dy_;
    std::vector<double> target_charge_dz_;
    
    std::vector<double> source_charge_;
    std::vector<double> source_charge_dx_;
    std::vector<double> source_charge_dy_;
    std::vector<double> source_charge_dz_;

    std::vector<std::size_t> order_;
    
    void generate_particles(Params::Mesh, double, double);
    void update_source_term_on_host() const;
    
public:
    Particles(const class Molecule&, const struct Params&, struct Timers_Particles&);
    ~Particles() = default;
    
    int partition_8(std::size_t, std::size_t, std::array<std::size_t, 16>&);
    void reorder();
    void unorder(std::vector<double>& potential);
    
    void compute_charges(const double* potential);
    
    const std::array<double, 6> bounds(std::size_t begin, std::size_t end) const;
    double compute_solvation_energy(std::vector<double>& potential) const;
    
    void output_VTK(const std::vector<double>& potential) const;
    
    std::size_t num() const { return num_; };
    double surface_area() const { return surface_area_; };
    
    const double* x_ptr() const { return x_.data(); };
    const double* y_ptr() const { return y_.data(); };
    const double* z_ptr() const { return z_.data(); };
    
    const double* nx_ptr() const { return nx_.data(); };
    const double* ny_ptr() const { return ny_.data(); };
    const double* nz_ptr() const { return nz_.data(); };
    
    const double* area_ptr() const { return area_.data(); };
    const double* source_term_ptr() const { return source_term_.data(); };
    
    const double* target_charge_ptr()    const { return target_charge_.data(); };
    const double* target_charge_dx_ptr() const { return target_charge_dx_.data(); };
    const double* target_charge_dy_ptr() const { return target_charge_dy_.data(); };
    const double* target_charge_dz_ptr() const { return target_charge_dz_.data(); };
    
    const double* source_charge_ptr()    const { return source_charge_.data(); };
    const double* source_charge_dx_ptr() const { return source_charge_dx_.data(); };
    const double* source_charge_dy_ptr() const { return source_charge_dy_.data(); };
    const double* source_charge_dz_ptr() const { return source_charge_dz_.data(); };
    
    void compute_source_term();
    void copyin_to_device() const;
    void delete_from_device() const;
};


struct Timers_Particles
{
    Timer ctor;
    Timer compute_source_term;
    Timer compute_charges;
    Timer compute_solvation_energy;
    Timer copyin_to_device;
    Timer delete_from_device;
    Timer output_VTK;
    
    void print() const;
    std::string get_durations() const;
    std::string get_headers() const;

    Timers_Particles() = default;
    ~Timers_Particles() = default;
};

#endif /* H_PARTICLE_STRUCT_H */
