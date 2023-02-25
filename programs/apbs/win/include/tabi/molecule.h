#ifndef H_TABIPB_MOLECULE_STRUCT_H
#define H_TABIPB_MOLECULE_STRUCT_H

#include <vector>
#include <string>
#include <fstream>
#include <cstddef>

#include "timer.h"
#include "params.h"

#ifdef TABIPB_APBS
    #include "generic/valist.h"
#endif

struct Timers_Molecule;

class Molecule
{
private:
    const struct Params& params_;
    struct Timers_Molecule& timers_;
    
    std::size_t num_atoms_;
    std::vector<double> coords_;
    std::vector<double> charge_;
    std::vector<double> radius_;

    double coulombic_energy_;

public:
    Molecule(struct Params&, struct Timers_Molecule&);
    ~Molecule() = default;
    
#ifdef TABIPB_APBS
    Molecule(Valist*, struct Params&, struct Timers_Molecule&);
#endif
    
    void build_xyzr_file() const;
    
    std::size_t num_atoms() const { return num_atoms_; };
    double coulombic_energy() const { return coulombic_energy_; };
    const double* coords_ptr() const { return coords_.data(); };
    const double* charge_ptr() const { return charge_.data(); };
    const double* radius_ptr() const { return radius_.data(); };
    
    void compute_coulombic_energy();
    void copyin_to_device() const;
    void delete_from_device() const;
};


struct Timers_Molecule
{
    Timer ctor;
    Timer compute_coulombic_energy;
    Timer build_xyzr_file;
    Timer copyin_to_device;
    Timer delete_from_device;
    
    void print() const;
    std::string get_durations() const;
    std::string get_headers() const;

    Timers_Molecule() = default;
    ~Timers_Molecule() = default;
};

#endif /* H_MOLECULE_STRUCT_H */
