#ifndef H_TABIPB_TIMERS_STRUCT_H
#define H_TABIPB_TIMERS_STRUCT_H

#include <iostream>
#include <iomanip>
#include <chrono>

#include "timer.h"
#include "molecule.h"
#include "particles.h"
#include "tree.h"
#include "clusters.h"
#include "interaction_list.h"
#include "boundary_element.h"

struct Timers
{
    Timers_Molecule molecule;
    Timers_Particles particles;
    Timers_Tree tree;
    Timers_Clusters clusters;
    Timers_InteractionList interaction_list;
    Timers_BoundaryElement boundary_element;

    Timer tabipb;

    Timers() = default;
    ~Timers() = default;


    void print() const
    {
        std::cout.setf(std::ios::fixed, std::ios::floatfield);
        std::cout.precision(5);
        
        std::cout << "|...Total TABIPB time (s)..........: ";
        std::cout << std::setw(12) << std::right << tabipb.elapsed_time() << std::endl;
        std::cout << "|" << std::endl;
        
        molecule         .print();
        particles        .print();
        tree             .print();
        clusters         .print();
        interaction_list .print();
        boundary_element .print();
    }


    std::string get_durations() const
    {
        std::string durations;
        durations.append(std::to_string(tabipb.elapsed_time())).append(", ");
        
        durations.append(molecule         .get_durations());
        durations.append(particles        .get_durations());
        durations.append(tree             .get_durations());
        durations.append(clusters         .get_durations());
        durations.append(interaction_list .get_durations());
        durations.append(boundary_element .get_durations());

        return durations;
    }


    std::string get_headers() const
    {
        std::string headers;
        headers.append("TABIPB Total time, ");
        
        headers.append(molecule         .get_headers());
        headers.append(particles        .get_headers());
        headers.append(tree             .get_headers());
        headers.append(clusters         .get_headers());
        headers.append(interaction_list .get_headers());
        headers.append(boundary_element .get_headers());

        return headers;
    }
};

#endif
