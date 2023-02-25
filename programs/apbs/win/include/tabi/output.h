#ifndef H_OUTPUT_H
#define H_OUTPUT_H

#include <array>

#include "boundary_element.h"
#include "tabipb_timers.h"

std::array<double, 3> Output(const BoundaryElement& bem, const Timers& timers);

#endif
