#ifndef H_TABIPB_CONSTANTS_H
#define H_TABIPB_CONSTANTS_H

namespace constants {
    constexpr double PI           = 3.14159265358979324;
    constexpr double ONE_OVER_4PI = 0.079577471545948;
    constexpr double KCAL_TO_KJ   = 4.184;
    constexpr double BULK_COEFF   = 2529.12179861515279;
    constexpr double UNITS_COEFF  = 1389.3875744; // 332.0716 * kcal2kj
    constexpr double UNITS_PARA   = 8729.779593448; // 2 * UNITS_COEFF * PI
}

#endif
