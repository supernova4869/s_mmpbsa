use std::collections::HashMap;

// get atom radius, returns 1.5 if not specified
pub fn get_radi(at_type: &str) -> f64 { // mBondi from AMBER20/parmed/tools/changeradii.py
    let rad_bondi: HashMap<&str, f64> = vec![("C", 1.7),
                                             ("H", 1.2),
                                             ("N", 1.55),
                                             ("HC", 1.3),
                                             ("O", 1.5),
                                             ("HN", 1.3),
                                             ("F", 1.5),
                                             ("HP", 1.3),
                                             ("SI", 2.1),
                                             ("HO", 0.8),
                                             ("P", 1.85),
                                             ("HS", 0.8),
                                             ("S", 1.8),
                                             ("CL", 1.7),
                                             ("BR", 1.85),
                                             ("I", 1.98)].into_iter().collect();
    let at_type = at_type.to_uppercase();
    let mut radius = 1.5;
    if at_type.len() >= 2 {
        let r = rad_bondi.get(&at_type[0..2]);
        if let Some(m) = r {
            radius = *m;
        } else {
            let r = rad_bondi.get(&at_type[0..1]);
            if let Some(m) = r {
                radius = *m;
            }
        }
    } else {
        let r = rad_bondi.get(&at_type[0..1]);
        if let Some(m) = r {
            radius = *m;
        }
    }
    return radius;
}