// APBS vmgrid.rs - Multi-resolution grid hierarchy
// Port of src/mg/vmgrid.h / vmgrid.c

use crate::vgrid::Vgrid;
use apbs_generic::error::{ApbsError, ApbsResult};

/// Maximum number of grid levels in hierarchy
pub const VMGRIDMAX: usize = 20;

/// Multi-resolution grid hierarchy
pub struct Vmgrid {
    /// Number of grids in hierarchy
    pub ngrids: usize,
    /// Array of grids (finest first, by convention)
    pub grids: Vec<Vgrid>,
}

impl Vmgrid {
    pub fn new() -> Self {
        Self {
            ngrids: 0,
            grids: Vec::with_capacity(VMGRIDMAX),
        }
    }

    /// Add a grid to the hierarchy
    pub fn add_grid(&mut self, grid: Vgrid) -> ApbsResult<()> {
        if self.grids.len() >= VMGRIDMAX {
            return Err(ApbsError::Grid("Maximum grid levels exceeded".to_string()));
        }
        self.grids.push(grid);
        self.ngrids = self.grids.len();
        Ok(())
    }

    /// Get grid by level number
    pub fn get_grid_by_num(&self, num: usize) -> Option<&Vgrid> {
        self.grids.get(num)
    }

    /// Get the finest grid containing a point
    pub fn get_grid_by_point(&self, pt: [f64; 3]) -> Option<&Vgrid> {
        for grid in &self.grids {
            if pt[0] >= grid.xmin && pt[0] <= grid.xmax
                && pt[1] >= grid.ymin && pt[1] <= grid.ymax
                && pt[2] >= grid.zmin && pt[2] <= grid.zmax
            {
                return Some(grid);
            }
        }
        None
    }

    /// Interpolate value at a point
    pub fn value(&self, pt: [f64; 3]) -> ApbsResult<f64> {
        if let Some(grid) = self.get_grid_by_point(pt) {
            grid.value(pt)
        } else {
            Err(ApbsError::Grid("Point outside all grids".to_string()))
        }
    }
}

impl Default for Vmgrid {
    fn default() -> Self {
        Self::new()
    }
}
