// APBS vmatrix.rs - 2D/3D array indexing on flat buffers
// Port of src/generic/vmatrix.h

/// 2D array storage
pub struct Mat2 {
    data: Vec<f64>,
    dx: usize,
    dy: usize,
}

impl Mat2 {
    pub fn new(dx: usize, dy: usize) -> Self {
        Self {
            data: vec![0.0; dx * dy],
            dx,
            dy,
        }
    }

    pub fn with_data(dx: usize, dy: usize, data: Vec<f64>) -> Self {
        assert_eq!(data.len(), dx * dy);
        Self { data, dx, dy }
    }

    /// 1-based indexing (matching C VAT2 macro)
    #[inline(always)]
    pub fn get(&self, x: usize, y: usize) -> f64 {
        debug_assert!(x >= 1 && x <= self.dx && y >= 1 && y <= self.dy);
        self.data[(y - 1) * self.dx + (x - 1)]
    }

    /// 1-based mutable indexing
    #[inline(always)]
    pub fn get_mut(&mut self, x: usize, y: usize) -> &mut f64 {
        debug_assert!(x >= 1 && x <= self.dx && y >= 1 && y <= self.dy);
        &mut self.data[(y - 1) * self.dx + (x - 1)]
    }

    /// 0-based indexing
    #[inline(always)]
    pub fn get0(&self, x: usize, y: usize) -> f64 {
        debug_assert!(x < self.dx && y < self.dy);
        self.data[y * self.dx + x]
    }

    #[inline(always)]
    pub fn set0(&mut self, x: usize, y: usize, val: f64) {
        debug_assert!(x < self.dx && y < self.dy);
        self.data[y * self.dx + x] = val;
    }

    pub fn data(&self) -> &[f64] {
        &self.data
    }

    pub fn data_mut(&mut self) -> &mut [f64] {
        &mut self.data
    }

    pub fn dx(&self) -> usize {
        self.dx
    }

    pub fn dy(&self) -> usize {
        self.dy
    }
}

/// 3D array storage
pub struct Mat3 {
    data: Vec<f64>,
    dx: usize,
    dy: usize,
    dz: usize,
}

impl Mat3 {
    pub fn new(dx: usize, dy: usize, dz: usize) -> Self {
        Self {
            data: vec![0.0; dx * dy * dz],
            dx,
            dy,
            dz,
        }
    }

    pub fn with_data(dx: usize, dy: usize, dz: usize, data: Vec<f64>) -> Self {
        assert_eq!(data.len(), dx * dy * dz);
        Self { data, dx, dy, dz }
    }

    /// 1-based indexing (matching C VAT3 macro)
    /// mat[(z-1)*dy*dx + (y-1)*dx + (x-1)]
    #[inline(always)]
    pub fn get(&self, x: usize, y: usize, z: usize) -> f64 {
        debug_assert!(x >= 1 && x <= self.dx && y >= 1 && y <= self.dy && z >= 1 && z <= self.dz);
        self.data[(z - 1) * self.dy * self.dx + (y - 1) * self.dx + (x - 1)]
    }

    /// 1-based mutable indexing
    #[inline(always)]
    pub fn get_mut(&mut self, x: usize, y: usize, z: usize) -> &mut f64 {
        debug_assert!(x >= 1 && x <= self.dx && y >= 1 && y <= self.dy && z >= 1 && z <= self.dz);
        &mut self.data[(z - 1) * self.dy * self.dx + (y - 1) * self.dx + (x - 1)]
    }

    /// 0-based indexing
    #[inline(always)]
    pub fn get0(&self, x: usize, y: usize, z: usize) -> f64 {
        debug_assert!(x < self.dx && y < self.dy && z < self.dz);
        self.data[z * self.dy * self.dx + y * self.dx + x]
    }

    #[inline(always)]
    pub fn set0(&mut self, x: usize, y: usize, z: usize, val: f64) {
        debug_assert!(x < self.dx && y < self.dy && z < self.dz);
        self.data[z * self.dy * self.dx + y * self.dx + x] = val;
    }

    pub fn data(&self) -> &[f64] {
        &self.data
    }

    pub fn data_mut(&mut self) -> &mut [f64] {
        &mut self.data
    }

    pub fn dx(&self) -> usize {
        self.dx
    }

    pub fn dy(&self) -> usize {
        self.dy
    }

    pub fn dz(&self) -> usize {
        self.dz
    }
}
