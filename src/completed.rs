use crate::fq::Fq;
use crate::extended::ExtendedPoint;

/// This is a "completed" point produced during a point doubling or
/// addition routine. These points exist in the `(U:Z, V:T)` model
/// of the curve. This is not exposed in the API because it is
/// an implementation detail.
pub (crate) struct CompletedPoint {
    pub(crate)u: Fq,
    pub(crate)v: Fq,
    pub(crate)z: Fq,
    pub(crate)t: Fq,
}

impl CompletedPoint {
    /// This converts a completed point into an extended point by
    /// homogenizing:
    ///
    /// (u/z, v/t) = (u/z * t/t, v/t * z/z) = (ut/zt, vz/zt)
    ///
    /// The resulting T coordinate is utvz/zt = uv, and so
    /// T1 = u, T2 = v, without loss of generality.
    pub (crate) fn into_extended(&self) -> ExtendedPoint {
        ExtendedPoint {
            u: &self.u * &self.t,
            v: &self.v * &self.z,
            z: &self.z * &self.t,
            t1: self.u,
            t2: self.v,
        }
    }
}