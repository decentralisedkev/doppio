use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};

use crate::affine::{AffineNielsPoint, AffinePoint};
use crate::completed::CompletedPoint;
use crate::curveconstants::{EDWARDS_D2, FR_MODULUS_BYTES};
use crate::fq::Fq;
use crate::fr::Fr;

/// This represents an extended point `(U, V, Z, T1, T2)`
/// with `Z` nonzero, corresponding to the affine point
/// `(U/Z, V/Z)`. We always have `T1 * T2 = UV/Z`.
///
/// You can do the following things with a point in this
/// form:
///
/// * Convert it into a point in the affine form.
/// * Add it to an `ExtendedPoint`, `AffineNielsPoint` or `ExtendedNielsPoint`.
/// * Double it using `double()`.
/// * Compare it with another extended point using `PartialEq` or `ct_eq()`.
#[derive(Clone, Copy, Debug)]
pub struct ExtendedPoint {
    pub(crate) u: Fq,
    pub(crate) v: Fq,
    pub(crate) z: Fq,
    pub(crate) t1: Fq,
    pub(crate) t2: Fq,
}

/// This is a pre-processed version of an extended point `(U, V, Z, T1, T2)`
/// in the form `(V + U, V - U, Z, T1 * T2 * 2d)`.
#[derive(Clone, Copy, Debug)]
pub struct ExtendedNielsPoint {
    pub(crate) v_plus_u: Fq,
    pub(crate) v_minus_u: Fq,
    pub(crate) z: Fq,
    pub(crate) t2d: Fq,
}

impl ConstantTimeEq for ExtendedPoint {
    fn ct_eq(&self, other: &Self) -> Choice {
        // (u/z, v/z) = (u'/z', v'/z') is implied by
        //      (uz'z = u'z'z) and
        //      (vz'z = v'z'z)
        // as z and z' are always nonzero.

        (&self.u * &other.z).ct_eq(&(&other.u * &self.z))
            & (&self.v * &other.z).ct_eq(&(&other.v * &self.z))
    }
}

impl ConditionallySelectable for ExtendedPoint {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        ExtendedPoint {
            u: Fq::conditional_select(&a.u, &b.u, choice),
            v: Fq::conditional_select(&a.v, &b.v, choice),
            z: Fq::conditional_select(&a.z, &b.z, choice),
            t1: Fq::conditional_select(&a.t1, &b.t1, choice),
            t2: Fq::conditional_select(&a.t2, &b.t2, choice),
        }
    }
}

impl PartialEq for ExtendedPoint {
    fn eq(&self, other: &Self) -> bool {
        self.ct_eq(other).unwrap_u8() == 1
    }
}

impl Default for ExtendedPoint {
    /// Returns the identity.
    fn default() -> ExtendedPoint {
        ExtendedPoint::identity()
    }
}

impl Neg for ExtendedPoint {
    type Output = ExtendedPoint;

    /// Computes the negation of a point `P = (U, V, Z, T)`
    /// as `-P = (-U, V, Z, -T1, T2)`. The choice of `T1`
    /// is made without loss of generality.
    #[inline]
    fn neg(self) -> ExtendedPoint {
        ExtendedPoint {
            u: -self.u,
            v: self.v,
            z: self.z,
            t1: -self.t1,
            t2: self.t2,
        }
    }
}

impl From<AffinePoint> for ExtendedPoint {
    /// Constructs an extended point (with `Z = 1`) from
    /// an affine point using the map `(u, v) => (u, v, 1, u, v)`.
    fn from(affine: AffinePoint) -> ExtendedPoint {
        ExtendedPoint {
            u: affine.u,
            v: affine.v,
            z: Fq::one(),
            t1: affine.u,
            t2: affine.v,
        }
    }
}

impl ExtendedPoint {
    /// Constructs an extended point from the neutral element `(0, 1)`.
    pub const fn identity() -> Self {
        ExtendedPoint {
            u: Fq::zero(),
            v: Fq::one(),
            z: Fq::one(),
            t1: Fq::zero(),
            t2: Fq::zero(),
        }
    }

    /// Determines if this point is the identity.
    pub fn is_identity(&self) -> Choice {
        // If this point is the identity, then
        //     u = 0 * z = 0
        // and v = 1 * z = z
        self.u.ct_eq(&Fq::zero()) & self.v.ct_eq(&self.z)
    }

    /// Determines if this point is of small order.
    pub fn is_small_order(&self) -> Choice {
        // We only need to perform two doublings, since the 2-torsion
        // points are (0, 1) and (0, -1), and so we only need to check
        // that the u-coordinate of the result is zero to see if the
        // point is small order.
        self.double().double().u.ct_eq(&Fq::zero())
    }

    /// Determines if this point is torsion free and so is contained
    /// in the prime order subgroup.
    pub fn is_torsion_free(&self) -> Choice {
        self.multiply(&FR_MODULUS_BYTES).is_identity()
    }

    /// Determines if this point is prime order, or in other words that
    /// the smallest scalar multiplied by this point that produces the
    /// identity is `r`. This is equivalent to checking that the point
    /// is both torsion free and not the identity.
    pub fn is_prime_order(&self) -> Choice {
        self.is_torsion_free() & (!self.is_identity())
    }

    /// Multiplies this element by the cofactor `8`.
    pub fn mul_by_cofactor(&self) -> ExtendedPoint {
        self.double().double().double()
    }

    /// Performs a pre-processing step that produces an `ExtendedNielsPoint`
    /// for use in multiple additions.
    pub fn to_niels(&self) -> ExtendedNielsPoint {
        ExtendedNielsPoint {
            v_plus_u: &self.v + &self.u,
            v_minus_u: &self.v - &self.u,
            z: self.z,
            t2d: &self.t1 * &self.t2 * EDWARDS_D2,
        }
    }

    /// Computes the doubling of a point more efficiently than a point can
    /// be added to itself.
    pub fn double(&self) -> ExtendedPoint {
        // Doubling is more efficient (three multiplications, four squarings)
        // when we work within the projective coordinate space (U:Z, V:Z). We
        // rely on the most efficient formula, "dbl-2008-bbjlp", as described
        // in Section 6 of "Twisted Edwards Curves" by Bernstein et al.
        //
        // See <https://hyperelliptic.org/EFD/g1p/auto-twisted-projective.html#doubling-dbl-2008-bbjlp>
        // for more information.
        //
        // We differ from the literature in that we use (u, v) rather than
        // (x, y) coordinates. We also have the constant `a = -1` implied. Let
        // us rewrite the procedure of doubling (u, v, z) to produce (U, V, Z)
        // as follows:
        //
        // B = (u + v)^2
        // C = u^2
        // D = v^2
        // F = D - C
        // H = 2 * z^2
        // J = F - H
        // U = (B - C - D) * J
        // V = F * (- C - D)
        // Z = F * J
        //
        // If we compute K = D + C, we can rewrite this:
        //
        // B = (u + v)^2
        // C = u^2
        // D = v^2
        // F = D - C
        // K = D + C
        // H = 2 * z^2
        // J = F - H
        // U = (B - K) * J
        // V = F * (-K)
        // Z = F * J
        //
        // In order to avoid the unnecessary negation of K,
        // we will negate J, transforming the result into
        // an equivalent point with a negated z-coordinate.
        //
        // B = (u + v)^2
        // C = u^2
        // D = v^2
        // F = D - C
        // K = D + C
        // H = 2 * z^2
        // J = H - F
        // U = (B - K) * J
        // V = F * K
        // Z = F * J
        //
        // Let us rename some variables to simplify:
        //
        // UV2 = (u + v)^2
        // UU = u^2
        // VV = v^2
        // VVmUU = VV - UU
        // VVpUU = VV + UU
        // ZZ2 = 2 * z^2
        // J = ZZ2 - VVmUU
        // U = (UV2 - VVpUU) * J
        // V = VVmUU * VVpUU
        // Z = VVmUU * J
        //
        // We wish to obtain two factors of T = UV/Z.
        //
        // UV/Z = (UV2 - VVpUU) * (ZZ2 - VVmUU) * VVmUU * VVpUU / VVmUU / (ZZ2 - VVmUU)
        //      = (UV2 - VVpUU) * VVmUU * VVpUU / VVmUU
        //      = (UV2 - VVpUU) * VVpUU
        //
        // and so we have that T1 = (UV2 - VVpUU) and T2 = VVpUU.

        let uu = self.u.square();
        let vv = self.v.square();
        let zz2 = self.z.square().double();
        let uv2 = (&self.u + &self.v).square();
        let vv_plus_uu = &vv + &uu;
        let vv_minus_uu = &vv - &uu;

        // The remaining arithmetic is exactly the process of converting
        // from a completed point to an extended point.
        CompletedPoint {
            u: &uv2 - &vv_plus_uu,
            v: vv_plus_uu,
            z: vv_minus_uu,
            t: &zz2 - &vv_minus_uu,
        }
        .into_extended()
    }

    #[inline]
    pub(crate) fn multiply(self, by: &[u8; 32]) -> Self {
        let zero = ExtendedPoint::identity().to_niels();
        let base = self.to_niels();

        let mut acc = ExtendedPoint::identity();

        // This is a simple double-and-add implementation of point
        // multiplication, moving from most significant to least
        // significant bit of the scalar.
        //
        // We skip the leading four bits because they're always
        // unset for Fr.
        for bit in by
            .iter()
            .rev()
            .flat_map(|byte| (0..8).rev().map(move |i| Choice::from((byte >> i) & 1u8)))
            .skip(4)
        {
            acc = acc.double();
            acc = acc + ExtendedNielsPoint::conditional_select(&zero, &base, bit);
        }

        acc
    }

    /// This is only for debugging purposes and not
    /// exposed in the public API. Checks that this
    /// point is on the curve.
    #[cfg(test)]
    pub(crate) fn is_on_curve_vartime(&self) -> bool {
        let affine = AffinePoint::from(*self);

        self.z != Fq::zero()
            && affine.is_on_curve_vartime()
            && (affine.u * affine.v * self.z == self.t1 * self.t2)
    }
}

impl<'a, 'b> Mul<&'b Fr> for &'a ExtendedPoint {
    type Output = ExtendedPoint;

    fn mul(self, other: &'b Fr) -> ExtendedPoint {
        self.multiply(&other.into_bytes())
    }
}

impl_binops_multiplicative!(ExtendedPoint, Fr);

impl<'a, 'b> Add<&'b ExtendedNielsPoint> for &'a ExtendedPoint {
    type Output = ExtendedPoint;

    fn add(self, other: &'b ExtendedNielsPoint) -> ExtendedPoint {
        // We perform addition in the extended coordinates. Here we use
        // a formula presented by Hisil, Wong, Carter and Dawson in
        // "Twisted Edward Curves Revisited" which only requires 8M.
        //
        // A = (V1 - U1) * (V2 - U2)
        // B = (V1 + U1) * (V2 + U2)
        // C = 2d * T1 * T2
        // D = 2 * Z1 * Z2
        // E = B - A
        // F = D - C
        // G = D + C
        // H = B + A
        // U3 = E * F
        // Y3 = G * H
        // Z3 = F * G
        // T3 = E * H

        let a = (&self.v - &self.u) * &other.v_minus_u;
        let b = (&self.v + &self.u) * &other.v_plus_u;
        let c = &self.t1 * &self.t2 * &other.t2d;
        let d = (&self.z * &other.z).double();

        // The remaining arithmetic is exactly the process of converting
        // from a completed point to an extended point.
        CompletedPoint {
            u: &b - &a,
            v: &b + &a,
            z: &d + &c,
            t: &d - &c,
        }
        .into_extended()
    }
}

impl<'a, 'b> Sub<&'b ExtendedNielsPoint> for &'a ExtendedPoint {
    type Output = ExtendedPoint;

    fn sub(self, other: &'b ExtendedNielsPoint) -> ExtendedPoint {
        let a = (&self.v - &self.u) * &other.v_plus_u;
        let b = (&self.v + &self.u) * &other.v_minus_u;
        let c = &self.t1 * &self.t2 * &other.t2d;
        let d = (&self.z * &other.z).double();

        CompletedPoint {
            u: &b - &a,
            v: &b + &a,
            z: &d - &c,
            t: &d + &c,
        }
        .into_extended()
    }
}

impl_binops_additive!(ExtendedPoint, ExtendedNielsPoint);

impl<'a, 'b> Add<&'b AffineNielsPoint> for &'a ExtendedPoint {
    type Output = ExtendedPoint;

    fn add(self, other: &'b AffineNielsPoint) -> ExtendedPoint {
        // This is identical to the addition formula for `ExtendedNielsPoint`,
        // except we can assume that `other.z` is one, so that we perform
        // 7 multiplications.

        let a = (&self.v - &self.u) * &other.v_minus_u;
        let b = (&self.v + &self.u) * &other.v_plus_u;
        let c = &self.t1 * &self.t2 * &other.t2d;
        let d = self.z.double();

        // The remaining arithmetic is exactly the process of converting
        // from a completed point to an extended point.
        CompletedPoint {
            u: &b - &a,
            v: &b + &a,
            z: &d + &c,
            t: &d - &c,
        }
        .into_extended()
    }
}

impl<'a, 'b> Sub<&'b AffineNielsPoint> for &'a ExtendedPoint {
    type Output = ExtendedPoint;

    fn sub(self, other: &'b AffineNielsPoint) -> ExtendedPoint {
        let a = (&self.v - &self.u) * &other.v_plus_u;
        let b = (&self.v + &self.u) * &other.v_minus_u;
        let c = &self.t1 * &self.t2 * &other.t2d;
        let d = self.z.double();

        CompletedPoint {
            u: &b - &a,
            v: &b + &a,
            z: &d - &c,
            t: &d + &c,
        }
        .into_extended()
    }
}

impl<'a, 'b> Add<&'b ExtendedPoint> for &'a ExtendedPoint {
    type Output = ExtendedPoint;

    #[inline]
    fn add(self, other: &'b ExtendedPoint) -> ExtendedPoint {
        self + other.to_niels()
    }
}

impl<'a, 'b> Sub<&'b ExtendedPoint> for &'a ExtendedPoint {
    type Output = ExtendedPoint;

    #[inline]
    fn sub(self, other: &'b ExtendedPoint) -> ExtendedPoint {
        self - other.to_niels()
    }
}

impl_binops_additive!(ExtendedPoint, ExtendedPoint);

impl<'a, 'b> Add<&'b AffinePoint> for &'a ExtendedPoint {
    type Output = ExtendedPoint;

    #[inline]
    fn add(self, other: &'b AffinePoint) -> ExtendedPoint {
        self + other.to_niels()
    }
}

impl<'a, 'b> Sub<&'b AffinePoint> for &'a ExtendedPoint {
    type Output = ExtendedPoint;

    #[inline]
    fn sub(self, other: &'b AffinePoint) -> ExtendedPoint {
        self - other.to_niels()
    }
}

impl ConditionallySelectable for ExtendedNielsPoint {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        ExtendedNielsPoint {
            v_plus_u: Fq::conditional_select(&a.v_plus_u, &b.v_plus_u, choice),
            v_minus_u: Fq::conditional_select(&a.v_minus_u, &b.v_minus_u, choice),
            z: Fq::conditional_select(&a.z, &b.z, choice),
            t2d: Fq::conditional_select(&a.t2d, &b.t2d, choice),
        }
    }
}

impl ExtendedNielsPoint {
    /// Constructs this point from the neutral element `(0, 1)`.
    pub const fn identity() -> Self {
        ExtendedNielsPoint {
            v_plus_u: Fq::one(),
            v_minus_u: Fq::one(),
            z: Fq::one(),
            t2d: Fq::zero(),
        }
    }
}
