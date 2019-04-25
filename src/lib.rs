//! This crate provides an implementation of the **Jubjub** elliptic curve and its associated
//! field arithmetic. See [`README.md`](https://github.com/zkcrypto/jubjub/blob/master/README.md) for more details about Jubjub.
//!
//! # API
//!
//! * `AffinePoint` / `ExtendedPoint` which are implementations of Jubjub group arithmetic
//! * `AffineNielsPoint` / `ExtendedNielsPoint` which are pre-processed Jubjub points
//! * `Fq`, which is the base field of Jubjub
//! * `Fr`, which is the scalar field of Jubjub
//! * `batch_normalize` for converting many `ExtendedPoint`s into `AffinePoint`s efficiently.
//!
//! # Constant Time
//!
//! All operations are constant time unless explicitly noted; these functions will contain
//! "vartime" in their name and they will be documented as variable time.
//!
//! This crate relies on the `subtle` crate for achieving constant time arithmetic. It is
//! recommended to enable the `nightly` feature on this crate (which enables the `nightly`
//! feature in the `subtle` crate) to defend against compiler optimizations that may
//! compromise constant time arithmetic. However, this requires use of the nightly version
//! of the Rust compiler.
//!
//! # Features
//!
//! * `nightly`: This enables `subtle/nightly` which attempts to prevent the compiler from
//! performing optimizations that could compromise constant time arithmetic. It is
//! recommended to enable this if you are able to use a nightly version of the Rust compiler.

#![no_std]
#![deny(missing_debug_implementations)]
#![deny(missing_docs)]
#![deny(unsafe_code)]

#[cfg(feature = "std")]
#[macro_use]
extern crate std;

use core::ops::{Add, AddAssign, Sub, SubAssign};

#[macro_use]
mod util;

mod ctoption;
pub use ctoption::CtOption;

mod fqconstants;
mod frconstants;
mod curveconstants;
pub use curveconstants::*;

mod fq;
mod fr;
pub use fq::Fq;
pub use fr::Fr;

mod affine;
pub use affine::{AffinePoint, AffineNielsPoint};
mod extended;
pub use extended::{ExtendedPoint, ExtendedNielsPoint};
mod completed;


impl_binops_additive!(ExtendedPoint, AffineNielsPoint);

impl_binops_additive!(ExtendedPoint, AffinePoint);

/// This takes a mutable slice of `ExtendedPoint`s and "normalizes" them using
/// only a single inversion for the entire batch. This normalization results in
/// all of the points having a Z-coordinate of one. Further, an iterator is
/// returned which can be used to obtain `AffinePoint`s for each element in the
/// slice.
///
/// This costs 5 multiplications per element, and a field inversion.
pub fn batch_normalize<'a>(v: &'a mut [ExtendedPoint]) -> impl Iterator<Item = AffinePoint> + 'a {
    let mut acc = Fq::one();
    for p in v.iter_mut() {
        // We use the `t1` field of `ExtendedPoint` to store the product
        // of previous z-coordinates seen.
        p.t1 = acc;
        acc *= &p.z;
    }

    // This is the inverse, as all z-coordinates are nonzero.
    acc = acc.invert().unwrap();

    for p in v.iter_mut().rev() {
        let mut q = *p;

        // Compute tmp = 1/z
        let tmp = q.t1 * acc;

        // Cancel out z-coordinate in denominator of `acc`
        acc *= &q.z;

        // Set the coordinates to the correct value
        q.u *= &tmp; // Multiply by 1/z
        q.v *= &tmp; // Multiply by 1/z
        q.z = Fq::one(); // z-coordinate is now one
        q.t1 = q.u;
        q.t2 = q.v;

        *p = q;
    }

    // All extended points are now normalized, but the type
    // doesn't encode this fact. Let us offer affine points
    // to the caller.

    v.iter().map(|p| AffinePoint { u: p.u, v: p.v })
}

#[test]
fn test_is_on_curve_var() {
    assert!(AffinePoint::identity().is_on_curve_vartime());
}

#[test]
fn test_d_is_non_quadratic_residue() {
    assert!(EDWARDS_D.sqrt().is_none().unwrap_u8() == 1);
    assert!((-EDWARDS_D).sqrt().is_none().unwrap_u8() == 1);
    assert!((-EDWARDS_D).invert().unwrap().sqrt().is_none().unwrap_u8() == 1);
}

#[test]
fn test_affine_niels_point_identity() {
    assert_eq!(
        AffineNielsPoint::identity().v_plus_u,
        AffinePoint::identity().to_niels().v_plus_u
    );
    assert_eq!(
        AffineNielsPoint::identity().v_minus_u,
        AffinePoint::identity().to_niels().v_minus_u
    );
    assert_eq!(
        AffineNielsPoint::identity().t2d,
        AffinePoint::identity().to_niels().t2d
    );
}

#[test]
fn test_extended_niels_point_identity() {
    assert_eq!(
        ExtendedNielsPoint::identity().v_plus_u,
        ExtendedPoint::identity().to_niels().v_plus_u
    );
    assert_eq!(
        ExtendedNielsPoint::identity().v_minus_u,
        ExtendedPoint::identity().to_niels().v_minus_u
    );
    assert_eq!(
        ExtendedNielsPoint::identity().z,
        ExtendedPoint::identity().to_niels().z
    );
    assert_eq!(
        ExtendedNielsPoint::identity().t2d,
        ExtendedPoint::identity().to_niels().t2d
    );
}

#[test]
fn test_assoc() {
    let p = ExtendedPoint::from(AffinePoint {
        u: Fq([
            0xc0115cb656ae4839,
            0x623dc3ff81d64c26,
            0x5868e739b5794f2c,
            0x23bd4fbb18d39c9c,
        ]),
        v: Fq([
            0x7588ee6d6dd40deb,
            0x9d6d7a23ebdb7c4c,
            0x46462e26d4edb8c7,
            0x10b4c1517ca82e9b,
        ]),
    }).mul_by_cofactor();
    assert!(p.is_on_curve_vartime());

    assert_eq!(
        (p * Fr::from(1000u64)) * Fr::from(3938u64),
        p * (Fr::from(1000u64) * Fr::from(3938u64)),
    );
}

#[cfg(feature = "std")]
#[test]
fn test_batch_normalize() {
    let mut p = ExtendedPoint::from(AffinePoint {
        u: Fq([
            0xc0115cb656ae4839,
            0x623dc3ff81d64c26,
            0x5868e739b5794f2c,
            0x23bd4fbb18d39c9c,
        ]),
        v: Fq([
            0x7588ee6d6dd40deb,
            0x9d6d7a23ebdb7c4c,
            0x46462e26d4edb8c7,
            0x10b4c1517ca82e9b,
        ]),
    }).mul_by_cofactor();

    let mut v = vec![];
    for _ in 0..10 {
        v.push(p);
        p = p.double();
    }

    for p in &v {
        assert!(p.is_on_curve_vartime());
    }

    let expected: std::vec::Vec<_> = v.iter().map(|p| AffinePoint::from(*p)).collect();
    let result1: std::vec::Vec<_> = batch_normalize(&mut v).collect();
    for i in 0..10 {
        assert!(expected[i] == result1[i]);
        assert!(v[i].is_on_curve_vartime());
        assert!(AffinePoint::from(v[i]) == expected[i]);
    }
    let result2: std::vec::Vec<_> = batch_normalize(&mut v).collect();
    for i in 0..10 {
        assert!(expected[i] == result2[i]);
        assert!(v[i].is_on_curve_vartime());
        assert!(AffinePoint::from(v[i]) == expected[i]);
    }
}

#[cfg(test)]
const FULL_GENERATOR: AffinePoint = AffinePoint::from_raw_unchecked(
    Fq::from_raw([
        0xe4b3d35df1a7adfe,
        0xcaf55d1b29bf81af,
        0x8b0f03ddd60a8187,
        0x62edcbb8bf3787c8,
    ]),
    Fq::from_raw([0xb, 0x0, 0x0, 0x0]),
);

#[cfg(test)]
const EIGHT_TORSION: [AffinePoint; 8] = [
    AffinePoint::from_raw_unchecked(
        Fq::from_raw([
            0xd92e6a7927200d43,
            0x7aa41ac43dae8582,
            0xeaaae086a16618d1,
            0x71d4df38ba9e7973,
        ]),
        Fq::from_raw([
            0xff0d2068eff496dd,
            0x9106ee90f384a4a1,
            0x16a13035ad4d7266,
            0x4958bdb21966982e,
        ]),
    ),
    AffinePoint::from_raw_unchecked(
        Fq::from_raw([
            0xfffeffff00000001,
            0x67baa40089fb5bfe,
            0xa5e80b39939ed334,
            0x73eda753299d7d47,
        ]),
        Fq::from_raw([0x0, 0x0, 0x0, 0x0]),
    ),
    AffinePoint::from_raw_unchecked(
        Fq::from_raw([
            0xd92e6a7927200d43,
            0x7aa41ac43dae8582,
            0xeaaae086a16618d1,
            0x71d4df38ba9e7973,
        ]),
        Fq::from_raw([
            0xf2df96100b6924,
            0xc2b6b5720c79b75d,
            0x1c98a7d25c54659e,
            0x2a94e9a11036e51a,
        ]),
    ),
    AffinePoint::from_raw_unchecked(
        Fq::from_raw([0x0, 0x0, 0x0, 0x0]),
        Fq::from_raw([
            0xffffffff00000000,
            0x53bda402fffe5bfe,
            0x3339d80809a1d805,
            0x73eda753299d7d48,
        ]),
    ),
    AffinePoint::from_raw_unchecked(
        Fq::from_raw([
            0x26d19585d8dff2be,
            0xd919893ec24fd67c,
            0x488ef781683bbf33,
            0x218c81a6eff03d4,
        ]),
        Fq::from_raw([
            0xf2df96100b6924,
            0xc2b6b5720c79b75d,
            0x1c98a7d25c54659e,
            0x2a94e9a11036e51a,
        ]),
    ),
    AffinePoint::from_raw_unchecked(
        Fq::from_raw([0x1000000000000, 0xec03000276030000, 0x8d51ccce760304d0, 0x0]),
        Fq::from_raw([0x0, 0x0, 0x0, 0x0]),
    ),
    AffinePoint::from_raw_unchecked(
        Fq::from_raw([
            0x26d19585d8dff2be,
            0xd919893ec24fd67c,
            0x488ef781683bbf33,
            0x218c81a6eff03d4,
        ]),
        Fq::from_raw([
            0xff0d2068eff496dd,
            0x9106ee90f384a4a1,
            0x16a13035ad4d7266,
            0x4958bdb21966982e,
        ]),
    ),
    AffinePoint::from_raw_unchecked(
        Fq::from_raw([0x0, 0x0, 0x0, 0x0]),
        Fq::from_raw([0x1, 0x0, 0x0, 0x0]),
    ),
];

#[test]
fn find_eight_torsion() {
    let g = ExtendedPoint::from(FULL_GENERATOR);
    assert!(g.is_small_order().unwrap_u8() == 0);
    let g = g.multiply(&FR_MODULUS_BYTES);
    assert!(g.is_small_order().unwrap_u8() == 1);

    let mut cur = g;

    for (i, point) in EIGHT_TORSION.iter().enumerate() {
        let tmp = AffinePoint::from(cur);
        if &tmp != point {
            panic!("{}th torsion point should be {:?}", i, tmp);
        }

        cur += &g;
    }
}

#[test]
fn find_curve_generator() {
    let mut trial_bytes = [0; 32];
    for _ in 0..255 {
        let a = AffinePoint::from_bytes(trial_bytes);
        if a.is_some().unwrap_u8() == 1 {
            let a = a.unwrap();
            assert!(a.is_on_curve_vartime());
            let b = ExtendedPoint::from(a);
            let b = b.multiply(&FR_MODULUS_BYTES);
            assert!(b.is_small_order().unwrap_u8() == 1);
            let b = b.double();
            assert!(b.is_small_order().unwrap_u8() == 1);
            let b = b.double();
            assert!(b.is_small_order().unwrap_u8() == 1);
            if b.is_identity().unwrap_u8() == 0 {
                let b = b.double();
                assert!(b.is_small_order().unwrap_u8() == 1);
                assert!(b.is_identity().unwrap_u8() == 1);
                assert_eq!(FULL_GENERATOR, a);
                assert!(a.mul_by_cofactor().is_torsion_free().unwrap_u8() == 1);
                return;
            }
        }

        trial_bytes[0] += 1;
    }

    panic!("should have found a generator of the curve");
}

#[test]
fn test_small_order() {
    for point in EIGHT_TORSION.iter() {
        assert!(point.is_small_order().unwrap_u8() == 1);
    }
}

#[test]
fn test_is_identity() {
    let a = EIGHT_TORSION[0].mul_by_cofactor();
    let b = EIGHT_TORSION[1].mul_by_cofactor();

    assert_eq!(a.u, b.u);
    assert_eq!(a.v, a.z);
    assert_eq!(b.v, b.z);
    assert!(a.v != b.v);
    assert!(a.z != b.z);

    assert!(a.is_identity().unwrap_u8() == 1);
    assert!(b.is_identity().unwrap_u8() == 1);

    for point in EIGHT_TORSION.iter() {
        assert!(point.mul_by_cofactor().is_identity().unwrap_u8() == 1);
    }
}

#[test]
fn test_mul_consistency() {
    let a = Fr([
        0x21e61211d9934f2e,
        0xa52c058a693c3e07,
        0x9ccb77bfb12d6360,
        0x07df2470ec94398e,
    ]);
    let b = Fr([
        0x03336d1cbe19dbe0,
        0x0153618f6156a536,
        0x2604c9e1fc3c6b15,
        0x04ae581ceb028720,
    ]);
    let c = Fr([
        0xd7abf5bb24683f4c,
        0x9d7712cc274b7c03,
        0x973293db9683789f,
        0x0b677e29380a97a7,
    ]);
    assert_eq!(a * b, c);
    let p = ExtendedPoint::from(AffinePoint {
        u: Fq([
            0xc0115cb656ae4839,
            0x623dc3ff81d64c26,
            0x5868e739b5794f2c,
            0x23bd4fbb18d39c9c,
        ]),
        v: Fq([
            0x7588ee6d6dd40deb,
            0x9d6d7a23ebdb7c4c,
            0x46462e26d4edb8c7,
            0x10b4c1517ca82e9b,
        ]),
    }).mul_by_cofactor();
    assert_eq!(p * c, (p * a) * b);
}

#[test]
fn test_serialization_consistency() {
    let gen = FULL_GENERATOR.mul_by_cofactor();
    let mut p = gen;

    let v = vec![
        [
            203, 85, 12, 213, 56, 234, 12, 193, 19, 132, 128, 64, 142, 110, 170, 185, 179, 108, 97,
            63, 13, 211, 247, 120, 79, 219, 110, 234, 131, 123, 19, 215,
        ],
        [
            113, 154, 240, 230, 224, 198, 208, 170, 104, 15, 59, 126, 151, 222, 233, 195, 203, 195,
            167, 129, 89, 121, 240, 142, 51, 166, 64, 250, 184, 202, 154, 177,
        ],
        [
            197, 41, 93, 209, 203, 55, 164, 174, 88, 0, 90, 199, 1, 156, 149, 141, 240, 29, 14, 82,
            86, 225, 126, 129, 186, 157, 148, 162, 219, 51, 156, 199,
        ],
        [
            182, 117, 250, 241, 81, 196, 199, 227, 151, 74, 243, 17, 221, 97, 200, 139, 192, 83,
            231, 35, 214, 14, 95, 69, 130, 201, 4, 116, 177, 19, 179, 0,
        ],
        [
            118, 41, 29, 200, 60, 189, 119, 252, 78, 40, 230, 18, 208, 221, 38, 214, 176, 250, 4,
            10, 77, 101, 26, 216, 193, 198, 226, 84, 25, 177, 230, 185,
        ],
        [
            226, 189, 227, 208, 112, 117, 136, 98, 72, 38, 211, 167, 254, 82, 174, 113, 112, 166,
            138, 171, 166, 113, 52, 251, 129, 197, 138, 45, 195, 7, 61, 140,
        ],
        [
            38, 198, 156, 196, 146, 225, 55, 163, 138, 178, 157, 128, 115, 135, 204, 215, 0, 33,
            171, 20, 60, 32, 142, 209, 33, 233, 125, 146, 207, 12, 16, 24,
        ],
        [
            17, 187, 231, 83, 165, 36, 232, 184, 140, 205, 195, 252, 166, 85, 59, 86, 3, 226, 211,
            67, 179, 29, 238, 181, 102, 142, 58, 63, 57, 89, 174, 138,
        ],
        [
            210, 159, 80, 16, 181, 39, 221, 204, 224, 144, 145, 79, 54, 231, 8, 140, 142, 216, 93,
            190, 183, 116, 174, 63, 33, 242, 177, 118, 148, 40, 241, 203,
        ],
        [
            0, 143, 107, 102, 149, 187, 27, 124, 18, 10, 98, 28, 113, 123, 121, 185, 29, 152, 14,
            130, 149, 28, 87, 35, 135, 135, 153, 54, 112, 53, 54, 68,
        ],
        [
            178, 131, 85, 160, 214, 51, 208, 157, 196, 152, 247, 93, 202, 56, 81, 239, 155, 122,
            59, 188, 237, 253, 11, 169, 208, 236, 12, 4, 163, 211, 88, 97,
        ],
        [
            246, 194, 231, 195, 159, 101, 180, 133, 80, 21, 185, 220, 195, 115, 144, 12, 90, 150,
            44, 117, 8, 156, 168, 248, 206, 41, 60, 82, 67, 75, 57, 67,
        ],
        [
            212, 205, 171, 153, 113, 16, 194, 241, 224, 43, 177, 110, 190, 248, 22, 201, 208, 166,
            2, 83, 134, 130, 85, 129, 166, 136, 185, 191, 163, 38, 54, 10,
        ],
        [
            8, 60, 190, 39, 153, 222, 119, 23, 142, 237, 12, 110, 146, 9, 19, 219, 143, 64, 161,
            99, 199, 77, 39, 148, 70, 213, 246, 227, 150, 178, 237, 178,
        ],
        [
            11, 114, 217, 160, 101, 37, 100, 220, 56, 114, 42, 31, 138, 33, 84, 157, 214, 167, 73,
            233, 115, 81, 124, 134, 15, 31, 181, 60, 184, 130, 175, 159,
        ],
        [
            141, 238, 235, 202, 241, 32, 210, 10, 127, 230, 54, 31, 146, 80, 247, 9, 107, 124, 0,
            26, 203, 16, 237, 34, 214, 147, 133, 15, 29, 236, 37, 88,
        ],
    ];

    for expected_serialized in v {
        assert!(p.is_on_curve_vartime());
        let affine = AffinePoint::from(p);
        let serialized = affine.into_bytes();
        let deserialized = AffinePoint::from_bytes(serialized).unwrap();
        assert_eq!(affine, deserialized);
        assert_eq!(expected_serialized, serialized);
        p = p + &gen;
    }
}
