pub use crate::fq::Fq;

/// Constant representing the modulus
/// q = 0x1000000000000000000000000000000014def9dea2f79cd65812631a5cf5d3ed
pub const MODULUS: Fq = Fq([
    0x5812631a5cf5d3ed,
    0x14def9dea2f79cd6,
    0x0000000000000000,
    0x1000000000000000,
]);

/// INV = -(q^{-1} mod 2^64) mod 2^64
pub const INV: u64 = 0xd2b51da312547e1b;

/// R = 2^256 mod q
pub const R: Fq = Fq([
    0xd6ec31748d98951d,
    0xc6ef5bf4737dcf70,
    0xfffffffffffffffe,
    0x0fffffffffffffff,
]);

/// R^2 = 2^512 mod q
pub const R2: Fq = Fq([
    0xa40611e3449c0f01,
    0xd00e1ba768859347,
    0xceec73d217f5be65,
    0x0399411b7c309a3d,
]);

/// R^3 = 2^768 mod q
pub const R3: Fq = Fq([
    0x2a9e49687b83a2db,
    0x278324e6aef7f3ec,
    0x8065dc6c04ec5b65,
    0x0e530b773599cec7,
]);

pub const S: u32 = 2;

/// GENERATOR^t where t * 2^s + 1 = q
/// with t odd. In other words, this
/// is a 2^s root of unity.
pub const ROOT_OF_UNITY: Fq = Fq([
    0xbe8775dfebbe07d4,
    0x0ef0565342ce83fe,
    0x7d3d6d60abc1c27a,
    0x094a7310e07981e7,
]);
