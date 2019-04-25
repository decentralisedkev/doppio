
pub use crate::fq::Fq;

/// Constant representing the modulus
/// q = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
pub const MODULUS: Fq = Fq([
    0xffffffff00000001,
    0x53bda402fffe5bfe,
    0x3339d80809a1d805,
    0x73eda753299d7d48,
]);

/// INV = -(q^{-1} mod 2^64) mod 2^64
pub const INV: u64 = 0xfffffffeffffffff;

/// R = 2^256 mod q
pub const R: Fq = Fq([
    0x00000001fffffffe,
    0x5884b7fa00034802,
    0x998c4fefecbc4ff5,
    0x1824b159acc5056f,
]);

/// R^2 = 2^512 mod q
pub const R2: Fq = Fq([
    0xc999e990f3f29c6d,
    0x2b6cedcb87925c23,
    0x05d314967254398f,
    0x0748d9d99f59ff11,
]);

/// R^3 = 2^768 mod q
pub const R3: Fq = Fq([
    0xc62c1807439b73af,
    0x1b3e0d188cf06990,
    0x73d13c71c7b5f418,
    0x6e2a5bb9c8db33e9,
]);

// /// 7*R mod q
// pub const GENERATOR: Fq = Fq([
//     0x0000000efffffff1,
//     0x17e363d300189c0f,
//     0xff9c57876f8457b0,
//     0x351332208fc5a8c4,
// ]);

pub const S: u32 = 32;

/// GENERATOR^t where t * 2^s + 1 = q
/// with t odd. In other words, this
/// is a 2^s root of unity.
pub const ROOT_OF_UNITY: Fq = Fq([
    0xb9b58d8c5f0e466a,
    0x5b1b4c801819d7ec,
    0x0af53ae352a31e64,
    0x5bf3adda19e9b27b,
]);