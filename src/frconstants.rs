pub use crate::fr::Fr;

/// Constant representing the modulus
/// r = 0x1fffffffffffffffffffffffffffffff49b2bf0e49f58d726a9d3de35b7a1e7
pub const MODULUS: Fr = Fr([
    0x26a9d3de35b7a1e7,
    0xf49b2bf0e49f58d7,
    0xffffffffffffffff,
    0x01ffffffffffffff,
]);

/// INV = -(r^{-1} mod 2^64) mod 2^64
pub const INV: u64 = 0x2a81f20882b21e29;

/// R = 2^256 mod r
pub const R: Fr = Fr([
    0xab1610e5242f0c80,
    0xb26a078db053946c,
    0x0000000000000005,
    0x0000000000000000,
]);

/// R^2 = 2^512 mod r
pub const R2: Fr = Fr([
    0x6921bd75f1e321aa,
    0x016f997a4e557d3f,
    0xfe677f26b8e821f2,
    0x007be9f42e0719ec,
]);

/// R^2 = 2^768 mod r
pub const R3: Fr = Fr([
    0xd26d015bfa48a645,
    0x6d5e224bce0e0539,
    0x9525cfe838da6011,
    0x003991e9df0f24dc,
]);
