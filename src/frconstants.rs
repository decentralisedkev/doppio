
pub use crate::fr::Fr;

/// Constant representing the modulus
/// r = 0x0e7db4ea6533afa906673b0101343b00a6682093ccc81082d0970e5ed6f72cb7
pub const MODULUS: Fr = Fr([
    0xd0970e5ed6f72cb7,
    0xa6682093ccc81082,
    0x06673b0101343b00,
    0x0e7db4ea6533afa9,
]);

/// INV = -(r^{-1} mod 2^64) mod 2^64
pub const INV: u64 = 0x1ba3a358ef788ef9;

/// R = 2^256 mod r
pub const R: Fr = Fr([
    0x25f80bb3b99607d9,
    0xf315d62f66b6e750,
    0x932514eeeb8814f4,
    0x09a6fc6f479155c6,
]);

/// R^2 = 2^512 mod r
pub const R2: Fr = Fr([
    0x67719aa495e57731,
    0x51b0cef09ce3fc26,
    0x69dab7fac026e9a5,
    0x04f6547b8d127688,
]);

/// R^2 = 2^768 mod r
pub const R3: Fr = Fr([
    0xe0d6c6563d830544,
    0x323e3883598d0f85,
    0xf0fea3004c2e2ba8,
    0x05874f84946737ec,
]);