use crate::fq::Fq;

/// `d = -(10240/10241)`
pub const EDWARDS_D: Fq = Fq::from_raw([
    0x01065fd6d6343eb1,
    0x292d7f6d37579d26,
    0xf5fd9207e6bd7fd4,
    0x2a9318e74bfa2b48,
]);

/// `2*d`
pub const EDWARDS_D2: Fq = Fq::from_raw([
    0x020cbfadac687d62,
    0x525afeda6eaf3a4c,
    0xebfb240fcd7affa8,
    0x552631ce97f45691,
]);

/// byte representation of the scalar modulus
pub const FR_MODULUS_BYTES: [u8; 32] = [
    183, 44, 247, 214, 94, 14, 151, 208, 130, 16, 200, 204, 147, 32, 104, 166, 0, 59, 52, 1, 1, 59,
    103, 6, 169, 175, 51, 101, 234, 180, 125, 14,
];