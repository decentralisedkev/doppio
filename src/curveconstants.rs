use crate::fq::Fq;

/// `d = -(86649/86650)`
pub const EDWARDS_D: Fq = Fq::from_raw([
    0x9403565f5a8d532d,
    0xf2f07621646802fb,
    0x1babc9915ffa2370,
    0x02a84588fbffedd9,
]);

/// `2*d`
pub const EDWARDS_D2: Fq = Fq::from_raw([
    0x16b34ccf69bf0be0,
    0x37e44f22cc2babfe,
    0x17529def406b3085,
    0x06263f0b3951c45b,
]);

/// byte representation of the scalar modulus
pub const FR_MODULUS_BYTES: [u8; 32] = [
    1, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 244, 155, 43,
    240, 228, 159, 88, 215, 38, 169, 211, 222, 53, 183, 161, 231,
];
