
This is a pure Rust implementation of the Doppio elliptic curve group and its associated fields.

* **This implementation has not been reviewed or audited. Use at your own risk.**
* This implementation targets Rust `1.33` or later.
* All operations are constant time unless explicitly noted.

## Features

* `std` (on by default): Enables APIs that leverage the Rust standard library.
* `nightly`: Enables `subtle/nightly` which prevents compiler optimizations that could jeopardize constant time operations.

## Curve Description

Jubjub is the [twisted Edwards curve](https://en.wikipedia.org/wiki/Twisted_Edwards_curve) `-u^2 + v^2 = 1 + d.u^2.v^2` of rational points over `GF(q)` with a subgroup of prime order `r` and cofactor `8`.

```
q = 0x1000000000000000000000000000000014def9dea2f79cd65812631a5cf5d3ed
r = 0x1fffffffffffffffffffffffffffffff49b2bf0e49f58d726a9d3de35b7a1e7
d = -(86649/86650)
```

The choice of `GF(q)` is made to be the scalar field of the Ristretto255 elliptic curve construction.

Jubjub is birationally equivalent to a [Montgomery curve](https://en.wikipedia.org/wiki/Montgomery_curve) `y^2 = x^3 + Ax^2 + x` over the same field with `A = 346598`. This value of `A` is the smallest integer such that `(A - 2) / 4` is a small integer, `A^2 - 4` is nonsquare in `GF(q)`, and the Montgomery curve and its quadratic twist have small cofactors `8` and `4`, respectively. This is identical to the relationship between Curve25519 and ed25519.


## Acknowledgements

This is a fork from Jubjub which was designed by Sean Bowe. Daira Hopwood is responsible for its name and specification.

Please see `Cargo.toml` for a list of primary authors of this codebase.

## License

Licensed under either of

 * Apache License, Version 2.0, ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

### Contribution

Unless you explicitly state otherwise, any contribution intentionally
submitted for inclusion in the work by you, as defined in the Apache-2.0
license, shall be dual licensed as above, without any additional terms or
conditions.
