# ddg-tools

Rust implementations of discrete differential geometry papers, with Blender add-ons built on them.

## Papers implemented

- Yousuf Soliman, Dejan Slepčev, Keenan Crane. *Optimal Cone Singularities for Conformal Flattening.*
  ACM Transactions on Graphics 37(4), 2018. https://doi.org/10.1145/3197517.3201367
- Rohan Sawhney, Keenan Crane. *Boundary First Flattening.*
  ACM Transactions on Graphics 37(1), 2017. https://doi.org/10.1145/3132705
- Nicholas Sharp, Keenan Crane. *Variational Surface Cutting.*
  ACM Transactions on Graphics 37(4), 2018. https://doi.org/10.1145/3197517.3201356
  A clean-room implementation from the text of the paper (`crates/ddg-variational-cuts`; its
  `SPEC.md` is the discrete method and lists every difference from the paper as printed).

See `CLAUDE.md` for the repository layout and the mathematical conventions shared by all crates.
`ROADMAP.md` lists what is scheduled next and the known gaps.

## License

MIT. Ported code keeps its original license notice in the crate that contains it.
