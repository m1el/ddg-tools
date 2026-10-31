# Datasets

Every dataset this workspace uses, with its source, citation and terms. The data itself lives in
`data/` (gitignored); each folder there has its own README with checksums or details. Only the meshes
under "Committed files" are part of the repository.

## In `data/` (not committed)

| folder | dataset | source | citation | terms |
|---|---|---|---|---|
| `models/` | Keenan Crane's 3D model repository: Spot, Bob, Blub, Crumpled Developable | https://www.cs.cmu.edu/~kmcrane/Projects/ModelRepository/ | Keenan Crane | CC0 1.0 (public domain) |
| `optimal-cones/` | optimal cone singularity examples (meshes, cones, UVs) | `optimal-cones-data.zip`, http://www.cs.cmu.edu/~kmcrane/Projects/OptimalCones/ | Soliman, Slepčev, Crane, "Optimal Cone Singularities for Conformal Flattening", ACM TOG 37(4), 2018 | no licence stated; research use with credit; `bunny` is the Stanford Bunny (Stanford Computer Graphics Laboratory) |
| `intrinsic-reference/` | outputs of the authors' demo, generated here | https://github.com/nmwsharp/navigating-intrinsic-triangulations-demo, run by `scripts/intrinsic-reference.sh` | Sharp, Soliman, Crane, "Navigating Intrinsic Triangulations", ACM TOG 38(4), 2019 | demo MIT; inputs from `models/`, `optimal-cones/` and our generator |
| `repulsive-curves/KnotDatasets/` | Knot128, Trefoil100 | `KnotDatasets.zip`, https://www.cs.cmu.edu/~kmcrane/Projects/RepulsiveCurves/ | Yu, Schumacher, Crane, "Repulsive Curves", ACM TOG 40(2), 2021; prototypes from KnotPlot (Robert Scharein) | no licence stated |
| `repulsive-curves/FreedmanHeWang/` | the Freedman–He–Wang unknot | `FreedmanHeWang.zip`, same page | Freedman, He, Wang, "Möbius Energy of Knots and Unknots", Annals of Mathematics 139, 1994 | no licence stated |
| `repulsive-curves/scenes/` | the authors' example scenes | https://github.com/icethrush/repulsive-curves (`scenes/`) | Yu, Schumacher, Crane 2021 | MIT, Copyright (c) 2019 Christopher Yu |
| `flipout/mpz/` | MPZ input models | https://vcg.isti.cnr.it/publication/2014/MPZ14/ (`inputmodels.zip`) | Myles, Pietroni, Zorin, "Robust Field-aligned Global Parametrization", ACM TOG 33(4), 2014 | no licence stated; research use with credit |
| `flipout/psb/` | Princeton Shape Benchmark v1 | https://shape.cs.princeton.edu/benchmark/ (`psb_v1.zip`) | Shilane, Min, Kazhdan, Funkhouser, "The Princeton Shape Benchmark", Shape Modeling International 2004 | academic use only, acknowledge the source; Copyright (c) 2003 Princeton University |
| `thingi10k/` | Thingi10K, 10 000 Thingiverse models | https://huggingface.co/datasets/Thingi10K/Thingi10K, revision `v1.5.0` | Zhou, Jacobson, "Thingi10K: A Dataset of 10,000 3D-Printing Models", arXiv:1605.04797, 2016 | per model, chosen by its author (mostly Creative Commons); every file's author and licence in `data/thingi10k/ATTRIBUTION.tsv`; models named in this repository below |
| `anatomical-boyer2011/` | teeth, metatarsal and radius surfaces (222 meshes) | arXiv ancillary files, https://arxiv.org/abs/1110.3649v2 | Boyer et al., "Algorithms to automatically quantify the geometric similarity of anatomical surfaces", PNAS 108(45), 2011 | no licence stated |
| `brain-yeo2011/` | Yeo 2011 cortical parcellation, fsaverage meshes, Schaefer 2018 and Yan 2023 parcellations, CBIG templates | https://github.com/ThomasYeoLab/CBIG at commit 634f676 (data only) | Yeo et al., J Neurophysiol 106(3), 2011; Schaefer et al., Cereb Cortex 28(9), 2018; Yan et al., NeuroImage 273, 2023 | CBIG MIT; fsaverage under the FreeSurfer licence; fs_LR/Conte69 from HCP / Van Essen lab; MNI152, Colin27 and CoRR_HNU under their own terms (`data/brain-yeo2011/README.md`) |

## Committed files

| path | content | source | terms |
|---|---|---|---|
| `crates/ddg-integer-coords/tests/assets/fox.ply`, `cat_head.obj` | geometry-central's test meshes | https://github.com/nmwsharp/geometry-central, `test/assets/` at ab116f1 | MIT, Copyright (c) 2017-2019 Nicholas Sharp and contributors (`crates/ddg-integer-coords/LICENSES/geometry-central-MIT.txt`) |
| `crates/ddg-repulsive-curves/scenes/` | five small scenes | https://github.com/icethrush/repulsive-curves, `scenes/` | MIT, Copyright (c) 2019 Christopher Yu (`crates/ddg-repulsive-curves/LICENSES/repulsive-curves-MIT.txt`) |

## Thingi10K models named in this repository

Each is a Thingiverse model by the author shown, under the licence shown (from the dataset's metadata).
The READMEs and contracts refer to them by Thingi10K file id; no mesh is redistributed. Any render or
figure of one must carry this credit and respect its licence (NC: no commercial use; ND: no modified
versions; SA: share alike).

| file id | Thingiverse thing | name | author | licence | named in |
|---|---|---|---|---|---|
| 44395 | [thing:13391](https://www.thingiverse.com/thing:13391) | Miniature Queen Anne Chair | PrettySmallThings | CC BY-NC-SA | `crates/ddg-integer-coords/README.md`, `crates/ddg-integer-coords/benches/representation.rs`, `scripts/ic-timings.sh` |
| 54621 | [thing:16636](https://www.thingiverse.com/thing:16636) | Akashi Head | bondi | CC BY-SA | `crates/ddg-integer-coords/README.md`, `crates/ddg-integer-coords/SPEC.md` |
| 54674 | [thing:16650](https://www.thingiverse.com/thing:16650) | James Tiberius Kirk #1 | MicahGanske | CC BY-NC-ND | `crates/ddg-integer-coords/README.md`, `crates/ddg-integer-coords/benches/representation.rs`, `scripts/ic-timings.sh` |
| 56527 | [thing:17309](https://www.thingiverse.com/thing:17309) | Pet Monster Valentine - Heartless | cyclone | CC BY-NC | `crates/ddg-integer-coords/README.md` |
| 65399 | [thing:20896](https://www.thingiverse.com/thing:20896) | Bowstring Equatorial Sundial | alany | CC BY | `crates/ddg-integer-coords/README.md` |
| 73178 | [thing:24155](https://www.thingiverse.com/thing:24155) | The Brute | megamanex | CC BY-NC | `crates/ddg-integer-coords/README.md`, `crates/ddg-integer-coords/SPEC.md` |
| 74463 | [thing:24670](https://www.thingiverse.com/thing:24670) | GOLF-IN-MINIATURE : The Desktop 18 Hole Miniature Golf Course | beekeeper | CC BY-SA | `crates/ddg-integer-coords/README.md` |
| 78735 | [thing:26334](https://www.thingiverse.com/thing:26334) | Almost Impossible Heart | GregFrost | CC BY | `crates/ddg-integer-coords/README.md` |
| 81313 | [thing:27317](https://www.thingiverse.com/thing:27317) | MakerBot Mixtape | MakerBot | CC BY-SA | `crates/ddg-integer-coords/README.md` |
| 82541 | [thing:27520](https://www.thingiverse.com/thing:27520) | Sleep No More Mask | gianteye | CC BY-NC | `crates/ddg-integer-coords/README.md` |
| 97660 | [thing:33913](https://www.thingiverse.com/thing:33913) | Big Ben | MakeALot | CC BY-SA | `crates/ddg-integer-coords/README.md`, `crates/ddg-integer-coords/SPEC.md` |
| 106830 | [thing:37525](https://www.thingiverse.com/thing:37525) | Blizzard of Unique Snowflakes | laird | CC BY-NC | `crates/ddg-integer-coords/README.md` |
| 106831 | [thing:37525](https://www.thingiverse.com/thing:37525) | Blizzard of Unique Snowflakes | laird | CC BY-NC | `crates/ddg-integer-coords/README.md` |
| 106839 | [thing:37525](https://www.thingiverse.com/thing:37525) | Blizzard of Unique Snowflakes | laird | CC BY-NC | `crates/ddg-integer-coords/README.md` |
| 106840 | [thing:37525](https://www.thingiverse.com/thing:37525) | Blizzard of Unique Snowflakes | laird | CC BY-NC | `crates/ddg-integer-coords/README.md` |
| 107202 | [thing:37525](https://www.thingiverse.com/thing:37525) | Blizzard of Unique Snowflakes | laird | CC BY-NC | `crates/ddg-integer-coords/README.md` |
| 176193 | [thing:74322](https://www.thingiverse.com/thing:74322) | Customizable Lithopane | MakerBot | CC BY | `crates/ddg-integer-coords/README.md`, `crates/ddg-integer-coords/SPEC.md`, `crates/ddg-integer-coords/tests/thingi_regressions.rs`, `crates/ddg-intrinsic/src/refine.rs` |
| 190253 | [thing:84142](https://www.thingiverse.com/thing:84142) | MaxGlide - Glidecam simple; efficient and cheap! | MaximSachs | CC BY | `crates/ddg-integer-coords/README.md`, `crates/ddg-integer-coords/SPEC.md`, `crates/ddg-integer-coords/tests/thingi_regressions.rs`, `crates/ddg-intrinsic/src/refine.rs` |
| 264708 | [thing:130890](https://www.thingiverse.com/thing:130890) | PLA Gear Kit | whpthomas | CC BY-SA | `crates/ddg-integer-coords/README.md` |
| 285440 | [thing:146293](https://www.thingiverse.com/thing:146293) | Verlan Dress from New Skins with Francis Bitonti Studio | MakerBot | CC BY | `crates/ddg-integer-coords/README.md` |
| 286163 | [thing:146400](https://www.thingiverse.com/thing:146400) | Spoonriders | JuliaTruchsess | CC BY-SA | `crates/ddg-integer-coords/README.md` |
| 298321 | [thing:131054](https://www.thingiverse.com/thing:131054) | Bunny Peel with meshmixer Support | meshmixer | CC BY | `crates/ddg-integer-coords/README.md` |
| 298322 | [thing:131054](https://www.thingiverse.com/thing:131054) | Bunny Peel with meshmixer Support | meshmixer | CC BY | `crates/ddg-integer-coords/README.md` |
| 298323 | [thing:131054](https://www.thingiverse.com/thing:131054) | Bunny Peel with meshmixer Support | meshmixer | CC BY | `crates/ddg-integer-coords/README.md` |
| 376253 | [thing:196038](https://www.thingiverse.com/thing:196038) | Winged Victory of Samothrace | CosmoWenman | CC BY | `crates/ddg-integer-coords/README.md` |
| 521601 | [thing:279098](https://www.thingiverse.com/thing:279098) | Calendar Blocks | Kraaij | CC BY-SA | `crates/ddg-integer-coords/README.md` |
| 521603 | [thing:279098](https://www.thingiverse.com/thing:279098) | Calendar Blocks | Kraaij | CC BY-SA | `crates/ddg-integer-coords/README.md` |
| 527631 | [thing:295142](https://www.thingiverse.com/thing:295142) | Holden Observatory -- 3D Scanning a Full Size Building! | ibudmen | CC BY-ND | `crates/ddg-integer-coords/README.md` |
| 551019 | [thing:311002](https://www.thingiverse.com/thing:311002) | Famous Paris buildings | LeFabShop | CC BY-SA | `crates/ddg-integer-coords/README.md`, `crates/ddg-integer-coords/SPEC.md` |
| 719791 | [thing:427024](https://www.thingiverse.com/thing:427024) | Cable Spool | cheewee2000 | CC BY | `crates/ddg-integer-coords/README.md`, `crates/ddg-integer-coords/SPEC.md`, `scripts/ic-timings.sh` |
