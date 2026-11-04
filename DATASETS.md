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
| `repulsive-surfaces/scenes/` | the authors example scenes (scene files + OBJs), copied by `scripts/rs-fetch-scenes.sh` | local checkout `../RepulsiveSurfaces/repulsive-surfaces` (`scenes/`) | Yu, Brakensiek, Schumacher, Crane, "Repulsive Surfaces", ACM TOG 40(6), 2021 | MIT, Copyright (c) 2019 Nick Sharp |
| (not in the repository) `ishii2012.pdf` | Ishii et al. 2012, the table of knotted genus-2 surfaces used for Fig. 16 | supplied by the user, kept at `/var/data/discrete-geometry/ishii2012.pdf`, outside the repository; no arXiv version exists | Ishii et al. 2012 | not redistributed |
| `flipout/mpz/` | MPZ input models | https://vcg.isti.cnr.it/publication/2014/MPZ14/ (`inputmodels.zip`) | Myles, Pietroni, Zorin, "Robust Field-aligned Global Parametrization", ACM TOG 33(4), 2014 | no licence stated; research use with credit |
| `flipout/psb/` | Princeton Shape Benchmark v1 | https://shape.cs.princeton.edu/benchmark/ (`psb_v1.zip`) | Shilane, Min, Kazhdan, Funkhouser, "The Princeton Shape Benchmark", Shape Modeling International 2004 | academic use only, acknowledge the source; Copyright (c) 2003 Princeton University |
| `thingi10k/` | Thingi10K, 10 000 Thingiverse models | https://huggingface.co/datasets/Thingi10K/Thingi10K, revision `v1.5.0` | Zhou, Jacobson, "Thingi10K: A Dataset of 10,000 3D-Printing Models", arXiv:1605.04797, 2016 | per model, chosen by its author (mostly Creative Commons); every file's author and licence in `data/thingi10k/ATTRIBUTION.tsv`; models named in this repository below |
| `points-as-tori/` and the authors' release (`PointsAsTori/`, outside the repository) | the Points-as-Tori release: code and weights (converted into `crates/ddg-points-as-tori/weights/`), `PAT_training_data` (ABC and blob training clouds `nice/`, validation `ABC17_validation_3.npz`; used by the retraining, the evaluation's ABC and blob sets and the golden files); `data/points-as-tori/analytic.json` is ours | https://github.com/nzfeng/points-as-tori at `253bf9f` | Feng, Gkioulekas, Crane, "Points as Tori: Fast Pointwise Signed Distance for Point Clouds", ACM TOG 45(4), 2026; ABC models: Koch et al., "ABC: A Big CAD Model Dataset for Geometric Deep Learning", CVPR 2019 | release MIT, Copyright (c) 2026 Nicole Feng (`crates/ddg-points-as-tori/LICENSES/points-as-tori-MIT.txt`); the ABC models keep their own terms (research use, credit) |
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
| 38741 | [thing:12034](https://www.thingiverse.com/thing:12034) | Snowflake | JonMonaghan | CC BY | `crates/ddg-points-as-tori/EVALUATION.md` |
| 40993 | [thing:12652](https://www.thingiverse.com/thing:12652) | Creatures of Thingiverse | hlbrosnihan | CC BY-SA | `crates/ddg-points-as-tori/EVALUATION.md` |
| 41984 | [thing:11836](https://www.thingiverse.com/thing:11836) | Functional Differential Gear System | Thing-O-Fun | CC BY-SA | `crates/ddg-points-as-tori/EVALUATION.md` |
| 44395 | [thing:13391](https://www.thingiverse.com/thing:13391) | Miniature Queen Anne Chair | PrettySmallThings | CC BY-NC-SA | `crates/ddg-integer-coords/README.md`, `crates/ddg-integer-coords/benches/representation.rs`, `scripts/ic-timings.sh` |
| 49167 | [thing:14599](https://www.thingiverse.com/thing:14599) | Tri-Track dual servo robot tank chassis | Ktronics | CC BY-SA | `crates/ddg-points-as-tori/EVALUATION.md` |
| 49544 | [thing:15007](https://www.thingiverse.com/thing:15007) | Kindle Fire Handle | dougc314 | CC BY-SA | `crates/ddg-points-as-tori/EVALUATION.md` |
| 54621 | [thing:16636](https://www.thingiverse.com/thing:16636) | Akashi Head | bondi | CC BY-SA | `crates/ddg-integer-coords/README.md`, `crates/ddg-integer-coords/SPEC.md` |
| 54674 | [thing:16650](https://www.thingiverse.com/thing:16650) | James Tiberius Kirk #1 | MicahGanske | CC BY-NC-ND | `crates/ddg-integer-coords/README.md`, `crates/ddg-integer-coords/benches/representation.rs`, `scripts/ic-timings.sh` |
| 56527 | [thing:17309](https://www.thingiverse.com/thing:17309) | Pet Monster Valentine - Heartless | cyclone | CC BY-NC | `crates/ddg-integer-coords/README.md` |
| 56954 | [thing:17455](https://www.thingiverse.com/thing:17455) | Kazoo | wes | CC BY-SA | `crates/ddg-points-as-tori/EVALUATION.md` |
| 58938 | [thing:18132](https://www.thingiverse.com/thing:18132) | MakerBot Wizard | wes | CC BY | `crates/ddg-points-as-tori/EVALUATION.md` |
| 59340 | [thing:18284](https://www.thingiverse.com/thing:18284) | GOYLE | fong182 | CC BY-SA | `crates/ddg-points-as-tori/EVALUATION.md` |
| 59936 | [thing:18392](https://www.thingiverse.com/thing:18392) | Mini Rover | Kepler | CC BY | `crates/ddg-points-as-tori/EVALUATION.md` |
| 61039 | [thing:18928](https://www.thingiverse.com/thing:18928) | Actin Protein Structure | destroyer2012 | CC BY-SA | `crates/ddg-points-as-tori/EVALUATION.md` |
| 62287 | [thing:19518](https://www.thingiverse.com/thing:19518) | 60s Sofa | PrettySmallThings | CC BY | `crates/ddg-points-as-tori/EVALUATION.md` |
| 63244 | [thing:19938](https://www.thingiverse.com/thing:19938) | GarabatoBOT (DoodleBOT) | MAF | CC BY | `crates/ddg-points-as-tori/EVALUATION.md` |
| 65399 | [thing:20896](https://www.thingiverse.com/thing:20896) | Bowstring Equatorial Sundial | alany | CC BY | `crates/ddg-integer-coords/README.md` |
| 72613 | [thing:22858](https://www.thingiverse.com/thing:22858) | Lightcycle Model Kit | cbabbage | CC BY | `crates/ddg-points-as-tori/EVALUATION.md` |
| 73178 | [thing:24155](https://www.thingiverse.com/thing:24155) | The Brute | megamanex | CC BY-NC | `crates/ddg-integer-coords/README.md`, `crates/ddg-integer-coords/SPEC.md` |
| 73736 | [thing:24308](https://www.thingiverse.com/thing:24308) | A Dodecahedron Speaker for Desktop Printers | seanmichaelragan | CC BY | `crates/ddg-points-as-tori/EVALUATION.md` |
| 74463 | [thing:24670](https://www.thingiverse.com/thing:24670) | GOLF-IN-MINIATURE : The Desktop 18 Hole Miniature Golf Course | beekeeper | CC BY-SA | `crates/ddg-integer-coords/README.md` |
| 75098 | [thing:24908](https://www.thingiverse.com/thing:24908) | Robots | Oblomobka | CC BY-SA | `crates/ddg-points-as-tori/EVALUATION.md` |
| 75282 | [thing:24999](https://www.thingiverse.com/thing:24999) | the Addalock | Mr_MegaTronic | CC BY | `crates/ddg-points-as-tori/EVALUATION.md` |
| 76291 | [thing:25344](https://www.thingiverse.com/thing:25344) | Magic Block Piercing Trick | dnewman | CC BY-SA | `crates/ddg-points-as-tori/EVALUATION.md` |
| 77014 | [thing:25680](https://www.thingiverse.com/thing:25680) | Lynx Building Toy | OSTycoon | CC BY | `crates/ddg-points-as-tori/EVALUATION.md` |
| 78735 | [thing:26334](https://www.thingiverse.com/thing:26334) | Almost Impossible Heart | GregFrost | CC BY | `crates/ddg-integer-coords/README.md` |
| 79193 | [thing:26519](https://www.thingiverse.com/thing:26519) | Automaton Hand | MakerBot | CC BY | `crates/ddg-points-as-tori/EVALUATION.md` |
| 80354 | [thing:26920](https://www.thingiverse.com/thing:26920) | PIXAR's Luxo Jr. | Mattraptions | CC BY-SA | `crates/ddg-points-as-tori/EVALUATION.md` |
| 81313 | [thing:27317](https://www.thingiverse.com/thing:27317) | MakerBot Mixtape | MakerBot | CC BY-SA | `crates/ddg-integer-coords/README.md` |
| 82541 | [thing:27520](https://www.thingiverse.com/thing:27520) | Sleep No More Mask | gianteye | CC BY-NC | `crates/ddg-integer-coords/README.md` |
| 87395 | [thing:29745](https://www.thingiverse.com/thing:29745) | Sugru Fixbot | mikevanis | CC BY-SA | `crates/ddg-points-as-tori/EVALUATION.md` |
| 91611 | [thing:31497](https://www.thingiverse.com/thing:31497) | Simple Adjacent Line Tool | Halodine | CC BY-SA | `crates/ddg-points-as-tori/EVALUATION.md` |
| 94153 | [thing:32510](https://www.thingiverse.com/thing:32510) | Shack Bot | MakerBot | CC BY | `crates/ddg-points-as-tori/EVALUATION.md` |
| 97660 | [thing:33913](https://www.thingiverse.com/thing:33913) | Big Ben | MakeALot | CC BY-SA | `crates/ddg-integer-coords/README.md`, `crates/ddg-integer-coords/SPEC.md` |
| 97729 | [thing:33940](https://www.thingiverse.com/thing:33940) | Universal bottle cap tripod | Camos | CC BY-SA | `crates/ddg-points-as-tori/EVALUATION.md` |
| 103538 | [thing:36200](https://www.thingiverse.com/thing:36200) | Gunny Sacks | BenitoSanduchi | CC BY-SA | `crates/ddg-points-as-tori/EVALUATION.md` |
| 106830 | [thing:37525](https://www.thingiverse.com/thing:37525) | Blizzard of Unique Snowflakes | laird | CC BY-NC | `crates/ddg-integer-coords/README.md` |
| 106831 | [thing:37525](https://www.thingiverse.com/thing:37525) | Blizzard of Unique Snowflakes | laird | CC BY-NC | `crates/ddg-integer-coords/README.md` |
| 106839 | [thing:37525](https://www.thingiverse.com/thing:37525) | Blizzard of Unique Snowflakes | laird | CC BY-NC | `crates/ddg-integer-coords/README.md` |
| 106840 | [thing:37525](https://www.thingiverse.com/thing:37525) | Blizzard of Unique Snowflakes | laird | CC BY-NC | `crates/ddg-integer-coords/README.md` |
| 107202 | [thing:37525](https://www.thingiverse.com/thing:37525) | Blizzard of Unique Snowflakes | laird | CC BY-NC | `crates/ddg-integer-coords/README.md` |
| 107587 | [thing:36554](https://www.thingiverse.com/thing:36554) | Advent Calendar | pleppik | CC BY | `crates/ddg-points-as-tori/EVALUATION.md` |
| 109772 | [thing:36554](https://www.thingiverse.com/thing:36554) | Advent Calendar | pleppik | CC BY | `crates/ddg-points-as-tori/EVALUATION.md` |
| 117647 | [thing:42094](https://www.thingiverse.com/thing:42094) | 3d Print - Avocado Pit Germination Holder - Cup Edition | sa7ina | CC0 | `crates/ddg-points-as-tori/EVALUATION.md` |
| 119918 | [thing:42681](https://www.thingiverse.com/thing:42681) | Truck V1 | reza | CC BY-SA | `crates/ddg-points-as-tori/EVALUATION.md` |
| 128482 | [thing:47654](https://www.thingiverse.com/thing:47654) | Spin Top | mariothemagician | CC BY-SA | `crates/ddg-points-as-tori/EVALUATION.md` |
| 153951 | [thing:61532](https://www.thingiverse.com/thing:61532) | RasPi Case | MakerBot | CC BY | `crates/ddg-points-as-tori/EVALUATION.md` |
| 162575 | [thing:20748](https://www.thingiverse.com/thing:20748) | Fennec Fox | Physics_Dude | CC BY | `crates/ddg-points-as-tori/EVALUATION.md` |
| 176193 | [thing:74322](https://www.thingiverse.com/thing:74322) | Customizable Lithopane | MakerBot | CC BY | `crates/ddg-integer-coords/README.md`, `crates/ddg-integer-coords/SPEC.md`, `crates/ddg-integer-coords/tests/thingi_regressions.rs`, `crates/ddg-intrinsic/src/refine.rs` |
| 178197 | [thing:76369](https://www.thingiverse.com/thing:76369) | Turbine Rotary Tool 60000 rpm | RichMac | CC BY | `crates/ddg-points-as-tori/EVALUATION.md` |
| 178201 | [thing:76369](https://www.thingiverse.com/thing:76369) | Turbine Rotary Tool 60000 rpm | RichMac | CC BY | `crates/ddg-points-as-tori/EVALUATION.md` |
| 190253 | [thing:84142](https://www.thingiverse.com/thing:84142) | MaxGlide - Glidecam simple; efficient and cheap! | MaximSachs | CC BY | `crates/ddg-integer-coords/README.md`, `crates/ddg-integer-coords/SPEC.md`, `crates/ddg-integer-coords/tests/thingi_regressions.rs`, `crates/ddg-intrinsic/src/refine.rs` |
| 238413 | [thing:114349](https://www.thingiverse.com/thing:114349) | Traxbot - an EZ-robot build | infinityplusplus | CC BY-SA | `crates/ddg-points-as-tori/EVALUATION.md` |
| 250042 | [thing:121456](https://www.thingiverse.com/thing:121456) | Wheeled Inchworm | someandy | CC BY-SA | `crates/ddg-points-as-tori/EVALUATION.md` |
| 264708 | [thing:130890](https://www.thingiverse.com/thing:130890) | PLA Gear Kit | whpthomas | CC BY-SA | `crates/ddg-integer-coords/README.md` |
| 285440 | [thing:146293](https://www.thingiverse.com/thing:146293) | Verlan Dress from New Skins with Francis Bitonti Studio | MakerBot | CC BY | `crates/ddg-integer-coords/README.md` |
| 286163 | [thing:146400](https://www.thingiverse.com/thing:146400) | Spoonriders | JuliaTruchsess | CC BY-SA | `crates/ddg-integer-coords/README.md` |
| 298321 | [thing:131054](https://www.thingiverse.com/thing:131054) | Bunny Peel with meshmixer Support | meshmixer | CC BY | `crates/ddg-integer-coords/README.md` |
| 298322 | [thing:131054](https://www.thingiverse.com/thing:131054) | Bunny Peel with meshmixer Support | meshmixer | CC BY | `crates/ddg-integer-coords/README.md` |
| 298323 | [thing:131054](https://www.thingiverse.com/thing:131054) | Bunny Peel with meshmixer Support | meshmixer | CC BY | `crates/ddg-integer-coords/README.md` |
| 376253 | [thing:196038](https://www.thingiverse.com/thing:196038) | Winged Victory of Samothrace | CosmoWenman | CC BY | `crates/ddg-integer-coords/README.md` |
| 462523 | [thing:252162](https://www.thingiverse.com/thing:252162) | 41 Piece Spirograph | Hercemer | CC BY-SA | `crates/ddg-points-as-tori/EVALUATION.md` |
| 472053 | [thing:259005](https://www.thingiverse.com/thing:259005) | 4-8-8-4 Big Boy Locomotive | MakerBot | CC BY | `crates/ddg-points-as-tori/EVALUATION.md` |
| 472093 | [thing:259005](https://www.thingiverse.com/thing:259005) | 4-8-8-4 Big Boy Locomotive | MakerBot | CC BY | `crates/ddg-points-as-tori/EVALUATION.md` |
| 472101 | [thing:259005](https://www.thingiverse.com/thing:259005) | 4-8-8-4 Big Boy Locomotive | MakerBot | CC BY | `crates/ddg-points-as-tori/EVALUATION.md` |
| 472174 | [thing:259005](https://www.thingiverse.com/thing:259005) | 4-8-8-4 Big Boy Locomotive | MakerBot | CC BY | `crates/ddg-points-as-tori/EVALUATION.md` |
| 509310 | [thing:282868](https://www.thingiverse.com/thing:282868) | Catalan Wireframe Polyhedra | mathgrrl | CC BY-SA | `crates/ddg-points-as-tori/EVALUATION.md` |
| 521601 | [thing:279098](https://www.thingiverse.com/thing:279098) | Calendar Blocks | Kraaij | CC BY-SA | `crates/ddg-integer-coords/README.md` |
| 521603 | [thing:279098](https://www.thingiverse.com/thing:279098) | Calendar Blocks | Kraaij | CC BY-SA | `crates/ddg-integer-coords/README.md` |
| 527631 | [thing:295142](https://www.thingiverse.com/thing:295142) | Holden Observatory -- 3D Scanning a Full Size Building! | ibudmen | CC BY-ND | `crates/ddg-integer-coords/README.md` |
| 550964 | [thing:311002](https://www.thingiverse.com/thing:311002) | Famous Paris buildings | LeFabShop | CC BY-SA | `crates/ddg-points-as-tori/EVALUATION.md` |
| 551019 | [thing:311002](https://www.thingiverse.com/thing:311002) | Famous Paris buildings | LeFabShop | CC BY-SA | `crates/ddg-integer-coords/README.md`, `crates/ddg-integer-coords/SPEC.md` |
| 700904 | [thing:33357](https://www.thingiverse.com/thing:33357) | A useful 3D-printed tripod (1/4-20 screw reqd) - UPDATED (July 31; 2014) | cmceezee | CC BY-SA | `crates/ddg-points-as-tori/EVALUATION.md` |
| 719791 | [thing:427024](https://www.thingiverse.com/thing:427024) | Cable Spool | cheewee2000 | CC BY | `crates/ddg-integer-coords/README.md`, `crates/ddg-integer-coords/SPEC.md`, `scripts/ic-timings.sh` |
| 1436278 | [thing:909118](https://www.thingiverse.com/thing:909118) | 3D-printable sand play set | CreativeTools | CC BY-SA | `crates/ddg-points-as-tori/EVALUATION.md` |
| 1514498 | [thing:957943](https://www.thingiverse.com/thing:957943) | Vienna Giant Wheel (Riesenrad) | stylesuxx | CC BY-SA | `crates/ddg-points-as-tori/EVALUATION.md` |
