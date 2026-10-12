import DdgProofs

/-!
# Axiom audit

`lake env lean DdgProofs/Audit.lean` prints, for each main theorem, the axioms its proof depends
on. The expected answer everywhere is the three standard axioms of Lean/Mathlib:
`propext`, `Classical.choice`, `Quot.sound`. (`scripts/check-formulas.sh --build` checks this.)

This file is not imported by the library.
-/

open DdgProofs

-- A.1 scalar convex facts
#print axioms isSubgradientWithinAt_abs_iff
#print axioms isSubgradientWithinAt_abs_sub_iff_eq_clampE
#print axioms isSubgradientWithinAt_abs_iff_eq_clampE
#print axioms eq_softThreshold_add_iff_isSubgradient
#print axioms eq_softThreshold_add_iff_cases
#print axioms eq_softThreshold_add_iff
#print axioms eq_softThreshold_add_iff_complementarity
#print axioms eq_softThreshold_div_iff_isSubgradient
#print axioms minimizes_reg_iff_eq_clampE
#print axioms Site.stationary_iff_eq_prox
-- A.2 optimality system
#print axioms ConeProblem.KKT.objective_add_le
#print axioms ConeProblem.KKT.isMinimizer
#print axioms ConeProblem.KKT.minimizes_reduced_of_posDef
#print axioms ConeProblem.KKT.minimizes_closed
#print axioms ConeProblem.IsMinimizer.unique
#print axioms ConeProblem.IsMinimizer.existsUnique_kkt
#print axioms ConeProblem.isMinimizer_iff_exists_kkt_of_posDef
-- A.3 guarantee
#print axioms ConeProblem.IsMinimizer.distortion_le
-- A.4 certificate
#print axioms ConeProblem.gap_identity
#print axioms ConeProblem.weak_duality
#print axioms ConeProblem.duality_gap_certificate
#print axioms ConeProblem.adjoint_gap_certificate
#print axioms ConeProblem.KKT.dualObjective_eq
-- A.5 Newton
#print axioms paper_stageI_stopping_rule
#print axioms paper_stageII_stopping_rule
#print axioms ConeProblem.residual_eq_zero_of_newton_step
#print axioms ConeProblem.isMinimizer_of_newton_step
#print axioms ConeProblem.residualReg_eq_zero_of_newton_step
#print axioms ConeProblem.isRegMinimizer_of_newton_step
#print axioms ConeProblem.IsRegMinimizer.objective_le
-- A.6, A.7 closed meshes, scaling
#print axioms ConeProblem.sum_cone_eq_sum_defect
#print axioms ConeProblem.distortion_add_const_eq
#print axioms ConeProblem.KKT.weighted_mean_eq_zero
#print axioms ConeProblem.IsMinimizer.weighted_mean_eq_zero
#print axioms ConeProblem.isMinimizer_scaleMass_iff
#print axioms ConeProblem.normalizeMass_scaleMass
-- B.8 triangles
#print axioms cos_cornerAngle
#print axioms sin_cornerAngle
#print axioms cornerCot_eq_cot
#print axioms cornerAngle_add_cornerAngle_add_cornerAngle
#print axioms cornerAngle_eq_two_arctan
#print axioms angle_eq_cornerAngle
#print axioms half_mul_sin_angle_eq_triArea
#print axioms dist_eq_sqrt_sum_sq
-- B.9, B.10 cotan Laplacian, mass
#print axioms cotan_formula
#print axioms cotEnergy_nonneg
#print axioms TriMesh.laplacian_isSymm
#print axioms TriMesh.laplacian_mulVec_one
#print axioms TriMesh.dotProduct_laplacian_mulVec_eq_sum_faceEnergy
#print axioms TriMesh.laplacian_posSemidef
#print axioms TriMesh.laplacian_mulVec_eq_zero_iff
#print axioms TriMesh.laplacian_interior_posDef
#print axioms TriMesh.sum_mass
#print axioms TriMesh.laplacian_scale
#print axioms TriMesh.mass_scale
-- B.11 Gauss–Bonnet
#print axioms TriMesh.gauss_bonnet
#print axioms TriMesh.sum_cone_eq_two_pi_eulerChar
-- Statements about the oracle's definitions
#print axioms Oracle.inSubdiff_iff_eq_boundedThreshold
#print axioms Oracle.inSubdiff_sub_div_iff_eq_regularizedMeasure
#print axioms Oracle.eq_regularizedMeasure_add_iff
#print axioms Oracle.regularizedMeasure_isPiecewiseAffine
#print axioms Oracle.cornerAngle_eq
#print axioms Ddg.Problem.objective_le_of_F_eq_zero
#print axioms Ddg.Problem.unique_of_minimizers
#print axioms Ddg.Problem.F_eq_zero_of_newton_step
#print axioms Ddg.Problem.Freg_eq_zero_of_newton_step
#print axioms Ddg.Problem.objective_le_add_gap
#print axioms Oracle.laplacian_symm
#print axioms Oracle.laplacian_quadratic_nonneg
#print axioms Oracle.gauss_bonnet
#print axioms Oracle.sum_lumpedMass
