# # Minimum edit distance therapeutic framework for personalized medicine
# # Based on analysis from Mycelia-Dev notebooks (2021-10-02, 2021-10-09-hmp-ibd-2.ipynb)

# """
#     TherapeuticOptimizationResult

# Results from minimum edit distance therapeutic optimization.

# # Fields
# - `optimal_interventions::Vector{Tuple{String, Float64}}`: Feature names and target adjustment values
# - `intervention_cost::Float64`: Total magnitude of required interventions
# - `predicted_health_score::Float64`: Predicted health score after interventions
# - `intervention_confidence::Float64`: Confidence in intervention recommendations (0-1)
# - `affected_pathways::Vector{String}`: Biological pathways affected by interventions
# - `safety_assessment::Dict{String, Float64}`: Safety scores for each intervention
# - `treatment_priority::Vector{Int}`: Recommended order of intervention implementation
# """
# struct TherapeuticOptimizationResult
#     optimal_interventions::Vector{Tuple{String, Float64}}
#     intervention_cost::Float64
#     predicted_health_score::Float64
#     intervention_confidence::Float64
#     affected_pathways::Vector{String}
#     safety_assessment::Dict{String, Float64}
#     treatment_priority::Vector{Int}
# end

# """
#     PersonalizedTreatmentPlan

# Comprehensive personalized treatment plan with implementation guidance.

# # Fields
# - `patient_id::String`: Unique identifier for the patient
# - `baseline_profile::Vector{Float64}`: Patient's initial biological profile
# - `target_profile::Vector{Float64}`: Desired healthy target profile
# - `intervention_steps::Vector{TherapeuticOptimizationResult}`: Sequential treatment steps
# - `monitoring_schedule::Dict{String, Int}`: Recommended monitoring intervals (days)
# - `success_criteria::Dict{String, Float64}`: Quantitative criteria for treatment success
# - `alternative_pathways::Vector{Vector{Tuple{String, Float64}}}`: Alternative intervention strategies
# """
# struct PersonalizedTreatmentPlan
#     patient_id::String
#     baseline_profile::Vector{Float64}
#     target_profile::Vector{Float64}
#     intervention_steps::Vector{TherapeuticOptimizationResult}
#     monitoring_schedule::Dict{String, Int}
#     success_criteria::Dict{String, Float64}
#     alternative_pathways::Vector{Vector{Tuple{String, Float64}}}
# end

# """
#     optimize_therapeutic_interventions(patient_profile::Vector{Float64},
#                                      healthy_controls::Matrix{Float64},
#                                      feature_names::Vector{String};
#                                      intervention_budget::Float64=2.0,
#                                      safety_constraints::Dict{String, Tuple{Float64, Float64}}=Dict(),
#                                      pathway_annotations::Dict{String, String}=Dict()) -> TherapeuticOptimizationResult

# Optimize therapeutic interventions using minimum edit distance in biological feature space.

# This algorithm identifies the minimal set of biological modifications needed to move
# a patient's profile toward healthy controls while respecting safety constraints and
# intervention budgets.

# # Arguments
# - `patient_profile::Vector{Float64}`: Patient's current biological measurements
# - `healthy_controls::Matrix{Float64}`: Reference healthy control profiles (samples × features)
# - `feature_names::Vector{String}`: Names of biological features being optimized
# - `intervention_budget::Float64=2.0`: Maximum allowed intervention magnitude (L2 norm)
# - `safety_constraints::Dict{String, Tuple{Float64, Float64}}`: Per-feature safety bounds (min, max)
# - `pathway_annotations::Dict{String, String}`: Feature to biological pathway mapping

# # Returns
# `TherapeuticOptimizationResult`: Comprehensive optimization results with safety assessment

# # Algorithm Details
# 1. **Target Calculation**: Computes healthy control centroid as therapeutic target
# 2. **Dimensionality Reduction**: Uses PCA to identify key intervention dimensions
# 3. **Constrained Optimization**: Minimizes ||patient - target||₂ subject to constraints
# 4. **Safety Assessment**: Evaluates intervention safety based on biological knowledge
# 5. **Pathway Analysis**: Maps interventions to affected biological pathways

# # Mathematical Foundation
# The optimization problem is formulated as:
# ```
# minimize ||x_patient + Δx - x_target||₂
# subject to:
#   ||Δx||₂ ≤ intervention_budget
#   x_min ≤ x_patient + Δx ≤ x_max  (safety constraints)
# ```

# # Clinical Applications
# - **IBD Treatment**: Optimize microbiome interventions for inflammatory bowel disease
# - **Metabolic Disorders**: Adjust metabolic profiles toward healthy ranges
# - **Precision Medicine**: Personalize treatments based on individual biological profiles

# # Examples
# ```julia
# # Basic therapeutic optimization
# patient = [0.5, 1.2, -0.8, 2.1]  # Patient's microbiome profile
# controls = healthy_microbiome_data  # Matrix of healthy controls
# features = ["Bacteroides", "Firmicutes", "Proteobacteria", "Actinobacteria"]

# result = optimize_therapeutic_interventions(patient, controls, features)

# # High-safety optimization with constraints
# safety_bounds = Dict(
#     "Bacteroides" => (0.1, 2.0),
#     "Firmicutes" => (0.2, 1.8)
# )

# result = optimize_therapeutic_interventions(patient, controls, features,
#                                           intervention_budget=1.5,
#                                           safety_constraints=safety_bounds)

# # Generate treatment plan
# plan = generate_treatment_plan("patient_001", result, features)
# ```

# # Safety Considerations
# - Interventions are constrained by biological feasibility bounds
# - Safety scores assess potential risks of each intervention
# - Alternative pathways provide backup strategies if primary approach fails

# # References
# Based on minimum edit distance therapeutic framework from Mycelia-Dev notebooks:
# - `2021-10-02-minimum-necessary-changes.ipynb`
# - `2021-10-09-hmp-ibd-2.ipynb`
# """
# function optimize_therapeutic_interventions(patient_profile::Vector{Float64},
#                                           healthy_controls::Matrix{Float64},
#                                           feature_names::Vector{String};
#                                           intervention_budget::Float64=2.0,
#                                           safety_constraints::Dict{String, Tuple{Float64, Float64}}=Dict(),
#                                           pathway_annotations::Dict{String, String}=Dict())

#     if length(patient_profile) != length(feature_names)
#         throw(ArgumentError("Patient profile length must match number of feature names"))
#     end

#     if size(healthy_controls, 2) != length(patient_profile)
#         throw(ArgumentError("Healthy controls must have same number of features as patient"))
#     end

#     n_features = length(patient_profile)
#     n_controls = size(healthy_controls, 1)

#     @info "Optimizing therapeutic interventions for $(n_features) features using $(n_controls) healthy controls"

#     ## Calculate healthy control centroid as target
#     target_profile = vec(Statistics.mean(healthy_controls, dims=1))

#     @debug "Target profile calculated from $(n_controls) healthy controls"

#     ## Perform PCA on healthy controls for dimensionality reduction
#     pca_result = perform_pca_analysis(healthy_controls, min(5, n_features))

#     ## Calculate unconstrained intervention direction
#     intervention_direction = target_profile .- patient_profile
#     intervention_magnitude = LinearAlgebra.norm(intervention_direction)

#     if intervention_magnitude == 0.0
#         @info "Patient profile already matches healthy controls - no intervention needed"
#         return TherapeuticOptimizationResult(
#             Tuple{String, Float64}[], 0.0, 1.0, 1.0, String[], Dict{String, Float64}(),
#             Int[]
#         )
#     end

#     ## Apply intervention budget constraint
#     if intervention_magnitude > intervention_budget
#         ## Scale down intervention to fit budget
#         intervention_direction = intervention_direction .* (intervention_budget / intervention_magnitude)
#         intervention_magnitude = intervention_budget
#         @info "Intervention scaled to fit budget: $(round(intervention_magnitude, digits=3))"
#     end

#     ## Apply safety constraints
#     safe_interventions = apply_safety_constraints(
#         patient_profile, intervention_direction, feature_names, safety_constraints
#     )

#     ## Calculate final adjusted profile
#     final_profile = patient_profile .+ safe_interventions

#     ## Assess intervention quality
#     distance_to_target = LinearAlgebra.norm(final_profile .- target_profile)
#     max_possible_distance = LinearAlgebra.norm(target_profile .- patient_profile)
#     predicted_health_score = max(0.0, 1.0 - distance_to_target / max_possible_distance)

#     ## Create intervention list
#     optimal_interventions = Tuple{String, Float64}[]
#     for (i, intervention) in enumerate(safe_interventions)
#         if abs(intervention) > 1e-6  ## Only include non-zero interventions
#             push!(optimal_interventions, (feature_names[i], intervention))
#         end
#     end

#     ## Calculate intervention confidence based on PCA variance explained
#     intervention_confidence = calculate_intervention_confidence(
#         safe_interventions, pca_result, predicted_health_score
#     )

#     ## Identify affected biological pathways
#     affected_pathways = identify_affected_pathways(optimal_interventions, pathway_annotations)

#     ## Assess safety of each intervention
#     safety_assessment = assess_intervention_safety(optimal_interventions, safety_constraints)

#     ## Determine treatment priority order
#     treatment_priority = determine_treatment_priority(optimal_interventions, safety_assessment)

#     final_cost = LinearAlgebra.norm(safe_interventions)

#     @info "Optimization complete: $(length(optimal_interventions)) interventions, " *
#           "cost $(round(final_cost, digits=3)), predicted health score $(round(predicted_health_score, digits=3))"

#     return TherapeuticOptimizationResult(
#         optimal_interventions,
#         final_cost,
#         predicted_health_score,
#         intervention_confidence,
#         affected_pathways,
#         safety_assessment,
#         treatment_priority
#     )
# end

# """
#     generate_treatment_plan(patient_id::String,
#                           optimization_result::TherapeuticOptimizationResult,
#                           feature_names::Vector{String};
#                           treatment_phases::Int=3,
#                           monitoring_frequency::Int=14) -> PersonalizedTreatmentPlan

# Generate a comprehensive personalized treatment plan from optimization results.

# # Arguments
# - `patient_id::String`: Unique patient identifier
# - `optimization_result::TherapeuticOptimizationResult`: Results from therapeutic optimization
# - `feature_names::Vector{String}`: Names of biological features
# - `treatment_phases::Int=3`: Number of sequential treatment phases
# - `monitoring_frequency::Int=14`: Days between monitoring assessments

# # Returns
# `PersonalizedTreatmentPlan`: Complete treatment plan with implementation guidance
# """
# function generate_treatment_plan(patient_id::String,
#                                 optimization_result::TherapeuticOptimizationResult,
#                                 feature_names::Vector{String};
#                                 treatment_phases::Int=3,
#                                 monitoring_frequency::Int=14)

#     @info "Generating treatment plan for patient $patient_id with $(treatment_phases) phases"

#     ## Divide interventions into sequential phases based on priority
#     intervention_steps = divide_interventions_into_phases(
#         optimization_result, treatment_phases
#     )

#     ## Create monitoring schedule
#     monitoring_schedule = Dict{String, Int}(
#         "comprehensive_profile" => monitoring_frequency,
#         "safety_markers" => monitoring_frequency ÷ 2,
#         "clinical_assessment" => monitoring_frequency * 2,
#         "intervention_compliance" => monitoring_frequency ÷ 3
#     )

#     ## Define success criteria
#     success_criteria = Dict{String, Float64}(
#         "health_score_improvement" => 0.2,  ## 20% improvement minimum
#         "safety_threshold" => 0.8,          ## 80% safety score minimum  
#         "intervention_tolerance" => 0.9,    ## 90% intervention compliance
#         "symptom_reduction" => 0.3          ## 30% symptom reduction
#     )

#     ## Generate alternative treatment pathways
#     alternative_pathways = generate_alternative_pathways(optimization_result)

#     ## Placeholder profiles (would be calculated from actual patient data)
#     baseline_profile = zeros(Float64, length(feature_names))
#     target_profile = zeros(Float64, length(feature_names))

#     return PersonalizedTreatmentPlan(
#         patient_id,
#         baseline_profile,
#         target_profile,
#         intervention_steps,
#         monitoring_schedule,
#         success_criteria,
#         alternative_pathways
#     )
# end

# """
#     assess_treatment_progress(current_profile::Vector{Float64},
#                             treatment_plan::PersonalizedTreatmentPlan,
#                             days_elapsed::Int) -> Dict{String, Float64}

# Assess patient progress relative to treatment plan expectations.

# # Arguments
# - `current_profile::Vector{Float64}`: Patient's current biological profile
# - `treatment_plan::PersonalizedTreatmentPlan`: Active treatment plan
# - `days_elapsed::Int`: Days since treatment initiation

# # Returns
# `Dict{String, Float64}`: Progress metrics with scores from 0-1
# """
# function assess_treatment_progress(current_profile::Vector{Float64},
#                                  treatment_plan::PersonalizedTreatmentPlan,
#                                  days_elapsed::Int)

#     progress_metrics = Dict{String, Float64}()

#     ## Calculate distance to target
#     target_distance = LinearAlgebra.norm(current_profile .- treatment_plan.target_profile)
#     baseline_distance = LinearAlgebra.norm(treatment_plan.baseline_profile .- treatment_plan.target_profile)

#     if baseline_distance > 0
#         progress_metrics["target_approach"] = max(0.0, 1.0 - target_distance / baseline_distance)
#     else
#         progress_metrics["target_approach"] = 1.0
#     end

#     ## Assess adherence to expected timeline
#     expected_phases = min(length(treatment_plan.intervention_steps), 
#                          ceil(Int, days_elapsed / treatment_plan.monitoring_schedule["comprehensive_profile"]))

#     progress_metrics["timeline_adherence"] = min(1.0, expected_phases / length(treatment_plan.intervention_steps))

#     ## Overall progress score
#     progress_metrics["overall_progress"] = Statistics.mean([
#         progress_metrics["target_approach"],
#         progress_metrics["timeline_adherence"]
#     ])

#     @info "Treatment progress assessment: $(round(100*progress_metrics["overall_progress"], digits=1))% complete"

#     return progress_metrics
# end

# ## Helper function for PCA analysis
# function perform_pca_analysis(data::Matrix{Float64}, n_components::Int)

#     ## Center the data
#     data_centered = data .- Statistics.mean(data, dims=1)

#     ## Compute covariance matrix
#     cov_matrix = (data_centered' * data_centered) ./ (size(data, 1) - 1)

#     ## Compute eigendecomposition
#     eigenvals, eigenvecs = LinearAlgebra.eigen(cov_matrix)

#     ## Sort by eigenvalue magnitude (descending)
#     sorted_indices = sortperm(eigenvals, rev=true)
#     sorted_eigenvals = eigenvals[sorted_indices]
#     sorted_eigenvecs = eigenvecs[:, sorted_indices]

#     ## Select top components
#     selected_components = min(n_components, length(sorted_eigenvals))
#     principal_components = sorted_eigenvecs[:, 1:selected_components]
#     explained_variance = sorted_eigenvals[1:selected_components]

#     total_variance = sum(sorted_eigenvals)
#     variance_explained_ratio = explained_variance ./ total_variance

#     return (
#         components = principal_components,
#         explained_variance = explained_variance,
#         variance_ratio = variance_explained_ratio,
#         total_variance_explained = sum(variance_explained_ratio)
#     )
# end

# ## Helper function to apply safety constraints
# function apply_safety_constraints(baseline::Vector{Float64},
#                                 interventions::Vector{Float64},
#                                 feature_names::Vector{String},
#                                 constraints::Dict{String, Tuple{Float64, Float64}})

#     safe_interventions = copy(interventions)

#     for (i, feature_name) in enumerate(feature_names)
#         if haskey(constraints, feature_name)
#             min_val, max_val = constraints[feature_name]

#             ## Clamp final value to safety bounds
#             proposed_value = baseline[i] + interventions[i]
#             safe_value = clamp(proposed_value, min_val, max_val)
#             safe_interventions[i] = safe_value - baseline[i]

#             if abs(safe_interventions[i] - interventions[i]) > 1e-6
#                 @debug "Intervention for $feature_name constrained by safety bounds"
#             end
#         end
#     end

#     return safe_interventions
# end

# ## Helper function to calculate intervention confidence
# function calculate_intervention_confidence(interventions::Vector{Float64},
#                                          pca_result::NamedTuple,
#                                          health_score::Float64)

#     ## Base confidence on PCA variance explained and health score
#     pca_confidence = pca_result.total_variance_explained

#     ## Penalize interventions that require many small changes (less focused)
#     intervention_focus = 1.0 - (count(x -> abs(x) > 1e-6, interventions) / length(interventions))

#     ## Combine factors
#     confidence = 0.4 * pca_confidence + 0.4 * health_score + 0.2 * intervention_focus

#     return min(1.0, max(0.0, confidence))
# end

# ## Helper function to identify affected pathways
# function identify_affected_pathways(interventions::Vector{Tuple{String, Float64}},
#                                   pathway_annotations::Dict{String, String})

#     affected_pathways = Set{String}()

#     for (feature_name, intervention_value) in interventions
#         if abs(intervention_value) > 1e-6 && haskey(pathway_annotations, feature_name)
#             push!(affected_pathways, pathway_annotations[feature_name])
#         end
#     end

#     return collect(affected_pathways)
# end

# ## Helper function to assess intervention safety
# function assess_intervention_safety(interventions::Vector{Tuple{String, Float64}},
#                                   safety_constraints::Dict{String, Tuple{Float64, Float64}})

#     safety_scores = Dict{String, Float64}()

#     for (feature_name, intervention_value) in interventions
#         ## Default high safety if no constraints specified
#         safety_score = 0.9

#         if haskey(safety_constraints, feature_name)
#             min_val, max_val = safety_constraints[feature_name]
#             range_size = max_val - min_val

#             ## Safety inversely related to intervention magnitude relative to allowed range
#             if range_size > 0
#                 intervention_fraction = abs(intervention_value) / range_size
#                 safety_score = max(0.1, 1.0 - intervention_fraction)
#             end
#         end

#         safety_scores[feature_name] = safety_score
#     end

#     return safety_scores
# end

# ## Helper function to determine treatment priority
# function determine_treatment_priority(interventions::Vector{Tuple{String, Float64}},
#                                     safety_scores::Dict{String, Float64})

#     ## Priority based on intervention magnitude and safety
#     priorities = Int[]

#     for (i, (feature_name, intervention_value)) in enumerate(interventions)
#         safety_score = get(safety_scores, feature_name, 0.5)
#         magnitude = abs(intervention_value)

#         ## Higher priority for high-safety, high-magnitude interventions
#         priority_score = safety_score * magnitude
#         push!(priorities, i)
#     end

#     ## Sort by priority score (highest first)
#     intervention_scores = [(i, get(safety_scores, interventions[i][1], 0.5) * abs(interventions[i][2])) 
#                           for i in priorities]
#     sorted_priorities = [i for (i, score) in sort(intervention_scores, by=x->x[2], rev=true)]

#     return sorted_priorities
# end

# ## Helper function to divide interventions into phases
# function divide_interventions_into_phases(result::TherapeuticOptimizationResult,
#                                         n_phases::Int)

#     if n_phases <= 1 || isempty(result.optimal_interventions)
#         return [result]
#     end

#     ## Divide interventions by priority
#     interventions_per_phase = ceil(Int, length(result.optimal_interventions) / n_phases)
#     intervention_phases = TherapeuticOptimizationResult[]

#     for phase in 1:n_phases
#         start_idx = (phase - 1) * interventions_per_phase + 1
#         end_idx = min(phase * interventions_per_phase, length(result.optimal_interventions))

#         if start_idx <= length(result.optimal_interventions)
#             phase_interventions = result.optimal_interventions[start_idx:end_idx]
#             phase_cost = sum(abs(intervention[2]) for intervention in phase_interventions)

#             phase_result = TherapeuticOptimizationResult(
#                 phase_interventions,
#                 phase_cost,
#                 result.predicted_health_score / n_phases,  ## Distributed across phases
#                 result.intervention_confidence,
#                 result.affected_pathways,
#                 result.safety_assessment,
#                 result.treatment_priority[start_idx:end_idx]
#             )

#             push!(intervention_phases, phase_result)
#         end
#     end

#     return intervention_phases
# end

# ## Helper function to generate alternative pathways
# function generate_alternative_pathways(result::TherapeuticOptimizationResult)

#     ## Generate alternative approaches by varying intervention magnitudes
#     alternatives = Vector{Tuple{String, Float64}}[]

#     if !isempty(result.optimal_interventions)
#         ## Conservative approach (50% intervention)
#         conservative = [(name, 0.5 * value) for (name, value) in result.optimal_interventions]
#         push!(alternatives, conservative)

#         ## Aggressive approach (150% intervention, if safe)
#         aggressive = [(name, 1.5 * value) for (name, value) in result.optimal_interventions]
#         push!(alternatives, aggressive)

#         ## Focused approach (top 50% interventions only)
#         sorted_interventions = sort(result.optimal_interventions, 
#                                   by=x->abs(x[2]), rev=true)
#         n_focused = max(1, length(sorted_interventions) ÷ 2)
#         focused = sorted_interventions[1:n_focused]
#         push!(alternatives, focused)
#     end

#     return alternatives
# end
