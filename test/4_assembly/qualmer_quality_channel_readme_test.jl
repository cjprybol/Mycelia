# From the Mycelia base directory, run the tests with:
#
# ```bash
# julia --project=. -e 'include("test/4_assembly/qualmer_quality_channel_readme_test.jl")'
# ```
#
# td-7k68: the quality-channel probe README must agree with the sibling
# Stage-C TSV and with the ADR errata. The committed probe TSVs must stay
# unchanged.

import Test
import SHA
import DelimitedFiles

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const PROBE_DIR = joinpath(REPO_ROOT, "benchmarking", "results", "qualmer_quality_channel_probe")
const README_PATH = joinpath(PROBE_DIR, "README.md")
const STAGE_C_PATH = joinpath(PROBE_DIR, "qualmer_corrector_stageC_end_to_end.tsv")

# sha256 of the committed probe TSVs at master 57cf00d1. The README fix must
# not rewrite these files.
const PROBE_TSV_SHA256 = Dict(
    "qualmer_corrector_quality_verdicts.tsv" =>
        "7c262a79352e864e443e7d84e99d16cb207febbe1552c5c8fe82aec3c9f6b2c9",
    "qualmer_corrector_stageA_unit.tsv" =>
        "dc38820803fdd9df2698aafebf80b9f08dec3e1db746dae8e96ddc2938d4f6dd",
    "qualmer_corrector_stageB_decisions.tsv" =>
        "e05ffc20e76c6f48b406df83d7268cc3448677e0794790c56cb0a7d6dde78254",
    "qualmer_corrector_stageC_end_to_end.tsv" =>
        "16d5ae3f93c68081ee50c9b7e3a977532e873b517c6c45d64c9b554f647d287e",
    "qualmer_joint_quality_calibration_bins.tsv" =>
        "b12d8471d7e61a1d552a80a4cf1f40a4ec1f7f003a91a6b54e68f92d28979ea8",
    "qualmer_joint_quality_saturation.tsv" =>
        "d62d929263e7609460fdb19fb48fb9ae7f159616c613a8c68ab40dacde1d6777",
    "qualmer_quality_channel_probe_runs.tsv" =>
        "8cce0075a3a06ebfb0278f2197065f8e8093b14ce37700a237aa3e69ae5e91d2",
    "qualmer_quality_channel_probe_verdicts.tsv" =>
        "9b7f48c098995f82a432f363081c6ed352b7e1240b847a80982121182b60915d",
)

function file_sha256(path)
    return bytes2hex(SHA.sha256(read(path)))
end

function parse_fraction_vector(s::AbstractString)
    inner = strip(s)
    inner = replace(inner, r"[\[\]]" => "")
    return [parse(Float64, strip(p)) for p in split(inner, ",") if !isempty(strip(p))]
end

function stage_c_by_chemistry()
    raw = DelimitedFiles.readdlm(STAGE_C_PATH, '\t', String)
    header = raw[1, :]
    col = Dict(header[i] => i for i in eachindex(header))
    out = Dict{String, Vector{Dict{String, String}}}()
    for r in 2:size(raw, 1)
        row = Dict(String(header[i]) => String(raw[r, i]) for i in eachindex(header))
        push!(get!(out, row["chemistry"], Dict{String, String}[]), row)
    end
    return out
end

function decoder_ever_ran(rows)
    return any(maximum(parse_fraction_vector(r["decode_fraction_per_pass"])) > 0 for r in rows)
end

Test.@testset "quality-channel README agrees with ADR errata (td-7k68)" begin
    Test.@test isfile(README_PATH)
    Test.@test isfile(STAGE_C_PATH)
    readme = read(README_PATH, String)

    Test.@testset "committed probe TSVs are unchanged" begin
        for (name, expected) in PROBE_TSV_SHA256
            path = joinpath(PROBE_DIR, name)
            Test.@test isfile(path)
            Test.@test file_sha256(path) == expected
        end
    end

    by_chem = stage_c_by_chemistry()
    Test.@test issetequal(keys(by_chem), Set(["illumina", "pacbio", "ont"]))
    cheap = Dict(c => only(unique(parse(Int, r["cheap_corrections_total"]) for r in rows))
                 for (c, rows) in by_chem)
    Test.@test cheap["illumina"] == 74
    Test.@test cheap["pacbio"] == 329
    Test.@test cheap["ont"] == 120
    Test.@test decoder_ever_ran(by_chem["illumina"])
    Test.@test !decoder_ever_ran(by_chem["pacbio"])
    Test.@test !decoder_ever_ran(by_chem["ont"])

    # Illumina oracle / constant40 fractions used by ADR erratum E2.
    function decode_for(chem, cond)
        rows = filter(r -> r["condition"] == cond, by_chem[chem])
        return parse_fraction_vector(only(rows)["decode_fraction_per_pass"])
    end
    illumina_oracle = decode_for("illumina", "oracle")
    illumina_flat = decode_for("illumina", "constant40")
    Test.@test illumina_oracle[1] == 0.0
    Test.@test illumina_oracle[2] > 0.0
    Test.@test maximum(illumina_flat) > maximum(illumina_oracle[2:end]) ||
              illumina_flat[end] > illumina_oracle[end]

    Test.@testset "Stage C decoder finding is not the retracted all-zero claim" begin
        # The retracted claim said every chemistry, every condition was [0,0,0,0]
        # and that the decoder never executed once. The TSV contradicts that for
        # Illumina.
        Test.@test !occursin(r"decode_fraction_per_pass was `\[0\.0, 0\.0, 0\.0, 0\.0\]` in every condition", readme)
        Test.@test !occursin(r"never executed once", readme)
        Test.@test !occursin(r"NEVER EXECUTED", readme)
        # The Stage-C table must not list Illumina as decoder-ever-ran = no.
        Test.@test !occursin(r"(?i)illumina\s*\|\s*yes\s*\|\s*yes\s*\|\s*\*\*no\*\*", readme)
        # The corrected finding: the decoder ran, on a quality-dependent fraction.
        Test.@test occursin(r"(?i)decoder ran", readme)
        Test.@test occursin("0.058", readme) || occursin("5.8%", readme) || occursin("0.0583", readme)
        # Quality-dependent rates from ADR E2 (flat ~3.6% vs oracle ~1.7%).
        Test.@test occursin("3.6%", readme) || occursin("0.036", readme)
        Test.@test occursin("1.7%", readme) || occursin("0.017", readme)
        # PacBio / ONT remain genuinely all-zero; assemblies still byte-identical.
        Test.@test occursin(r"(?i)(pacbio|PacBio).{0,80}(all-zero|all zero|did not run)", readme) ||
                   occursin(r"(?i)genuinely all-zero", readme)
        # Section 6 must not still say Stage C never ran the decoder.
        Test.@test !occursin(r"Stage C never ran the decoder", readme)
    end

    Test.@testset "120 cheap corrections are not presented as all-chemistry" begin
        Test.@test !occursin(r"only corrections applied were 120", readme)
        # Per-chemistry totals from the sibling TSV must appear.
        Test.@test occursin(r"74", readme)
        Test.@test occursin(r"329", readme)
        Test.@test occursin(r"120", readme)
        # 120, if mentioned as a correction count, must be tied to ONT.
        if occursin(r"120", readme)
            Test.@test occursin(r"(?i)ont.{0,40}120|120.{0,40}ont", readme)
        end
    end

    Test.@testset "withdrawn repeat-collapse claim is not stated as fact" begin
        # ADR §5 / §7 withdraws the claim that quality distinguishes
        # repeat-collapsed k-mers. The README must not still make that claim.
        Test.@test !occursin(r"repeat-collapsed k-mers \*\*also\*\* recur", readme)
        Test.@test !occursin(r"they recur at _low_ quality, so\s+their joint score stays low", readme)
    end

    Test.@testset "255 clamp is not presented as the feature working" begin
        Test.@test !occursin(r"400\s*→\s*clamped 255", readme)
        Test.@test !occursin(r"400 -> clamped 255", readme)
        # ADR §3 / §6: saturation at ordinary Illumina depth; Q255 is not a
        # meaningful claim.
        Test.@test occursin("97.6%", readme) || occursin(r"97\.6", readme)
        Test.@test occursin(r"Q255 is not a meaningful claim", readme) ||
                   occursin(r"not a meaningful claim", readme)
    end

    Test.@testset "production-caller list is not stale" begin
        # combine_phred_scores is called by get_vertex_joint_quality;
        # get_vertex_min_quality is called by filter_by_quality.
        Test.@test !occursin(r"Zero non-test, non-docstring call sites anywhere in `src/`", readme)
        Test.@test occursin("get_vertex_joint_quality", readme)
        Test.@test occursin("combine_phred_scores", readme)
    end

    Test.@testset "path:LINE citations that remain still resolve" begin
        # Prefer file+symbol. Any leftover path:LINE cite must still point at
        # a line that exists. Known-broken cites from the pre-errata README
        # must be gone.
        Test.@test !occursin("src/rhizomorph/assembly.jl:4057", readme)
        Test.@test !occursin("src/rhizomorph/assembly.jl:43,47", readme)
        Test.@test !occursin("src/rhizomorph/assembly.jl:2461", readme)
        Test.@test !occursin("src/iterative-assembly.jl:6190", readme)
        Test.@test !occursin("quality-functions.jl:271", readme)
        rx = r"`((?:src|docs|planning-docs)/[\w./-]+\.\w+):(\d+)(?:-(\d+))?(?:,(\d+))?`"
        repo = REPO_ROOT
        for m in eachmatch(rx, readme)
            path = joinpath(repo, m.captures[1])
            Test.@test isfile(path)
            nlines = count(==('\n'), read(path, String)) + 1
            for cap in m.captures[2:end]
                cap === nothing && continue
                ln = parse(Int, cap)
                Test.@test 1 <= ln <= nlines
            end
        end
        # The previously rot-prone cites are now file+symbol.
        Test.@test occursin("_qualmer_vertex_quality_scores", readme)
        Test.@test occursin("src/rhizomorph/assembly.jl", readme)
        Test.@test occursin("_qualmer_vertex_score", readme)
        Test.@test occursin("_rhizomorph_joint_probability", readme)
        Test.@test occursin("src/iterative-assembly.jl", readme)
    end
end
