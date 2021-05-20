using Coverage
# process '*.cov' files
coverage = process_folder() # defaults to src/; alternatively, supply the folder name as argument
# process '*.info' files
coverage = merge_coverage_counts(coverage, filter!(
    let prefixes = (joinpath(pwd(), "src", ""))
        c -> any(p -> startswith(c.filename, p), prefixes)
    end,
    LCOV.readfolder("test")))
# Get total coverage for all Julia files
covered_lines, total_lines = get_summary(coverage)
percent_coverage = round((covered_lines / total_lines) * 100, digits = 3)
println("total_coverage = $(percent_coverage)")
# Or process a single file
# @show get_summary(process_file(joinpath("src", "Mycelia.jl")))