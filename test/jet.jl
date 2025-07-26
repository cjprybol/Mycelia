# run from the root of the Mycelia package directory:
# julia --project=. --color=yes test/jet.jl

import JET

# Run analysis and capture the result
result = JET.report_package("Mycelia")

# Print summary
if isempty(JET.get_reports(result))
    println("\n✅ No type instabilities or errors detected by JET!")
else
    println("\n⚠️  Found $(length(JET.get_reports(result))) issues")
    display(JET.get_reports(result))
end