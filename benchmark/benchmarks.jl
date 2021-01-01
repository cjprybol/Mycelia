# all of these must be added as a test dependency!!
using BenchmarkTools
using Random
using Eisenia

const SUITE = BenchmarkGroup()

SUITE["utf8"] = BenchmarkGroup(["string", "unicode"])
teststr = String(join(rand(MersenneTwister(1), 'a':'d', 10^4)))
SUITE["utf8"]["replace"] = @benchmarkable replace($teststr, "a" => "b")
SUITE["utf8"]["join"] = @benchmarkable join($teststr, $teststr)
SUITE["utf8"]["plots"] = BenchmarkGroup()

SUITE["trigonometry"] = BenchmarkGroup(["math", "triangles"])
SUITE["trigonometry"]["circular"] = BenchmarkGroup()
for f in (sin, cos, tan)
    for x in (0.0, pi)
        SUITE["trigonometry"]["circular"][string(f), x] = @benchmarkable ($f)($x)
    end
end

SUITE["trigonometry"]["hyperbolic"] = BenchmarkGroup()
for f in (sin, cos, tan)
    for x in (0.0, pi)
        SUITE["trigonometry"]["hyperbolic"][string(f), x] = @benchmarkable ($f)($x)
    end
end

# any time that I update the benchmarks or want to modify the tunings,
# just delete the tune.json file and then re run this script
# https://github.com/JuliaCI/BenchmarkTools.jl/blob/master/doc/manual.md#caching-parameters
tuning_file = "$(@__DIR__)/tune.json"
if !isfile(tuning_file)
    
    # tune the suite to configure benchmark parameters
    tune!(SUITE);

    # save the suite's parameters using a thin wrapper
    # over JSON (this wrapper maintains compatibility
    # across BenchmarkTools versions)
    BenchmarkTools.save("tune.json", params(SUITE));
end
