using Eisenia
using PkgBenchmark

# compares this current state against master branch
results = judge(Eisenia, "master")
export_markdown(stdout, results)