import Test
import Mycelia

@eval Mycelia begin
    function add_bioconda_env(pkg::AbstractString; force = false, quiet = false)
        return nothing
    end
end

function write_fake_conda_runner_additional(path::AbstractString)
    script = """
#!/usr/bin/env bash
set -euo pipefail

if [[ "\${1:-}" == "run" ]]; then
    shift
fi

if [[ "\${1:-}" == "--live-stream" ]]; then
    shift
fi

if [[ "\${1:-}" == "-n" ]]; then
    shift 2
fi

tool="\${1:-}"
if [[ -n "\$tool" ]]; then
    shift
fi

if [[ "\$tool" == "blastdbcmd" ]]; then
    db=""
    metadata=0
    tax_info=0
    info=0
    entry_batch=""
    taxidlist=""

    while [[ \$# -gt 0 ]]; do
        case "\$1" in
            -db)
                db="\${2:-}"
                shift 2
                ;;
            -metadata)
                metadata=1
                shift
                ;;
            -tax_info)
                tax_info=1
                shift
                ;;
            -info)
                info=1
                shift
                ;;
            -entry)
                shift 2
                ;;
            -entry_batch)
                entry_batch="\${2:-}"
                shift 2
                ;;
            -taxidlist)
                taxidlist="\${2:-}"
                shift 2
                ;;
            -outfmt)
                shift 2
                ;;
            *)
                shift
                ;;
        esac
    done

    if [[ \$metadata -eq 1 ]]; then
        case "\$db" in
            nt)
                printf '{"dbtype":"Nucleotide","last-updated":"2024-01-02T03:04:05"}'
                ;;
            nr)
                printf '{"dbtype":"Protein","last-updated":"2024-02-03T04:05:06"}'
                ;;
            weird)
                printf '{"dbtype":"Unsupported","last-updated":"2024-03-04T05:06:07"}'
                ;;
        esac
        exit 0
    fi

    if [[ \$tax_info -eq 1 ]]; then
        printf '# comment\\n'
        printf '2\\tBacillus subtilis\\tbacterium\\tBacteria\\tBacteria\\t1\\n'
        exit 0
    fi

    if [[ \$info -eq 1 ]]; then
        printf 'Database: %s\\n' "\$db"
        exit 0
    fi

    if [[ -n "\$entry_batch" ]]; then
        printf '>entries\\n'
        tr '\\n' ',' < "\$entry_batch"
        printf '\\n'
        exit 0
    fi

    if [[ -n "\$taxidlist" ]]; then
        printf '>taxids\\n'
        tr '\\n' ',' < "\$taxidlist"
        printf '\\n'
        exit 0
    fi

    printf '>all\\nSEQUENCE\\n'
    exit 0
fi

if [[ "\$tool" == "pigz" ]]; then
    cat
    exit 0
fi

printf 'unexpected tool: %s\\n' "\$tool" >&2
exit 1
"""

    write(path, script)
    chmod(path, 0o755)
    return path
end

function with_fake_conda_runner_additional(f::Function)
    mktempdir() do dir
        runner_path = Mycelia.CONDA_RUNNER
        backup_path = joinpath(dir, "conda-runner-backup")
        runner_existed = isfile(runner_path)
        if runner_existed
            cp(runner_path, backup_path; force = true)
        else
            mkpath(dirname(runner_path))
        end
        write_fake_conda_runner_additional(runner_path)
        try
            return f()
        finally
            if runner_existed
                mv(backup_path, runner_path; force = true)
            else
                rm(runner_path; force = true)
            end
        end
    end
end

Test.@testset "Reference Database Additional Coverage" begin
    with_fake_conda_runner_additional() do
        mktempdir() do dir
            cd(dir) do
                metadata = Mycelia.get_blastdb_metadata(blastdb = "nt")
                Test.@test metadata["dbtype"] == "Nucleotide"
                Test.@test metadata["last-updated"] == "2024-01-02T03:04:05"

                blast_info_process = Mycelia.get_blastdb_info(blastdb = "nr")
                Test.@test blast_info_process isa Base.Process
                Test.@test success(blast_info_process)

                write("existing.fna.gz", "preserved")
                cached_out = Mycelia.blastdb_to_fasta(
                    blastdb = "nt",
                    outfile = "existing.fna.gz",
                    force = false
                )
                Test.@test cached_out == "existing.fna.gz"
                Test.@test read(cached_out, String) == "preserved"

                nucleotide_out = Mycelia.blastdb_to_fasta(
                    blastdb = "nt",
                    force = true
                )
                Test.@test nucleotide_out == "nt.2024-01-02.fna.gz"
                Test.@test read(nucleotide_out, String) == ">all\nSEQUENCE\n"

                entries_out = Mycelia.blastdb_to_fasta(
                    blastdb = "nt",
                    entries = ["seqA", "seqB"],
                    outfile = joinpath("nested", "entries.fna.gz"),
                    force = true
                )
                Test.@test entries_out == joinpath("nested", "entries.fna.gz")
                Test.@test read(entries_out, String) == ">entries\nseqA,seqB,\n"

                taxids_out = Mycelia.blastdb_to_fasta(
                    blastdb = "nr",
                    taxids = [2, 2157],
                    outfile = joinpath("nested", "taxids.faa.gz"),
                    force = true
                )
                Test.@test taxids_out == joinpath("nested", "taxids.faa.gz")
                Test.@test read(taxids_out, String) == ">taxids\n2,2157,\n"

                write("forced.faa.gz", "stale")
                protein_out = Mycelia.blastdb_to_fasta(
                    blastdb = "nr",
                    outfile = "forced.faa.gz",
                    force = true,
                    max_cores = 2
                )
                Test.@test protein_out == "forced.faa.gz"
                Test.@test read(protein_out, String) == ">all\nSEQUENCE\n"

                Test.@test_throws ErrorException Mycelia.blastdb_to_fasta(
                    blastdb = "weird",
                    force = true
                )
            end
        end

        default_tax_info = Mycelia.get_blastdb_tax_info(blastdb = "nr")
        Test.@test Mycelia.DataFrames.nrow(default_tax_info) == 1
        Test.@test default_tax_info[1, "taxid"] == 2

        entries_tax_info = Mycelia.get_blastdb_tax_info(
            blastdb = "nr",
            entries = ["seqA", "seqB"]
        )
        Test.@test Mycelia.DataFrames.nrow(entries_tax_info) == 1
        Test.@test entries_tax_info[1, "BLAST name"] == "Bacteria"

        taxid_tax_info = Mycelia.get_blastdb_tax_info(
            blastdb = "nr",
            taxids = [2, 2157]
        )
        Test.@test Mycelia.DataFrames.nrow(taxid_tax_info) == 1
        Test.@test taxid_tax_info[1, "scientific name"] == "Bacillus subtilis"

        Test.@test_throws ErrorException Mycelia.get_blastdb_tax_info(
            blastdb = "nr",
            entries = ["seqA"],
            taxids = [2]
        )
    end
end
