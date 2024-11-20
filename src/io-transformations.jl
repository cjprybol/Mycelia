# Save the distance matrix to a JLD2 file
function save_matrix_jld2(;matrix, filename)
    if !isfile(filename) || (filesize(filename) == 0)
        JLD2.@save filename matrix
    else
        @warn "$(filename) already exists and is non-empty, skipping..."
    end
    return filename
end

# Load the distance matrix from a JLD2 file
function load_matrix_jld2(filename)
    return JLD2.load(filename, "matrix")
end

function load_jld2(filename)
    return JLD2.load(filename)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function genbank_to_fasta(;genbank, fasta=genbank * ".fna", force=false)
    add_bioconda_env("emboss")
    if !isfile(fasta) || force
        run(`$(Mycelia.CONDA_RUNNER) run -n emboss --live-stream seqret $(genbank) fasta:$(fasta)`)
    end
end


"""
    function ncbi_genome_download_accession(;
            accession,
            outdir = pwd(),
            outpath = joinpath(outdir, accession * ".zip"),
            include_string = "genome"
        )

Download an accession using NCBI datasets command line tool

the .zip download output to outpath will be unzipped

returns the outfolder

ncbi's default include string is 
include_string = "gff3,rna,cds,protein,genome,seq-report"
"""
function ncbi_genome_download_accession(;
        accession,
        outdir = pwd(),
        outpath = joinpath(outdir, accession * ".zip"),
        include_string = "genome"
    )
    outfolder = joinpath(outdir, accession)
    if !isdir(outfolder)
        add_bioconda_env("ncbi-datasets-cli")
        if isfile(outpath)
            @info "$(outpath) already exists, skipping download..."
        else
            mkpath(outdir)
            run(`$(Mycelia.CONDA_RUNNER) run --live-stream -n ncbi-datasets-cli datasets download genome accession $(accession) --include $(include_string) --filename $(outpath) --no-progressbar`)
        end
        run(`unzip -q -d $(outfolder) $(outpath)`)
    end
    final_outfolder = joinpath(outfolder, "ncbi_dataset", "data", accession)
    isfile(outpath) && rm(outpath)
    return final_outfolder
end


"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function load_graph(file)
    return FileIO.load(file)["graph"]
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function parse_transterm_output(transterm_output)
    
   #     3. FORMAT OF THE TRANSTERM OUTPUT

#     The organism's genes are listed sorted by their end coordinate and terminators
#     are output between them. A terminator entry looks like this:

#         TERM 19  15310 - 15327  -      F     99      -12.7 -4.0 |bidir
#         (name)   (start - end)  (sense)(loc) (conf) (hp) (tail) (notes)

#     where 'conf' is the overall confidence score, 'hp' is the hairpin score, and
#     'tail' is the tail score. 'Conf' (which ranges from 0 to 100) is what you
#     probably want to use to assess the quality of a terminator. Higher is better.
#     The confidence, hp score, and tail scores are described in the paper cited
#     above.  'Loc' gives type of region the terminator is in:

#         'G' = in the interior of a gene (at least 50bp from an end),
#         'F' = between two +strand genes,
#         'R' = between two -strand genes,
#         'T' = between the ends of a +strand gene and a -strand gene,
#         'H' = between the starts of a +strand gene and a -strand gene,
#         'N' = none of the above (for the start and end of the DNA)

#     Because of how overlapping genes are handled, these designations are not
#     exclusive. 'G', 'F', or 'R' can also be given in lowercase, indicating that
#     the terminator is on the opposite strand as the region.  Unless the
#     --all-context option is given, only candidate terminators that appear to be in
#     an appropriate genome context (e.g. T, F, R) are output. 

#     Following the TERM line is the sequence of the hairpin and the 5' and 3'
#     tails, always written 5' to 3'.
    
    transterm_table = DataFrames.DataFrame()
    chromosome = ""
    for line in Iterators.filter(x -> occursin(r"^\s*(SEQUENCE|TERM)", x), eachline(transterm_output))
        line = strip(line)
        if occursin(r"^SEQUENCE", line)
            chromosome = split(line)[2]
        else
            transterm_regex = r"(TERM \d+)\s+(\d+) - (\d+)\s+(\S)\s+(\w+)\s+(\d+)\s+(\S+)\s+(\S+)\s+\|(.*)"
            term_id, start, stop, strand, location, confidence, hairpin_score, tail_score, notes = match(transterm_regex, line).captures
            notes = strip(notes)
            row = (;chromosome, term_id, start, stop, strand, location, confidence, hairpin_score, tail_score, notes)
            push!(transterm_table, row)
        end
    end
    return transterm_table
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function transterm_output_to_gff(transterm_output)
    transterm_table = parse_transterm_output(transterm_output)
    transterm_table[!, "source"] .= "transterm"
    transterm_table[!, "type"] .= "terminator"
    transterm_table[!, "phase"] .= "."
    transterm_table[!, "attributes"] = map(x -> "label=" * replace(x, " " => "_"), transterm_table[!, "term_id"])
    DataFrames.rename!(transterm_table,
        ["chromosome" => "#seqid",
            "stop" => "end",
            "confidence" => "score",
        ]
    )
    transterm_table = transterm_table[!, ["#seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]]
    transterm_gff = write_gff(gff=transterm_table, outfile=transterm_output * ".gff")
    uCSV.write(transterm_gff, transterm_table, delim='\t')
    return transterm_gff
end



# TODO: switch to using GenomicAnnotations if GFF3 package isn't updated
# PR -> https://github.com/BioJulia/GFF3.jl/pull/12
"""
$(DocStringExtensions.TYPEDSIGNATURES)

Get dna (db = "nuccore") or protein (db = "protein") sequences from NCBI
or get fasta directly from FTP site
"""
function get_gff(;db=""::String, accession=""::String, ftp=""::String)
    if !isempty(db) && !isempty(accession)
        # API will block if we request more than 3 times per second, so set a 1/2 second sleep to set max of 2 requests per second when looping
        sleep(0.5)
        url = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=$(db)&report=gff3&id=$(accession)"
        return IOBuffer(HTTP.get(url).body)
    elseif !isempty(ftp)
        return CodecZlib.GzipDecompressorStream(IOBuffer(HTTP.get(ftp).body))
    else
        @error "invalid call"
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Get dna (db = "nuccore") or protein (db = "protein") sequences from NCBI
or get fasta directly from FTP site
"""
function get_genbank(;db=""::String, accession=""::String, ftp=""::String)
    if !isempty(db) && !isempty(accession)
        # API will block if we request more than 3 times per second, so set a 1/2 second sleep to set max of 2 requests per second when looping
        sleep(0.5)
        # url = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=$(db)&report=genbank&id=$(accession)&rettype=text"
        url = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=$(db)&report=genbank&id=$(accession)&retmode=text"
        # readgbk can't read from an io buffer, so need to download to a temp file
        # outfile = tempname()
        # open(outfile, "w") do io
        #     write(io, HTTP.get(url).body)
        # end
        # genbank_data = GenomicAnnotations.readgbk(outfile)
        # rm(outfile)
        # return genbank_data
        return GenomicAnnotations.GenBank.Reader(IOBuffer(HTTP.get(url).body))
    elseif !isempty(ftp)
        return GenomicAnnotations.GenBank.Reader(CodecZlib.GzipDecompressorStream(IOBuffer(HTTP.get(ftp).body)))
    else
        @error "invalid call"
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function fasta_and_gff_to_genbank(;fasta, gff, genbank)
    add_bioconda_env("emboss")
    # https://www.insdc.org/submitting-standards/feature-table/
    genbank_directory = dirname(genbank)
    genbank_basename = basename(genbank)
    genbank_prefix = replace(genbank, r"\.(genbank|gb|gbk|gbff)$" => "")
    # genbank_extension = "genbank"
    # https://bioinformatics.stackexchange.com/a/11140
    # https://www.biostars.org/p/72220/#72272
    # seqret -sequence aj242600.fasta -feature -fformat gff -fopenfile aj242600.gff -osformat genbank -auto
#     -osname
    # seqret -sequence {genome file} -feature -fformat gff -fopenfile {gff file} -osformat genbank -osname_outseq {output prefix} -ofdirectory_outseq gbk_file -auto
    run(`$(Mycelia.CONDA_RUNNER) run -n emboss --live-stream seqret -sequence $(fasta) -feature -fformat gff -fopenfile $(gff) -osformat genbank -osname_outseq $(genbank_prefix) -ofdirectory_outseq gbk_file -auto`)
    # return genbank
end

# function gff_to_genbank(gff, genbank)
#     # https://www.insdc.org/submitting-standards/feature-table/
#     genbank_directory = dirname(genbank)
#     genbank_basename = basename(genbank)
#     genbank_prefix = replace(genbank, r"\.(genbank|gb|gbk|gbff)$" => "")
#     # genbank_extension = "genbank"
#     # https://bioinformatics.stackexchange.com/a/11140
#     # https://www.biostars.org/p/72220/#72272
#     # seqret -sequence aj242600.fasta -feature -fformat gff -fopenfile aj242600.gff -osformat genbank -auto
# #     -osname
#     # seqret -sequence {genome file} -feature -fformat gff -fopenfile {gff file} -osformat genbank -osname_outseq {output prefix} -ofdirectory_outseq gbk_file -auto
#     run(`seqret -sequence $(fasta) -feature -fformat gff -fopenfile $(gff) -osformat genbank -osname_outseq $(genbank_prefix) -ofdirectory_outseq gbk_file -auto`)
#     # return genbank
# end

# function genbank_to_fasta_and_gff(genbank_file)

# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function open_genbank(genbank_file)
    return GenomicAnnotations.readgbk(genbank_file)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function parse_virsorter_score_tsv(virsorter_score_tsv)
    data, header = uCSV.read(virsorter_score_tsv, delim='\t', header=1)
    if length(data) == 0
        data = [[] for i in 1:length(header)]
    end
    return DataFrames.DataFrame(data, header)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function parse_mmseqs_tophit_aln(tophit_aln)
    data, header = uCSV.read(tophit_aln, delim='\t')
    # (1,2) identifiers for query and target sequences/profiles, (3) sequence identity, (4) alignment length, (5) number of mismatches, (6) number of gap openings, (7-8, 9-10) domain start and end-position in query and in target, (11) E-value, and (12) bit score.
    # query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits
    header = [
        "query",
        "target",
        "percent identity",
        "alignment length",
        "number of mismatches",
        "number of gaps",
        "query start",
        "query end",
        "target start",
        "target end",
        "evalue",
        "bit score"
    ]
    DataFrames.DataFrame(data, header)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function parse_mmseqs_easy_taxonomy_tophit_report(tophit_report)
    data, header = uCSV.read(tophit_report, delim='\t')
    # tophit_report
    # (1) Target identifier 
    # (2) Number of sequences aligning to target
    # (3) Unique coverage of target uniqueAlignedResidues / targetLength
    # (4) Target coverage alignedResidues / targetLength
    # (5) Average sequence identity
    # (6) Taxonomical information identifier, species, lineage
    header = [
        "target_id",
        "number of sequences aligning to target",
        "unique coverage of target (uniqueAlignedResidues / targetLength)",
        "Target coverage (alignedResidues / targetLength)",
        "Average sequence identity",
        "taxon_id",
        "taxon_rank",
        "taxon_name"
    ]
    DataFrames.DataFrame(data, header)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function parse_mmseqs_easy_taxonomy_lca_tsv(lca_tsv)
    data, header = uCSV.read(lca_tsv, delim='\t')
    # contig
    # (1) a single taxonomy numeric identifier
    # (2) a taxonomic rank column
    # (3) taxonomic name column
    # fragments retained
    # fragments taxonomically assigned
    # fragments in agreement with the contig label (i.e. same taxid or have it as an ancestor)
    # the support received -log(E-value)
    header = [
        "contig_id",
        "taxon_id",
        "taxon_rank",
        "taxon_name",
        "fragments_retained",
        "fragments_taxonomically_assigned",
        "fragments_in_agreement_with_assignment",
        "support -log(E-value)"
    ]
    return DataFrames.DataFrame(data, header)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function read_mmseqs_easy_search(mmseqs_file; top_hit_only=false)
    mmseqs_results = DataFrames.DataFrame(uCSV.read(mmseqs_file, header=1, delim='\t')...)
    if top_hit_only
        gdf = DataFrames.groupby(mmseqs_results, "query")
        for g in gdf
            @assert issorted(g[!, "evalue"])
        end
        top_hits = DataFrames.combine(gdf, first)
        mmseqs_results = top_hits
    end
    return mmseqs_results
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function update_gff_with_mmseqs(gff_file, mmseqs_file)
    top_hits = read_mmseqs_easy_search(mmseqs_file, top_hit_only=true)

    id_to_product = Dict{String, String}()
    for row in DataFrames.eachrow(top_hits)
        id = Dict(a => b for (a, b) in split.(split(last(split(row["qheader"], " # ")), ';'), '='))["ID"]
        product = replace(row["theader"], " " => "__")
        id_to_product[id] = product
    end

    gff_table = Mycelia.read_gff(gff_file)
    for (i, row) in enumerate(DataFrames.eachrow(gff_table))
        id = Dict(a => b for (a,b) in split.(filter(!isempty, split(row["attributes"], ';')), '='))["ID"]
        product = get(id_to_product, id, "")
        gff_table[i, "attributes"] = "label=\"$(product)\";product=\"$(product)\";" * row["attributes"]
    end
    return gff_table
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function open_gff(path::String)
    if isfile(path)
        io = open(path)
    elseif occursin(r"^ftp", path) || occursin(r"^http", path)
        path = replace(path, r"^ftp:" => "http:")
        io = IOBuffer(HTTP.get(path).body)
    else
        error("unable to locate file $path")
    end
    path_base = basename(path)
    if occursin(r"\.gz$", path_base)
        io = CodecZlib.GzipDecompressorStream(io)
    end
    return io
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function read_gff(gff::AbstractString)
    return read_gff(open_gff(gff))
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function read_gff(gff_io)
    data, header = uCSV.read(gff_io, delim='\t', header=0, comment='#')
    header = [
        "#seqid",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes"
    ]
    types=[
        String,
        String,
        String,
        Int,
        Int,
        Int,
        String,
        Int,
        String
    ]
    DataFrames.DataFrame(data, header)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function write_gff(;gff, outfile)
    uCSV.write(outfile, gff, delim='\t')
    return outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function read_kraken_report(kraken_report)
    kraken_report_header = [
        "percentage_of_fragments_at_or_below_taxon",
        "number_of_fragments_at_or_below_taxon",
        "number_of_fragments_assigned_directly_to_taxon",
        "rank",
        "ncbi_taxonid",
        "scientific_name"
    ]

    data, header = uCSV.read(kraken_report, delim='\t')
    kraken_report_table = DataFrames.DataFrame(data, kraken_report_header)
    kraken_report_table[!, "scientific_name"] = string.(strip.(kraken_report_table[!, "scientific_name"]))
    return kraken_report_table
end

# function diamond_line_to_named_tuple(diamond_line)
#     sline = split(line)
#     values_named_tuple = (
#         qseqid = sline[1],
#         sseqid = sline[2],
#         pident = parse(Float64, sline[3]),
#         length = parse(Int, sline[4]),
#         mismatch = parse(Int, sline[5]),
#         gapopen = parse(Int, sline[6]),
#         qlen = parse(Int, sline[7]),
#         qstart = parse(Int, sline[8]),
#         qend = parse(Int, sline[9]),
#         slen = parse(Int, sline[10]),
#         sstart = parse(Int, sline[11]),
#         send = parse(Int, sline[12]),
#         evalue = parse(Float64, sline[13]),
#         bitscore = parse(Float64, sline[14])
#         )
#     return values_named_tuple
# end

# function read_diamond_alignments_file(diamond_file)
#     column_names_to_types = [
#         "qseqid" => String,
#         "sseqid" => String,
#         "pident" => Float64,
#         "length" => Int,
#         "mismatch" => Int,
#         "gapopen" => Int,
#         "qlen" => Int,
#         "qstart" => Int,
#         "qend" => Int,
#         "slen" => Int,
#         "sstart" => Int,
#         "send" => Int,
#         "evalue" => Float64,
#         "bitscore" => Float64,
#     ]
#     types = Dict(i => t for (i, t) in enumerate(last.(column_names_to_types)))
    
#     data, header = uCSV.read(diamond_file, header=0, delim='\t', types = types)
#     header = first.(column_names_to_types)    
    
#     # data, header = uCSV.read(diamond_file, header=1, delim='\t', types = types)
#     # @assert header == first.(column_names_to_types)
    
#     table = DataFrames.DataFrame(data, header)
#     return table
# end

# function add_header_to_diamond_file(infile, outfile=replace(infile, ".tsv" => ".with-header.tsv"))
#     column_names = [
#         "qseqid",
#         "sseqid",
#         "pident",
#         "length",
#         "mismatch",
#         "gapopen",
#         "qlen",
#         "qstart",
#         "qend",
#         "slen",
#         "sstart",
#         "send",
#         "evalue",
#         "bitscore"
#     ]
#     # dangerous but fast
#     # try
#     #     inserted_text = join(columns_names, '\t') * '\n'
#     #     sed_cmd = "1s/^/$(inserted_text)/"
#     #     full_cmd = `sed -i $sed_cmd $infile`
#     # catch
#     open(outfile, "w") do io
#         println(io, join(column_names, "\t"))
#         for line in eachline(infile)
#             println(io, line)
#         end
#     end
#     return outfile
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Parse the contig coverage information from qualimap bamqc text report, which looks like the following:

```
# this is spades
>>>>>>> Coverage per contig

	NODE_1_length_107478_cov_9.051896	107478	21606903	201.0355886786133	60.39424208607496
	NODE_2_length_5444_cov_1.351945	5444	153263	28.152645113886848	5.954250612823136
	NODE_3_length_1062_cov_0.154390	1062	4294	4.043314500941619	1.6655384692688975
	NODE_4_length_776_cov_0.191489	776	3210	4.13659793814433	2.252009588980858

# below is megahit
>>>>>>> Coverage per contig

	k79_175	235	3862	16.43404255319149	8.437436249612457
	k79_89	303	3803	12.551155115511552	5.709975376279777
	k79_262	394	6671	16.931472081218274	7.579217802849293
	k79_90	379	1539	4.060686015831134	1.2929729111266581
	k79_91	211	3749	17.767772511848342	11.899185693011933
	k79_0	2042	90867	44.49902056807052	18.356525483516613
```

To make this more robust, consider reading in the names of the contigs from the assembled fasta
"""
function parse_qualimap_contig_coverage(qualimap_report_txt)
    coverage_line_regex = r"\t.*?\t\d+\t\d+\t[\d\.]+\t[\d\.]+$"
    lines = filter(x -> occursin(coverage_line_regex, x), readlines("$(qualimap_report_txt)"))
    io = IOBuffer(join(map(x -> join(split(x, '\t')[2:end], '\t'), lines), '\n'))
    header = ["Contig", "Length", "Mapped bases", "Mean coverage", "Standard Deviation"]
    types = [String, Int, Int, Float64, Float64]
    data, _ = uCSV.read(io, delim='\t', types=types)
    qualimap_results = DataFrames.DataFrame(data, header)
    qualimap_results[!, "% Mapped bases"] = qualimap_results[!, "Mapped bases"] ./ sum(qualimap_results[!, "Mapped bases"]) .* 100
    return qualimap_results
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Imports results of fastani
"""
function read_fastani(path::String)
    data, header = uCSV.read(path, delim='\t', typedetectrows=100)
    header = [
        "query",
        "reference",
        "%_identity",
        "fragments_mapped",
        "total_query_fragments"
    ]
    ani_table = DataFrames.DataFrame(data, header)
    ani_table[!, "query_identifier"] = replace.(basename.(ani_table[!, "query"]), r"\.(fasta|fna|fa)$" => "")
    ani_table[!, "reference_identifier"] = replace.(basename.(ani_table[!, "reference"]), r"\.(fasta|fna|fa)$" => "")
    columns = [
        "query",
        "query_identifier",
        "reference",
        "reference_identifier",
        "%_identity",
        "fragments_mapped",
        "total_query_fragments"
    ]
    ani_table = ani_table[!, columns]    
    return ani_table
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Imports results of Diamond (or blast) in outfmt 6 as a DataFrame
# """
# function read_diamond(path::String)
#   diamond_colnames = [ "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore" ]
  
#   diamond_coltypes = Dict(
#      1 => String, 
#      2 => String, 
#      3 => Float64,
#      4 => Int,
#      5 => Int,
#      6 => Int,
#      7 => Int,
#      8 => Int,
#      9 => Int,
#      10 => Int,
#      11 => String,
#      12 => String
#   )
#     return DataFrames.DataFrame(uCSV.read(open(path), delim = '\t', header = diamond_colnames, types = diamond_coltypes)...)
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Expects output type 7 from BLAST, default output type 6 doesn't have the header comments and won't auto-parse
"""
function parse_blast_report(blast_report)
    # example header line 
    # "# Fields: query id, subject id, subject acc., subject acc.ver, subject title, query length, subject length, q. start, q. end, s. start, s. end, evalue, bit score, score, alignment length, % identity, identical, mismatches, subject tax id"
    header_lines = collect(Iterators.filter(x -> occursin(r"# Fields:", x), eachline(blast_report)))
    if isempty(header_lines)
        @info "not hits found, returning empty table"
        return DataFrames.DataFrame()
    end
    header_line = first(header_lines)
    header = split(last(split(header_line, ": ")), ", ")
    blast_col_types = Dict(
        "query id" => String,
        "query title" => String,
        "subject id" => String,
        "subject gi" => String,
        "subject acc." => String,
        "subject acc.ver" => String,
        "subject title" => String,
        "query length" => Int,
        "subject length" => Int,
        "q. start" => Int,
        "q. end" => Int,
        "s. start" => Int,
        "s. end" => Int,
        "evalue" => Float64,
        "bit score" => Float64,
        "score" => Float64,
        "alignment length" => Int,
        "% identity" => Float64,
        "identical" => Int,
        "mismatches" => Int,
        "subject tax id" => Int,
        "subject sci name" => String,
        "subject com names" => String,
        "subject blast name" => String,
        "subject super kingdom" => String,
        "subject tax ids" => String,
        "subject sci names" => String,
        "subject com names" => String,
        "subject blast names" => String,
        "subject super kingdoms" => String,
        "subject title" => String,
        "subject titles" => String
    )
    data, _ = uCSV.read(
        blast_report,
        delim='\t',
        comment='#',
        skipmalformed=true,
        allowmissing=true,
        encodings=Dict("N/A" => missing),
        types=[blast_col_types[h] for h in header])
    return DataFrames.DataFrame(data, header, makeunique=true)
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# Parse a GFA file into a genome graph - need to finish implementation and assert contig normalization (i.e. is canonical) before using with my code
# """
# function parse_gfa(gfa)
    
#     # collect(GraphicalFragmentAssembly.Reader(open(primary_contig_gfa)))
    
#     gfa_record_types = Dict(
#         '#' => "Comment",
#         'H' => "Header",
#         'S' => "Segment",
#         'L' => "Link",
#         'J' => "Jump",
#         'C' => "Containment",
#         'P' => "Path",
#         'W' => "Walk"
#     )

#     gfa_graph = MetaGraphs.MetaDiGraph()
#     MetaGraphs.set_prop!(gfa_graph, :paths, Dict{String, Any}())
#     for line in eachline(gfa)
#         record_type = gfa_record_types[line[1]]
#         if record_type == "Header"
#             # metadata
#             sline = split(line)
#             # add me later
#         elseif record_type == "Comment"
#             # metadata
#             # add me later
#         elseif record_type == "Segment"
#             # node
#             record_type, record_name, sequence = split(line, '\t')
#             Graphs.add_vertex!(gfa_graph)
#             node_index = Graphs.nv(gfa_graph)
#             MetaGraphs.set_prop!(gfa_graph, node_index, :identifier, record_name)
#             MetaGraphs.set_indexing_prop!(gfa_graph, :identifier)
#             MetaGraphs.set_prop!(gfa_graph, node_index, :sequence, sequence)
#         elseif record_type == "Link"
#             record_type, source_identifier, source_orientation, destination_identifier, destination_orientation, overlap_CIGAR = split(line, '\t')
#             source_index = gfa_graph[source_identifier, :identifier]
#             destination_index = gfa_graph[destination_identifier, :identifier]
#             edge = Graphs.Edge(source_index, destination_index)
#             Graphs.add_edge!(gfa_graph, edge)
#             MetaGraphs.set_prop!(gfa_graph, edge, :source_identifier, source_identifier)
#             MetaGraphs.set_prop!(gfa_graph, edge, :source_orientation, source_orientation)
#             MetaGraphs.set_prop!(gfa_graph, edge, :destination_identifier, destination_identifier)
#             MetaGraphs.set_prop!(gfa_graph, edge, :destination_orientation, destination_orientation)
#             MetaGraphs.set_prop!(gfa_graph, edge, :overlap_CIGAR, overlap_CIGAR)
#         elseif record_type == "Path"
#             record_type, path_identifier, segments, overlaps = split(line, '\t')
#             gfa_graph.gprops[:paths][path_identifier] = Dict("segments" => segments, "overlaps" => overlaps)
#         else
#             @warn "GFA line type $(record_type) not currently handled by the import - please add"
#         end
#     end
#     return gfa_graph
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function graph_to_gfa(;graph, outfile)
    # kmer_vertices = collect(MetaGraphs.filter_vertices(graph, :TYPE, Kmers.kmertype(Kmers.DNAKmer{kmer_size})))
    # add fastq here too???
    # record_vertices = collect(MetaGraphs.filter_vertices(graph, :TYPE, FASTX.FASTA.Record))
    # edgemer_edges = collect(MetaGraphs.filter_edges(graph, :TYPE, Kmers.kmertype(Kmers.DNAKmer{kmer_size+1})))
    open(outfile, "w") do io
        println(io, "H\tVN:Z:1.0")
        for vertex in Graphs.vertices(graph)
            if haskey(graph.vprops[vertex], :kmer)
                sequence = graph.vprops[vertex][:kmer]
            else
                sequence = graph.vprops[vertex][:sequence]
            end
            depth = graph.vprops[vertex][:count]
            # depth = length(graph.vprops[vertex][:evidence])
            # total_count = 0
            # vertex_outneighbors = Graphs.outneighbors(graph, vertex)
            # for connected_record in vertex_outneighbors
            #     edge = Graphs.Edge(vertex, connected_record)
            #     total_count += MetaGraphs.get_prop(graph, edge, :count)
            # end

            fields = ["S", "$vertex", sequence, "DP:f:$(depth)"]
            line = join(fields, '\t')
            println(io, line)
        end
        overlap = MetaGraphs.get_prop(graph, :k) - 1
        for edge in collect(Graphs.edges(graph))
            # @show edge
            # @show MetaGraphs.props(graph, edge)
            unique_orientations = unique(observation.orientation for observation in MetaGraphs.get_prop(graph, edge, :observations))
            # @show unique_orientations
            source_kmer = edge.src
            destination_kmer = edge.dst
            for unique_orientation in unique_orientations
                # source_kmer, dest_kmer = Mycelia.edgemer_to_vertex_kmers(MetaGraphs.get_prop(graph, edgemer_edge, :sequence))
                link = ["L",
                            edge.src,
                            first(unique_orientation) ? '+' : '-',
                            edge.dst,
                            last(unique_orientation) ? '+' : '-',
                            "$(overlap)M"]
                line = join(link, '\t')
                println(io, line)
            end
                
            # source_kmer, dest_kmer = Mycelia.edgemer_to_vertex_kmers(MetaGraphs.get_prop(graph, edgemer_edge, :sequence))
            # link = ["L",
            #             edgemer_edge.src,
            #             BioSequences.iscanonical(source_kmer) ? '+' : '-',
            #             edgemer_edge.dst,
            #             BioSequences.iscanonical(dest_kmer) ? '+' : '-',
            #             "$(overlap)M"]
            # line = join(link, '\t')
            # println(io, line)
        end
    end
    return outfile
end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function graph_to_gfa(graph, kmer_size, outfile="$(kmer_size).gfa")
#     kmer_vertices = collect(MetaGraphs.filter_vertices(graph, :TYPE, Kmers.kmertype(Kmers.DNAKmer{kmer_size})))
#     # add fastq here too???
#     record_vertices = collect(MetaGraphs.filter_vertices(graph, :TYPE, FASTX.FASTA.Record))
#     edgemer_edges = collect(MetaGraphs.filter_edges(graph, :TYPE, Kmers.kmertype(Kmers.DNAKmer{kmer_size+1})))
#     open(outfile, "w") do io
#         println(io, "H\tVN:Z:1.0")
#         for vertex in kmer_vertices
#             if haskey(graph.vprops[vertex], :kmer)
#                 sequence = graph.vprops[vertex][:kmer]
#             else
#                 sequence = graph.vprops[vertex][:sequence]
#             end
#             # depth = graph.vprops[vertex][:count]
#             # depth = length(graph.vprops[vertex][:evidence])
#             total_count = 0
#             vertex_outneighbors = Graphs.outneighbors(graph, vertex)
#             for connected_record in intersect(vertex_outneighbors, record_vertices)
#                 edge = Graphs.Edge(vertex, connected_record)
#                 total_count += MetaGraphs.get_prop(graph, edge, :count)
#             end

#             fields = ["S", "$vertex", sequence, "DP:f:$(total_count)"]
#             line = join(fields, '\t')
#             println(io, line)
#         end
#         overlap = kmer_size - 1
#         for edgemer_edge in edgemer_edges
#             source_kmer, dest_kmer = Mycelia.edgemer_to_vertex_kmers(MetaGraphs.get_prop(graph, edgemer_edge, :sequence))
#             link = ["L",
#                         edgemer_edge.src,
#                         BioSequences.iscanonical(source_kmer) ? '+' : '-',
#                         edgemer_edge.dst,
#                         BioSequences.iscanonical(dest_kmer) ? '+' : '-',
#                         "$(overlap)M"]
#             line = join(link, '\t')
#             println(io, line)
#         end
#     end
#     return
# end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function open_fastx(path::AbstractString)
    if isfile(path)
        io = open(path)
    elseif occursin(r"^ftp", path) || occursin(r"^http", path)
        path = replace(path, r"^ftp:" => "http:")
        io = IOBuffer(HTTP.get(path).body)
    else
        error("unable to locate file $path")
    end
    path_base = basename(path)
    if occursin(r"\.gz$", path_base)
        io = CodecZlib.GzipDecompressorStream(io)
        path_base = replace(path_base, ".gz" => "")
    end
    if occursin(r"\.(fasta|fna|faa|fa)$", path_base)
        fastx_io = FASTX.FASTA.Reader(io)
    elseif occursin(r"\.(fastq|fq)$", path_base)
        fastx_io = FASTX.FASTQ.Reader(io)
    else
        @show path_base
        error()
    end
    return fastx_io
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
Writes FASTA records to a file, optionally gzipped.

# Arguments
- `outfile::AbstractString`: Path to the output FASTA file.  Will append ".gz" if `gzip` is true.
- `records::Vector{FASTX.FASTA.Record}`: A vector of FASTA records.
- `gzip::Bool=false`: Whether to compress the output with gzip.

# Returns
- `outfile::String`: The path to the output FASTA file (including ".gz" if applicable).
"""
function write_fasta(;outfile::AbstractString=tempname()*".fna", records::Vector{FASTX.FASTA.Record}, gzip::Bool=false)
    # Determine if gzip compression should be used based on both the filename and the gzip argument
    gzip = occursin(r"\.gz$", outfile) || gzip

    # Append ".gz" to the filename if gzip is true and the extension isn't already present
    outfile = gzip && !occursin(r"\.gz$", outfile) ? outfile * ".gz" : outfile  # More concise way to handle filename modification

    # Use open with do block for automatic resource management (closing the file)
    open(outfile, "w") do io
        if gzip
            io = CodecZlib.GzipCompressorStream(io)  # Wrap the io stream for gzip compression
        end

        FASTX.FASTA.Writer(io) do fastx_io  # Use do block for automatic closing of the FASTA writer
            for record in records
                write(fastx_io, record)
            end
        end # fastx_io automatically closed here
    end # io automatically closed here

    return outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)

Plots a histogram of kmer counts against # of kmers with those counts

Returns the plot object for adding additional layers and saving
"""
function plot_kmer_frequency_spectra(counts; log_scale = log2, kwargs...)
    kmer_counts_hist = StatsBase.countmap(c for c in counts)
    xs = collect(keys(kmer_counts_hist))
    ys = collect(values(kmer_counts_hist))
    if isa(log_scale, Function)
        xs = log_scale.(xs)
        ys = log_scale.(ys)
    end
    
    p = StatsPlots.plot(
        xs,
        ys,
        xlims = (0, maximum(xs) + max(1, ceil(0.1 * maximum(xs)))),
        ylims = (0, maximum(ys) + max(1, ceil(0.1 * maximum(ys)))),
        seriestype = :scatter,
        legend = false,
        xlabel = isa(log_scale, Function) ? "$(log_scale)(observed frequency)" : "observed frequency",
        ylabel = isa(log_scale, Function) ? "$(log_scale)(# of kmers)" : "observed frequency",
        ;kwargs...
    )
    return p
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function convert(args)
    @show args
    in_type = missing
    out_type = missing
    if occursin(r"\.jld2$", args["in"])
        in_type = :jld2
    elseif occursin(r"\.gfa$", args["in"])
        in_type = :gfa
    elseif occursin(r"\.neo4j$", args["in"])
        in_type = :neo4j
    end

    if occursin(r"\.jld2$", args["out"])
        out_type = :jld2
    elseif occursin(r"\.gfa$", args["out"])
        out_type = :gfa
    elseif occursin(r"\.neo4j$", args["out"])
        out_type = :neo4j
    end
    
    if ismissing(in_type) || ismissing(out_type)
        error("unable to determine in and out types")
    end
    
    if (in_type == :jld2) && (out_type == :jld2)
        # done
    elseif (in_type == :jld2) && (out_type != :jld2)
        # convert out of jld2
        if out_type == :gfa
            loaded = FileIO.load(args["in"])
            @assert haskey(loaded, "graph")
            graph = loaded["graph"]
            graph_to_gfa(graph, args["out"])
        end
    elseif (in_type != :jld2) && (out_type == :jld2)
        # convert into jld2
    else
        # need to convert from input to jld2 as an intermediate first
    end
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function plot_graph(graph)
    
#     kmer_counts = MetaGraphs.get_prop(graph, :kmer_counts)
    kmers = [MetaGraphs.get_prop(graph, v, :kmer) for v in Graphs.vertices(graph)]
    counts = [MetaGraphs.get_prop(graph, v, :count) for v in Graphs.vertices(graph)]
    scale = 150
    
    n = Graphs.nv(graph)
    p = GraphRecipes.graphplot(
        graph,
#         markersize = 1/log2(n),
        markersize = 1/2^2,
        size = (2 * scale * log(n), scale * log(n)),
        node_weights = counts,
        names = kmers)
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function save_graph(graph::Graphs.AbstractGraph, outfile::String)
    if !occursin(r"\.jld2$", outfile)
        outfile *= ".jld2"
    end
    FileIO.save(outfile, Dict("graph" => graph))
    return outfile
end

"""
$(DocStringExtensions.TYPEDSIGNATURES)
"""
function load_graph(file::String)
    return FileIO.load(file)["graph"]
end

# OLD FOR SIMPLE KMER GRAPHS
# # """
# # $(DocStringExtensions.TYPEDSIGNATURES)

# # Description

# # ```jldoctest
# # julia> 1 + 1
# # 2
# # ```
# # """
# function graph_to_gfa(graph, outfile)
#     open(outfile, "w") do io
#         println(io, "H\tVN:Z:1.0")
#         for vertex in Graphs.vertices(graph)
#             if haskey(graph.vprops[vertex], :kmer)
#                 sequence = graph.vprops[vertex][:kmer]
#             else
#                 sequence = graph.vprops[vertex][:sequence]
#             end
# #             depth = graph.vprops[vertex][:weight]
#             depth = graph.vprops[vertex][:count]
#             fields = ["S", "$vertex", sequence, "DP:f:$(depth)"]
#             line = join(fields, '\t')
#             println(io, line)
#         end
#         for edge in Graphs.edges(graph)
#             overlap = graph.gprops[:k] - 1
#             for o in graph.eprops[edge][:orientations]
# #                 if !(!o.source_orientation && !o.destination_orientation)
#                 link = ["L",
#                             edge.src,
#                             o.source_orientation ? '+' : '-',
#                             edge.dst,
#                             o.destination_orientation ? '+' : '-',
#                             "$(overlap)M"]
#                 line = join(link, '\t')
#                 println(io, line)
# #                 end
#             end
#         end
#     end
#     return outfile
# end

# """
# $(DocStringExtensions.TYPEDSIGNATURES)

# A short description of the function

# ```jldoctest
# julia> 1 + 1
# 2
# ```
# """
# function my_plot(graph::KmerGraph)
#     graph_hash = hash(sort(graph.graph.fadjlist), hash(graph.graph.ne))
#     filename = "/assets/images/$(graph_hash).svg"
#     p = plot_graph(graph)
#     Plots.savefig(p, dirname(pwd()) * filename)
#     display(p)
#     display("text/markdown", "![]($filename)")
# end
