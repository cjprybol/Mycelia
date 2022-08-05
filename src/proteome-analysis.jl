function merge_fasta_files(list_of_fasta_files, joint_fasta_file)
    open(joint_fasta_file, "w") do io
        for f in list_of_fasta_files
            write(io, read(f))
        end
    end
    return joint_fasta_file
end

function make_diamond_db(fasta_file, db_file=fasta_file)
    @time run(`diamond makedb --in $(fasta_file) -d $(db_file)`)
end