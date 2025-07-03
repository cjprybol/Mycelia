"""
    parallel_pyrodigal(normalized_fastas::Vector{String})

Runs Mycelia.run_pyrodigal on a list of FASTA files in parallel using Threads.

Args:
    normalized_fastas: A vector of strings, where each string is a path to a FASTA file.

Returns:
    A tuple containing two elements:
    1. successes (Vector{Tuple{String, Any}}): A vector of tuples, where each tuple contains the
       filename and the result returned by a successful Mycelia.run_pyrodigal call.
    2. failures (Vector{Tuple{String, String}}): A vector of tuples, where each tuple contains the
       filename and the error message string for a failed Mycelia.run_pyrodigal call.
"""
function parallel_pyrodigal(normalized_fastas::Vector{String})
    num_files = Base.length(normalized_fastas)
    Base.println("Processing $(num_files) FASTA files using $(Threads.nthreads()) threads...")

    # Create a Progress object for manual updates
    p = ProgressMeter.Progress(num_files, 1, "Running Pyrodigal: ", 50)

    # Use Channels to collect results and failures thread-safely
    # Channel{Tuple{Filename, ResultType}} - adjust ResultType if known
    successes = Base.Channel{Tuple{String, Any}}(num_files)
    failures = Base.Channel{Tuple{String, String}}(num_files)

    # Use Threads.@threads for parallel execution
    Threads.@threads for fasta_file in normalized_fastas
        result = nothing # Initialize result variable in the loop's scope
        try
            # --- Execute the function ---
            # Base.println("Thread $(Threads.threadid()) processing: $(fasta_file)") # Optional: for debugging
            result = Mycelia.run_pyrodigal(fasta_file = fasta_file) # Capture the result

            # --- Store success ---
            Base.put!(successes, (fasta_file, result))

        catch e
            # --- Store failure ---
            err_msg = Base.sprint(Base.showerror, e) # Get the error message as a string
            Base.println(Base.stderr, "ERROR processing $(fasta_file) on thread $(Threads.threadid()): $(err_msg)")
            Base.put!(failures, (fasta_file, err_msg))
        finally
            # --- Always update progress ---
            ProgressMeter.next!(p)
        end
    end

    # Close channels now that all threads are done writing
    Base.close(successes)
    Base.close(failures)

    # Collect results and failures from the channels
    successful_results = Base.collect(successes)
    failed_files = Base.collect(failures)

    # --- Report Summary ---
    Base.println("\n--- Pyrodigal Processing Summary ---")
    num_success = Base.length(successful_results)
    num_failed = Base.length(failed_files)
    Base.println("Successfully processed: $(num_success)")
    Base.println("Failed: $(num_failed)")

    if !Base.isempty(failed_files)
        Base.println("\nFailures:")
        for (file, err) in failed_files
            Base.println("- File: $(file)\n  Error: $(err)")
        end
    end
    Base.println("------------------------------------")

    return successful_results, failed_files # Return both successes and failures
end