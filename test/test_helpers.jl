import Test

const COMMON_ERROR_MESSAGE_FRAGMENTS = (
    "Invalid",
    "invalid",
    "Missing",
    "missing",
    "not found",
    "does not exist",
    "empty",
    "Empty",
    "must",
    "required",
    "Unsupported",
    "unsupported",
    "not supported",
    "only use SingleStrand",
    "Unknown",
    "unknown",
    "failed",
    "Failed",
    "No ",
    "cannot",
    "Cannot",
    "does not contain",
    "will not fit",
    "exceeds requested length",
    "positive",
    "range",
)

function test_throws_message(f, expected_type::Type{<:Exception}, expected_fragments)
    try
        f()
        Test.@test false
    catch err
        Test.@test err isa expected_type
        message = sprint(showerror, err)
        fragments = expected_fragments isa AbstractString ? (expected_fragments,) : expected_fragments
        Test.@test any(fragment -> occursin(fragment, message), fragments)
    end
end
