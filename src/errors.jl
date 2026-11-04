# Exceptions of Aeolis. Every exception the package defines lives in this file, with its
# `showerror` method, and subtypes the root `AeolisError`.
#
# - Carry the offending names (parameter, value, variable) in the fields and the message.
#   Library code never calls a bare `error("...")`.
# - Invalid arguments and keywords raise `ArgumentError`, not a type from this file.

"""
    AeolisError

Abstract supertype of every exception that Aeolis defines. Invalid arguments and keywords
raise Base's `ArgumentError` instead, which is outside this root.
"""
abstract type AeolisError <: Exception end

"""
    ConfigurationError(message)

Raised when the parameters do not describe a runnable model, for example when the grid or bed
is missing or an input series has too few columns.
"""
struct ConfigurationError <: AeolisError
    message::String
end
Base.showerror(io::IO, e::ConfigurationError) = print(io, "ConfigurationError: ", e.message)

"""
    UnknownOptionError(parameter, value, allowed)

Raised when a method or boundary parameter names an option that does not exist; `allowed`
lists the options that do.
"""
struct UnknownOptionError <: AeolisError
    parameter::Symbol
    value::Symbol
    allowed::Vector{Symbol}
end
function Base.showerror(io::IO, e::UnknownOptionError)
    return print(io, "UnknownOptionError: parameter :", e.parameter, " has no option :",
                 e.value, "; options are ",
                 join((":" * String(a) for a in e.allowed), ", "))
end

"""
    UnsupportedOptionError(parameter, value, reason)

Raised when an option exists in AeoLiS but has no implementation there either (for example a
cross-shore circular boundary or the dry-layer threshold), or is not available in this port.
"""
struct UnsupportedOptionError <: AeolisError
    parameter::Symbol
    value::Any
    reason::String
end
function Base.showerror(io::IO, e::UnsupportedOptionError)
    return print(io, "UnsupportedOptionError: :", e.parameter, " = ", repr(e.value), ": ",
                 e.reason)
end

function _check_option(name::Symbol, value::Symbol, allowed)
    value in allowed || throw(UnknownOptionError(name, value, collect(Symbol, allowed)))
    return value
end
