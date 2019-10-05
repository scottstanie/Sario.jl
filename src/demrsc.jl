using Parameters

import OrderedCollections: OrderedDict
import Base: show, convert

# DemRsc = Dict{String, Any}

@with_kw mutable struct DemRsc
    width::Int
    file_length::Int
    x_first::Float64
    y_first::Float64
    x_step::Float64
    y_step::Float64
    x_unit::String = "degrees"
    y_unit::String = "degrees"
    z_offset::Int = 0
    z_scale::Int = 1
    projection::String = "LL"
end

# E.g. DemRsc(width=10, file_length=10, x_step=1., y_step=1., x_first=1., y_first=1.)


# If we just provide the basics as positional:
DemRsc(width, file_length, x_first, y_first, 
       x_step, y_step) = DemRsc(width=width, file_length=file_length, 
                                x_first=x_first, y_first=y_first, 
                                x_step=x_step, y_step=y_step)

# For use in unpacking
_symdict(d::AbstractDict{String, Any}) = Dict(Symbol(k) => v for (k, v) in d)
DemRsc(d::AbstractDict{String, Any}) = DemRsc(; _symdict(d)...)


# Note: fn is a Symbol below
function stringdict(x::DemRsc)::OrderedDict{String, Any}
    d = OrderedDict{String, Any}()
    for fn in fieldnames(typeof(x))
        d[string(fn)] = getfield(x, fn) 
    end
    return d
end
# Do I need a symbol dict ever?
# typedict(x)::Dict{Symbol, Any} = Dict(fn => getfield(x, fn) for fn in fieldnames(typeof(x)))


# For converting to/from Dicts for python/HDF5 saving
Base.convert(::Type{DemRsc}, x::OrderedDict{String, Any}) = DemRsc(x)
Base.convert(::Type{OrderedDict{String, Any}}, x::DemRsc) = stringdict(x)
Base.convert(::Type{Dict{String, Any}}, x::DemRsc) = stringdict(x)

# Now we can either do the 6 positional, or all as keywords

# function format_dem_rsc(rsc_dict)
#     outstring = ""
#     rsc_dict = Dict(String(k) => v for (k, v) in rsc_dict)
#     # for field, value in rsc_dict.items():
#     for field in RSC_KEYS
#         # Make sure to skip extra keys that might be in the dict
#         if !(field in RSC_KEYS)
#             continue
#         end
# 
#         value = rsc_dict.get(field, DEFAULT_KEYS.get(field))
#         if value is nothing
#             error("$field is necessary for .rsc file: missing from dict")
#         end
# 
#         # Files seemed to be left justified with 14 spaces? Not sure why 14
#         # Apparently it was an old fortran format, where they use "read(15)"
#         if field in ("x_step", "y_step"):
#             # give step floats proper sig figs to not output scientific notation
#             outstring += "{field:<14s}{val:0.12f}\n".format(field=field.upper(), val=value)
#         else:
#             outstring += "{field:<14s}{val}\n".format(field=field.upper(), val=value)
#         end
#     end
# 
#     return outstring
# end
