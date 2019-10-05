
import Printf: @sprintf
import OrderedCollections: OrderedDict
import Base: show, convert, print

using Parameters

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


# If we just provide the basics as positional (to do only the 6 positional)
DemRsc(width, file_length, x_first, y_first, 
       x_step, y_step) = DemRsc(width=width, file_length=file_length, 
                                x_first=x_first, y_first=y_first, 
                                x_step=x_step, y_step=y_step)

# For use in unpacking
_symdict(d::AbstractDict{String, Any}) = Dict(Symbol(k) => v for (k, v) in d)
DemRsc(d::AbstractDict{String, Any}) = DemRsc(; _symdict(d)...)


# Note: fn is a Symbol below
function stringdict(x::DemRsc)  #::OrderedDict{String, Any}
    d = OrderedDict{String, Any}()
    for fn in fieldnames(typeof(x))
        d[string(fn)] = getfield(x, fn) 
    end
    return d
end

# For converting to/from Dicts for python/HDF5 saving
Base.convert(::Type{DemRsc}, x::OrderedDict{String, Any}) = DemRsc(x)
Base.convert(::Type{OrderedDict{String, Any}}, x::DemRsc) = stringdict(x)
# Base.convert(::Type{Dict{String, Any}}, x::DemRsc) = stringdict(x)

# TODO: iterate
# julia> iterate(dd)
# ("width" => 2, 2)
# Base.iterate(demrsc::DemRsc)

function format_dem_rsc(demrsc::DemRsc)
    outstring = ""

    # rsc_dict = Dict(String(k) => v for (k, v) in rsc_dict)
    # for field, value in rsc_dict.items():
    for field in fieldnames(typeof(demrsc))
        value = getproperty(demrsc, field)

        # Files seemed to be left justified with 14 spaces? Not sure why 14
        # Apparently it was an old fortran format, where they use "read(15)"
        if field in (:x_step, :y_step)
            # give step floats proper sig figs to not output scientific notation
            outstring *= @sprintf("%-14s%0.12f\n", uppercase(string(field)), value)
            # outstring *= @sprint"{field:<14s}{val:0.12f}\n".format(field=field.upper(), val=value)
        else
            outstring *= @sprintf("%-14s%s\n", uppercase(string(field)), value)
            # outstring *= "{field:<14s}{val}\n".format(field=field.upper(), val=value)
        end
    end

    return outstring
end
Base.print(io::IO, x::DemRsc) = print(format_dem_rsc(x))
