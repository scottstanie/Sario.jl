using Parameters

DemRsc = Dict{String, Any}

@with_kw mutable struct DemRsc2
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
# E.g. DemRsc2(width=10, file_length=10, x_step=1., y_step=1., x_first=1., y_first=1.)
# TODO
