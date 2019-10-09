"""Module Sario
"""


__precompile__(true)

module Sario

include("./demrsc.jl")

export load, save, DemRsc, load_stack, load_hdf5_stack, load_mask, load_geolist_from_h5,
       load_intlist_from_h5, load_dem_from_h5, get_file_ext, find_files, save_hdf5_stack

# TODO: python transfers:
# const LOAD_IN_PYTHON = vcat(IMAGE_EXTS, [".geojson", ".npy"])
# GeoJSON.geo2dict(GeoJSON.parsefile("../../bbox2.geojson"))["coordinates"]

import Glob
import JSON
using HDF5
using Dates
using Parameters

const DATE_FMT = "yyyymmdd"
const SENTINEL_EXTS = [".geo", ".cc", ".int", ".amp", ".unw", ".unwflat"]
const COMPLEX_EXTS = [".int", ".slc", ".geo", ".cc", ".unw", ".unwflat", ".mlc", ".grd"]
const REAL_EXTS = [".amp", ".cor", ".mlc", ".grd"]
const ELEVATION_EXTS = [".dem", ".hgt"]
# These file types are not simple complex matrices: see load_stacked_img for detail
# .unwflat are same as .unw, but with a linear ramp removed
const STACKED_FILES = [".cc", ".unw", ".unwflat"]
const BOOL_EXTS = [".mask"]

const IMAGE_EXTS = [".png", ".tif", ".tiff", ".jpg"]
const NOT_IMPLEMENTED = vcat(IMAGE_EXTS, [".geojson", ".npy"])

# For HDF5 saving
const DEM_RSC_DSET = "dem_rsc"
const STACK_DSET = "stack"
const GEOLIST_DSET = "geo_dates"
const INTLIST_DSET = "int_dates"

RSC_KEY_TYPES = [
    ("width", Int),
    ("file_length", Int),
    ("x_first", Float64),
    ("y_first", Float64),
    ("x_step", Float64),
    ("y_step", Float64),
    ("x_unit", String),
    ("y_unit", String),
    ("z_offset", Int),
    ("z_scale", Int),
    ("projection", String),
]
RSC_KEYS = [tup[1] for tup in RSC_KEY_TYPES]


find_files(ext, directory=".") = sort(Glob.glob("*"*ext, directory))

"""Extracts the file extension, including the "." (e.g.: .slc)"""
get_file_ext(filename::AbstractString) = splitext(filename)[end]

_is_h5(filename) = Sario.get_file_ext(filename) in (".h5", ".hdf5")

import Base.size
function Base.size(h5file::String, dset::String)
    h5open(h5file) do f
        return size(f[dset])
    end
end

"""Examines file type for real/complex and runs appropriate load

Raises:
    ValueError: if sentinel files loaded without a .rsc file in same path
        to give the file width
"""
function load(filename::AbstractString; 
              rsc_file::Union{AbstractString, Nothing}=nothing,
              looks::Tuple{Int, Int}=(1, 1), do_permute=true, 
              dset_name::AbstractString="", return_amp::Bool=false)
    ext = get_file_ext(filename)

    # For now, just pass through unimplemented extensions to Python
    if ext in NOT_IMPLEMENTED
        error("$ext is not yet implemented")
    elseif ext == ".rsc"
        return _get_rsc_data(filename, filename)
    elseif ext in ELEVATION_EXTS
        return take_looks(load_elevation(filename), looks...)
    elseif _is_h5(filename)
        return take_looks(_load_hdf5(filename, dset_name; do_permute=do_permute), looks...)
    end

    # Sentinel files should have .rsc file: check for dem.rsc, or elevation.rsc
    demrsc = _get_rsc_data(filename, rsc_file)

    if ext in STACKED_FILES
        return take_looks(load_stacked_img(filename, demrsc, do_permute=do_permute, return_amp=return_amp),
                          looks...)
    elseif ext in BOOL_EXTS
        return take_looks(load_bool(filename, demrsc, do_permute=do_permute), looks...)
    else
        return take_looks(load_complex(filename, demrsc, do_permute=do_permute), looks...)
    # having demrsc implies that this is not a UAVSAR file, so is complex
    # TODO: haven"t transferred over UAVSAR functions, so no load_real yet
    end

end

"""Load one element of a file on disk (avoid reading in all of huge file"""
# TODO: Load a chunk of a file now?
function load(filename::AbstractString, row_col::Tuple{Int, Int}; rsc_file::Union{AbstractString, Nothing}=nothing)
    data_type, demrsc, num_rows, num_cols = _file_info(filename, rsc_file)

    row, col = row_col
    _check_bounds(row, col, num_rows, num_cols)

    seek_pos = _get_seek_position(row, col, num_cols, data_type)
    
    open(filename) do f
        seek(f, seek_pos)
        # This read syntax loads single `data_type`
        return read(f, data_type)
    end
end

"""For single element reading in binary files, seek to the right row, col"""
_get_seek_position(row, col, num_cols, data_type) = sizeof(data_type) * ((col - 1) + (num_cols * (row - 1)) )

function _check_bounds(row, col, num_rows, num_cols)
    if row < 1 || col < 1 || row > num_rows || col > num_cols
        throw(DomainError((row, col), " out of bounds for $filename of size ($num_rows, $num_cols)"))
    end
end

function _file_info(filename, rsc_file)
    data_type = _get_data_type(filename)
    demrsc = _get_rsc_data(filename, rsc_file)
    @unpack rows, cols = demrsc
    return data_type, demrsc, rows, cols
end

# NOTE: the subset will only work for sequential data (complex), not the stacked filetypes
"""Load subset of a file on disk using range"""
RangeTuple = Tuple{T, T} where { T <: OrdinalRange }
function load(filename::AbstractString, idxs::RangeTuple; 
              rsc_file::Union{AbstractString, Nothing}=nothing,  do_permute=true)
    data_type, demrsc, num_rows, num_cols = _file_info(filename, rsc_file)

    # Note: reversing these since julia uses column-major
    cols, rows = idxs
    _check_bounds(rows.start, cols.start, num_rows, num_cols)
    
    outrows = rows.stop  - rows.start + 1
    outcols = cols.stop  - cols.start + 1
    out = zeros(data_type, (outrows, outcols))
    buf = zeros(data_type, outcols)

    open(filename) do f
        for (idx, r) in enumerate(rows)
            seek_pos = _get_seek_position(r, cols.start, num_cols, data_type)
            seek(f, seek_pos)
            # This read syntax loads single `data_type`
            # read!(out[r, :], f, data_type)
            read!(f, buf)
            out[:, idx] .= buf  
        end
    end
    return do_permute ? permutedims(out) : out
end

function _load_hdf5(h5file::AbstractString, dset_name::AbstractString=""; do_permute=true)
    h5open(h5file, "r") do f
        if isempty(dset_name)
            if length(names(f)) == 1 && typeof(f[names(f)[1]]) == HDF5Dataset
                dset_name = names(f)[1]
                println("Loading \"$dset_name\" from $h5file")
            else
                error("More than 1 dset exists in $h5file, provide name to `load`")
            end
        end
        order = ndims(f[dset_name]) == 3 ? (2, 1, 3) : (2, 1)  # Handle 3D permuting
        return do_permute ? permutedims(read(f, dset_name), order) : read(f, dset_name)
        k
    end
end


"""Handle getting the .rsc data from either an image filename, or .rsc filename"""
function _get_rsc_data(filename, rsc_file::Union{AbstractString, Nothing})
    ext = get_file_ext(filename)

    demrsc = nothing
    if !isnothing(rsc_file)
        demrsc = load_dem_rsc(rsc_file)
    elseif ext in vcat(SENTINEL_EXTS, ELEVATION_EXTS, BOOL_EXTS)
        rsc_file = find_rsc_file(filename)
        demrsc = load_dem_rsc(rsc_file)
    end

    return demrsc
end

function find_rsc_file(filename=nothing; directory=nothing, verbose=false)
    if !isnothing(filename)
        dir_parts = splitpath(filename)[1:end-1]
        directory = isempty(dir_parts) ? "." : joinpath(dir_parts...)
    end
    # Should be just elevation.dem.rsc (for .geo folder) or dem.rsc (for igrams)
    possible_rscs = find_files("*.rsc", directory)
    if verbose
        println("Searching $directory for rsc files")
        println("Possible rsc files:")
        println(possible_rscs)
    end
    if length(possible_rscs) < 1
        logger.info("No .rsc file found in $directory")
        return nothing
    elseif length(possible_rscs) > 1
        error("$filename has multiple .rsc files in its directory: $possible_rscs")
    end
    return abspath(possible_rscs[1])
end

"""Convert the text file into the DemRsc struct
Starts with a dict to gather all fiields, then unpacks using keywords args"""
function load_dem_rsc(filename, kwargs...)::DemRsc
    output_data = Dict{Symbol, Any}()

    for line in readlines(filename)
        for (field, valtype) in RSC_KEY_TYPES
            if startswith(line, uppercase(field))
                val = valtype <: AbstractString ? split(line)[2] : parse(valtype, split(line)[2]) 
                output_data[Symbol(lowercase(field))] = val
            end
        end
    end
    return DemRsc(;output_data...)
end

function _get_data_type(filename)
    ext = get_file_ext(filename)
    if ext in ELEVATION_EXTS
        return Int16
    elseif ext in STACKED_FILES
        return Float32
    else
        return ComplexF32
    end
end


"""Loads a digital elevation map from either .hgt file or .dem

.hgt is the NASA SRTM files given. Documentation on format here:
https://dds.cr.usgs.gov/srtm/version2_1/Documentation/SRTM_Topo.pdf
Key point: Big-endian 2 byte (16-bit) integers

.dem is format used by Zebker geo-coded and ROI-PAC SAR software
Only difference is data is stored little-endian (like other SAR data)

Note on both formats: gaps in coverage are given by INT_MIN -32768,
so either manually set data(data == np.min(data)) = 0,
or something like 
    data = clamp(data, -10000, Inf)
"""
function load_elevation(filename; do_permute=true)
    ext = get_file_ext(filename)
    data_type = Int16 

    if ext == ".dem"
        rsc_file = find_rsc_file(filename)
        dem_rsc = load(rsc_file)
        rows, cols = (dem_rsc["file_length"], dem_rsc["width"])
        data = Array{data_type, 2}(undef, (cols, rows))

        read!(filename, data)

        # # TODO: Verify that the min real value will be above -1000
        # min_valid = -10000
        # # Set NaN values to 0
        # @. data[data < min_valid] = 0
        return do_permute ? permutedims(data) : data
    else
        # swap_bytes = (ext == ".hgt")
        throw("$ext not implemented")
    end
end


function load_complex(filename::AbstractString, demrsc::DemRsc; do_permute=true)
    return _load_bin_matrix(filename, demrsc, ComplexF32, do_permute)
end

function load_bool(filename::AbstractString, demrsc::DemRsc; do_permute=true)
    return _load_bin_matrix(filename, demrsc, Bool, do_permute)
end

function _load_bin_matrix(filename, demrsc::DemRsc, dtype, do_permute::Bool)
    @unpack rows, cols = demrsc
    # Note: must be saved in c/row-major order, so loading needs a transpose
    out = Array{dtype, 2}(undef, (cols, rows))

    read!(filename, out)
    return do_permute ? permutedims(out) : out
end


"""load_stacked_img is for the weird "ALT_LINE_DATA" formatted images:
https://web.stanford.edu/group/radar/softwareandlinks/sw/snaphu/snaphu_man1.html#FILE%20FORMATS

Format is two stacked matrices:
    [[first], [second]] where the first "cols" number of floats
    are the first matrix, next "cols" are second, etc.
For .unw height files, the first is amplitude, second is phase (unwrapped)
For .cc correlation files, first is amp, second is correlation (0 to 1)
"""
function load_stacked_img(filename::AbstractString, demrsc::DemRsc; do_permute=true, return_amp=false)
    @unpack rows, cols = demrsc
    # Note: must be saved in c/row-major order, so loading needs a transpose
    # out_left usually is amplitude data
    # Usually we are interested in out_right
    #
    # First make container for all of data
    out = Array{Float32, 2}(undef, (2cols, rows))
    read!(filename, out)

    out_amp = @view out[1:cols, :]
    out_data = @view out[cols+1:end, :]

    if return_amp
        return do_permute ? permutedims(cat(out_amp, out_data, dims=3), (2, 1, 3)) : cat(out_amp, out_data, dims=3)
    else
        return do_permute ? permutedims(out_data) : out_data
    end
end

function load_geolist_from_h5(h5file::AbstractString)
    geo_strings = h5read(h5file, GEOLIST_DSET)
    return parse_geolist_strings(geo_strings)
end

function load_geolist_from_h5(h5file::AbstractString, dset::AbstractString)
    geo_strings = h5readattr(h5file, dset)[GEOLIST_DSET]
    return parse_geolist_strings(geo_strings)
end


function load_intlist_from_h5(h5file)
    # Note transpose, since it's stored as a N x 2 array
    # (which is flipped to 2 x N for julia
    int_strings = h5read(h5file, INTLIST_DSET)
    return parse_intlist_strings([Tuple(int_strings[:, ii]) for ii in 1:size(int_strings, 2)])
end

function load_intlist_from_h5(h5file, dset::AbstractString)
    int_strings = h5readattr(h5file, dset)[INTLIST_DSET]
    return parse_intlist_strings([Tuple(int_strings[:, ii]) for ii in 1:size(int_strings, 2)])
end

_parse(datestr) = Date(datestr, DATE_FMT)

function parse_intlist_strings(date_pairs::AbstractArray{<:AbstractString, 1})
    # For passing list of filenames ["20150101_20160101.int", ...]
    # `collect` used to make into an array of chars for strip
    date_pairs = [split(strip(d, collect(".int")), "_")[1:2] for d in date_pairs]
    date_tups = [Tuple(d) for d in date_pairs]
    return parse_intlist_strings(date_tups)
end

parse_intlist_strings(date_pairs) = [_parse.(pair) for pair in date_pairs]

parse_geolist_strings(geolist_str::AbstractArray{String}) = _parse.(geolist_str)


import Base.string
string(arr::AbstractArray{Date}, fmt="yyyymmdd") = Dates.format.(arr, fmt)
string(arr::AbstractArray{Tuple{Date, Date}}, fmt="yyyymmdd") = [Dates.format.(tup, fmt) for tup in arr]

"""Save the geolist as a list of strings to an dataset `h5file`"""
function save_geolist_to_h5(h5file::String, geolist::AbstractArray{Date}; overwrite=false)
    !check_dset(h5file, GEOLIST_DSET, overwrite) && return
    h5write(h5file, GEOLIST_DSET, string(geolist))
end

"""In this version, save the geolist to an attribute of `object` already in `h5file`"""
function save_geolist_to_h5(h5file::String, object::String, geolist::AbstractArray{Date}; overwrite=false)
    h5open(h5file, "cw") do f
        obj = f[object]
        # Delete if exists
        # TODO: wrap this attr one into check_dset
        overwrite && exists(attrs(obj), GEOLIST_DSET) && a_delete(obj, GEOLIST_DSET)

        attrs(obj)[GEOLIST_DSET] = string(geolist)
    end
end

"""Save the geolist as a list of strings to an dataset `h5file`"""
function save_intlist_to_h5(h5file::String, intlist::AbstractArray{Tuple{Date, Date}}; overwrite=false)
    !check_dset(h5file, INTLIST_DSET, overwrite) && return
    h5write(h5file, INTLIST_DSET, string(intlist))
end

# In HDF5, DemRscs are stored as JSON dicts
load_dem_from_h5(h5file, dset=DEM_RSC_DSET) = DemRsc(JSON.parse(h5read(h5file, dset)))

function save_dem_to_h5(h5file, dem_rsc::DemRsc, dset_name=DEM_RSC_DSET; overwrite=true)
    !check_dset(h5file, dset_name, overwrite) && return
    h5write(h5file, dset_name, JSON.json(dem_rsc))
end
# TODO:
# save_dem_to_h5(h5file, dem_rsc::DemRsc, dset_name=DEM_RSC_DSET; overwrite=true)

function save_hdf5_stack(h5file::AbstractString, dset_name::AbstractString, stack; overwrite::Bool=false, do_permute=true)
    # TODO: is there a way to combine this with normal saving of files?
    mode = overwrite ? "w" : "cw"  # Append mode is 'cw'  (why different than "a"?) 
    h5open(h5file, mode) do f 
        do_permute ? write(f, dset_name, permutedims(stack, (2, 1, 3))) :
                     write(f, dset_name, stack)
    end
end

"""Wrapper around h5read to account for the transpose
necessary when reading from python-written stacks
"""
function load_hdf5_stack(h5file::AbstractString, dset_name::AbstractString; do_permute=true)
    h5open(h5file) do f
        return do_permute ? permutedims(read(f[dset_name]), (2, 1, 3)) : read(f[dset_name])
    end
end

# If loading only certain layers, don't read all into memory
function load_hdf5_stack(h5file::AbstractString, dset_name::AbstractString, valid_layer_idxs; do_permute=true)
    nrows, ncols, _ = size(h5file, dset_name)
    h5open(h5file) do f
        dset = f[dset_name]
        out = Array{eltype(dset), ndims(dset)}(undef, (nrows, ncols, length(valid_layer_idxs)))
        for (tot_idx, v_idx) in enumerate(valid_layer_idxs)
            out[:, :, tot_idx] = dset[:, :, v_idx]
        end
        return do_permute ? permutedims(out, (2, 1, 3)) : out
    end
end

"""Sum the 3rd dim (layers) of a stack"""
function sum_hdf5_stack(h5file::AbstractString, dset_name::AbstractString, valid_layer_idxs)
    rows, cols, _ = size(h5file, dset_name)
    # For summing, using int instead of bool
    dtype = eltype(h5file, dset_name) <: Bool ? Int : eltype(h5file, dset_name)
    out = zeros(dtype, (rows, cols))

    h5open(h5file) do f
        d = f[dset_name]
        for ii in valid_layer_idxs
            out .+= @view d[:, :, ii][:, :, 1]   
        end
    end
    return out
end

"""Get the composite mask from the stack, true only where ALL pixels are masked"""
function load_mask(geolist::AbstractArray{Date, 1}; do_permute::Bool=true, fname="masks.h5", dset_name="geo")
    geolist_full = load_geolist_from_h5(fname)
    idxs = indexin(geolist, geolist_full)
    out = convert(Array{Bool}, sum_hdf5_stack(fname, dset_name, idxs) .> 0)
    return do_permute ? permutedims(out) : out
end

# TODO: this sucks and doesnt work well... if one day masks out half the map,
# the geo_sum will be 0 on half, 1 on half. if we ignore that date... this still
# says "it equals the maximum of the sum, so it's masked
"""Get the composite mask from the stack, true only where ALL pixels are masked"""
function load_mask(;do_permute::Bool=true, fname="masks.h5", dset="geo_sum")
    mask = h5read(fname, dset)
    mval = max(1, maximum(mask))  # Make sure we dont mask all 0s
    return do_permute ? permutedims(mask .== mval) : mask .== mval
end

"""Returns false if the dataset exists and overwrite is false

If overwrite is set to true, will delete the dataset to make
sure a new one can be created"""
function check_dset(h5file, dset_name, overwrite)
    h5open(h5file, "cw") do f
        if dset_name in names(f)
            println("$dset_name already exists in $h5file")
            if overwrite
                println("Overwrite true: Deleting.")
                o_delete(f, dset_name)
            else
                println("Skipping $dset_name: not deleting")
                return false
            end
        end

        return true
    end
end


# TODO: probably a better way to do this.. but can't figure out how to
# without allocating so many arrays that it's slow as Python
function load_stack(; file_list::Union{AbstractArray{AbstractString}, Nothing}=nothing, 
                    directory::Union{AbstractString, Nothing}=nothing,
                    file_ext::Union{AbstractString, Nothing}=nothing)
    if isnothing(file_list)
        file_list = find_files(file_ext, directory)
    end

    # demrsc = load(find_rsc_file(basepath=directory))
    test_arr = load(file_list[1])
    rows, cols = size(test_arr)
    T = eltype(test_arr)

    stack, buffer = _return_array(T, rows, cols, length(file_list))
    for (idx, f) in enumerate(file_list)
         read!(f, buffer)
            stack[:, :, idx] = buffer
    end

    return _permute(stack, cols)
end

# Using multiple dispatch to avoid if statements for 2 types of images
function _return_array(T::Union{Type{ComplexF32}, Type{Bool}}, rows::Int, cols::Int, file_len::Int)
    stack = Array{T, 3}(undef, (cols, rows, file_len))
    buffer = Array{T, 2}(undef, (cols, rows))
    stack, buffer
end

function _return_array(T::Type{Float32}, rows::Int, cols::Int, file_len::Int)
    stack = Array{T, 3}(undef, (2cols, rows, file_len))
    buffer = Array{T, 2}(undef, (2cols, rows))
    stack, buffer
end

# For normal .int, .geo complex/ masks, just transpose stack
function _permute(stack::AbstractArray{T, 3}, cols::Int) where {T <: Number}
    return permutedims(stack, (2, 1, 3))
end

# For normal weird stacked types, pick just the right half
function _permute(stack::AbstractArray{Float32, 3}, cols::Int)
    return permutedims(stack[cols+1:end, :, :], (2, 1, 3))
end



# function _load_stack_complex(file_list::AbstractArray{AbstractString}, rows::Int, cols::Int)
#     stack = Array{ComplexF32, 3}(undef, (cols, rows, length(file_list)))
#     buffer = Array{ComplexF32, 2}(undef, (cols, rows))
# 
#     for (idx, f) in enumerate(file_list)
#          read!(f, buffer)
#             stack[:, :, idx] = buffer
#     end
# 
#     return permutedims(stack, (2, 1, 3))
# end
# 
# 
# function _load_stack_stacked(file_list::AbstractArray{AbstractString}, rows::Int, cols::Int)
#     stack = Array{Float32, 3}(undef, (2cols, rows, length(file_list)))
#     buffer = Array{Float32, 2}(undef, (2cols, rows))
# 
#     for (idx, f) in enumerate(file_list)
#          read!(f, buffer)
#             stack[:, :, idx] = buffer
#     end
#     return permutedims(stack[cols+1:end, :, :], (2, 1, 3))
# end


function save(filename::AbstractString, array; do_permute=true, kwargs...)
    ext = get_file_ext(filename)

    if ext in BOOL_EXTS
        tofile(filename, array, do_permute=do_permute)
    elseif (ext in vcat(COMPLEX_EXTS, REAL_EXTS, ELEVATION_EXTS)) && (!(ext in STACKED_FILES))
        tofile(filename, _force_float32(array), do_permute=do_permute)
    elseif ext in STACKED_FILES
        # ndims(array) != 3 && throw(DimensionMismatch("array must be 3D [amp; data] to save as $filename"))
        if ndims(array) == 3
            amp = view(array, :, :, 1)
            data = view(array, :, :, 2)
        else
            println("!!! Warning: saving $filename with 1s for amplitude")
            data = array
            amp = ones(size(array))
        end
        # Handle the permuting for stacked here outside of function
        out = do_permute ? transpose(hcat(amp, data)) : vcat(amp, data)
        tofile(filename, _force_float32(out), do_permute=false)
    else
        error("Filetype not implemented: $filename")
    end
end

function tofile(filename, array; do_permute=true)  #, overwrite=true)
    # mode = overwrite ? "w" : "a"
    mode = "w"
    open(filename, mode) do f 
        do_permute ? write(f, transpose(array)) : write(f, array)
    end
end


_iscomplex(array) = eltype(array) <: Complex
# _isfloat(array) = eltype(array) <: AbstractFloat

function _force_float32(array) 
    if eltype(array) <: Complex
        return ComplexF32.(array)
    elseif eltype(array) <: AbstractFloat
        return Float32.(array)
    else
        return array
    end
end

"""Downsample a matrix by summing blocks of (row_looks, col_looks)

Cuts off values if the size isn't divisible by num looks
size = floor(rows / row_looks, cols / col_looks)
"""
function take_looks(image::AbstractArray{T}, row_looks, col_looks) where {T <: Number}
    (row_looks == 1 && col_looks == 1) && return image

    if ndims(image) > 2
        return cat((take_looks(image[:, :, i], row_looks, col_looks)
                    for i in 1:size(image, 3))..., dims=3)
    end

    nrows, ncols = size(image)
    nr = div(nrows, row_looks)
    nc = div(ncols, col_looks)

    # Can't sum into a Bool array
    T == Bool ? outtype = Float32 : outtype = T
    out = zeros(outtype, nr, nc)

    @inbounds Threads.@threads for j = 1:nc
        @inbounds for i = 1:nr
            indx_i = 1+(i-1)*row_looks:i*row_looks
            indx_j = 1+(j-1)*col_looks:j*col_looks
            @views out[i, j] = sum(image[indx_i, indx_j])
        end
    end

    out ./= (row_looks * col_looks)
    return T == Bool ? Bool.(out .> 0) : out
end

# If we pass something else like a dem_rsc
take_looks(other, row_looks, col_looks) = other

# For future:
# cc_patch(a, b) = real(abs(sum(a .* conj.(b))) / sqrt(sum(a .* conj.(a)) * sum(b .* conj.(b))))


"""Reads the list of igrams to returns Array of Igrams

Args:
    directory (str): path to the igram directory
    parse (bool): output as parsed datetime tuples. False returns the filenames
    filename (string): name of a file with .geo filenames

Returns:
    tuple(date, date) of (early, late) dates for all igrams (if parse=True)
        if parse=False: returns list[str], filenames of the igrams
"""
function find_igrams(;directory::AbstractString=".", parse::Bool=true, filename::AbstractString="")
    igram_file_list = isempty(filename) ? find_files("*.int", directory) : readlines(filename)
    intlist_strings = [splitpath(fname)[end] for fname in igram_file_list]
    isempty(intlist_strings) && return intlist_strings
    return parse ? parse_intlist_strings(intlist_strings) : intlist_strings
end


"""Reads in the list of .geo files used, in time order

Can also pass a filename containing .geo files as lines.

Args:
    directory (str): path to the geolist file or directory
    parse (bool): output as parsed datetime tuples. False returns the filenames
    filename (string): name of a file with .geo filenames

Returns:
    list[date]: the parse dates of each .geo used, in date order

"""
# TODO: parsers.Sentinel...
function find_geos(;directory::AbstractString=".", parse::Bool=true, filename::AbstractString="")
    geo_file_list = isempty(filename) ? find_files("*.geo", directory) : readlines(filename)

    !parse && return geo_file_list

    # Stripped of path for parser
    geolist_strings = [splitpath(fname)[end] for fname in geo_file_list]

    return sort(parse_geolist_strings(_strip_geoname.(geolist_strings)))

    # return sorted(parse_geolist_strings([_strip_geoname(geo) for geo in geolist_strings]))
    # if re.match(r'S1[AB]_\d{8}\.geo', geolist_strings[0]):  # S1A_YYYYmmdd.geo
    #     return sorted([_parse(_strip_geoname(geo)) for geo in geolist_strings])
    # elseif re.match(r'\d{8}', geolist_strings[0]):  # YYYYmmdd , just a date string
    #     return sorted([_parse(geo) for geo in geolist_strings])
    # else  # Full sentinel product name
    #     return sorted([apertools.parsers.Sentinel(geo).start_time.date() for geo in geolist_strings])
    # end
end


"""Leaves just date from format S1A_YYYYmmdd.geo"""
_strip_geoname(name) = reduce(replace, ["S1A_" => "", "S1B_" => "", ".geo" => ""], init=name)


end # module
