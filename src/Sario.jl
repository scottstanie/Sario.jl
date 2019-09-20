"""Module Sario
"""


__precompile__(true)

module Sario

export load, save, load_stack, load_hdf5_stack, load_mask, load_geolist_from_h5,
       load_intlist_from_h5, get_file_ext, find_files, save_hdf5_stack

# Include apertools Python modeules here to make available to all
using PyCall
const pysario = PyNULL()

function __init__()
    # TODO: remove all the backups so we don't need this...
    copy!(pysario, pyimport("apertools.sario"))
end

import Glob
using HDF5
using Dates

const DATE_FMT = "yyyymmdd"
const SENTINEL_EXTS = [".geo", ".cc", ".int", ".amp", ".unw", ".unwflat"]
const COMPLEX_EXTS = [".int", ".slc", ".geo", ".cc", ".unw", ".unwflat", ".mlc", ".grd"]
const REAL_EXTS = [".amp", ".cor", ".mlc", ".grd"]
const ELEVATION_EXTS = [".dem", ".hgt"]
# These file types are not simple complex matrices: see load_stacked_img for detail
# .unwflat are same as .unw, but with a linear ramp removed
const STACKED_FILES = [".cc", ".unw", ".unwflat"]
const IMAGE_EXTS = [".png", ".tif", ".tiff", ".jpg"]
const BOOL_EXTS = [".mask"]

const LOAD_IN_PYTHON = vcat(IMAGE_EXTS, [".geojson", ".npy"])

# For HDF5 saving
const DEM_RSC_DSET = "dem_rsc"
const STACK_DSET = "stack"
const GEOLIST_DSET = "geo_dates"
const INTLIST_DSET = "int_dates"

find_files(ext, directory=".") = sort(Glob.glob("*"*ext, directory))

"""Extracts the file extension, including the "." (e.g.: .slc)"""
get_file_ext(filename::AbstractString) = splitext(filename)[end]


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
function load(filename::AbstractString; rsc_file::Union{AbstractString, Nothing}=nothing,
              looks::Tuple{Int, Int}=(1, 1), do_permute=true, return_amp=false)
    ext = get_file_ext(filename)

    # For now, just pass through unimplemented extensions to Python
    if ext in LOAD_IN_PYTHON
        return take_looks(pysario.load(filename), looks...)
    elseif ext == ".rsc"
        return _get_rsc_data(filename, filename)
    elseif ext in ELEVATION_EXTS
        return take_looks(load_elevation(filename), looks...)
    end

    # Sentinel files should have .rsc file: check for dem.rsc, or elevation.rsc
    rsc_data = _get_rsc_data(filename, rsc_file)

    if ext in STACKED_FILES
        return take_looks(load_stacked_img(filename, rsc_data, do_permute=do_permute, return_amp=return_amp), looks...)
    elseif ext in BOOL_EXTS
        return take_looks(load_bool(filename, rsc_data, do_permute=do_permute), looks...)
    else
        return take_looks(load_complex(filename, rsc_data, do_permute=do_permute), looks...)
    # having rsc_data implies that this is not a UAVSAR file, so is complex
    # TODO: haven"t transferred over UAVSAR functions, so no load_real yet
    end

end

"""Load one element of a file on disk (avoid reading in all of huge file"""
# TODO: Load a chunk of a file now?
function load(filename::AbstractString, row_col::Tuple{Int, Int}; rsc_file::Union{AbstractString, Nothing}=nothing)
    data_type = _get_data_type(filename)

    rsc_data = _get_rsc_data(filename, rsc_file)
    num_rows, num_cols = rsc_data["file_length"], rsc_data["width"]

    row, col = row_col
    if row < 1 || col < 1 || row > num_rows || col > num_cols
        throw(DomainError((row, col), " out of bounds for $filename of size ($num_rows, $num_cols)"))
    end
    seek_pos = _get_seek_position(row, col, num_cols, data_type)
    
    open(filename) do f
        seek(f, seek_pos)
        # This read syntax loads single `data_type`
        return read(f, data_type)
    end
end

"""For single element reading in binary files, seek to the right row, col"""
_get_seek_position(row, col, num_cols, data_type) = sizeof(data_type) * ((col - 1) + (num_cols * (row - 1)) )


function _get_rsc_data(filename, rsc_file)
    ext = get_file_ext(filename)

    rsc_data = nothing
    if !isnothing(rsc_file)
        rsc_data = pysario.load(rsc_file)
    elseif ext in vcat(SENTINEL_EXTS, ELEVATION_EXTS, BOOL_EXTS)
        rsc_file = pysario.find_rsc_file(filename)
        rsc_data = pysario.load(rsc_file)
    end

    if !isnothing(rsc_data)
        rsc_data = convert(Dict{String, Any}, rsc_data)
    end
    return rsc_data
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
        rsc_file = pysario.find_rsc_file(filename)
        dem_rsc = pysario.load(rsc_file)
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


function load_complex(filename::AbstractString, rsc_data::Dict{<:AbstractString, Any}; do_permute=true)
    return _load_bin_matrix(filename, rsc_data, ComplexF32, do_permute)
end

function load_bool(filename::AbstractString, rsc_data::Dict{<:AbstractString, Any}; do_permute=true)
    return _load_bin_matrix(filename, rsc_data, Bool, do_permute)
end

function _load_bin_matrix(filename, rsc_data, dtype, do_permute::Bool)
    rows = rsc_data["file_length"]
    cols = rsc_data["width"]
    # Note: must be saved in c/row-major order, so loading needs a transpose
    out = Array{dtype, 2}(undef, (cols, rows))

    read!(filename, out)
    return do_permute ? permutedims(out) : out
end

"""Get the composite mask from the stack, true only where ALL pixels are masked"""
function load_mask(do_permute=true, fname="masks.h5", dset="geo_sum")
    mask = h5read(fname, dset)
    return do_permute ? permutedims(mask .== maximum(mask)) : mask .== maximum(mask)
end


"""load_stacked_img is for weirdly formatted images:

Format is two stacked matrices:
    [[first], [second]] where the first "cols" number of floats
    are the first matrix, next "cols" are second, etc.
For .unw height files, the first is amplitude, second is phase (unwrapped)
For .cc correlation files, first is amp, second is correlation (0 to 1)
"""
function load_stacked_img(filename::AbstractString, rsc_data::Dict{<:AbstractString, Any}; do_permute=true, return_amp=false)
    rows = rsc_data["file_length"]
    cols = rsc_data["width"]
    # Note: must be saved in c/row-major order, so loading needs a transpose
    # out_left usually is amplitude data
    # Usually we are interested in out_right
    #
    # First make container for all of data
    out = Array{Float32, 2}(undef, (2cols, rows))
    read!(filename, out)

    # TODO: port over rest of code for handling amp (if we care about that)
    out_amp = @view out[1:cols, :]
    out_data = @view out[cols+1:end, :]
    if return_amp
        return do_permute ? permutedims(cat(out_amp, out_data, dims=3), (2, 1, 3)) : cat(out_amp, out_data, dims=3)
    else
        return do_permute ? permutedims(out_data) : out_data
    end
end

function load_geolist_from_h5(h5file::AbstractString)
    h5open(h5file) do f
        geo_strings = read(f, GEOLIST_DSET)
        return pysario.parse_geolist_strings(geo_strings)
    end
end


function load_intlist_from_h5(h5file)
    h5open(h5file) do f
        int_strings = read(f, INTLIST_DSET)
        # Note transpose, since it's stored as a N x 2 array
        # (which is flipped to 2 x N for julia
        return pysario.parse_intlist_strings(int_strings')
    end
end

function save_geolist_to_h5(h5file, geolist; overwrite=false)
    h5open(h5file, "cw") do f
        # Delete if exists
        overwrite && GEOLIST_DSET in names(f) && o_delete(f, GEOLIST_DSET)

        geo_strings = Dates.format.(geolist, DATE_FMT)
        write(f, GEOLIST_DSET, geo_strings)
    end
end

function save_hdf5_stack(h5file::AbstractString, dset_name::AbstractString, stack; overwrite::Bool=true, do_permute=true)
    # TODO: is there a way to combine this with normal saving of files?
    if overwrite
        mode = "w"
    else
        mode = "cw"  # Append mode  (why different than "a"?)
    end
    h5open(h5file, mode) do f 
        do_permute ? write(f, dset_name, permutedims(stack, (2, 1, 3))) :
                     write(f, dset_name, stack)
    end
end

"""Wrapper around h5read to account for the transpose
necessary when reading from python-written stacks

"""
function load_hdf5_stack(h5file::AbstractString, dset_name::AbstractString)
    h5open(h5file) do f
        return permutedims(read(f[dset_name]), (2, 1, 3))
    end
end

# If loading only certain layers, don't read all into memory
function load_hdf5_stack(h5file::AbstractString, dset_name::AbstractString, valid_layer_idxs)
    nrows, ncols, _ = size(h5file, dset_name)
    h5open(h5file) do f
        dset = f[dset_name]
        out = Array{eltype(dset), ndims(dset)}(undef, (nrows, ncols, length(valid_layer_idxs)))
        for (tot_idx, v_idx) in enumerate(valid_layer_idxs)
            out[:, :, tot_idx] = dset[:, :, v_idx]
        end
        return permutedims(out, (2, 1, 3))
    end
end


# TODO: probably a better way to do this.. but can't figure out how to
# without allocating so many arrays that it's slow as Python
function load_stack(; file_list::Union{Array{AbstractString}, Nothing}=nothing, 
                    directory::Union{AbstractString, Nothing}=nothing,
                    file_ext::Union{AbstractString, Nothing}=nothing)
    if isnothing(file_list)
        file_list = find_files(file_ext, directory)
    end

    # rsc_data = pysario.load(pysario.find_rsc_file(basepath=directory))
    test_arr = load(file_list[1])
    rows, cols = size(test_arr)
    T = eltype(test_arr)

    stack, buffer = _return_array(T, rows, cols, length(file_list))
    for (idx, f) in enumerate(file_list)
         read!(f, buffer)
            stack[:, :, idx] = buffer
    end

    return _permute(stack, cols)
    # throw(ArgumentError(file_ext, "cant make stack from $file_ext"))
    # end
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
function _permute(stack::Array{T, 3}, cols::Int) where {T <: Number}
    return permutedims(stack, (2, 1, 3))
end

# For normal weird stacked types, pick just the right half
function _permute(stack::Array{Float32, 3}, cols::Int)
    return permutedims(stack[cols+1:end, :, :], (2, 1, 3))
end



# function _load_stack_complex(file_list::Array{AbstractString}, rows::Int, cols::Int)
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
# function _load_stack_stacked(file_list::Array{AbstractString}, rows::Int, cols::Int)
#     stack = Array{Float32, 3}(undef, (2cols, rows, length(file_list)))
#     buffer = Array{Float32, 2}(undef, (2cols, rows))
# 
#     for (idx, f) in enumerate(file_list)
#          read!(f, buffer)
#             stack[:, :, idx] = buffer
#     end
#     return permutedims(stack[cols+1:end, :, :], (2, 1, 3))
# end


function save(filename::AbstractString, array ; kwargs...)
    ext = get_file_ext(filename)

    if ext in BOOL_EXTS
        tofile(filename, array, kwargs...)
    elseif (ext in vcat(COMPLEX_EXTS, REAL_EXTS, ELEVATION_EXTS)) && (!(ext in STACKED_FILES))
        tofile(filename, _force_float32(array))
    elseif ext in STACKED_FILES
        # ndims(array) != 3 && throw(DimensionMismatch("array must be 3D [amp; data] to save as $filename"))
        if ndims(array) == 3
            amp = view(array, :, :, 1)
            data = view(array, :, :, 2)
        else
            println("Warning: saving $filename with 1s for amplitude")
            data = array
            amp = ones(size(array))
        end
        tofile(filename, _force_float32(hcat(amp, data)))
    else
        pysario.save(filename, array, kwargs...)
    end
end

function tofile(filename, array, do_permute=true)  #, overwrite=true)
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
function take_looks(image::Array{T}, row_looks, col_looks) where {T <: Number}
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

end # module
