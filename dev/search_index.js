var documenterSearchIndex = {"docs":
[{"location":"#Sario.jl-1","page":"Home","title":"Sario.jl","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Modules = [Sario]","category":"page"},{"location":"#Sario.get_file_ext-Tuple{AbstractString}","page":"Home","title":"Sario.get_file_ext","text":"Extracts the file extension, including the \".\" (e.g.: .slc)\n\n\n\n\n\n","category":"method"},{"location":"#Sario.load-Tuple{AbstractString}","page":"Home","title":"Sario.load","text":"load(filename::AbstractString; \n     rsc_file::Union{AbstractString, Nothing}=nothing,\n     looks::Tuple{Int, Int}=(1, 1), do_permute=true, \n     dset_name::AbstractString=\"\", return_amp::Bool=false)\n\nMain entry point to loading files.\n\nExamines file type of filename and runs appropriate load function.\n\n\n\n\n\n","category":"method"},{"location":"#Sario.load_hdf5_stack-Tuple{AbstractString,AbstractString}","page":"Home","title":"Sario.load_hdf5_stack","text":"Wrapper around h5read to account for the transpose necessary when reading from python-written stacks\n\n\n\n\n\n","category":"method"},{"location":"#Sario.load_mask-Tuple{AbstractArray{Dates.Date,1}}","page":"Home","title":"Sario.load_mask","text":"Get the composite mask from the stack, true only where ALL pixels are masked\n\n\n\n\n\n","category":"method"},{"location":"#Sario.load_mask-Tuple{}","page":"Home","title":"Sario.load_mask","text":"Get the composite mask from the stack, true only where ALL pixels are masked\n\n\n\n\n\n","category":"method"},{"location":"#Sario.save-Tuple{AbstractString,Any}","page":"Home","title":"Sario.save","text":"save(filename::AbstractString, array; do_permute=true, kwargs...)\n\nSaves a 2D/3D array as filename. If the file was loaded with do_permute=false, it should also be saved as do_permute=false \n\nFor stacked files, the expected format of array is a 3D stack with array[:,:,1] being the amplitude image, array[:,:,2] being the data. If only a 2D image is passed, an dummy amplitude image of all 1s is made.\n\n\n\n\n\n","category":"method"},{"location":"#Sario.RangeTuple","page":"Home","title":"Sario.RangeTuple","text":"Load subset of a file on disk using range\n\n\n\n\n\n","category":"type"},{"location":"#Sario._get_rsc_data-Tuple{Any,Union{Nothing, AbstractString}}","page":"Home","title":"Sario._get_rsc_data","text":"Handle getting the .rsc data from either an image filename, or .rsc filename\n\n\n\n\n\n","category":"method"},{"location":"#Sario._get_seek_position-NTuple{4,Any}","page":"Home","title":"Sario._get_seek_position","text":"For single element reading in binary files, seek to the right row, col\n\n\n\n\n\n","category":"method"},{"location":"#Sario._strip_geoname-Tuple{Any}","page":"Home","title":"Sario._strip_geoname","text":"Leaves just date from format S1A_YYYYmmdd.geo\n\n\n\n\n\n","category":"method"},{"location":"#Sario.check_dset-Tuple{Any,Any,Any}","page":"Home","title":"Sario.check_dset","text":"Returns false if the dataset exists and overwrite is false\n\nIf overwrite is set to true, will delete the dataset to make sure a new one can be created\n\n\n\n\n\n","category":"method"},{"location":"#Sario.find_igrams-Tuple{}","page":"Home","title":"Sario.find_igrams","text":"find_igrams(;directory::AbstractString=\".\", parse::Bool=true, filename::AbstractString=\"\")\n\nReads the list of igrams to returns Array of Igrams\n\nArgs:     directory (str): path to the igram directory     parse (bool): output as parsed datetime tuples. False returns the filenames     filename (string): name of a file with .geo filenames\n\nReturns:     tuple(date, date) of (early, late) dates for all igrams (if parse=True)         if parse=False: returns list[str], filenames of the igrams\n\n\n\n\n\n","category":"method"},{"location":"#Sario.load_dem_rsc-Tuple{Any,Vararg{Any,N} where N}","page":"Home","title":"Sario.load_dem_rsc","text":"Convert the text file into the DemRsc struct Starts with a dict to gather all fiields, then unpacks using keywords args\n\n\n\n\n\n","category":"method"},{"location":"#Sario.load_elevation-Tuple{Any,DemRsc}","page":"Home","title":"Sario.load_elevation","text":"Loads a digital elevation map from either .hgt file or .dem\n\n.hgt is the NASA SRTM files given. Documentation on format here: https://dds.cr.usgs.gov/srtm/version21/Documentation/SRTMTopo.pdf Key point: Big-endian 2 byte (16-bit) integers\n\n.dem is format used by Zebker geo-coded and ROI-PAC SAR software Only difference is data is stored little-endian (like other SAR data)\n\nNote on both formats: gaps in coverage are given by INT_MIN -32768, so either manually set data(data == np.min(data)) = 0, or something like      data = clamp(data, -10000, Inf)\n\n\n\n\n\n","category":"method"},{"location":"#Sario.load_stacked_img-Tuple{AbstractString,DemRsc}","page":"Home","title":"Sario.load_stacked_img","text":"loadstackedimg is for the weird \"ALTLINEDATA\" formatted images: https://web.stanford.edu/group/radar/softwareandlinks/sw/snaphu/snaphu_man1.html#FILE%20FORMATS\n\nFormat is two stacked matrices: [[first], [second]] where the first \"cols\" number of floats are the first matrix, next \"cols\" are second, etc.\n\nFor .unw height files, the first is amplitude, second is phase (unwrapped) For .cc correlation files, first is amp, second is correlation (0 to 1)\n\n\n\n\n\n","category":"method"},{"location":"#Sario.save_geolist_to_h5-Tuple{String,AbstractArray{Dates.Date,N} where N}","page":"Home","title":"Sario.save_geolist_to_h5","text":"save_geolist_to_h5(h5file::String, geolist::AbstractArray{Date}; overwrite=false)\n\nSave the geolist as a list of strings to an dataset h5file\n\n\n\n\n\n","category":"method"},{"location":"#Sario.save_geolist_to_h5-Tuple{String,String,AbstractArray{Dates.Date,N} where N}","page":"Home","title":"Sario.save_geolist_to_h5","text":"In this version, save the geolist to an attribute of object already in h5file\n\n\n\n\n\n","category":"method"},{"location":"#Sario.save_intlist_to_h5-Tuple{String,AbstractArray{Tuple{Dates.Date,Dates.Date},N} where N}","page":"Home","title":"Sario.save_intlist_to_h5","text":"Save the geolist as a list of strings to an dataset h5file\n\n\n\n\n\n","category":"method"},{"location":"#Sario.sum_hdf5_stack-Tuple{AbstractString,AbstractString,Any}","page":"Home","title":"Sario.sum_hdf5_stack","text":"Sum the 3rd dim (layers) of a stack without loading all into memory\n\n\n\n\n\n","category":"method"},{"location":"#Sario.take_looks-Union{Tuple{T}, Tuple{AbstractArray{T,N} where N,Any,Any}} where T<:Number","page":"Home","title":"Sario.take_looks","text":"Downsample a matrix by summing blocks of (rowlooks, collooks)\n\nCuts off values if the size isn't divisible by num looks size = floor(rows / rowlooks, cols / collooks)\n\n\n\n\n\n","category":"method"}]
}
