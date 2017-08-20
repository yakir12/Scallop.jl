using Rayden, Scallop


function format_radii{T <: Number}(x::Tuple{T})
    xx1, = x
    x1 = format(xx1)
    return Vec(x1, x1, x1)
end
function format_radii{T <: Number}(x::Tuple{T, T})
    xx1, xx2 = x
    x1 = format(xx1)
    x2 = format(xx2)
    return Vec(x1, x1, x2)
end
function format_radii{T <: Number}(xx::T) 
    x = format(xx)
    return Vec(x, x, x)
end

_format(x::T) where T <: Number = uconvert(rad, x)
_format(x::T) where T <: Unitful.Length = uconvert(μm, x)
format(x::T) where T <: Number = Float64(ustrip(_format(x)))
format(x::Int) = Float64(x)
format(x::Float64) = x

morphing_factor = format(morphing_factor)
source_distance = format(source_distance)

min_aperture = format(min_aperture)
max_aperture = format(max_aperture)
min_distance = format(min_distance)
max_distance = format(max_distance)


