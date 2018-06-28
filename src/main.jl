using Unitful
import Unitful: nm, μm, mm, cm, m, km, inch, ft, mi, rad, °
# include("prepare.jl")
mutable struct RI 
    water
    cornea
    lens
    lens2distal_retina
    distal_retina2proximal_retina
    gap
    behind_mirror
end

mutable struct Thick
    water
    cornea
    lens
    lens2distal_retina
    distal_retina2proximal_retina
    gap
    behind_mirror
end

mutable struct Radii
    cornea_distal
    lens_distal
    lens_proximal
    distal_retina
    proximal_retina
    mirror
end

ri = RI(zeros(7)...)
thick = Thick(μm*zeros(7)...)
radii = Radii(μm*zeros(6)...)

##################################################
############## EDIT BELOW HERE ###################
##################################################
include("variables.jl")
#= Light properties
source_distance = 1m
source_angle    = 0°

# Eye morphology
aperture        = 251μm
morphing_factor = 1.12

## Volumes
### Refractive indices

ri.water                         = 1.334
ri.cornea                        = 1.37
ri.lens                          = 1.42
ri.lens2distal_retina            = 1.35
ri.distal_retina2proximal_retina = 1.35
ri.gap                           = 1.34

### Thicknesses

thick.cornea                        = 23μm
thick.lens                          = 215μm
thick.lens2distal_retina            = 6μm
thick.distal_retina2proximal_retina = 81μm
thick.gap                           = 49μm

## Interfaces
### Radii of the ellipsoids that form the interface

radii.cornea_distal   = (244μm, 227μm)
radii.lens_distal     = (195μm, 213μm)
radii.lens_proximal   = 337μm
radii.distal_retina   = 337μm
radii.proximal_retina = 337μm
radii.mirror          = 417μm

# Figure ranges
min_aperture = 100μm
max_aperture = 300μm
min_distance = 1cm
max_distance = 1m=#

##################################################
############## EDIT ABOVE HERE ###################
##################################################

# include("plotit.jl")

using Rayden, Scallop

thick.behind_mirror = Inf*μm
ri.behind_mirror = 0.0
thick.water = 0μm

function format_radii(x::Tuple{T}) where T <: Number
    xx1, = x
    x1 = format(xx1)
    return Vec(x1, x1, x1)
end
function format_radii(x::Tuple{T, T}) where T <: Number
    xx1, xx2 = x
    x1 = format(xx1)
    x2 = format(xx2)
    return Vec(x1, x1, x2)
end
function format_radii(xx::T) where T <: Number 
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

volumes = Dict(String(k) => Volume(Float64(getfield(ri, k)), format(getfield(thick, k))) for k in fieldnames(ri))
planes = Dict(String(k) => Plane(format_radii(getfield(radii, k)), !r"distal$"(String(k))) for k in fieldnames(radii))
aperture = format(aperture)

aperture, ellipsoids = getmembranes(volumes, planes, morphing_factor, aperture)
opticunits = scallop(ellipsoids, volumes)
l = getLight(source_distance, aperture, ellipsoids["cornea_distal"], format(source_angle))

min_aperture = format(min_aperture)
max_aperture = format(max_aperture)
min_distance = format(min_distance)
max_distance = format(max_distance)


# include("plottingFunctions.jl")

using Plots
plotlyjs()

# function ploteye(ellipsoids::Dict{String,Ellipsoid}, l::Light, opticunits::Vector{OpticUnit}; n_rays::Int = 15)
function ploteye(;n_rays::Int = 15)
    include("main.jl")
    x, y = coordinates2d(ellipsoids)
    plot(x, y, aspect_ratio = :equal, leg=false, color = :black, reuse=false)
    r = Ray()
    for b in linspace(0, 1, n_rays), b2 in [0, 0.5]
        getRay!(r, l, b, b2)
        x = Float64[]
        y = Float64[]
        push!(x, r.orig[1])
        push!(y, r.orig[3])
        complete = true
        for i in opticunits
            if raytrace!(r, i)
                complete = false
                break
            end
            push!(x, r.orig[1])
            push!(y, r.orig[3])
        end
        complete && plot!(x, y, color = :grey)
    end
    plot!(ylims = (-500,100), xlims = (-300,300), xlabel = "μm", ylabel = "μm")
    gui()
end


# function plotpixels(opticunits::Vector{OpticUnit}, l::Light; n_rays::Int = 100000)
function plotpixels(;n_rays::Int = 100000)
    include("main.jl")
    psf = getpsf(opticunits, l, n_rays)
    ph = []
    for (k, v) in psf
        push!(ph, heatmap(collect(v), title = k, aspect_ratio = :equal, reuse = false))
    end
    plot(reuse = false, ph...)
    gui()
end

# function plotfwhm_aperture(opticunits::Vector{OpticUnit}, source_distance::Number, min_aperture::Number, max_aperture::Number; n_data::Int = 10, n_rays::Int = 10000, angle_step = 0.25°)
function plotfwhm_aperture(;n_data::Int = 10, n_rays::Int = 10000, angle_step = 0.25°)
    include("main.jl")
    apertures = linspace(min_aperture, max_aperture, n_data)
    fwhms = Dict(k => zeros(n_data) for k in retinas)
    for (i, aperture) in enumerate(apertures), (k, v) in getfwhm(opticunits, source_distance, aperture, n_rays, format(angle_step))
        fwhms[k][i] = get(v)
    end
    plot(reuse=false)
    [plot!(apertures, v, label = k) for (k, v) in fwhms]
    plot!(xlabel = "Aperture (μm)", ylabel = "FWHM (°)")
    gui()
end

# function plotfwhm_distance(opticunits::Vector{OpticUnit}, aperture::Number, min_distance::Number, max_distance::Number; n_data::Int = 10, n_rays::Int = 10000, angle_step = 0.25°, u = cm)
function plotfwhm_distance(;n_data::Int = 10, n_rays::Int = 10000, angle_step = 0.25°)
    include("main.jl")
    viewdistances = logspace(log10(min_distance), log10(max_distance), n_data)
    fwhms = Dict(k => zeros(n_data) for k in retinas)
    for (i, viewdistance) in enumerate(viewdistances), (k, v) in getfwhm(opticunits, viewdistance, aperture, n_rays, format(angle_step))
        fwhms[k][i] = get(v)
    end
    plot(xscale = :log10, xlabel = "Source distance (cm)", ylabel = "FWHM (°)", reuse=false)
    [plot!(ustrip.(uconvert.(cm, μm.*viewdistances)), v, label = k) for (k, v) in fwhms]
    gui()
end

# function plotfwhm_aperture_distance(opticunits::Vector{OpticUnit}, min_aperture::Number, max_aperture::Number, min_distance::Number, max_distance::Number; n_data::Int = 10, n_rays::Int = 10000, angle_step = 0.25°)
function plotfwhm_aperture_distance(;n_data::Int = 10, n_rays::Int = 10000, angle_step = 0.25°)
    include("main.jl")
    viewdistances = logspace(log10(min_distance), log10(max_distance), n_data)
    apertures = linspace(min_aperture, max_aperture, n_data)
    fwhms = Dict(k => zeros(n_data,n_data) for k in retinas)
    for (i, viewdistance) in enumerate(viewdistances), (j, aperture) in enumerate(apertures), (k, v) in getfwhm(opticunits, viewdistance, aperture, n_rays, format(angle_step))
        fwhms[k][i,j] = get(v)
    end
    ph = []
    for (k, v) in fwhms
        push!(ph, heatmap(apertures, ustrip.(uconvert.(cm, μm.*viewdistances)), v, title = k, xlabel = "Aperture (μm)", ylabel = "Source distance (cm)", reuse = false))
    end
    plot(ph..., yscale = :log10, zlims = (0,20), reuse=false)
    gui()
end

#=l = getLight(1e9, aperture, ellipsoids["cornea_distal"], .05)
psf = get2dpsf(opticunits, l, 100000)
plot()
for (k, v) in psf
    plot!(collect(v), label = k)
end
gui()=#


#=nt = 20
thetas = linspace(0, .5, nt)
ratios = Dict(k => zeros(nt) for k in retinas)
for (i, θ) in enumerate(thetas)
    l = getLight(1e9, aperture, ellipsoids["cornea_distal"], θ)
    psf = get2dpsf(opticunits, l, 100000)
    for k in keys(ratios)
        ratios[k][i] = maximum(psf[k])/psf[k][0]
    end
end
plot(ylims = (0,5))
for (k, v) in ratios
    plot!(thetas, v, label = k)
end
gui()=#


