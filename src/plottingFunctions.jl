using Plots
plotlyjs()

function ploteye(ellipsoids::Dict{String,Ellipsoid}, l::Light, opticunits::Vector{OpticUnit}; n_rays::Int = 15)
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

function plotpixels(opticunits::Vector{OpticUnit}, l::Light; n_rays::Int = 100000)
    psf = getpsf(opticunits, l, n_rays)
    ph = []
    for (k, v) in psf
        push!(ph, heatmap(collect(v), title = k, aspect_ratio = :equal, reuse = false))
    end
    plot(reuse = false, ph...)
    gui()
end

function plotfwhm_aperture(opticunits::Vector{OpticUnit}, source_distance::Number, min_aperture::Number, max_aperture::Number; n_data::Int = 10, n_rays::Int = 10000, angle_step = 0.25°)
    apertures = linspace(min_aperture, max_aperture, n_data)
    fwhms = Dict(k => zeros(n_data) for k in retinas)
    for (i, aperture) in enumerate(apertures), (k, v) in getfwhm(opticunits, source_distance, aperture, n_rays, angle_conv(angle_step))
        fwhms[k][i] = get(v)
    end
    plot(reuse=false)
    [plot!(apertures, v, label = k) for (k, v) in fwhms]
    plot!(xlabel = "Aperture (μm)", ylabel = "FWHM (°)")
    gui()
end

function plotfwhm_distance(opticunits::Vector{OpticUnit}, aperture::Number, min_distance::Number, max_distance::Number; n_data::Int = 10, n_rays::Int = 10000, angle_step = 0.25°, u = cm)
    viewdistances = logspace(log10(min_distance), log10(max_distance), n_data)
    fwhms = Dict(k => zeros(n_data) for k in retinas)
    for (i, viewdistance) in enumerate(viewdistances), (k, v) in getfwhm(opticunits, viewdistance, aperture, n_rays, angle_conv(angle_step))
        fwhms[k][i] = get(v)
    end
    plot(xscale = :log10, xlabel = "Source distance ($u)", ylabel = "FWHM (°)", reuse=false)
    [plot!(ustrip.(uconvert.(u, μm.*viewdistances)), v, label = k) for (k, v) in fwhms]
    gui()
end

function plotfwhm_aperture_distance(opticunits::Vector{OpticUnit}, min_aperture::Number, max_aperture::Number, min_distance::Number, max_distance::Number; n_data::Int = 10, n_rays::Int = 10000, angle_step = 0.25°)
    viewdistances = logspace(log10(min_distance), log10(max_distance), n_data)
    apertures = linspace(min_aperture, max_aperture, n_data)
    fwhms = Dict(k => zeros(n_data,n_data) for k in retinas)
    for (i, viewdistance) in enumerate(viewdistances), (j, aperture) in enumerate(apertures), (k, v) in getfwhm(opticunits, viewdistance, aperture, n_rays, angle_conv(angle_step))
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


