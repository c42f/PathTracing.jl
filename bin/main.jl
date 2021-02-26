# For use as a JuliaHub Application

using PathTracing
using PNGFiles

w = parse(Int64, get(ENV, "width", "640"))
h = parse(Int64, get(ENV, "height", "480"))
spp = parse(Int64, get(ENV, "samples_per_pixel", "4"))

@info "Rendering $w×$h image (spp $spp)"

scene = PathTracing.cornell_box_scene()
c = PathTracing.render(scene, w, h, spp)
img = PathTracing.to_sRGB.(permutedims(reverse(c, dims=2)))
PNGFiles.save("image.png", img, filters=0)

ENV["RESULTS_FILE"] = "image.png"
