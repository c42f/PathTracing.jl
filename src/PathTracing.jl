module PathTracing

# The following started out as a port of smallpt
# https://www.kevinbeason.com/smallpt
#
# Unlike the original, we're not restricting ourselves to 99 lines - we're
# permitting ourselves to expand and generalize a bit for the sake of clarity.

using StaticArrays, LinearAlgebra, PNGFiles, Random, Colors

const Vec = SVector{3,Float64}
const V = SA{Float64}
const C = SA{Float32}
const CVec = SVector{3,Float32}

struct Ray
    o::Vec # Origin
    d::Vec # Direction (normalized)
end

@enum Refl_t DIFF SPEC REFR

struct Sphere
    rad::Float64
    p::Vec  # Position
    e::CVec # Emission
    c::CVec # Color
    refl::Refl_t
end

# Returns distance, 0 if nohit
function intersect(s::Sphere, r::Ray)::Float64
    # Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
    op = s.p - r.o
    b = dot(op, r.d)
    det = b^2 - dot(op,op) + s.rad^2
    det < 0 && return 0
    eps = 1e-4
    det = sqrt(det)
    t = b - det
    t > eps && return t
    t = b + det
    t > eps && return t
    return 0
end

function intersect(scene, r::Ray)
    inf = 1e20
    t = inf
    id = 0
    for i = 1:length(scene)
        d = intersect(scene[i], r)
        if d != 0 && d < t
            t = d
            id = i
        end
    end
    (t < inf, t, id)
end

clamp01(x) = clamp(x, zero(x), one(x))

function radiance(scene, r::Ray, depth::Int, Xi)
    (did_hit, t, id) = intersect(scene, r)
    did_hit || return C[0,0,0]      # if miss, return black
    obj = scene[id]                 # the hit object
    x = r.o + r.d*t
    n = normalize(x-obj.p)
    nl = dot(n, r.d) < 0 ? n : -n
    f = obj.c
    p = maximum(f) # max refl
    depth += 1
    if depth > 5
        if rand(Xi) < p
            f = f*(1/p)
        else
            return obj.e # R.R.
        end
    end
    if obj.refl == DIFF     # Ideal DIFFUSE reflection
        r1 = 2*π*rand(Xi)
        r2 = rand(Xi)
        r2s = sqrt(r2)
        w = nl
        u = normalize((abs(w[1]) > 0.1 ? V[0,1,0] : V[1,0,0]) × w)
        v = w × u
        d = normalize(u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2))
        return obj.e + f .* radiance(scene, Ray(x,d), depth, Xi)
    elseif obj.refl == SPEC # Ideal SPECULAR reflection
        return obj.e + f .* radiance(scene, Ray(x, r.d - n*2*dot(n,r.d)), depth, Xi)
    else # obj.refl == REFR # Ideal dielectric REFRACTION
        reflRay = Ray(x, r.d - n*2*dot(n,r.d))
        into = dot(n, nl) > 0  # Ray from outside going in?
        nc = 1
        nt = 1.5
        nnt = into ? nc/nt : nt/nc
        ddn = dot(r.d, nl)
        cos2t = 1 - nnt^2 * (1 - ddn^2)
        if cos2t < 0 # Total internal reflection
            return obj.e + f .* radiance(scene, reflRay, depth, Xi)
        end
        tdir = normalize(r.d*nnt - n*((into ? 1 : -1)*(ddn*nnt+sqrt(cos2t))))
        a = nt - nc
        b = nt + nc
        R0 = a^2/b^2
        c = 1 - (into ? -ddn : dot(tdir, n))
        Re = R0 + (1-R0)*c^5
        Tr = 1-Re
        P = 0.25 + 0.5*Re
        RP = Re/P
        TP = Tr/(1-P)
        return obj.e + f .* (depth>2 ?
            (rand(Xi)<P ? radiance(scene, reflRay,depth,Xi)*RP : radiance(scene, Ray(x,tdir),depth,Xi)*TP) :
            radiance(scene, reflRay,depth,Xi)*Re + radiance(scene, Ray(x,tdir),depth,Xi)*Tr
        )
    end
end

function render(scene, w, h, samps)
    subpix_samps = samps ÷ 4
    cam = Ray(V[50,52,295.6], normalize(V[0,-0.042612,-1])) # cam pos, dir
    cx = 0.5135 * V[w/h, 0, 0]
    cy = 0.5135 * normalize(cx × cam.d)
    c = zeros(CVec, w, h)
    Xi = Random.GLOBAL_RNG # TODO: row-based seeding for multithreading
    rowcount = Threads.Atomic{Int}(0)
    for y = 1:h # Loop over image rows
        rowcount[] += 1
        #@info "Rendering ($samps spp)" progress=(rowcount[])/h
        for x = 1:w # Loop cols
            # Xi = Random.MersenneTwister(y^3)
            for sy = 0:1      # 2x2 subpixel rows
                for sx = 0:1  # 2x2 subpixel cols
                    r = V[0,0,0]
                    for _ = 1:subpix_samps
                        # Bilinear (tent) filter for subpixels
                        r1 = 2*rand(Xi)
                        r2 = 2*rand(Xi)
                        dx = r1 < 1 ? sqrt(r1)-1 : 1-sqrt(2-r1)
                        dy = r2 < 1 ? sqrt(r2)-1 : 1-sqrt(2-r2)
                        d = normalize(cx*(((sx + 0.5 + dx)/2 + x-1)/w - 0.5) +
                                      cy*(((sy + 0.5 + dy)/2 + y-1)/h - 0.5) + cam.d)
                        r += radiance(scene, Ray(cam.o + d*140, d), 0, Xi)
                        # Camera rays are pushed ^^^^^ forward to start in interior
                    end
                    c[x,y] += 0.25 * clamp01.((1.00/subpix_samps) .* r)
                end 
            end
        end
    end
    c
end

to_sRGB(v) = RGB(clamp01.(v).^(1/2.2)...)

function main()
    scene = cornell_box_scene()
    c = render(scene, 128, 96, 128)
    #c = render(scene, 512, 384, 256)
    #c = render(scene, 1024, 768, 4)
    #c = render(scene, 2, 2, 4)
    img = to_sRGB.(permutedims(reverse(c, dims=2)))
    # PNGFiles.save("image.png", img, filters=0)
    img
end

function cornell_box_scene()
    # radius, position, emission, color, material
    [
      Sphere(1e5,  V[ 1e5+1,40.8,81.6],  C[0,0,0],    C[.75,.25,.25], DIFF), #Left
      Sphere(1e5,  V[-1e5+99,40.8,81.6], C[0,0,0],    C[.25,.25,.75], DIFF), #Rght
      Sphere(1e5,  V[50,40.8, 1e5],      C[0,0,0],    C[.75,.75,.75], DIFF), #Back
      Sphere(1e5,  V[50,40.8,-1e5+170],  C[0,0,0],    C[0,0,0],       DIFF), #Frnt
      Sphere(1e5,  V[50, 1e5, 81.6],     C[0,0,0],    C[.75,.75,.75], DIFF), #Botm
      Sphere(1e5,  V[50,-1e5+81.6,81.6], C[0,0,0],    C[.75,.75,.75], DIFF), #Top
      Sphere(16.5, V[27,16.5,47],        C[0,0,0],    C[1,1,1]*.999,  SPEC), #Mirr
      Sphere(16.5, V[73,16.5,78],        C[0,0,0],    C[1,1,1]*.999,  REFR), #Glas
      Sphere(600,  V[50,681.6-.27,81.6], C[12,12,12], C[0,0,0],       DIFF), #Lite
      # Sphere( 4.5, V[50,30,62],          C[12,12,12],    C[0,0,0],       DIFF), #Light
    ]
end

end
