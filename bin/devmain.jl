using Revise

using PathTracing

using ImageView, Gtk, GtkReactive

# Create a canvas `c`. There are other approaches, like stealing one from a previous call
# to `imshow`, or using GtkReactive directly.
guidict = imshow_gui((512, 384), name="Path Tracer")
frame = guidict["frame"]
canvas = guidict["canvas"]
# To see anything you have to call `showall` on the window (once)
Gtk.showall(guidict["window"])

# This dance with the frame and zoom region signal `zrsig` is to preserve the
# aspect ratio. Why is this *soo* annoying!? :-( :-(
img = PathTracing.main()
zrsig = Signal(ZoomRegion(img))
imgsig = Signal(img)
imshow(frame, canvas, imgsig, zrsig)

@async Revise.entr([], all=true, postpone=true) do
    global img = PathTracing.main()
    push!(imgsig, img)
    push!(zrsig, ZoomRegion(img))
end

