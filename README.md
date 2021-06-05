# Terminal 3D Renderer
A 3D renderer powered by [Notcurses](https://github.com/dankamongmen/notcurses).
As opposed to many other terminal 3D rendering programs, this uses the ability
of some terminals to accept and draw pixels directly to the screen, making it
possible to draw without compromising on resolution.

![Demo 1 (gif)](https://masflam.com/static/rend3d-demo-1.gif)

![Demo 2 (video)](https://masflam.com/static/rend3d-demo-2.mp4)

# Building
Make sure Notcurses is installed (on Debian it's `libnotcurses-dev`). Then run `make`
to compile and link the `rend3d-demo` executable.
