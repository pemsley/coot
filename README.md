# *Coot*

![Ubuntu binary build badge](https://github.com/pemsley/coot/actions/workflows/build-coot-ubuntu.yml/badge.svg)

![chapi binary build badge](https://github.com/pemsley/coot/actions/workflows/build-libcootapi-ubuntu.yml/badge.svg)

![macOS binary build badge](https://github.com/pemsley/coot/actions/workflows/build-coot-macos.yml/badge.svg)

*Coot* is a toolkit for Macromolecular Crystallography and
model-building.  *Coot* uses GTK widgets 
mmdb, clipper, and OpenGL to provide interactive tools for model-building,
refinement and validation.

## Installing using package managers

The simplest way to install *Coot* is using package managers.

### Homebrew (for Mac) / Linuxbrew (for Linux)

After [installing Homebrew or Linuxbrew](https://brew.sh/) run the commands below.

```shell
brew install brewsci/bio/coot
coot
```

### Flatpak (for Linux)

<a href="https://flathub.org/apps/io.github.pemsley.coot"><img width='240' alt='Download on Flathub' src='https://dl.flathub.org/assets/badges/flathub-badge-en.svg' align="right"/></a>

After [installing Flatpak and registering Flathub](https://flatpak.org/setup/) run the commands below.

```shell
flatpak install flathub io.github.pemsley.coot
flatpak run flathub io.github.pemsley.coot  # or simply click Coot's icon in the menu
```

## Building from source

See [this](https://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/web/build-install-coot-from-scratch.html).

## Blog

[Coot Development Blog](https://pemsley.github.io/coot/ "Coot Development Blog")

## Badges

[![Powered by RDKit](https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQBAMAAADt3eJSAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAFVBMVEXc3NwUFP8UPP9kZP+MjP+0tP////9ZXZotAAAAAXRSTlMAQObYZgAAAAFiS0dEBmFmuH0AAAAHdElNRQfmAwsPGi+MyC9RAAAAQElEQVQI12NgQABGQUEBMENISUkRLKBsbGwEEhIyBgJFsICLC0iIUdnExcUZwnANQWfApKCK4doRBsKtQFgKAQC5Ww1JEHSEkAAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAyMi0wMy0xMVQxNToyNjo0NyswMDowMDzr2J4AAAAldEVYdGRhdGU6bW9kaWZ5ADIwMjItMDMtMTFUMTU6MjY6NDcrMDA6MDBNtmAiAAAAAElFTkSuQmCC)](https://www.rdkit.org/)

[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/downloads/release/python-3114/)
