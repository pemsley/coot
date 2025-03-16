# *Coot*

![Ubuntu binary build badge](https://github.com/pemsley/coot/actions/workflows/build-coot-ubuntu.yml/badge.svg)
![chapi binary build badge](https://github.com/pemsley/coot/actions/workflows/build-libcootapi-ubuntu.yml/badge.svg)
![macOS binary build badge](https://github.com/pemsley/coot/actions/workflows/build-coot-macos.yml/badge.svg)
[![Powered by RDKit](https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQBAMAAADt3eJSAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAFVBMVEXc3NwUFP8UPP9kZP+MjP+0tP////9ZXZotAAAAAXRSTlMAQObYZgAAAAFiS0dEBmFmuH0AAAAHdElNRQfmAwsPGi+MyC9RAAAAQElEQVQI12NgQABGQUEBMENISUkRLKBsbGwEEhIyBgJFsICLC0iIUdnExcUZwnANQWfApKCK4doRBsKtQFgKAQC5Ww1JEHSEkAAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAyMi0wMy0xMVQxNToyNjo0NyswMDowMDzr2J4AAAAldEVYdGRhdGU6bW9kaWZ5ADIwMjItMDMtMTFUMTU6MjY6NDcrMDA6MDBNtmAiAAAAAElFTkSuQmCC)](https://www.rdkit.org/)
[![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/downloads/release/python-3129/)

*Coot* is a toolkit for Macromolecular Crystallography and
model-building.  *Coot* uses GTK widgets
mmdb, clipper, and OpenGL to provide interactive tools for model-building,
refinement and validation.

## Installation

### Homebrew (MacOS / Linux)

After [installing Homebrew](https://brew.sh/), run the commands below.

```shell
brew install brewsci/bio/coot
coot
```

To install the latest development version (bleeding edge), run the commands below.

```shell
brew install --HEAD brewsci/bio/coot
coot
```

### APT (Debian)

Install to run using commands below.

```shell
# Add sid (Debian's unstable package repo)
echo "deb http://deb.debian.org/debian unstable main" | sudo tee -a /etc/apt/sources.list
sudo apt update
sudo apt install coot
coot
```

### Flatpak (Linux)

![Flathub version](https://img.shields.io/flathub/v/io.github.pemsley.coot.svg?logo=flatpak&logoColor=white&color=blue&style=flat)

<div>
  <a href="https://flathub.org/apps/io.github.pemsley.coot">
    <img width='240' alt='Download on Flathub' src='https://dl.flathub.org/assets/badges/flathub-badge-en.svg' align="right"/>
  </a>
</div>

After [installing Flatpak and registering Flathub](https://flatpak.org/setup/), run the commands below.

```shell
flatpak install flathub io.github.pemsley.coot

# Simply click Coot's icon in the menu, or
flatpak run io.github.pemsley.coot
```

To install the beta build (more recent but possibly less stable), [register flahub-beta repository](https://github.com/flathub/io.github.pemsley.coot?tab=readme-ov-file#beta-build) and then run the commands below.

```shell
flatpak install flathub-beta io.github.pemsley.coot

# Simply click Coot's icon in the menu, or
flatpak run io.github.pemsley.coot//beta
```

You can see the flatpak manifest [here](https://github.com/flathub/io.github.pemsley.coot).

### Windows

Download the installer for WinCoot from [here](https://bernhardcl.github.io/coot/wincoot-download.html).

### As a component of CCP4 Suite (Linux / MacOS / Windows)

Coot-1 can be installed by using the package manager included in the CCP4 suite.

CCP4 Suite is can be downloaded from [CCP4 Download pages](https://www.ccp4.ac.uk/download/).

## Building from source

See [this](https://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/web/build-install-coot-from-scratch.html).

## Blog

[Coot Development Blog](https://pemsley.github.io/coot/ "Coot Development Blog")
