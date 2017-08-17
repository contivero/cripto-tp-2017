Kuang-Shyr Tsung-Ming Secret Sharing Scheme
====
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

Simple implementation of a secret image sharing scheme, as described in
[Kuang-Shyr Wu, Tsung-Ming Lo. An Efficient Secret Image Sharing Scheme. Applied Mechanics and Materials Vols. 284-287 (2013) pp 3025-3029]

The program hides 8-bit BMP images inside others. The 40 byte BITMAPINFOHEADER
format is assumed.

To build simply use `make`, the different flags can be found in `config.mk`

usage:

```
bmpsss (-d|-r) -secret <image> -k <number> -w <width> -h <height> [-s <seed>] [-n <number>] [-dir <directory>]

-d                  distribute image by hiding it on others
-r                  recover image hidden in others
-secret <image>     if -d was specified, image is the file name of the BMP file
                    to hide. Otherwise (if -r was specified), output file name
                    with the revealed  image.
-w <width>          width of the image to recover
-h <height>         height of the image to recover
-s <seed>           seed for the permutation. If non specified, uses 691.
-n <number>         amount of files in which to distribute the image. If not
                    specified, uses the total amount of files in the directory
-dir <directory>    directory in which to search for the images. If not
                    specified, use the current directory.
```

The paper was given to be used for the implementation project of the 2017
Cryptography and Security course at ITBA.
