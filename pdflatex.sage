#!/usr/bin/env python
# encoding: utf-8
import os


def pp(texCode, mode=1):
    """@todo: Docstring for pp.
    pp is 'pretty print'. It is a wrapper to the latex
    command, that additionally writes it to a fixed tex file
    and compiles it. That pdf file can remain open in a pdf
    viewer which support automatically reload. This gives
    sort of a "live view" of the whatever sage objects you
    want, from the shell.
    :returns: Exit status of the pdflatex command.
    """

    texFileOut = "output/sage_output.tex"

    if mode == 1:
        head = '\\documentclass[a4paper]{article}\n\\usepackage{breqn}\n\\begin{document}\n\\begin{dmath*}\n'
        foot = '\n\\end{dmath*}\n\\end{document}'
    if mode == 2:
        head = '\\documentclass{standalone}\n\\begin{document}\n\\begin{equation}\n'
        foot = '\n\\end{equation}\n\\end{document}'

    f = open(texFileOut, "w")
    f.write(head + latex(texCode) + foot)
    f.close()

    r = os.system('pdflatex ' + texFileOut + ' > /dev/null ')

    return r
