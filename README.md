# Gaia-microlensing-errors

last update: 24/05/2021

First part of my Master's Thesis

This script was designed as a workhorse to determine if the problem I was researching was even worth pursuing so a lot of the code here needs to be rewritten to be more efficient and clear. All of it were done before Nov 2020 and last update was launched only to draw graphs for my Mater Thesis. The thesis right now is available only in Polish, and you can get it by contacting me, but hopefully will be available in the form of a research paper soon.

It has 4 main purposes which you can choose when launching the code:
1. Draw positions of source, lens and blend of those two (blend based on lever rule) for blending parameter taken from standard input. It generates positions to a .dat file, uses the data to draw a graph and then saves everything to a file with name based on the blending parameter.
2. Second function tries to determine the parallax of source, lens, blend and prediction of blends parallax based on a relation I was looking for. It was designed to help me check how well my equation predicts parallax of the blend.
3. This option does basically the same thing as option 1 for a static set of parameters.
4. Last option draws a relation of parallax and blending parameter.

Whole script requires numpy and matplotlib to run apart from standard set of libraries.
