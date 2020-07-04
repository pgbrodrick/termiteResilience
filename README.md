# termiteResilience
This repository stores code written to perform the analyses in 
 
>[A.B. Davies, P.G. Brodrick, C.L. Parr, and G.P. Asner. Resistance of mound-buildingtermites to anthropogenic land-use change.Environmental Research Letters, 2020. doi:10.1088/1748-9326/aba0ff.](10.1088/1748-9326/aba0ff)

Code here requires mapped termite mounds, which were generated using
the code package [bfgn](https://github.com/pgbrodrick/bfg-nets/issues). 
Resulting mapped mound locations are provided as shapefiles in the directoy 
full_landscape_mound_predictions.  Use of these identified mound locations
should cite the manuscript above.

To reproduce the figures from the manuscript, from pre-extracted data
and saved bootstrap subsets, use the following two calls:

> python main_text_plots.py
>
> python si_plots.py.py
