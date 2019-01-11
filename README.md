# CGI-flux-ratio-plot
This contains the python code to plot the flux ratio vs separation plot. It also contains the text files with the contrast curves and planet flux ratio values. Each datafile has a full description of its comment section; a brief explanation (auto_caption.txt) is auto-generated by the script from these headers, including only those data you choose plot.

## python requirements
* python 2
* matplotlib
* astropy (tables, ascii, units)
* numpy
* yaml


## quickstart

See the `documentation` folder for an example plot and a complete description of the data and calculations.

If you would like to contribute to this repo, please [fork it](https://help.github.com/articles/fork-a-repo/). Otherwise, simply download or clone and you're ready to go.

Copy `example_files/default.yaml`  to the same directory as `plot_flux_ratio.py`.
Run `python plot_flux_ratio.py`
Your figure and an auto-generated description will be created in the "output" folder.


## Customizing your plot

### YAML file
To customize what is plotted, place a YAML configuration file in the same directory as the script plot_flux_ratio.py. See the "example_files" directory for an example.  The name of the config file will be appended to the output plot & caption filenames, so it's recommended to give your config file a descriptive name.

### Plot options

Available WFIRST modes:
* narrow FOV imaging requirement
* wide FOV imaging requirement
* spectroscopy requirement
* imaging & spectroscopy perfomance predictions

Available instruments:
* Gemini GPI IFS
* VLT SPHERE: IFS & IRDIS
* Magellan VisAO
* HST ACS
* HST NICMOS
* HST STIS
* JWST NIRCam

Available planets:
* known self-luminious directly imaged planets: measured H-band contrasts & predictions for CGI bandpasses.
* known RV planets: predicted reflected light brightness
* Earth & Jupiter at 10pc
* Tau Ceti e & f

Other available plotting options:
* color coding by bandpass: none, minimal color, simple color, full color

### Data file format

Datafiles should be formatted in the following way to allow for auto import by astropy.io.ascii & auto generation of caption.txt:
```
#col_name1  col_name2
#Short caption: Write a short explanation, suitable for a figure caption. This line is automatically fed into "auto_caption.txt" when contrast_plot.py is run. This line must begin with the phrase "short caption:" & be all on one line (ie: no "return" characters).
#Full explanation: Optional. Any other info that you don't want to appear in auto_caption.
#References
#URLs or paper references or other source information
data11  data12
data21	data22
...
```

### Submitting new data or features

Thanks for your help! Before you submit a pull request please:
* Ensure that the script still runs with doc.yml and default.yml
* If your changes affect flux_ratio_doc.png, place your new version in the Documentation folder.
* Update description.md, if necessary.
* Update the short caption header in any changed data files.
