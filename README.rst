Pytheas: an open-source software solution for local shear-wave splitting studies
================================================================================

.. image:: https://zenodo.org/badge/181635686.svg
   :target: https://zenodo.org/badge/latestdoi/181635686

**IMPORTANT**
While version [0.3.0+d14] offers new functionality and a much much stabler Catalogue CA scheme, it is not yet fully test/validated and some features (mainly GUI-related), as well as proper documentation, is yet to be implemented. However, providing a way to analyze large datasets reliably is in the core of the program, so we decided to release this version (see also the relevant commit comments). Let us know if you have any questions/inquiries in the e-mail below! Please make sure you read comments of the latest commits in the devel branch, to be informed on known bugs.

Pytheas is a tool that aims to introduce a new mentality in shear-wave splitting analysis from local recordings, incorporating manual, semi- and fully- automatic methods under one Graphical User Interface. Integrating databases and offering compatibility with popular data and metadata file formats, Pytheas is designed with the simplification of analysis in mind, while, at the same time, enhanching the effectiveness of processing and quality control of results.

Along with the program itself, you can find the following files:

* *docs:* Contains a detailed User's Guide as well as other relevant documentation.
* *acquisition_Kscripts:* Scripts for acquiring event data and metadata, as well as station metadata from EIDA and IRIS FDSN services.
* *example:* Sample waveforms, catalogues, station file and velocity model for trying out the program.

Pytheas is released under the GNU GPLv3 license.

Authors: Spingos I. & Kaviris G. (c) 2019-2024

Special thanks to Millas C. for testing the software and providing valuable feedback from the very early stages of this endeavor!

For any issues, comments or suggestions please contact us at pytheas.splitting@gmail.com or through `GitHub <https://www.github.com/ispingos/pytheas-splitting>`_.


**Demonstrations**

Pytheas on `YouTube <https://www.youtube.com/channel/UC7USfZT9PfnNTNqMiY1AgTg>`_

Short hands-on lecture about Pytheas in `CRL School 2020 <https://www.youtube.com/watch?v=cUB5qNdUFh0>`_.


**Installation**

There are two scripts (`install_reqs_linux.sh` and `install_reqs_windows.bat`) that install requirements through PIP (`requirements.txt`). This means that you *need* to have Python installed (latest code has been tested with Python 3.10).

If you want to create an environment (e.g. through `venv`, `virtualenv` or conda), you will have to **activate** it, install the requirements through the script or requirements file and **then** start the software through the terminal with the activated environment.

**How to cite**

   Spingos, I., Kaviris, G., Millas, C., Papadimitriou, P., Voulgaris, N., 2020. 
   Pytheas: An open-source software solution for local shear-wave splitting studies. Comput. Geosci. 134, 104346. 
   doi: 10.1016/j.cageo.2019.104346

The published article can be accessed through `Elsevier <https://www.sciencedirect.com/science/article/pii/S0098300419303784>`_.
